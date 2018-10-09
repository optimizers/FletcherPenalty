%FLETCHER_SOLVER Solver entry point
classdef fletcher_solver
% Copyright (C) 2018  Ron Estrin, Michael P. Friedlander, 
% Dominique Orban, and Michael A. Saunders.
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
   
   properties
      exit = 0       % exit flag
      sigma_damp = 1 % Damping parameter for penalty parameter
      fid            % file ID for printing output
      log_level       % verbosity level
      mf             % merit function object
      nlp            % original nlp object
      nlnFlag        % indicate if linear constraints treated explicitly
      sigma_max      % penalty parameter max
      sigma_min      % penalty parameter min
      delta_min      % regularization parameter
      exit_msg       % string with exit message
      du_feas        % dual feasibility error
      pr_feas        % primal feasibility error
      du_feas0       % dual feasibility error (initial)
      pr_feas0       % primal feasibility error (initial)
      mf_grad        % projected gradient of the merit function
      sol            % solution structure
      subsolver_info % subsolver's information structure
      last_x         % latest iterate
      iteration      % current iteration number
      max_iterations % maximum number of allowed iterations
      tol            % tolerance structure
      check_grad     % flag if gradient check requested
      sigma_strategy % 'adaptive' or 'fixed'
      delta_dec      % Function controlling the decrease in delta
   end
   
   properties (Access=public, Constant)
   % Options for solving augmented linear system
      LSQ_SNE = 1;      % Seminormal equations
      LSQ_LDL = 2;      % LDL Factorization
      LSQ_MINRES = 3;   % MINRES
      LSQ_LNLQ = 4;     % LNLQ
   end
   
   properties (Access=public, Constant)
   % Options for subsolver to minimize penalty function
      BCFLASH = 1;
      KNITRO  = 2;
      IPOPT   = 3;
      SNOPT   = 4;
      TRPCG   = 5;
   end
   
   properties (SetAccess=private)
      time_prox_point = 0;  % Currently unused
      time_total = 0;
   end
   
   properties (Hidden=true, Constant)
       SUBSOLVER_NONE            = 1; % Do nothing
       SUBSOLVER_USER            = 2; % Subsolver exit because of user
       SUBSOLVER_MODIFY          = 3; % Modify the problem
   end
   
   properties (Hidden=true, Constant)
       SIGMA_ADAPTIVE            = 1;
       SIGMA_FIXED               = 2;
   end
   
   properties (Hidden=true, Constant)
      EXIT_OPTIMAL               = 1;
      EXIT_MERIT_OPTIMAL         = 2;
      EXIT_PENALTY_TOO_LARGE     = 3;
      EXIT_ITERATIONS            = 4;
      EXIT_SUBSOLVER             = 5;
      EXIT_UNKNOWN               = 6;
      EXIT_MSG = {'Optimal solution found'
                  'Merit function sufficiently optimal'
                  'Penalty parameter too large'
                  'Too many iterations'
                  'Early subsolver exit'
                  'Unknown termination condition'};

      % Log header and body formats.
      logH = '\n%5s  %13s  %13s  %9s  %9s  %9s  %9s  %9s  %10s  %7s  %9s  ';
      logB = '%5i  %13.6e  %13.6e  %9.3e  %9.3e  %9.3e  %9.3e  %9.3e  %10s  %7s  %9.3e  ';
      logT = {'iter','merit','objective','opt Error','du  Feas','nln Feas',...
              'lin Feas','||y||','penalty   ','delta  ', '||s||'};
   end
   
   methods
      
      function self = fletcher_solver(nlp, varargin)
         %FLETCHER_SOLVER Entrypoint for solver
         %
         % Required inputs:
         %  nlp      Optimization problem to be minimized
         %
         % Optional inputs:
         %  fid               Save iteration log?
         %  log_level         Show output?
         %  sigma             Initial penalty parameter
         %  sigma_max         Maximum allowable penalty parameter
         %  sigma_min         Minimum allowable penalty parameter
         %  sigma_strategy    SIGMA_ADAPTIVE or SIGMA_FIXED
         %  delta             Regularization parameter
         %  delta_min         Maximum allowable penalty parameter
         %  delta_dec         Function to compute next delta
         %  lsq_options       Struct for options for augmented system solver
         %  lsq_method        Method for solving augmented system
         %  subsolver         Subsolver for minimizing penalty function
         %  subsolver_options Options struct for subsolver
         %  optTolAbs         Absolute tolerance for grad phi or Lagrangian grad
         %  optTolRes         Relative tolerance for grad phi or Lagrangian grad
         %  feaTolAbs         Absolute tolerance for primal feasibility
         %  optTolRel         Relative tolerance for primal feasibility
         %  max_iteraitons    Maximum number of outer iterations
         %  check_grad        Call gradient checker every iteration
         %  merit_options     Options struct for penalty function
         %  x0                Initial point
         
         % ---------------------------------------------------------------------
         % Parse input parameters and initialize local variables.
         % ---------------------------------------------------------------------
         p = inputParser;
         p.KeepUnmatched = false;
         p.addParameter('log_file', '');
         p.addParameter('log_level', 1);

         p.addParameter('sigma', 0);
         p.addParameter('sigma_max', 1e+6);
         p.addParameter('sigma_min', 1e-6);
         p.addParameter('sigma_strategy', fletcher_solver.SIGMA_ADAPTIVE);

         p.addParameter('delta', 1e-8);
         p.addParameter('delta_min', 1e-8);
         p.addParameter('delta_dec', @(d) d/10);
         
         p.addParameter('lsq_options', struct(), @isstruct);
         p.addParameter('lsq_method', fletcher_solver.LSQ_SNE);
         
         p.addParameter('subsolver', fletcher_solver.BCFLASH);
         
         p.addParameter('optTolAbs', 1e-6);
         p.addParameter('feaTolAbs', 1e-6);
         p.addParameter('optTolRel', 1e-6);
         p.addParameter('feaTolRel', 1e-6);
         
         p.addParameter('max_iterations', 500);
         p.addParameter('check_grad', false);
         p.addParameter('subsolver_options', struct(), @isstruct);
         p.addParameter('merit_options', struct(), @isstruct);
         p.addParameter('x0', nlp.x0);
         p.parse(varargin{:});
         
         x0 = p.Results.x0;
         delta = p.Results.delta;
         sigma = p.Results.sigma;
         subsolver_options = p.Results.subsolver_options;
         lsq_options = p.Results.lsq_options;
         lsq_method = p.Results.lsq_method;
         subsolver = p.Results.subsolver;
         self.sigma_strategy = p.Results.sigma_strategy;
         
         % ---------------------------------------------------------------------
         % Ensure that the penalty parameter isn't larger than allowed.
         % ---------------------------------------------------------------------
         if sigma > p.Results.sigma_max
            warning('Clipping sigma (%10.1e) to sigma_max (%10.1e)',...
               sigma, p.Results.sigma_max);
            sigma = p.Results.sigma_max;
         end
         
         % ---------------------------------------------------------------------
         % Subsolver iteration log.
         % ---------------------------------------------------------------------
         self.fid = 1; % If no output log, print to screen
         if ~strcmp(p.Results.log_file, '')
            self.fid = fopen(p.Results.log_file, 'w');
         end
         
         % ---------------------------------------------------------------------
         % Create the least-squares factorizer.
         % ---------------------------------------------------------------------
         switch lsq_method
             case fletcher_solver.LSQ_SNE
                 lsq = least_squares.lssne(nlp.Jpattern', lsq_options);
             case fletcher_solver.LSQ_LDL
                 lsq = least_squares.lsldl(nlp.Jpattern', lsq_options);
             case fletcher_solver.LSQ_MINRES
                 lsq = least_squares.lsminres(nlp.m, nlp.n, lsq_options);
             case fletcher_solver.LSQ_LNLQ
                 lsq = least_squares.lslnlq(nlp.m, nlp.n, lsq_options);
             otherwise
                 ME = MException('FLETCHER_SOLVER:LSQ_METHOD',...
                     'Invalid least-squares method');
                 throw(ME);
         end
         
         % ---------------------------------------------------------------------
         % Create the merit function object, and initialize mult estimate.
         % This estimate is based on the initial sigma -- likely 0, which
         % gives the least-squares multiplier estimate.
         % ---------------------------------------------------------------------
         mfl = fletcher_merit_function(...
            nlp, lsq,...
            'x0', x0,...
            'delta', delta, ...
            'sigma', sigma,...
            p.Results.merit_options);
         y = mfl.val.y;
         
         % ---------------------------------------------------------------------
         % Store various objects and parameters.
         % ---------------------------------------------------------------------
         self.mf = mfl;
         self.nlnFlag = mfl.nlnFlag;
         self.log_level = p.Results.log_level;
         self.nlp = nlp;
         self.sigma_max = p.Results.sigma_max;
         self.sigma_min = p.Results.sigma_min;
         self.delta_min = p.Results.delta_min;
         self.tol.optAbs = p.Results.optTolAbs;
         self.tol.optRel = p.Results.optTolRel;
         self.tol.feaAbs = p.Results.feaTolAbs;
         self.tol.feaRel = p.Results.feaTolRel;
         self.max_iterations = p.Results.max_iterations;
         self.check_grad = p.Results.check_grad;
         self.iteration = 0;         
         self.delta_dec = p.Results.delta_dec;
         self.last_x = mfl.x0;
                  
         % ---------------------------------------------------------------------
         % Start the clock!
         % ---------------------------------------------------------------------
         self.time_total = tic;
         
         % ---------------------------------------------------------------------
         % Determine stopping tests. Set these before penalty parameter is
         % increased; otherwise, the stopping tests can be artificially
         % weakened.
         % ---------------------------------------------------------------------
         gPhi = self.mf.gobj(x0);
         c = self.mf.val.c_nlp;
         g = self.mf.val.g;
         Jtprod = self.mf.val.Jtprod;
                  
         [nlnFeas, linFeas, bndFeas] = self.nlp.prResidual(x0, c);
         
         self.mf_grad  = self.merit_optimality(x0, gPhi);
         self.du_feas  = self.nlp.duResidual(x0, c, g, Jtprod, y);
         self.du_feas0 = self.du_feas;
         self.pr_feas0 = max([nlnFeas  linFeas  bndFeas]);
         
         % The subsolver needs the same dual stopping tolerance.
         stop_tol = self.tol.optAbs + self.tol.optRel*(norm(y,inf)+self.du_feas0);
         
         % ---------------------------------------------------------------------
         % Set penalty parameter a least as large as 2|PHP| or input value.
         % The constant-objective case (g==0) seems to come up frequently.
         % In that case, any positive value of sigma will do, so set it to 1.
         % ---------------------------------------------------------------------
         % MPF (Jul 23, 14): seems to interfere with aircrfta; disable for now.
%          if false % norm(gPhi) == 0
%             sigma = 1; %#ok<UNRCH>
%             
%          else
%             normPHP = self.mf.norm_proj_hess();
%             % normAYt = self.mf.norm_AYt();
%             % sigma = max([1 sigma, normPHP, 2*normAYt]);
%             sigma = max([sigma, self.sigma_min, 2*normPHP]);
%             
%          end
%         sigma = min([sigma, self.sigma_max]);
%         self.mf = set_sigma(self.mf, sigma);

         fPhi = self.mf.fobj(x0); % Evaluate merit func with the new sigma.
         f = self.mf.val.f;       % Grab the corresponding objective value.
         
         % ---------------------------------------------------------------------
         % Check gradients, if requested.
         % ---------------------------------------------------------------------
         if self.check_grad
            self.check_gradient(x0);
         end
                 
         % ----------------------------------------------------------------
         % Construct the subsolver
         % ----------------------------------------------------------------
         switch subsolver
             case fletcher_solver.BCFLASH
                 subsolver = solvers.bcflash_solver(self, ...
                     'subsolver_options', subsolver_options);
             case fletcher_solver.KNITRO
                 subsolver = solvers.knitro_solver(self, ...
                     subsolver_options);                 
             case fletcher_solver.IPOPT
                 subsolver = solvers.ipopt_solver(self, ...
                     subsolver_options);
             case fletcher_solver.SNOPT
                 subsolver = solvers.snopt_solver(self, ...
                     subsolver_options);
             case fletcher_solver.TRPCG
                 subsolver = solvers.trpcg_solver(self, ...
                     subsolver_options);
             otherwise
                 ME = MException('FLETCHER_SOLVER:SUBSOLVER',...
                     'Invalid subsolver selection');
                 throw(ME);
         end
         
         % ----------------------------------------------------------------
         % Print headers.
         % ----------------------------------------------------------------
         if self.log_level
            self.printf('\n');
            self.printf('%s\n',repmat('=',1,80));
            self.printf('FLASH \n');
            self.printf('%s\n\n',repmat('=',1,80));
            self.printf(self.nlp.formatting())
            self.printf('\nParameters\n----------\n')
            self.printf('%-22s: %4s %7.1e'  ,'sigma max','',self.sigma_max);
            self.printf('%10s','');
            self.printf('%-22s: %4s %7.1e\n','cond limit','',self.mf.lsq.cond_max);
            self.printf('%-22s: %4s %7.1e'  ,'delta max','',self.mf.lsq.delta_max);
            self.printf('%10s','');
            self.printf('%-22s: %4s %7i\n'  ,'linear explicit','',self.nlnFlag);
            self.printf('%-22s: %4s %7.1e'  ,'delta min','',self.delta_min);
            self.printf('%10s','');
            self.printf('%-22s: %4s %7i\n'  ,'','','');
            self.printf('%-22s: %4s %7i'    ,'nnz(A)','',nnz(nlp.Jpattern));
            self.printf('%10s','');
            self.printf('%-22s: %4s %7s\n'  ,'','','');
            if lsq_method == fletcher_solver.LSQ_SNE
                self.printf('%-22s: %4s %7i'    ,'nnz(R)','',nnz(self.mf.lsq.Rl));
                self.printf('%10s','');
                self.printf('%-22s: %4s %7s\n'  ,'','','');
            elseif lsq_method == fletcher_solver.LSQ_LDL && self.mf.lsq.static_p
                self.printf('%-22s: %4s %7i'    ,'nnz(L)','',nnz(self.mf.lsq.Ll));
                self.printf('%10s','');
                self.printf('%-22s: %4s %7s\n'  ,'','','');
            end
            self.logHeader(subsolver.logHeader());
            sigstr = sprintf('%7.1e%1s%1s%1s',self.mf.sigma,'','','');
            delstr = sprintf('%7.1e',self.mf.delta);
            self.printf(self.logB, ...
               self.iteration, fPhi, f, self.mf_grad, self.du_feas, ...
               nlnFeas, linFeas, norm(y), sigstr, delstr);
            self.printf('\n');
         end

         % ----------------------------------------------------------------
         % Fire up the solver.
         % ----------------------------------------------------------------
         [self, info] = subsolver.solve();
         
         self.sol.x = info.sol.x;
         self.sol.y = info.sol.y;
         self.sol.f = info.sol.f;
         self.sol.info = info;
         self.exit = info.exit;
         self.exit_msg = info.exit_msg;
         
         % ----------------------------------------------------------------
         % Stop the clock!
         % ----------------------------------------------------------------
         self.time_total = toc(self.time_total);
         
         % ----------------------------------------------------------------
         % Print exit log.
         % ----------------------------------------------------------------
         if self.log_level
            self.printf('\n EXIT: %s\n', self.exit_msg);
            self.printf('\n')
            self.printf(' %-27s  %10i     %-22s  %15.8e\n',...
               'No. of iterations', self.iteration,...
               'Objective value', self.sol.f);
            self.printf(' %-27s  %10i     %-27s  %10i\n',...
               'No. of calls to objective' , self.nlp.ncalls_fobj,...
               'No. of calls to constraint', self.nlp.ncalls_fcon);
            self.printf(' %-27s  %10i     %-27s  %10i\n',...
               'No. of calls to gradient', self.nlp.ncalls_gobj, ...
               'No. of calls to Jacobian', self.nlp.ncalls_gcon);
            self.printf(' %-27s  %10i     %-27s  %10i\n',...
               'No. of Jac-vector prods', self.nlp.ncalls_jprod, ...
               'No. of Adj Jac-vector prods', self.nlp.ncalls_jtprod);
            self.printf(' %-27s  %10i     %-27s  %10.2e\n',...
               'No. of Hessian-vector prods', self.nlp.ncalls_hvp,...
               'Final penalty parameter', self.mf.sigma);
            self.printf(' %-27s  %10.2e     %-27s  %10.2e\n',...
               'Feasibility error (abs)', self.pr_feas,...
               'Optimality error (abs)', self.du_feas);
            self.printf(' %-27s  %10.2e     %-27s  %10.2e\n',...
               'Feasibility error (rel)',self.pr_feas/(1+norm(self.sol.x,inf)+self.pr_feas0),...
               'Optimality error (rel)',self.du_feas/(1+norm(self.sol.y,inf)+self.du_feas0));
            self.printf('\n');
            tt = self.time_total;
            t1 = self.nlp.time_fobj+self.nlp.time_fcon; t1t = round(100*t1/tt);
            t2 = self.nlp.time_gcon+self.nlp.time_gcon; t2t = round(100*t2/tt);            
            self.printf(' %-24s %10.2f (%3d%%)  %-25s %6.2f (%3d%%)\n',...
               'Time for function evals' , t1, t1t,...
               'Time for gradient evals', t2, t2t);
            t1 = self.nlp.time_hvp; t1t = round(100*t1/tt);
            t2 = self.mf.lsq.time;  t2t = round(100*t2/tt);
            self.printf(' %-24s %5.2f (%3d%%)  %-25s %6.2f (%3d%%)\n',...
               'Time for Hessian-vector prods', t1, t1t,...
               'Time for factorize', t2, t2t);
            t1 = self.time_prox_point; t1t = round(100*t1/tt);
            t2 = self.time_total;      t2t = round(100*t2/tt);
            self.printf(' %-24s %10.2f (%3d%%)  %-25s %6.2f (%3d%%)\n',...
               'Time for proximal point', t1, t1t,...
               'Time for solve', t2, t2t);
            t1 = self.nlp.time_ghiv; t1t = round(100*t1/tt);
            t2 = t2/self.nlp.ncalls_ghiv;
            self.printf(' %-24s %10.2f (%3d%%)  %-25s %6.2f\n',...
               'Time for ghiv prods', t1, t1t,...
               'Time/eval for ghiv', t2);
         end
         
         % ---------------------------------------------------------------------
         % Close subsolver iteration log.
         % ---------------------------------------------------------------------
         if self.fid > 1
            fclose(self.fid);
         end
         
      end % constructor

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function [self, info] = post_iteration(self, x, successful, logHeader, log)
         % POST_ITERATION Called at end of every iteration.
         %
         % This gets called at the end of each major iteration of the
         % subsolver being used. The subsolver implements a callback using
         % an interface expected by the subsolver's implementation, and
         % the subsolver's callback acts as a wrapper which calls this
         % function in a standardized form. Thus this function does not
         % care what the underlying subsolver is as long as the subsolver's
         % callback can provide the inputs:
         %
         % Inputs:
         %    x           current iterate
         %
         %    successful  whether the update made to the current iterate
         %    is successful. e.g.: trust-region methods may fail to update
         %    x at a given iteration due to a poor trust-region radius.
         %
         %    log         formatted string of additional values desired to
         %    be logged.
         %
         % The function returns:
         %    self   mutation is weird in Matlab
         %
         %    info   structure containing instructions for the subsolver on
         %    how to proceed after this iteration. Most commonly it will
         %    contain the fields:
         %        - flag   information on convergence, or whether to
         %        increase/decrease sigma/delta.
         %
         %        - val    new value for sigma/delta if they are to be modified
         
         % Check gradient, if requested.
         if self.check_grad
            self.check_gradient(x)
         end

         % Grab various statistics: iteration, merit, and objective vals.
         self.iteration = self.iteration + 1;
         fPhi = self.mf.fobj(x);
         gPhi = self.mf.gobj(x);
         f = self.mf.val.f;
         y = self.mf.val.y;
         c = self.mf.val.c_nlp;
         g = self.mf.val.g;
         Jtprod = self.mf.val.Jtprod;
         sigma = self.mf.sigma;
         c_exp = self.mf.fcon(x);

         % ------------------------------------------------------------
         % Compute primal-dual infeasibility. The subsolver had the right
         % dual infeasibility, and so we'll just use that.
         % Could use gPhiNrm = norm(self.mf.val.gPhi, inf);
         % ------------------------------------------------------------
         pr_feas_prev = self.pr_feas;
         [nlnFeas, linFeas, bndFeas] = self.nlp.prResidual(x, c);
         assert(bndFeas < eps);
         self.mf_grad = self.merit_optimality(x,gPhi);
         self.du_feas = self.nlp.duResidual(x, c, g, Jtprod, y);
         self.pr_feas = max(nlnFeas, linFeas);
         pr_exp_feas  = norm(c_exp - self.mf.cL, 'inf');

         % ------------------------------------------------------------
         % Check various exit flags. First initialize defaults.
         % ------------------------------------------------------------
         info.flag = self.SUBSOLVER_NONE;
         info.delta = -1;
         info.sigma = -1;
         decreased_sigma = false;
         increased_sigma = false;
         decreased_delta = false;

         % Check primal-dual optimality.
         stop_d = self.tol.optAbs + self.tol.optRel*(norm(y,inf)+self.du_feas0);
         stop_p = self.tol.feaAbs + self.tol.feaRel*(norm(x,inf)+self.pr_feas0);

         du_optimal = self.du_feas < stop_d;
         pr_optimal = self.pr_feas < stop_p;
         pr_exp_opt = pr_exp_feas  < stop_p;
         if ~self.exit && (du_optimal && pr_optimal)
            self.exit = self.EXIT_OPTIMAL;
         end

         %-----------------------------------------------------------------
         % If problem is regularized, check if we've minimized enough to
         % decrease delta
         % TODO: Might cause early exit if merit function is optimal but
         % primal/dual feasibility not sufficiently tight
         %-----------------------------------------------------------------
         delta_next = self.delta_dec(self.mf.delta);
         delta_small_enough = self.mf.delta <= self.delta_min;
         if delta_small_enough
             stop_d_del = stop_d;
         else
             stop_d_del = delta_next;
         end
         
         du_optimal_del = self.mf_grad < stop_d_del;
         
         if delta_small_enough && du_optimal_del && pr_exp_opt
             self.exit = self.EXIT_MERIT_OPTIMAL;
         end

         % Check too many iterations.
         if ~self.exit && (self.iteration >= self.max_iterations)
            self.exit = self.EXIT_ITERATIONS;
         end

         % Flag if sigma should change. These tests are mutually exclusive.
         too_infeasible = (100*self.mf_grad < self.pr_feas) ...
             && self.sigma_strategy == self.SIGMA_ADAPTIVE;
         too_feasible   = 100*self.pr_feas < self.mf_grad;

         if ~self.exit && too_feasible && successful ...
                 && self.sigma_strategy == self.SIGMA_ADAPTIVE
            % -----------------------------------------------------------------
            % Decrease the penalty parameter. Try a value of sigma that is
            % the geometric mean of the previous sigma and a damped value
            % of sig_min. Increase the damping term, to ensure that this
            % doesn't happen infinitely many times. If the trial value is
            % significantly smaller (by half, say), than take it, but don't
            % let it get too smaller than the smallest allowed sigma_min.
            % -----------------------------------------------------------------
            sig_dmp = self.sigma_damp;
            sig_min = self.sigma_min;
            sig_nrm = 0.5*self.mf.norm_proj_hess();
            sig_new = sqrt( sigma*(sig_dmp + sig_nrm) );
            sig_new = max([sig_new, sig_nrm, sig_min]);

            if sig_new < 0.5*sigma
               sig_dmp = min(2*sig_dmp, 1/eps^(4/5));
               decreased_sigma = true;
               self.sigma_damp = sig_dmp;
               % Subsolver is responsible for updating sigma
               info.flag = self.SUBSOLVER_MODIFY;
               info.sigma = sig_new;
            end

         end

         if ~self.exit && too_infeasible && successful

            % -----------------------------------------------------------------
            % Too infeasible. Two options: increase the penalty, or
            % decrease the regularization parameter.
            % -----------------------------------------------------------------
            sig_old = self.mf.sigma;            
            normPHP = self.mf.norm_proj_hess();
            sigma = min([max([2*sig_old, 2*normPHP]), self.sigma_max]);

            if sigma >= self.sigma_max
               % No room to increase the sigma parameter.
               % Only option is to decrease delta, if it's not already 0.
               delta = self.mf.delta;
               if delta > self.delta_min
                  % Subsolver is responsible for updating delta
                  % self.mf.delta_max = self.mf.lsq.delta_max / 10;
                  % self.mf.delta = delta / 10;
                  decreased_delta = true;
                  info.flag = self.SUBSOLVER_MODIFY;
                  info.delta = delta / 10;
               elseif self.sigma_strategy == self.SIGMA_ADAPTIVE
                  if sig_old < self.sigma_max
                      % Clip sigma at sigma_max and hope for the best
                      sigma = min([sigma self.sigma_max]);
                      % Subsolver is responsible for updating sigma
                      % self.mf = set_sigma(self.mf, sigma);
                      info.flag = self.SUBSOLVER_MODIFY;
                      info.sigma = sigma;
                      increased_sigma = true;
                  elseif self.pr_feas >= pr_feas_prev
                      % No progress in feasibility with largest possible
                      % penalty, so fail at this point
                      self.exit = self.EXIT_PENALTY_TOO_LARGE;
                  end
                  % If progress was made in feasibility since the previous
                  % iteration, then keep going and hope that eventually
                  % it's not 'too infeasible'
               end

            elseif self.sigma_strategy == self.SIGMA_ADAPTIVE
               % Increase the penalty parameter.
               % Subsolver is responsible for updating sigma
               % self.mf = set_sigma(self.mf, sigma);
               info.flag = self.SUBSOLVER_MODIFY;
               info.sigma = sigma;
               increased_sigma = true;
            end

         end

         if ~self.exit && du_optimal_del && ~too_infeasible
            % -----------------------------------------------------------------
            % Solved problem sufficiently for given delta. Need to decrease
            % delta to continue.
            %
            % TODO: Make sure this plays nice with the too_infeasible case
            % -----------------------------------------------------------------
            delta = max([delta_next self.delta_min]);
            decreased_delta = true;

            if pr_optimal
                % Check if we can just quit now
                % TODO: Make this more efficient
                % Right now Matlab is making a copy of the merit function
                % so that we can test the new problem for convergence
                % without it being passed to the subsolver
                mfl = self.mf;
                mfl.delta = delta;
                fPhi = mfl.fobj(x);
                yy = mfl.val.y;
                self.du_feas = self.nlp.duResidual(x, c, g, Jtprod, yy);
                if self.du_feas < stop_d
                    self.exit = self.EXIT_OPTIMAL;
                    info.delta = delta;
                else
                    info.flag = self.SUBSOLVER_MODIFY;
                    info.delta = delta;
                end
            else
                info.flag = self.SUBSOLVER_MODIFY;
                info.delta = delta;
            end
         end

         % ------------------------------------------------------------
         % Print log.
         % ------------------------------------------------------------
         if self.log_level
            if mod(self.iteration, 20) == 0
               self.logHeader(logHeader)
            end

            if increased_sigma || decreased_sigma
               flag1 = ''; flag2 = ''; flag3 = '';
               if too_infeasible, flag1 = 'i'; end
%               if du_optimal, flag2 = 'd'; end
               if too_feasible, flag3 = char(8600); end
               sigstr = sprintf('%7.1e%1s%1s%1s',info.sigma,flag1,flag2,flag3);
            else
               sigstr = '';
            end

            if decreased_delta
               delstr = sprintf('%7.1e',info.delta);
            else
               delstr = '';
            end

            default_log = sprintf(self.logB, ...
               self.iteration, fPhi, f, self.mf_grad, self.du_feas, ...
               nlnFeas, linFeas, norm(y), sigstr, delstr, norm(self.last_x - x));
           
            self.printf([default_log log]);
           
            self.printf('\n');   
         end

         % Exit flag determines if subsolver exits prematurely.
         self.last_x = x;

         if self.exit
            info.flag = self.SUBSOLVER_USER;
         end

      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function check_gradient(self, x)
         %CHECK_GRADIENT
         [result, eabs, erel] = helpers.gradient_check(self.mf, x);
         if ~result
            self.printf('Gradient check failed : eAbs = %1.3e, eRel = %1.3e\n', eabs, erel);
         end
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function mf_opt = merit_optimality(self, x, gPhi)
         %MERIT_OPTIMALITY  Estimate dual optimality to penalty problem
         %
         % Let B be Jacobian of explicit constraints
         %    y = argmin_y 0.5*|gPhi - B'*y|^2_Q^2
         %    mf_opt = norm(gPhi - B'*y, 'inf')
         c = self.mf.fcon(x);
         B = self.mf.gcon(x);
         if self.mf.m == 0
             y = zeros(0,1);
         else
             % TODO: Can we make this more efficient?
             % For equality only case, can do factorization once and for all
             Q = self.mf.val.Q;
             y = spqr_solve(Q*B', Q*gPhi);
         end
         mf_opt = self.mf.duResidual(x, c, gPhi, B, y);
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function printf(self, varargin)
         fprintf(self.fid, varargin{:});
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function logHeader(self, subsolver_header)
         self.printf(self.logH, self.logT{:});
         self.printf(subsolver_header);
         self.printf('\n');
      end

   end % methods

end % classdef
