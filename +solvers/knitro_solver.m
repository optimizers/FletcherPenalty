% KNITRO_SOLVER Minimizes the penalty function with Knitro
classdef knitro_solver < matlab.mixin.Copyable
% Minimizes the penalty function using KNITRO as subsolver.
%
% Knitro can be obtained from here:
%      https://www.artelys.com/en/optimization-tools/knitro
%
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
      fletcher_solver % handle to fletcher_solver object
      x0              % Initial point
      opt_file        % File for options

      options = struct();  % Options struct for knitro

      new_sigma       % Updated value for sigma
      new_delta       % Updated value for delta

      prev_x          % Previous iterate

      logH = ''
      logB = ''
      logT = {}; 
   end
    
   methods
     
      function self = knitro_solver(fletcher_solver, varargin)
         % KNITRO_SOLVER Constructor
         %
         % Required Inputs:
         %    fletcher_solver   pointer to supersolver
         %
         % Optional Inputs:
         %    x0                initial point
         %    opt_file          filename for KNITRO options file
         p = inputParser;
         p.KeepUnmatched = false;
         p.addParameter('x0', fletcher_solver.mf.x0);
         p.addParameter('opt_file', []);
         p.parse(varargin{:});
         
         self.x0 = p.Results.x0;
         self.opt_file = p.Results.opt_file;
         n = length(self.x0);
         
         self.fletcher_solver = fletcher_solver;
         
         self.options = optimset('Algorithm', 'interior-point', ...
             'GradObj', 'on', 'GradConstr', 'on', ...
             'AlwaysHonorConstraints', 'bounds', ...
             'TolCon', 1e-16, 'TolFun', 1e-16, 'TolX', 1e-16, ... % Supersolver controls exit
             'SubproblemAlgorithm', 'cg', ...
             'JacobPattern', zeros(0,n), ...
             'Hessian', 'user-supplied', ...
             'HessMult', @(x,y,s) self.fletcher_solver.mf.hlagprod(x, y.eqlin, s), ...
             'MaxIter', 1e6, ... % Supersolver controls number of iterations
             'OutputFcn', @(x, optimValues, state) self.post_iteration(x, optimValues, state),...
             'Display', 'off');
      end
     
      function [fletcher, info] = solve(self)
         
         info = struct();
         x = self.x0;
         self.prev_x = x;
         n = length(x);
         
         fun = @(x) self.objgrad(x);
         nonlcon = @(x) self.nonlcon(x);
         A = [];
         b = [];
         Aeq = self.fletcher_solver.mf.gcon(x);           
         beq = self.fletcher_solver.mf.cL ...
               - self.fletcher_solver.mf.fcon(zeros(n,1));
         lb = self.fletcher_solver.mf.bL;
         ub = self.fletcher_solver.mf.bU;
         
         extendedFeatures = [];
         
         while true
             
            [x,fval,exitflag,output,lambda,grad,hessian] = ...
                knitromatlab(fun,self.x0,A,b,Aeq,beq,lb,ub,...
                    nonlcon,extendedFeatures,self.options,self.opt_file);
            
            if exitflag == -504 && self.fletcher_solver.exit == 0
               % User requested stop (i.e. supersolver)
               % need to update values
               self.fletcher_solver.mf = ...
                  set_sigma(self.fletcher_solver.mf, self.new_sigma);
               self.fletcher_solver.mf.delta = self.new_delta;
               self.x0 = x;
            elseif self.fletcher_solver.exit == 0
               % Knitro quit for some reason
               info.exit = self.fletcher_solver.EXIT_SUBSOLVER;
               info.exit_msg = self.getMessage(exitflag);
               break;
            else
               info.exit = self.fletcher_solver.exit;
               info.exit_msg = ...
                  self.fletcher_solver.EXIT_MSG{self.fletcher_solver.exit};
               break;
            end
         end
         
         info.sol.x = x;
         self.fletcher_solver.mf.fobj(info.sol.x); % Evaluate one final time
         info.sol.y = self.fletcher_solver.mf.val.y;
         info.sol.f = self.fletcher_solver.mf.val.f;
         
         fletcher = self.fletcher_solver;
      end
     
      function b = post_iteration(self, x, optimValues, state)
         % Prepare the additional logging
         logstr = '';
         
         % Check whether the iterate has changed since the previous one
         if norm(self.prev_x - x) == 0
            successful = false;
         else
            successful = true;
            self.prev_x = x;
         end
         
         % Call the generic post_iteration hook
         [self.fletcher_solver, info] = ...
            self.fletcher_solver.post_iteration(x, successful, ...
            self.logHeader(), logstr);
         
         % Adjust problem according to info
         if info.flag == self.fletcher_solver.SUBSOLVER_NONE
            % Do nothing, just keep going
            b = false;
            return;
         end
         
         if info.flag == self.fletcher_solver.SUBSOLVER_USER
            % Supersolver wants to quit for some reason
            b = true;
            return;
         end
         
         if info.flag == self.fletcher_solver.SUBSOLVER_MODIFY
            if info.sigma >= 0
               self.new_sigma = info.sigma;
            else
               self.new_sigma = self.fletcher_solver.mf.sigma;
            end

            if info.delta >= 0
               % Need to update delta
               self.new_delta = info.delta;
            else
               self.new_delta = self.fletcher_solver.mf.delta;
            end
            b = true;
            return;
         end
         
         % No known flag was triggered, throw error?
         ME = MException('KNITRO_SOLVER:POSTITERATION',...
              'Unknown flag encountered');
         throw(ME);
      end
     
      function [f,g] = objgrad(self, x)
         f = self.fletcher_solver.mf.fobj(x);
         g = self.fletcher_solver.mf.gobj(x);
      end
     
      function [c,ceq,J,Jeq] = nonlcon(self, x)
         n = length(x);
         c = zeros(0,n);
         ceq = zeros(0,n);
         J = zeros(0,n);
         Jeq = zeros(0,n);
      end
     
      function str = logHeader(self)
         str = sprintf(self.logH, self.logT{:});
      end 
     
      function msg = getMessage(self, flag)
         switch flag
            case 0,    msg = 'solved';
            case -100, msg = 'Current feasible solution estimate cannot be improved. Nearly optimal';
            case -101, msg = 'Relative change in feasible solution estimate < xtol';
            case -102, msg = 'Current feasible solution estimate cannot be improved';
            case -103, msg = 'Relative change in feasible objective < ftol for ftol_iters';
            case -200, msg = 'Convergence to infeasible point. Problem may be locally infeasible';
            case -201, msg = 'Relative change in infeasible solution estimate < xtol';
            case -202, msg = 'Current infeasible solution estimate cannot be improved';
            case -203, msg = 'Multistart: No primal feasible point found';
            case -204, msg = 'Problem determined to be infeasible with respect to constraint bounds';
            case -205, msg = 'Problem determined to be infeasible with respect to variable bounds';
            case -300, msg = 'Problem appears to be unbounded';
            case -400, msg = 'Iteration limit reached. Current point is feasible';
            case -401, msg = 'Time limit reached. Current point is feasible';
            case -402, msg = 'Function evaluation limit reached. Current point is feasible';
            case -403, msg = 'MIP: All nodes have been explored. Integer feasible point found';
            case -404, msg = 'MIP: Integer feasible point found';
            case -405, msg = 'MIP: Subproblem solve limit reached. Integer feasible point found';
            case -406, msg = 'MIP: Node limit reached. Integer feasible point found';
            case -410, msg = 'Iteration limit reached. Current point is infeasible';
            case -411, msg = 'Time limit reached. Current point is infeasible';
            case -412, msg = 'Function evaluation limit reached. Current point is infeasible';
            case -413, msg = 'MIP: All nodes have been explored. No integer feasible point found';
            case -415, msg = 'MIP: Subproblem solve limit reached. No integer feasible point found';
            case -416, msg = 'MIP: Node limit reached. No integer feasible point found';
            case -501, msg = 'LP solver error';
            case -502, msg = 'Evaluation error';
            case -503, msg = 'Not enough memory';
            case -504, msg = 'Terminated by user';
            case -523, msg = 'Derivative check failed';
            case -524, msg = 'Derivative check finished';
            case -600, msg = 'Internal Knitro error';
            otherwise, msg = 'unknown error';
         end
      end
     
   end % methods

end % classdef