% IPOPT_SOLVER Interface for calling IPOPT to optimize Fletcher penalty
classdef ipopt_solver < matlab.mixin.Copyable
% Minimizes the penalty function using IPOPT as subsolver.
%
% Not particularly advised for minimizing the penalty function due to the use
% of linesearches, and lack of second-derivative support because the Hessian
% is not explicit.
%
% IPOPT can be obtained from here:
%      https://projects.coin-or.org/Ipopt
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
      x0              % initial point
  
      funcs = struct();    % callback functions for ipopt
      options = struct();  % options struct for ipopt
  
      new_sigma       % updated value of sigma
      new_delta       % updated value of delta

      prev_x          % previous iterate

      logH = ''
      logB = ''
      logT = {}; 
   end
    
   methods
        
      function self = ipopt_solver(fletcher_solver, varargin)

         p = inputParser;
         p.KeepUnmatched = false;
         p.addParameter('mu_strategy', 'adaptive');
         p.addParameter('x0', fletcher_solver.mf.x0);
         p.addParameter('start_feasible', true);
         p.parse(varargin{:});
         
         mu_strategy = p.Results.mu_strategy;
         self.x0 = p.Results.x0;
         start_feasible = p.Results.start_feasible;
         
         if start_feasible && size(fletcher_solver.mf.linear,1) > 0
            % Redefine x0 as the least-squares projection onto the
            % linear constraints
             
            % Assume no bound constraints for now, will need actual
            % solver for that
            if sum(fletcher_solver.mf.jFre) ~= fletcher_solver.mf.n
               ME = MException('FLETCHER_SOLVER:IPOPT_SOLVER',...
                    'Cannot project onto linear constraints with bounds');
               throw(ME);
            end
            mfl = fletcher_solver.mf;
            % Merit function problem contains only linear constraints
            A = mfl.gcon(self.x0);
            r = mfl.fcon(self.x0) - mfl.cL;
            z = spqr_solve(A, -r, struct('solution', 'min2norm'));
            self.x0 = self.x0 + z;
         end
         
         % The callback functions.
         self.funcs.objective         = @(x) fletcher_solver.mf.fobj(x);
         self.funcs.gradient          = @(x) fletcher_solver.mf.gobj(x);
         self.funcs.constraints       = @(x) fletcher_solver.mf.fcon(x);
         self.funcs.jacobian          = @(x) fletcher_solver.mf.gcon(x);
         self.funcs.jacobianstructure = @( ) fletcher_solver.mf.Jpattern;
         self.funcs.iterfunc          = @(t, f, s) self.post_iteration(t, f, s);

         % Set the IPOPT options.
         self.options.cl = fletcher_solver.mf.cL;
         self.options.cu = fletcher_solver.mf.cU;
         self.options.lb = fletcher_solver.mf.bL;
         self.options.ub = fletcher_solver.mf.bU;
         self.options.ipopt.print_level = 0;
         self.options.ipopt.mu_strategy = mu_strategy;
         self.options.ipopt.max_iter    = 1e6; % Supersolver controls max_iter
         self.options.ipopt.tol         = 1e-32; % ipopt won't exit on its own
         self.options.ipopt.hessian_approximation = 'limited-memory';
         
         self.fletcher_solver = fletcher_solver;
      end

      function [fletcher, info] = solve(self)
         
         x = self.x0;
         self.prev_x = x;
         
         while true
             
            [x, info] = ipopt(x,self.funcs,self.options);
            
            if info.status == 5 && self.fletcher_solver.exit == 0
               % User requested stop (i.e. supersolver)
               % need to update values
               self.fletcher_solver.mf = ...
                  set_sigma(self.fletcher_solver.mf, self.new_sigma);
               self.fletcher_solver.mf.delta = self.new_delta;
            elseif self.fletcher_solver.exit == 0
               % Ipopt quit for some reason
               info.exit = self.fletcher_solver.EXIT_SUBSOLVER;
               info.exit_msg = self.getMessage(info.status);
               break;
            else
               info.exit = self.fletcher_solver.exit;
               info.exit_msg = ...
                  self.fletcher_solver.EXIT_MSG{self.fletcher_solver.exit};
               break;
            end
         end

         info.sol.x = self.fletcher_solver.last_x;
         self.fletcher_solver.mf.fobj(info.sol.x);
         info.sol.y = self.fletcher_solver.mf.val.y;
         info.sol.f = self.fletcher_solver.mf.val.f;
         
         fletcher = self.fletcher_solver;
      end
           
      function b = post_iteration(self, t, f, s)
         % Prepare the additional logging
         logstr = '';
         x = s.x;

         % Sometimes the current iterate has size zero
         % This might be due to Ipopt going into the restoration phase
         % (I think I read that somewhere).
         if size(x) == 0
            b = true;
            return;
         end
         
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
            b = true;
            return;
         end
         
         if info.flag == self.fletcher_solver.SUBSOLVER_USER
            % Supersolver wants to quit for some reason
            b = false;
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
            b = false;
            return;
         end
         
         % No known flag was triggered, throw error?
         ME = MException('IPOPT_SOLVER:POSTITERATION', ...
             'Unknown flag encountered');
         throw(ME);
      end

      function str = logHeader(self)
         str = sprintf(self.logH, self.logT{:});
      end 

      function msg = getMessage(self, flag)
         switch flag
            case 0,    msg = 'solved';
            case 1,    msg = 'solved to acceptable level';
            case 2,    msg = 'infeasible problem detected';
            case 3,    msg = 'search direction becomes too small';
            case 4,    msg = 'diverging iterates';
            case 5,    msg = 'user requested stop';
            case -1,   msg = 'maximum number of iterations exceeded';
            case -2,   msg = 'restoration phase failed';
            case -3,   msg = 'error in step computation';
            case -10,  msg = 'not enough degrees of freedom';
            case -11,  msg = 'invalid problem definition';
            case -12,  msg = 'invalid option';
            case -13,  msg = 'invalid number detected';
            case -100, msg = 'unrecoverable exception';
            case -101, msg = 'non-ipopt exception thrown';
            case -102, msg = 'insufficient memory';
            case -103, msg = 'internal error';
            otherwise, msg = 'unknown error';
         end
      end
        
   end % methods
    
end % classdef