% SNOPT_SOLVER Minimizes the penalty function with SNOPT
classdef trpcg_solver < matlab.mixin.Copyable
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
      fid             % File for log output

      options = struct();  % Options struct for ipopt

      new_sigma       % Updated value for sigma
      new_delta       % Updated value for delta

      logH = '%5s  %4s'
      logB = '%5i  %4s'
      logT = {'CGits', 'stat'};
   end
    
   methods
     
      function self = trpcg_solver(fletcher_solver, varargin)
         % SNOPT_SOLVER Constructor
         %
         % Required Inputs:
         %    fletcher_solver   pointer to supersolver
         %
         % Optional Inputs:
         %    x0                initial point
         %    fid               filename for log

         if ~solvers.trpcg_solver.mf_valid(fletcher_solver.mf)
             ME = MException('FLETCHER_SOLVER:TRPCG_SOLVER',...
                  'Problem cannot have bounds.');
             throw(ME);
         end
         
         p = inputParser;
         p.KeepUnmatched = false;
         p.addParameter('fid', '');
         p.parse(varargin{:});
         
         self.fid = p.Results.fid;
         
         self.fletcher_solver = fletcher_solver;
         
         % Set parameters so that TRPCG only successfully exits on user request
         self.options.atol     = 1e-16;
         self.options.rtol     = 1e-16;
         self.options.maxfev   = 1e6;
         self.options.callback = @(x,s,cgits) self.post_iteration(x,s,cgits);
      end
        
      function [fletcher, info] = solve(self)
         
         while true
             
            [x, ~, exitflag, ~] = trpcg(self.fletcher_solver.mf, self.options);
            
            if exitflag == 2 && self.fletcher_solver.exit == 0 
               % User requested stop (i.e. supersolver)
               % need to update values
               self.fletcher_solver.mf = ...
                  set_sigma(self.fletcher_solver.mf, self.new_sigma);
               self.fletcher_solver.mf.delta = self.new_delta;
            elseif self.fletcher_solver.exit == 0
               % Ipopt quit for some reason
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
         self.fletcher_solver.mf.fobj(info.sol.x);
         info.sol.y = self.fletcher_solver.mf.val.y;
         info.sol.f = self.fletcher_solver.mf.val.f;
         
         fletcher = self.fletcher_solver;
      end
        
      function b = post_iteration(self, x, succ, cgits)
         % Prepare the additional logging
         if succ
            succstr = '';
         else
            succstr = 'rej';
         end
         
         logstr = sprintf(self.logB, cgits, succstr);
         
         % Call the generic post_iteration hook
         [self.fletcher_solver, info] = ...
             self.fletcher_solver.post_iteration(x, succ, ...
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
         ME = MException('BCFLASH_SOLVER:POSTITERATION',...
             'Unknown flag encountered');
         throw(ME);
      end

      function str = logHeader(self)
         str = sprintf(self.logH, self.logT{:});
      end 
        
   end % methods
    
   methods(Static)
      function valid = mf_valid(mfl)
         % merit function must have only equality linear constraints
         n = mfl.n;
         m = mfl.m;
         
         lin = mfl.linear == true(m,1);
         
         % Check for lower bounds only with zero bounds
         bnds = (mfl.jUpp == false(n,1)) ...
              & (mfl.jLow == false(n,1));

         valid = (sum(bnds) == n) && (sum(lin) == m);
     end
   end
    
end % classdef