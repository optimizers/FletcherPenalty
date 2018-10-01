% BCFLASH_SOLVER Interface for calling BCFLASH to optimize Fletcher penalty
classdef bcflash_solver < bcflash
% Minimizes the penalty function using BCFLASH as subsolver, and must therefore
% be on the Matlab path.
% BCFLASH does not accept constraints on the penalty except for bounds.
%
% BCFLASH can be downloaded from:
%      https://github.com/restrin/bcflash
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
        
      logH = '%5s  %4s'
      logB = '%5i  %4s'
      logT = {'CGits', 'stat'};
   end

   methods

      function self = bcflash_solver(fletcher_solver, varargin)
         % BCFLASH_SOLVER Constructor
         %
         % Required Inputs:
         %    fletcher_solver   pointer to supersolver
         %
         % Optional Inputs:
         %    subsolver_options see BCFLASH documentation for options
         %    x0                initial point

         % Check that merit function is a valid formulation
         if ~solvers.bcflash_solver.mf_valid(fletcher_solver.mf)
            ME = MException('FLETCHER_SOLVER:BCFLASH_SOLVER',...
                 'Fletcher penalty must only have bound constraints.');
            throw(ME);
         end 
         p = inputParser;
         p.KeepUnmatched = false;
         p.addParameter('subsolver_options', struct(), @isstruct);
         p.addParameter('x0', fletcher_solver.mf.x0);
         p.parse(varargin{:}); 
         subsolver_options = p.Results.subsolver_options;
         self = self@bcflash(fletcher_solver.mf, subsolver_options);
         % Hack for setting callback hook
         % Makes post_iteration a class method instead of static
         self.callback = @post_iteration;
         
         self.x0 = p.Results.x0;
         self.fletcher_solver = fletcher_solver;
      end
    
      function [fletcher, info] = solve(self)
         % Call bcflash
         [~,~,self] = solve@bcflash(self, self.x0);

         % Retrieve an exit code. The super-solver's exit code takes
         % precedence. If it's zero, than it must be the subsolver that
         % requested exit, and then we report that exit code.
         if self.fletcher_solver.exit == 0
            % Subsolver requested the exit.
            info.exit = self.fletcher_solver.EXIT_SUBSOLVER;
            info.exit_msg = sprintf('%s (%s)', ...
               self.fletcher_solver.EXIT_MSG{info.exit}, self.exit_msg);
         else
            info.exit = self.fletcher_solver.exit;
            info.exit_msg = ...
               self.fletcher_solver.EXIT_MSG{self.fletcher_solver.exit};
         end
                     
         % Set values in return struct
         info.sol.x = self.fletcher_solver.last_x;
         info.sol.y = self.fletcher_solver.mf.val.y;
         info.sol.f = self.fletcher_solver.mf.val.f;
         
         fletcher = self.fletcher_solver;
      end
        
      function [self, exit_flag] = post_iteration(self, x, cgits, successful)
      % POST_ITERATION Callback function called at end of each iteration.
      
         exit_flag = self.getprop('EXIT_NONE');
         
         % Prepare additional logging
         if successful
            succstr = '';
         else
            succstr = 'rej';
         end
         
         logstr = sprintf(self.logB, cgits, succstr);
         
         % Call the generic post_iteration hook
         [self.fletcher_solver, info] = ...
            self.fletcher_solver.post_iteration(x, successful, ...
               self.logHeader(), logstr);
         
         % Adjust problem according to info
         if info.flag == self.fletcher_solver.SUBSOLVER_NONE
            % Do nothing, just keep going
            return;
         end
         
         if info.flag == self.fletcher_solver.SUBSOLVER_USER
            % Supersolver wants to quit for some reason
            exit_flag = self.getprop('EXIT_USER_REQUEST');
            return;
         end
         
         if info.flag == self.fletcher_solver.SUBSOLVER_MODIFY
            % Penalty function is being modified
            if info.sigma >= 0
               % Need to update sigma
               self.fletcher_solver.mf = ...
                  set_sigma(self.fletcher_solver.mf, info.sigma);
               exit_flag = self.getprop('EXIT_RESTART');
            end
            if info.delta >= 0
               % Need to update delta
               self.fletcher_solver.mf.delta = info.delta;
               exit_flag = self.getprop('EXIT_RESTART');
            end
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
         % MF_VALID Ensure penalty function is valid
         %
         % BCFLASH does not accept explicit non-bound constraints.
         valid = (mfl.m == 0);
      end
   end
    
end % classdef