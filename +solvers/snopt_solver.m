% SNOPT_SOLVER Minimizes the penalty function with SNOPT
classdef snopt_solver < matlab.mixin.Copyable
% Minimizes the penalty function using SNOPT as subsolver.
%
% Not particularly advised for minimizing the penalty function due to the use
% of linesearches, and lack of second-derivative support because the Hessian
% is not explicit.
%
% SNOPT can be obtained from here:
%      http://www.sbsi-sol-optimize.com/asp/sol_product_snopt.htm
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
      fid             % File for log output
      snfid           % File snopt output

      new_sigma       % Updated value for sigma
      new_delta       % Updated value for delta

      options = struct();

      prev_x;         % previous iterate

      logH = ''
      logB = ''
      logT = {}; 
   end
    
   methods
        
      function self = snopt_solver(fletcher_solver, varargin)
         % SNOPT_SOLVER Constructor
         %
         % Required Inputs:
         %    fletcher_solver   pointer to supersolver
         %
         % Optional Inputs:
         %    x0                initial point
         %    fid               filename for log
         %    snfid             filename for SNOPT log
         p = inputParser;
         p.KeepUnmatched = false;
         p.addParameter('x0', fletcher_solver.mf.x0);
         p.addParameter('fid', '');
         p.addParameter('snfid', '');
         p.parse(varargin{:});
         
         self.x0 = p.Results.x0;
         self.fid = p.Results.fid;
         self.snfid = p.Results.snfid;
         
         self.options.name = fletcher_solver.nlp.name;
         self.options.stop = @(itn, nMajor, nMinor, condZHZ, ...
                   obj, merit, step,  primalInf, dualInf, maxViol, ...
                   maxViolRel, x, xlow, xupp, xmul, xstate, ...
                   F, Flow, Fupp, Fmul, Fstate) ...
              self.post_iteration(itn, nMajor, nMinor, condZHZ, ...
                   obj, merit, step,  primalInf, dualInf, maxViol, ...
                   maxViolRel, x, xlow, xupp, xmul, xstate, ...
                   F, Flow, Fupp, Fmul, Fstate);
         
         self.fletcher_solver = fletcher_solver;
      end
        
      function [fletcher, info] = solve(self)
         
         x = self.x0;
         self.prev_x = x;
         
         mfl = self.fletcher_solver.mf;
         n = mfl.n;
         m = mfl.m;
         linear = mfl.linear;
         xlow = mfl.bL;
         xupp = mfl.bU;
         xmul = zeros(n,1);
         xstate = 2*ones(n,1);
         
         % F is structured as
         % F = [objective            ]
         %     [linear constraints   ]
         %     [nonlinear constraints]
         % Generally, there should only be linear constraints, but we'll
         % accept the general form for now
         
         % Need to undo equality slacks since SNOPT expects linear
         % constraints in the form B*x = d
         d = -mfl.fcon(zeros(n,1));
         d = d(linear);
         
         Flow = [-Inf;...
                 mfl.cL( linear) + d;...
                 mfl.cL(~linear)];
         Fupp = [ Inf;...
                 mfl.cU( linear) + d;...
                 mfl.cU(~linear)];
         Fmul = zeros(m+1,1);
         Fstate = zeros(m+1,1);
         ObjAdd = 0;
         ObjRow = 1;
         userfun = @(x) self.objgrad(x);
         [iAfun,jAvar,~] = find(mfl.Jpattern(linear,:));
         iAfun = iAfun + 1; % Since objective is in first row
         J = mfl.gcon(self.x0);
         [~,~,A] = find(J(linear,:));
         [iGfun,jGvar,~] = find(mfl.Jpattern(~linear,:));
         iGfun = iGfun + sum(linear) + 1; % Nonlinear constraints go last
         iGfun = [ones(n,1);iGfun];
         jGvar = [(1:n)';jGvar];
         
         % Set SNOPT Options
         snseti ('Superbasic limit', 2*n);
         snsetr ('Major feasibility tolerance', 1e-16);
         snsetr ('Major optimality  tolerance', 1e-16);
         snsetr ('Minor feasibility tolerance', 1e-8);
         
         ctr = 1;
         while true
             
            if ~strcmp(self.snfid, '')
                snprint([self.snfid '_' int2str(ctr) '.out']);
            end
             
            [x,~,~,xmul,Fmul,xstate,Fstate,info] = ...
                snopt( self.x0, xlow, xupp, xmul, xstate,  ...
                       Flow, Fupp, Fmul, Fstate,     ...
                       userfun, ObjAdd, ObjRow, ...
                       A, iAfun, jAvar, iGfun, jGvar, ...
                       self.options );
            
            snset ('Sticky parameters Yes');
                   
            snprint off;
                   
            info_code = 10*floor(info.info/10);
            
            if info_code == 70 && self.fletcher_solver.exit == 0 % User requested stop (i.e. supersolver)
                % need to update values
                self.fletcher_solver.mf = set_sigma(self.fletcher_solver.mf, self.new_sigma);
                self.fletcher_solver.mf.delta = self.new_delta;
            elseif self.fletcher_solver.exit == 0
                % snopt quit for some reason
                info.exit = self.fletcher_solver.EXIT_SUBSOLVER;
                info.exit_msg = self.getMessage(info.info);
                break;
            else
                info.exit = self.fletcher_solver.exit;
                info.exit_msg = self.fletcher_solver.EXIT_MSG{self.fletcher_solver.exit};
                break;
            end
            
            ctr = ctr+1;
         end
         
         snend;
         
         info.sol.x = x;
         self.fletcher_solver.mf.fobj(info.sol.x);
         info.sol.y = self.fletcher_solver.mf.val.y;
         info.sol.f = self.fletcher_solver.mf.val.f;
         
         fletcher = self.fletcher_solver;
      end
        
      function iAbort = post_iteration(self, itn, nMajor, nMinor, condZHZ, ...
                           obj, merit, step,  primalInf, dualInf, maxViol, ...
                           maxViolRel, x, xlow, xupp, xmul, xstate, ...
                           F, Flow, Fupp, Fmul, Fstate)
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
             iAbort = 0;
             return;
         end

         if info.flag == self.fletcher_solver.SUBSOLVER_USER
             % Supersolver wants to quit for some reason
             iAbort = -1;
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
             iAbort = -1;
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

      function [F,G] = objgrad(self,x)
         linear = self.fletcher_solver.mf.linear;
         
         f = self.fletcher_solver.mf.fobj(x);
         c = self.fletcher_solver.mf.fcon(x);
         F = [f;zeros(size(c(linear)));c(~linear)];
         
         g = self.fletcher_solver.mf.gobj(x);
         A = self.fletcher_solver.mf.gcon(x);

         % TODO: Check if this is most efficient way to get nonzero
         % values
         [~,~,s] = find(A(~linear,:));
         G = [g;s];
      end
        
      function msg = getMessage(self, flag)
         switch flag
            % Finished successfully
            case 1,    msg = 'optimality conditions satisfied';
            case 2,    msg = 'feasible point found';
            case 3,    msg = 'requested accuracy could not be achieved';
            case 4,    msg = 'elastic objective minimized';
            case 5,    msg = 'elastic infeasibilities minimized';
            % The problem appears to be infeasible
            case 11,   msg = 'infeasible linear constraints';
            case 12,   msg = 'infeasible linear equality constraints';
            case 13,   msg = 'nonlinear infeasibilties minimized';
            case 14,   msg = 'linear infeasibilities minimized';
            case 15,   msg = 'infeasible linear constraints in QP subproblem';
            case 16,   msg = 'infeasible nonelastic constraints';
            % The problem apppears to be unbounded
            case 21,   msg = 'unbounded objective';
            case 22,   msg = 'constraint violation limit reached';
            % Resource limit error
            case 31,   msg = 'iteration limit reached';
            case 32,   msg = 'major iteration limit reached';
            case 33,   msg = 'the superbasics limit is too small';
            case 34,   msg = 'time limit reached';
            % Terminated after numerical difficulties
            case 41,   msg = 'current point cannot be improved';
            case 42,   msg = 'singular basis';
            case 43,   msg = 'cannot satisfy the general constraints';
            case 44,   msg = 'ill-conditioned null-space basis';
            case 45,   msg = 'unable to compute acceptable LU factors';
            % Error in the user-supplied function
            case 51,   msg = 'incorrect objective derivatives';
            case 52,   msg = 'incorrect constraint derivatives';
            case 56,   msg = 'irregular or badly scaled problem functions';
            % Undefined user-supplied functions
            case 61,   msg = 'undefined function at the first feasible point';
            case 62,   msg = 'undefined function at initial point';
            case 63,   msg = 'unable to proceed into undefined region';
            % User requested termination
            case 71,   msg = 'terminated during function evaluation';
            case 74,   msg = 'terminated from monitor routine';
            % Insufficient storage allocated
            case 81,   msg = 'work arrays must have at least 500 elements';
            case 82,   msg = 'not enough character storage';
            case 83,   msg = 'not enough integer storage';
            case 84,   msg = 'not enough real storage';
            % Input arguments out of range
            case 91,   msg = 'invalid input argument';
            case 92,   msg = 'basis file dimension do not match the problem';
            % System error
            case 141,  msg = 'wrong number of basic variables';
            case 142,  msg = 'error in basis package';
            otherwise, msg = 'unknown error';
         end
      end
        
    end % methods

end % classdef