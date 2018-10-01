% LSLNLQ  LNLQ iterative solver for least squares.
classdef lslnlq < least_squares.least_squares
% Solve systems of the form
%
%    [I      Q*A ][r] = [b]
%    [Q*A' -d^2*I][x]   [c]
%
% via the iterative method LNLQ.
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
      Mprod     % Handle to solves with preconditioner
      Aprod     % Handle to product A*v
      Atprod    % Handle to product A'*v
                
      sigma     % underestimate to sigma_{min}(P\A)
      tol_cond  % termination condition (residual or error)
      
      tol       % Termination tolerance
      maxit     % Maximum number of iterations
   end
   
   methods

      function self = lslnlq(m, n, varargin)
         %LSLNLQ  Constructor
         %
         % Sets internal variables and constructs operator for linear system solve.
         %
         % Required Inputs:
         %    m            number of rows in A
         %    n            number of cols in A
         %
         % Optional Inputs:
         %    regularized  boolean to indicate if d>0
         %    delta_max    maximum regularization value
         %    termination_condition     least_squares.TOL_RESIDUAL or 
         %                              least_squares.TOL_ERROR
         %    tol          tolerance for linear solve

         % Construct superclass
         self = self@least_squares.least_squares();
          
         % Parse input parameters and initialize object properties.
         ip = inputParser;
         ip.KeepUnmatched = true;
         ip.addParameter('regularized',true);
         ip.addParameter('delta_max',eps^(1/4));
         ip.addParameter('termination_condition', self.TOL_RESIDUAL);
         ip.addParameter('tol', 1e-10);
         parse(ip, varargin{:});
         
         self.tol_cond = ip.Results.termination_condition;
         self.tol = ip.Results.tol;
         
         self.delta_max = ip.Results.delta_max;
         self.regularized = ip.Results.regularized;
         self.m = m;
         self.n = n;
         self.maxit = 20*(n+m);
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function self = preprocess_local(self, ~, Aprod, Atprod, Q, solve_opts)    
         %PREPROCESS_LOCAL  Called by preprocess in least_squares
         %
         % Constructs linear operator for iterative solve.
         
         ip = inputParser;
         ip.KeepUnmatched = true;
         ip.addParameter('preconditioner', @() (@(x) x));
         ip.addParameter('min_singular_val', 0);
         parse(ip,solve_opts);
         
         if isempty(Q) 
             Q = speye(self.n);
         end
         
         if ~self.regularized
             self.Aprod = @(v) Q*Aprod(v);
             self.Atprod = @(v) Atprod(Q*v);
             d = 0;
         else
             % If regularized, solve problem
             % [I             Q*A ][r] = [Q*b]
             % [        I    del I][z] = [ 0 ] 
             % [A'*Q  del*I       ][x]   [ c ]
             d = self.delta;
             self.Aprod = @(v) [Q*Aprod(v); d*v];
             self.Atprod = @(v) Atprod(v(1:self.n)) + d*v(self.n+1:end);
         end
         
         self.Mprod = ip.Results.preconditioner();
         self.sigma = ip.Results.min_singular_val;
         self.delta = d;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function [x, r] = lsq_local(self, b, c)
         %LSQ_LOCAL  Called by lsq3 in least_squares.
         
         if ~self.regularized
             bc = c - self.Atprod(b);
         else
             bc = c - self.Atprod([b;zeros(self.m,1)]);
         end
         
         atol = 0;
         if self.tol_cond == self.TOL_RESIDUAL || self.sigma == 0
             btol = self.tol;
             etol = -Inf;
         else
             btol = -Inf;
             etol = self.tol;
         end
         
         [r,x,flag,iter,normr,resvec,errvec] = lnlq(...
             @(x,t) self.Aop(x,t), bc, btol, atol, etol, [], ...
             self.maxit, self.Mprod, self.delta, self.sigma);
         r = r(1:self.n) + b; % If was regularized, remove z part
         
         if( flag ~= 0 && flag ~= 4)
             warning('LNLQ failed to converge to tolerance %1.1e\n', self.tol);
             warning('       achieved tolerance %1.1e', normr); 
         end
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function u = Aop(self, v, t)
          if t==1
              u = self.Atprod(v);
          else
              u = self.Aprod(v);
          end
      end
      
   end % methods
      
end % classdef