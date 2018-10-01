% LSMINRES  Minres solver for least squares.
classdef lsminres < least_squares.least_squares
% Solve systems of the form
%
%    [I      Q*A ][r] = [b]
%    [Q*A' -d^2*I][x]   [c]
%
% via the iterative method minres.
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
      Kprod     % handle to product with augmented system
      Mprod     % handle to solves with preconditioner
      Aprod     % handle to product A*v
      Atprod    % handle to product A'*v
                
      tol       % tolerance for iterative method
      maxit     % maximum number of iterations
   end
   
   methods

      function self = lsminres(m, n, varargin)
         %LSMINRES  Constructor
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
         %    tol          tolerance for linear solve

         % Construct superclass
         self = self@least_squares.least_squares();
          
         % Parse input parameters and initialize object properties.
         ip = inputParser;
         ip.KeepUnmatched = true;
         ip.addParameter('regularized',true);
         ip.addParameter('delta_max',eps^(1/4));
         ip.addParameter('tol', 1e-8);
         parse(ip, varargin{:});
         
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
         parse(ip,solve_opts);
         
         preconditioner = ip.Results.preconditioner();
         
         if isempty(Q)
             Q = speye(self.n);
         end
         
         if ~self.regularized; d=0; else d=self.delta; end
         
         self.Kprod = @(v) [v(1:self.n) + Q*Aprod(v(self.n+1:end)); ...
                            Atprod(Q*v(1:self.n)) - d^2*v(self.n+1:end)];
         self.Mprod = @(v) [v(1:self.n); preconditioner(v(self.n+1:end))];

         self.Aprod = Aprod;
         self.Atprod = Atprod;
         self.delta = d;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function [x, r] = lsq_local(self, b, c)
         %LSQ_LOCAL  Called by lsq3 in least_squares.
         
         [x, flag, relres, iter] = minres(self.Kprod, [b;c], self.tol, self.maxit, self.Mprod);
         r = x(1:self.n);
         x = x(self.n+1:end);
         
         if( flag ~= 0 )
             warning('MINRES failed to converge to tolerance %1.1e\n', self.tol);
             warning('       achieved tolerance %1.1e', relres); 
         end
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
   end % methods
      
end % classdef