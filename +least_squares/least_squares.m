% LEAST_SQUARES  Abstract factorization class for least squares.
classdef least_squares < matlab.mixin.Copyable
% This interface is responsible for solving systems of the form
%
%    [I       Q*A ][r] = [b]
%    [Q*A' -d^2*I ][x]   [c]
%
% which is equivalent to the problem
%
%    min_{x} 0.5*|A*x - b|^2_Q^2 + 0.5*d^2|x|^2 + c'*x, r = b - Q*A*x
%
% where A is an n-by-m matrix (n >= m), Q is a diagonal matrix with non-negative
% entries, and b, c are n- and m-vectors respectively.

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
   
   % The following properties are for iterative methods specifically
   properties (Access=public, Constant)
      TOL_RESIDUAL = 1; % Terminate based on residual norm
      TOL_ERROR = 2; % Terminate based on error norm
   end
    
   properties       
      time = 0  % time spent in factorization
      delta = 0 % current regularization parameter
      delta_max % maximum allowed regularization parameter
      increased_delta = false; % has delta been increased
      cond_max  % maximum allowed condition number of QR factor
      regularized % Boolean flag for regularization
      cond_num  % Condition number of matrix to be factorized
      
      n % Number of rows of A
      m % Number of columns of A
   end
   
   methods

      function self = preprocess(self, A, Aprod, Atprod, Q, solve_opts)
         %PREPROCESS  Prepare necessary quantities for future solves
         %
         % Prepares the underlying solver for solving the augmented systems.
         % For example, if a direct solver is used, factorization would occur here;
         % if iterative methods are used, the linear operators are be constructed.
         % 
         % Inputs:
         %    A            n-by-m matrix (n >= m)
         %    Aprod        function handle to compute products A*v
         %    Atprod       function handle to compute products A'*v
         %    Q            sparse diagonal matrix
         %    solve_opts   structure for optional arguments. Solver dependent.
         if nargin < 3 || isempty(Q)
             Q = [];
         end
         t = tic;
         self = self.preprocess_local(A, Aprod, Atprod, Q, solve_opts);
         self.time = self.time + toc(t);
      end

      function [x,r] = lsq3(self, b, c)
         %LSQ3  Solve the augmented system
         %
         % Solves the above augmented system for given right-hand-side.
         % Preprocess must be called beforehand.
      
         % LSQ2: b = 0.
         if isempty(b) || (isscalar(b) && b == 0)
            b = zeros(self.n,1);
         end
         
         % LSQ1: c = 0.
         if nargin < 3 || isempty(c) || (isscalar(c) && c == 0)
            c = zeros(self.m,1);
         end
          
         [x,r] = self.lsq_local(b, c);
      end
   end % methods
      
end % classdef