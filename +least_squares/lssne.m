% LSSNE  Seminormal equation solver for least squares.
classdef lssne < least_squares.least_squares
% Solve systems of the form
%
%    [I      Q*A ][r] = [b]
%    [Q*A' -d^2*I][x]   [c]
%
% using the seminormal equations using the Q-less QR factorization of Q*A.
% We obtain the R factor from (Q*A), and solve
%
%     (R'*R + d^2*I) x = A'*b - c
%
% with one step of iterative refinement.
%
% This implementation uses SuiteSparse's SPQR:
%     http://faculty.cse.tamu.edu/davis/suitesparse.html
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
      % These properties are post-fixed by "l" to distinguish them from the
      % local copies used in the methods.
      pl        % fill-reducing permutation
      Al        % copy of matrix that is being factorized (Q*A)
      Rl        % triangular factor
   end
   
   methods

      function self = lssne(A, varargin)
         %LSSNE  Constructor
         %
         % Sets internal variables and prepares for QR factorization.
         %
         % Required Inputs:
         %    A            n-by-m matrix
         %
         % Optional Inputs:
         %    regularized  boolean to indicate if d>0
         %    delta_max    maximum regularization value
         %    cond_max     maximum allowable condition number

         % Construct superclass
         self = self@least_squares.least_squares();
          
         % Parse input parameters and initialize object properties.
         ip = inputParser;
         ip.KeepUnmatched = true;
         ip.addParameter('regularized',true);
         ip.addParameter('delta_max',eps^(1/4));
         ip.addParameter('cond_max',(1/eps)^(1/2));
         parse(ip, varargin{:});

         self.delta_max = ip.Results.delta_max;
         self.cond_max = ip.Results.cond_max;
         self.regularized = ip.Results.regularized;
         
         % Fill-reducing column ordering.
         n = size(A,2);
         if self.regularized
            p = colamd([A; speye(n)]);
         else
            p = colamd(A);
         end
         self.pl = p;
         
         [n,m] = size(A);
         self.m = m;
         self.n = n;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function self = preprocess_local(self, A, ~, ~, Q, ~)
         %PREPROCESS_LOCAL  Called by preprocess in least_squares
         %
         % Computes R factor of QR factorization in preparation for
         % upcoming linear solves.

         n = size(A, 2);
         p = self.pl;
         if ~self.regularized; d=0; else d=self.delta; end

         if ~isempty(Q)
             A = Q*A;
         end
         
         while true
            if d == 0
               AA = A;
            else
               AA = [ A; d*speye(n) ];
            end
            
            spqr_opts = struct('Q','discard',...
                               'econ', n,...
                               'ordering','fixed');

            [~, R] = spqr(AA(:,p), spqr_opts);
            
            % Estimate the condition number of R.
            r = abs(spdiags(R,0));
            condR = max(r) / min(r);
            
            % If too ill conditioned throw exception
            % Assuming that this is a trial point being checked, and if it
            % is too singular we move onto the next point instead. A more
            % TODO: More sophisticated approach to allow for increasing the
            % regularization parameter.
            if n > 0 && self.regularized && condR > self.cond_max
%                d = max(r) / self.cond_max;
% %              d = max([ 2*d, 1e-8, 0.1*min(r) ]);
%                if d > self.delta_max
%                   ME = MException('FLETCHER_MERIT:regularization',...
%                      'regularization parameter exceeded maximum value');
%                   throw(ME);
%                else
%                   self.delta = d;
%                   warning('increased delta to %10.1e', d);
%                end
               self.increased_delta = true;
               ME = MException('FLETCHER_MERIT:regularization',...
                  'regularization parameters needs to be increased');
               throw(ME);
            else
               self.increased_delta = false;
               break
            end
         end
         
         self.delta = d;
         self.Al = A;
         self.Rl = R;
         
         self.cond_num = condR;
         
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function [x, r] = lsq_local(self, b, c)
         %LSQ_LOCAL  Called by lsq3 in least_squares.
         %
         % Apply the corrected semi-normal equations (SNE):
         %
         %     (A'*Q^2*A + del^2) x = A'Qb - c   (SNE)
         % ie, (R'R             ) x = A'Qb - c,
         %
         % plus a step of iterative refinement. The iterative refinement
         % step is easier to derive using the augmented system formulation:
         %
         % 1. Solve for initial solution:
         %
         % [I        Q A][r] = [b]  <==>  SNE above
         % [A'Q -del^2 I][x]   [c]
         %
         % 2. Compute residuals:
         %
         % [I        Q A][dr] = [b] - [I        Q A][r] = b - r - Ax
         % [A'Q -del^2 I][dx]   [c]   [A'Q -del^2 I][x] = c - A'Qr + del^2 x
         %
         % 3. Correct initial solutions:
         %
         % r = r + dr
         % x = x + dx
         
         R = self.Rl;
         A = self.Al;
         p = self.pl;
         if ~self.regularized; d=0; else d=self.delta; end
         [m, n] = size(A);
         
         % LSQ2: b = 0.
         if isempty(b) || (isscalar(b) && b == 0)
            b = zeros(m,1);
         end
         
         % LSQ1: c = 0.
         if nargin < 3 || isempty(c) || (isscalar(c) && c == 0)
            c = zeros(n,1);
         end
         
         % Step 1. Solve for initial solution
         s = A'*b - c;
         z = cs_utsolve(R, s(p)); % = R' \ s(p);
         x = cs_usolve(R, z); % = R \ z;
         x(p) = x;
         Ax = A*x;
         r = b - Ax;
         
         % Step 2. Compute residuals
         w = b - r - Ax;
         s = c - A'*r + d^2*x;
         
         % Step 3. Solve for corrections
         s = A'*w - s;
         z = cs_utsolve(R, s(p)); % = R' \ s(p);
         dx = cs_usolve(R, z); % = R \ z;
         dx(p) = dx;
         dr = w - A*dx;
         
         % Step 4. Update solution
         x = x + dx;
         r = r + dr;
         
      end

   end % methods
      
end % classdef