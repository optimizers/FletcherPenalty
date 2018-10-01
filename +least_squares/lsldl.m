% LSLDL  LDL Factorization class for least squares.
classdef lsldl < least_squares.least_squares
% Solve systems of the form
%
%    [a*I       Q*A  ][ r ] = [a*b]
%    [Q*A' -(d/a)^2*I][a*x]   [ c ]
%
% via an LDL factorization. The scalar a is used to improve
% the conditioning of the linear system.
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
      static_p % Use static permutation?

      % These properties are post-fixed by "l" to distinguish them from the
      % local copies used in the methods.
      Al        % copy of matrix that is being factorized
      Ll        % lower triangular factor
      Dl        % diagonal factor
      al = 1e-4 % KKT scaling factor
      pl        % permutation vector (if static permutation used)
      Kl        % Full saddle-point matrix
   end
   
   methods

      function self = lsldl(A, varargin)
         %LSLDL  Constructor
         % 
         % Sets internal variables and if static_p=true, then computes
         % matrix reordering for LDL factorization. Factorizations are
         % faster if static_p=true, but may be more unstable if d=0.
         % 
         % Required Inputs:
         %    A            n-by-m matrix (n>=m)
         %
         % Optional Inputs:
         %    regularized  boolean to indicate if d>0
         %    delta_max    maximum regularization value
         %    cond_max     maximum allowable condition number
         %    static_p     use static re-ordering?

         % Construct superclass
         self = self@least_squares.least_squares();
          
         % Parse input parameters and initialize object properties.
         ip = inputParser;
         ip.KeepUnmatched = true;
         ip.addParameter('regularized',true);
         ip.addParameter('delta_max',eps^(1/4));
         ip.addParameter('cond_max',(1/eps)^(1/2));
         ip.addParameter('static_p', false);
         parse(ip, varargin{:});

         self.delta_max = ip.Results.delta_max;
         self.cond_max = ip.Results.cond_max;
         self.regularized = ip.Results.regularized;
         self.static_p = ip.Results.static_p;

         if self.static_p
            % If we permute based on sparsity alone (rather than based on 
            % stability) then we should do so now
            
            % Fill-reducing column ordering.
            [n, m] = size(A);
            if self.regularized
               p = colamd([speye(n) A; A' speye(m)]);
            else
               p = colamd([speye(n) A; A' sparse(m,m)]);
            end
            self.pl = p;
         end

         [n,m] = size(A);
         self.m = m;
         self.n = n;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function self = preprocess_local(self, A, ~, ~, Q, ~)    
         %PREPROCESS_LOCAL  Called by preprocess in least_squares
         %
         % Factorizes saddle-point matrix with LDL factorization and
         % stores it for later use.
         
         [n, m] = size(A);
         if ~self.regularized; d=0; else d=self.delta; end
         
         if isempty(Q)
             Q = speye(n);
         end
         
         if self.static_p
            a = 1;
            thresh = 0;
         else
            % TODO: How to choose?
            a = self.al;
            thresh = 0.01;
         end
         
         while true
            % Q is diagonal, so inv should be fine
            K = [a*eye(n) Q*A; A'*Q -d^2*speye(m)/a];

            % TODO: Use scaling matrix?
            [L,D,p] = ldl(K, thresh, 'lower', 'vector');
            break;

            % TODO: Estimate condition number
%             cond_A = condest(L)*sqrt(condest(D));
%              
%             % If too ill conditioned, increase regularization parameter.
%             if m > 0 && self.regularized && cond_A > self.cond_max
%                % TODO: Come up with better way to get regularization
%                % parameter
%                d = max([1e-3, 2*d]);
%                if d > self.delta_max
%                   ME = MException('FLETCHER_MERIT:regularization',...
%                      'regularization parameter exceeded maximum value');
%                   throw(ME);
%                else
%                   self.delta = d;
%                   warning('increased delta to %10.1e', d);
%                end
%                self.increased_delta = true;
%             else
%                self.increased_delta = false;
%                break
%             end
         end
         
         self.Al = A;
         self.pl = p;
         self.Ll = L;
         self.Dl = D;
         self.delta = d;
         self.Kl = K;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function [x, r] = lsq_local(self, b, c)
         %LSQ_LOCAL  Called by lsq3 in least_squares.
         
         [n, m] = size(self.Al);
         
         K = self.Kl;
         a = self.al;
         bc = [a*b;c];
         
         x(self.pl,:) = self.Ll'\(self.Dl\(self.Ll\bc(self.pl)));
         
         resid = bc - K*x;
         dx(self.pl,:) = self.Ll'\(self.Dl\(self.Ll\resid(self.pl)));
         
         x = x + dx;
         
         r = x(1:n);
         x = x(n+1:end)/a;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
   end % methods
      
end % classdef