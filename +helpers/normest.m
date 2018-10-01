function [e,cnt] = normest(S, St, m, tol)
%NORMEST Estimate the matrix 2-norm.
%
%   normest(S, St, m, tol) is an estimate of the 2-norm of the matrix S with m rows.
%
%   [nrm,cnt] = normest(..) also gives the number of iterations used.
%
%   This function is a minor adaptation of Matlab's built-in NORMEST.
%
%   See also NORM, COND, RCOND, CONDEST.

if nargin < 4, tol = 1.0e-6; end
maxiter = 100;
cnt = 0;

% Compute an "estimate" of the ab-val column sums.
v = ones(m,1);
v(randn(m,1) < 0) = -1;
x = abs( St(v) );

e = norm(x);
if e == 0, return, end
x = x/e;
e0 = 0;
while abs(e-e0) > tol*e
   e0 = e;
   Sx = S(x);
   if nnz(Sx) == 0
      Sx = rand(size(Sx));
   end
   x = St(Sx);
   normx = norm(x);
   e = normx/norm(Sx);
   x = x/normx;
   cnt = cnt+1;
   if cnt > maxiter
      warning('normest:notconverge', '%s%d%s%g', ...
              'NORMEST did not converge for ', maxiter, ' iterations with tolerance ', tol);
      break;
   end
end
