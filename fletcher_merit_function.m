% FLETCHER_MERIT_FUNCTION Penalty function implementation
classdef fletcher_merit_function < model.nlpmodel & handle
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

   properties (SetAccess = private, Hidden=true)
      new_sigma = true   % flag if sigma has been updated before an obj eval
      new_delta = true   % flag if delta has been updated before an obj eval
      last_fx_hash = ''  % MD5 of iterate x at which obj was last evaluated
      last_gx_hash = ''  % MD5 of iterate x at which obj was last evaluated
      last_cx_hash = ''  % MD5 of iterate x at which obj was last evaluated
      last_Jx_hash = ''  % MD5 of iterate x at which obj was last evaluated
   end

   properties (Dependent)
      delta        % dual regularization parameter
      delta_max    % ...and its max value
   end
   
   properties
      nlnFlag      % linear constraints explicit?
   end

   properties (SetAccess = private, Hidden=false)
      nlp          % Original NLP object
      val          % structure to hold last iteration info
      sigma        % penalty parameter
      rho          % quadratic penalty parameter
      lsq          % least-squares solver object
      hprod        % Hessian-vector product option (1, 2, or 3)
      iExp         % Index set of explicit constraints
      cwidth       % width of complementarity approximation function
   end
   
   methods
      
      function self = fletcher_merit_function(nlp, lsq, varargin)

         % Parse input parameters and initialize object properties.
         p = inputParser;
         addParameter(p,'delta', 0);
         addParameter(p,'sigma', 0);
         addParameter(p,'rho', 0);
         addParameter(p,'x0', nlp.x0);
         addParameter(p,'hprod', 1, ...
            @(x)validateattributes(x,{'numeric'},{'scalar','>=',1,'<=',4}));
         addParameter(p,'lin_explicit', false);
         parse(p, varargin{:});

         name = sprintf('Fletcher merit (%s)',nlp.name);
         lin_explicit = p.Results.lin_explicit;

         if ~fletcher_merit_function.formulated_correctly(nlp)
             ME = MException('FLETCHER_SOLVER:MERIT_FCN',...
                 'Problem must be in slack formulation.');
             throw(ME);
         end
         
         if lin_explicit
            % Note: we're assuming nlp has been given in
            % slack formulation already based on above check
            cL = nlp.cL(nlp.linear);
            cU = nlp.cU(nlp.linear);
         else
            cL = [];
            cU = [];
         end
         
         self = self@model.nlpmodel(name, nlp.x0, cL, cU, nlp.bL, nlp.bU);
         
         if lin_explicit
             % Classify constraints as linear
             self.linear = true(size(cL));
             self.iExp = nlp.linear;
             self.Jpattern = nlp.Jpattern(self.iExp, :);
         else
             self.iExp = false(nlp.m,1);
         end
         
         self.nlp     = nlp;
         self.lsq     = lsq;
         self.delta   = p.Results.delta;
         self.sigma   = p.Results.sigma;
         self.rho     = p.Results.rho;
         self.hprod   = p.Results.hprod;
         self.cwidth  = self.jTwo .* min(1,(self.bU - self.bL)/2);
         x0           = p.Results.x0;
         self.nlnFlag = lin_explicit;

         % Evaluate the objective function and initialize multipliers.
         [~] = self.fobj(x0);
         [~] = self.gobj(x0);

      end % fletcher_merit_function

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function fPhi = fobj_local(self, x)
         %FOBJ_LOCAL  Objective of the merit function.

         % Check if this x was the same as the last. Quick exit if so.
         hash = utils.DataHash(x);
         if strcmp(hash, self.last_fx_hash) ...
            && ~self.new_sigma && ~self.new_delta

            fPhi = self.val.fPhi;
            return
         end
         self.last_fx_hash = hash;

         % NLP objective and gradient.
         f = self.nlp.fobj(x);
         g = self.nlp.gobj(x);
         
         c = self.nlp.fcon(x);
         self.val.c_nlp = c; % unshifted constraints

         % NLP constraint and Jacobian. Transform into c(x) = 0 by
         % substracting off any RHS.
         J = self.nlp.gcon(x);
         c = c - self.nlp.cL;
         [Jprod, Jtprod] = self.nlp.gconprod(x);
         
         % Scaling matrix.
         jLow = self.nlp.jLow;
         nLow = sum(jLow);
         jUpp = self.nlp.jUpp;
         nUpp = sum(jUpp);
         jTwo = self.nlp.jTwo;
         nTwo = sum(jTwo);
         bL = self.bL;
         bU = self.bU;
         
         % Sometimes the bounds aren't exactly satisfied (they're violated
         % by a tiny bit), which needs to be handled.
         assert( sum(x(jLow) - bL(jLow) >= 0) == nLow ...
            & sum(bU(jUpp) - x(jUpp) >= 0) == nUpp );
      
         Q = self.Qop(x);

         self.val.x = x;
         self.val.f = f;
         self.val.g = g;
         self.val.c = c;
         self.val.J = J;
         self.val.Jprod = Jprod;
         self.val.Jtprod = Jtprod;
         self.val.Q = Q;

         % Compute quantities needed for LS multipliers.
         try
            % Do not build preconditioner unless required
            P = @() self.nlp.preconditioner(x);
            min_singular_val = self.nlp.gcon_min_singular_value(x);
            solve_opts = struct('preconditioner', P, ...
                'min_singular_val', min_singular_val);
            self.lsq = preprocess(self.lsq, J', Jtprod, Jprod, Q, solve_opts);
         catch ME
             % If the Jacobian is too ill-conditioned, we treat the point
             % as a singularity
             if strcmp(ME.identifier, 'FLETCHER_MERIT:regularization')
                 self.val.fPhi = Inf;
                 fPhi = Inf;
                 self.val.gL = Inf*ones(size(g));
                 return
             else
                 rethrow(ME)
             end
         end
         
         [y,d] = self.multiplier(g, self.sigma*c); % multiplier estimate
         self.val.y = y;                           % ... store it
         
         if nLow == 0 && nUpp == 0
            % This should be more accurate in the equality only case
            gL = d;                    % gradient of Lagrangian
         else
            gL = g - Jtprod(y);
         end
         self.val.gL = gL;             % ... store it
         self.val.Qd = d;              % Q^{1/2}*gL
         
         yExp = y;
         yExp(~self.iExp) = 0;         % Get explicit constraints
         y = y(~self.iExp);            % Get multipliers for
         c = c(~self.iExp);            % ... non-explicit constraints

         % Matlab seems to have a bug where if a scalar has its entry
         % removed then it becomes a 0x0 matrix instead of a 0x1 matrix.
         if size(y,2) == 0
             y = zeros(size(y,1),1);
             c = zeros(size(c,1),1);
         end
         
         % gradient of 'reduced' Lagrangian (ignore linear constraints)
         gRL = gL + Jtprod(yExp);      
         self.val.gRL = gRL;           % ... store it
         
         fPhi = f - c'*y + 0.5*self.rho*(c'*c);    % merit function value
         self.val.fPhi = fPhi;                     % ... store it

         % This tells hprod that it's safe to evaluate Hessian.
         self.new_sigma = false;
         self.new_delta = false;
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function gPhi = gobj_local(self, x)
         %GOBJ_LOCAL  Gradient of the merit function.

         % Make sure that the merit function has been evaluated.
         self.fobj(x);
 
         % Check if this x was the same as the last. Quick exit if so.
         hash = utils.DataHash(x);
         if strcmp(hash, self.last_gx_hash) ...
            && ~self.new_sigma && ~self.new_delta
            gPhi = self.val.gPhi;  %#ok<UNRCH>
            return
         end
         self.last_gx_hash = hash;
                  
         c  = self.val.c;
         y  = self.val.y;
         gRL = self.val.gRL;
         Jtprod  = self.val.Jtprod;
         
         c(self.iExp) = 0; % Ignore Jacobian of explicit constraints
         
         % Merit function gradient.
         Yc = self.Yprod(self.sigma, y, c);
         gPhi = gRL - Yc;
         if self.rho ~= 0
            gPhi = gPhi + self.rho*Jtprod(c);
         end
         self.val.gPhi = gPhi;      
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function Hx = hlagprod_local(self, ~, ~, v)
         %HLAGPROD_LOCAL  Products with approximate Hessian of phi
         %
         % This implementation assumes only linear constraints are explicit.

         switch self.hprod
            case 1
               Hx = self.hprod1(v);
            case 2
               Hx = self.hprod2(v);
            case 3
               Hx = self.hprod3(v);
            case 4
               Hx = self.hprod4(v);  
         end
         
         x = self.val.x;
         c = self.val.c;
         c(self.iExp) = 0; % Ignore explicit constraints

         % Add Hessian of quadratic penalty term
         if self.rho ~= 0
            Jtprod = self.val.Jtprod;
            Jprod = self.val.Jprod;
            Jv = Jprod(v);
            Hx = Hx + self.rho*self.nlp.hconprod(x,c,v) + self.rho*Jtprod(Jv);
         end
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function [y,r] = multiplier(self, g, w)
         %MULTIPLIER  Compute first-order multiplier estimates
         %
         % y = argmin_y 0.5*|Ay - g|^2_Q^2 + sigma*c'*w, r = Q*g - Q*A*y
         
         if nargin < 3 || isempty(w)
            w = zeros(self.nlp.m, 1);
         end
         Q = self.val.Q;
         [y,r] = self.lsq.lsq3(Q*g, w);
         
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function z = Yprod(self, sigma, y, u)
         %YPROD  Compute Y*v, where Y is the Jacobian of y(x)

         nLow = sum(self.nlp.jLow);
         nUpp = sum(self.nlp.jUpp);
         x  = self.val.x;
         % A  = self.val.J';
         Jtprod = self.val.Jtprod;
         gL = self.val.gL;
         Qd = self.val.Qd; % Q^{1/2}*gL
         
         Rop = @(v)self.Rop(x, gL, v);
         Qprod = @(v)self.Qprod(x, v);
         
         % v solves the LS problem
         % minimize  1/2 |QA*v| - <v, u>,
         % ie, A'Q^2 A v = u.
         [v,t] = self.lsq.lsq3(0, -u);

         if nLow == 0 && nUpp == 0
            Av = -t;       % In equality case, this should be more accurate
         else
            Av = Jtprod(v);
         end
         Qgl = Qprod(Qd);
         Stv = self.nlp.hconprod(x, v, Qgl);
         HLQAv = self.nlp.hlagprod(x, y, Qprod(-t));
         z = Rop(Av) + HLQAv - sigma*Av + Stv;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function z = Ytprod(self, sigma, y, v, Hv)
         %YTPROD  Compute Y'v, where Y is the Jacobian of y(x)
         % Y'v solves the linear least-squares problem
         %
         %  minimize 1/2 |A*u - (HL - sigma I)v|^2 + <-SLv, u>.
         %
         % where SLv = S(x, gL)*v.
         
         x  = self.val.x;
         
         if nargin < 5 || isempty(Hv)
             Hv = self.nlp.hlagprod(x, y, v);
         end
         
         % A  = self.val.J';
         Jprod = self.val.Jprod;
         gL = self.val.gL;
         Qd = self.val.Qd;

         Rop = @(v)self.Rop(x, gL, v);
         Qprod = @(v)self.Qprod(x, v);
         SLv = self.nlp.ghivprod(x, Qprod(Qd), v); % S(x, gL)*v
         
         w = SLv + Jprod(Rop(v) - sigma*v);
         
         z = self.lsq.lsq3(Qprod(Hv), -w);
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function z = hprod1(self, v)
         %HPROD1  Form a matrix-vector product with Hessian of merit function.
         %
         % Computes expression:
         %     (Hess phi)*v = HL*v - A*Y'*v - Y*A'*v
         %
         % Requires the use of ghivprod operator.
         
         if self.new_sigma || self.new_delta
            error('Need to evaluate objective/gradient before Hessian.');
         end
         
         A = self.val.Jtprod;
         At = self.val.Jprod;
         x = self.val.x;
         y = self.val.y;
         sig = self.sigma;

         Hv = self.nlp.hlagprod(x, y, v);
         
         Ytv = self.Ytprod(sig, y, v, Hv);
         Ytv(self.iExp) = 0;
         
         Atv = At(v);
         Atv(self.iExp) = 0;
         YAtv = self.Yprod(sig, y, Atv);
         
         AYtv = A(Ytv);
         
         z = Hv - AYtv - YAtv;
       end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function z = hprod2(self, v)
         %HPROD2  Form a matrix-vector product with Hessian of merit function.
         %
         % Computes expression:
         %     (Hess \phi)*v = H*v - P*H*v - H*P*v
         %
         % P = projector only Range(A)
         %
         % Does not use ghivprod.
         % Cannot be used when explicit or bound constraints present.
         
         nLow = sum(self.nlp.jLow);
         nUpp = sum(self.nlp.jUpp);
         
         if self.new_sigma || self.new_delta
            error('Need to evaluate objective/gradient before Hessian.');
         end
         
         if sum(self.linear) > 0 || nLow > 0 || nUpp > 0
             error('Cannot use hprod2 with explicit linear constraints or bounds.');
         end             
         
         A = self.val.Jtprod;
         x = self.val.x;
         y = self.val.y;
         
         Hv = self.nlp.hlagprod(x, y, v);
         Pv = A(self.lsq.lsq3(v));

         PHv = A(self.lsq.lsq3(Hv));
         HPv = self.nlp.hlagprod(x, y, Pv);
         
         z = Hv - PHv - HPv + (2*self.sigma)*Pv;
       end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function z = hprod3(self, v)
         %HPROD  Form a matrix-vector product with Hessian of merit function.
         %
         % Hess \phi = P2 H P2 - P1 H P1 + 2 sigma P1, where
         %
         % P1 = projection onto Range(A);
         % P2 = projection onto complement of Range(A), ie, P2 = I - P1.
         %
         % This is algebraically equivalent to hprod2, which can be seen
         % by simply making the substitution P2 = I-P1. However, this
         % routine requires an additional `lsq3` call.
         % Therefore, hprod3 is not advised, use hprod2 instead.
         % Does not require ghivprod.
         %
         % This is the Hessian approx that is suggested by Fletcher, 1972.
         

         nLow = sum(self.nlp.jLow);
         nUpp = sum(self.nlp.jUpp);
         
         if self.new_sigma || self.new_delta
            error('Need to evaluate objective/gradient before Hessian.');
         end
         
         if sum(self.linear) > 0 || nLow > 0 || nUpp > 0
             error('Must use hprod1 with explicit linear constraints or bounds.');
         end             
                  
         A = self.val.Jtprod;
         x = self.val.x;
         y = self.val.y;
         
         % P1v and P2v
         [w, r] = self.lsq.lsq3(v);
         P1v = A(w);
         P2v = r;
         
         % HP1, HP2
         HP1v = self.nlp.hlagprod(x, y, P1v);
         HP2v = self.nlp.hlagprod(x, y, P2v);
         
         % P1 H P1   and   P2 H P2
         P1HP1v = A(self.lsq.lsq3(HP1v));
         
         [~, r] = self.lsq.lsq3(HP2v);
         P2HP2v = r;
         
         z = P2HP2v - P1HP1v + (2*self.sigma)*P1v;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function z = hprod4(self, v)
         %HPROD  Form a matrix-vector product with Hessian of merit function.
         %
         % Does not require ghivprod.
         % Similar to hprod2, but can handle bound constraints.

         if self.new_sigma || self.new_delta
            error('Need to evaluate objective/gradient before Hessian.');
         end
         
         A = self.val.Jtprod;
         At = self.val.Jprod;
         x = self.val.x;
         y = self.val.y;
         gL = self.val.gL;

         Rop = @(v)self.Rop(x, gL, v);
         Qprod = @(v)self.Qprod(x, v);
         
         Hv = self.nlp.hlagprod(x, y, v);
         
         Atv = At(v);
         Atv(self.iExp) = 0;
         [u,t] = self.lsq.lsq3(0, -Atv);
         
         HQt = self.nlp.hlagprod(x, y, Qprod(t));
         Au = A(u);
         
         YAtv = -HQt + Rop(Au) - self.sigma*Au;
         
         w = At(Rop(v) - self.sigma*v);
         Ytv = self.lsq.lsq3(Qprod(Hv), -w);
         Ytv(self.iExp) = 0;
         AYtv = A(Ytv);
         
         z = Hv - YAtv - AYtv;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function sigma = norm_proj_hess(self)
         %NORM_PROJ_HESS  Norm of the projected Hessian.
         %
         % See PROJECTED_HESS_PROD.
                  
         % Least-squares multiplier estimates (ie, sigma=0).
         y = self.multiplier(self.val.g);

         opts = struct('isreal',true,'issym',true,'tol',1e-1,'disp',0);
         PHPprod = @(v)self.projected_hess_prod(v, y);
         sigma = abs(eigs(PHPprod, self.nlp.n, 1, 'LM', opts));
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function PHPv = projected_hess_prod(self, v, y)
         %PROJECTED_HESS_PROD  Form product with projected Hessian Lag.

         x = self.val.x;
         A = self.val.Jtprod;
         
         % Compute (P*H(x,y)*P)*v
         Pv = A(self.lsq.lsq3(v));
         HPv = self.nlp.hlagprod(x, y, Pv);
         PHPv = A(self.lsq.lsq3(HPv));            
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function nrm = norm_AYt(self)
         %NORM_AYT  Estimate the norm of AYt.
         % Rather than use products with A, we can use the QR factors:
         % |AY'| = |QRY'| = |RY'|
         % because of orthogonality of Q. This allows us to
         % 1. avoid products with A;
         % 2. implicitly use any regularization, ie, delta > 0.
         
         R = self.lsq.R;
         m = size(R,1);
         y = self.multiplier(self.val.g);
         AYtprod = @(v) R*self.Ytprod(0, y, v);
         YAtprod = @(v) self.Yprod(0, y, R'*v);
         nrm = helpers.normest(AYtprod, YAtprod, m, 1e-3);
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function z = Qprod(self, x, v)
         %QPROD  Form product with the operator QOP(x).

         z = self.val.Q * v;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function z = Rop(self, x, g, v)
         %ROP  Form product between the operator R(g) and v.
         %
         % Product with operator R(g)*v = ((d/dx Q(x))*g).*v

         jLow = self.nlp.jLow;
         jUpp = self.nlp.jUpp;
         jTwo = self.nlp.jTwo;
         bL = self.bL;
         bU = self.bU;
         w = self.cwidth;
         n = length(g);
         z = zeros(n,1);
         z(jLow) =  g(jLow) .* v(jLow);
         z(jUpp) = -g(jUpp) .* v(jUpp);
         
         xx = x(jTwo);
         bl = bL(jTwo);
         bu = bU(jTwo);
         w  = w(jTwo);
         
         zt = sign((bu-xx) - (xx-bl));
         ix = abs(2*xx - bl - bu) < w;
         zt(ix) = 1 - (1./w(ix)).*(2*xx(ix) - bl(ix) - bu(ix) + w(ix));
         
         z(jTwo) = zt .* g(jTwo) .* v(jTwo);
      end
         
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function Q = Qop(self, x)
         %QOP Form square-root of Q operator
         %
         % Q(x)(i,i) = 1               ,   if bL = -inf, bU = inf
         %           = sqrt(x_i - bL_i),   if bU = inf
         %           = sqrt(bU_i - x_i),   if bL = -inf
         %           ~ sqrt( min( x_i - bL_i, bU_i - x_i ) ),  else
         %
         % We actually use a smooth version of the last case above

         bL = self.bL;
         bU = self.bU;
         jLow = self.jLow;
         jUpp = self.jUpp;
         jTwo = self.jTwo;
         nLow = sum(jLow);
         nUpp = sum(jUpp);
         nTwo = sum(jTwo);
         w = self.cwidth(jTwo);
         
         Q = speye(self.nlp.n);
         Q(jLow, jLow) = spdiags(sqrt(x(jLow) - bL(jLow)), 0, nLow, nLow);
         Q(jUpp, jUpp) = spdiags(sqrt(bU(jUpp) - x(jUpp)), 0, nUpp, nUpp);
         
         xx = x(jTwo);
         bl = bL(jTwo);
         bu = bU(jTwo);
         
         q = min([(xx-bl) (bu-xx)], [], 2);
         r = (xx-bl) - (0.25./w).*((2*xx - bl - bu + w).^2);
         ix = abs(2*xx - bl - bu) < w;
         q(ix) = r(ix);

         Q(jTwo, jTwo) = spdiags(sqrt(q), 0, nTwo, nTwo);
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function self = set_sigma(self, val)
         self.sigma = val;
         self.new_sigma = true;
      end
      
      function set.delta(self, val)
         self.new_delta = true;
         self.lsq.delta = val;
      end
      
      function set.delta_max(self, val)
         self.lsq.delta_max = val;
      end
      
      function val = get.delta_max(self)
         val = self.lsq.delta_max;
      end
      
      function val = get.delta(self)
         val = self.lsq.delta;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function c = fcon_local(self, x)
          hash = utils.DataHash(x);
          if strcmp(hash, self.last_cx_hash) ...
            && ~self.new_sigma && ~self.new_delta
             c = self.val.c_nlp(self.iExp);
             return
          end
          
          self.last_cx_hash = hash;

          c = self.nlp.fcon(x);
          c = c(self.iExp);
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function J = gcon_local(self, x)
          hash = utils.DataHash(x);
          if strcmp(hash, self.last_Jx_hash) ...
            && ~self.new_sigma && ~self.new_delta
             J = self.val.J(self.iExp,:);
             return
          end
          
          self.last_Jx_hash = hash;
          
          % TODO: Just check if it's the same x as the current value?
          J = self.nlp.gcon(x);
          J = J(self.iExp,:);
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
   end
   
   methods(Static)
      function result = formulated_correctly(nlp)
         % FORMULATED_CORRECTLY
         %
         % Check that nlp is only has equality constraints and bounds,
         % that constraints are equality only.

         result = (sum(nlp.iFix) == nlp.m);
      end
   end
   
end % classdef
