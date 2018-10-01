function [result, eabs, erel] = gradient_check(nlp, x, varargin)
    % Checks gradient against centered finite difference approximation.
    % If check_all=true, then check full gradient, otherwise pick its
    % random vectors and check the directional derivative against those
    % vectors.
    % result = true if gradient is sufficiently accurate.
    % If check_all=true, then errs is a vector of the relative errors.
    % Otherwise it is the largest error encountered.
    
    p = inputParser;
    p.KeepUnmatched = false;
    p.addParameter('atol', 1e-6); % Absolute tolerance
    p.addParameter('rtol', 1e-6); % Relative tolerance
    p.addParameter('check_all', false); % Check full gradient?
    p.addParameter('its', 5); % How many random vectors if not full grad
    p.parse(varargin{:});

    atol = p.Results.atol;
    rtol = p.Results.rtol;
    check_all = p.Results.check_all;
    its = p.Results.its;
    n = nlp.n;
    
    result = true;
    g = nlp.gobj(x);
    h = zeros(n,1);
    
    % Optimal-ish step for second-order centered finite differences.
    step = (eps / 3)^(1/3);
    if check_all
        erel = zeros(n,1);
        eabs = zeros(n,1);
        for i=1:n
            h(i) = step;
            dfdxi = (nlp.fobj(x+h) - nlp.fobj(x-h))/(2*step);
            err = abs(dfdxi - g(i));
            if err > atol + rtol*abs(dfdxi)
                result = false;
            end
            erel(i) = err/abs(dfdxi);
            eabs(i) = err;
            h(i) = 0;
        end
    else
        erel = 0;
        eabs = 0;
        for i=1:its
            h = 2*rand(n,1)-1;
            h = h/norm(h);
            
%             while ( sum((x+step*h - nlp.bL< 0) | (nlp.bU - (x+step*h) < 0) | ...
%                (x-step*h - nlp.bL< 0) | (nlp.bU - (x-step*h) < 0)) > 0)
%                 h( (x+step*h - nlp.bL< 0) | (nlp.bU - (x+step*h) < 0) | ...
%                    (x-step*h - nlp.bL< 0) | (nlp.bU - (x-step*h) < 0) ) = 0;
% 
%                 h = h/norm(h);
%             end
            
            dfdx = (nlp.fobj(x+step*h) - nlp.fobj(x-step*h))/(2*step);
            err = abs(dfdx - g'*h);
            if err > atol + rtol*abs(dfdx)
                result = false;
            end
            erel = max(erel, err/abs(dfdx));
            eabs = max(eabs, err);
        end
    end
end