%FLETCHER_CUTEST_DRIVER Driver for running Fletcher on CUTEst test problems
%
% PREREQUISITES:
%    Requires that cutest-mirror be installed from here:
%       https://github.com/optimizers/cutest-mirror
%    and that CUTEst problems are compiled
%
%    Requires that directory containing +model/ is added to Matlab's path

clear all
import model.*

% Set path ================================================================

% Path to directory containing compiled CUTEst problems
% For example, if problem 3pk is to be run, then the compiled problem
% should be in directory: cutest_path/3pk
problem_path = '/home/restrin/Libraries/cutest/problems/compiled';
% Problem name
pname = 'hs78';

fname = [problem_path '/' pname];

% Options for Fletcher solver =============================================

% Penalty parameter
sigma          = 1e2;
sigma_min      = min(sigma, 1e-6);
sigma_max      = 1e7;
rho            = 0;
sigma_strategy = fletcher_solver.SIGMA_FIXED;

% Regularization
regularized = false;
delta       = 0;
delta_max   = 1e-1;
delta_min   = 1e-8;
delta_dec   = @(d) d/10;

% Tolerances
optTolAbs = 1e-8;
feaTolAbs = 1e-8;
optTolRel = 1e-8;
feaTolRel = 1e-8;

% Maximum iterations
max_iterations = 1000;

% Augmented system solver
lsq_method = fletcher_solver.LSQ_SNE;

% Subsolver
subsolver = fletcher_solver.BCFLASH;

% Hessian approximation
hprod = 4;

% Log file (set to '' to log to console)
log_file = '';

% Subsolver specific options
switch subsolver
    case fletcher_solver.BCFLASH
        subsolver_options = struct('exit_user_only', true, ...
                           'frtol', 0);
        lin_explicit = false;
    case fletcher_solver.KNITRO
        subsolver_options = struct();
        lin_explicit = true;
    case fletcher_solver.IPOPT
        start_feasible = false;
        mu_strategy    = 'adaptive';
        subsolver_options = struct('mu_strategy', mu_strategy, ...
                           'start_feasible', start_feasible);
        lin_explicit = true;
    case fletcher_solver.SNOPT
        subsolver_options = struct();
        lin_explicit = true;
    case fletcher_solver.TRPCG
        subsolver_options = struct();
        lin_explicit = true;
end

% Begin solve procedure ===================================================

% Ensure reproducibility
s = rng;
rng('default');

% Build problem
nlp = cutestmodel(fname, true);
p = slackmodel(nlp);

% Set initial point between bounds
n = p.n;
bL = p.bL;
bU = p.bU;
jLow = p.jLow;
jUpp = p.jUpp;
if ~(sum(p.x0 > bL) == n) || ~(sum(p.x0 < bU) == n)
    fprintf('Changing initial point\n');
    p.x0(p.x0 <= bL) = bL(p.x0 <= bL) + 1;
    p.x0(p.x0 >= bU) = bU(p.x0 >= bU) - 1;
    p.x0(p.jTwo & ((p.x0 <= bL) | (p.x0 >= bU))) = ...
        0.5*(bL(p.jTwo & ((p.x0 <= bL) | (p.x0 >= bU))) ...
        + bU(p.jTwo & ((p.x0 <= bL) | (p.x0 >= bU))));
end

% Construct solver
fs = fletcher_solver(p, ...
                     'subsolver', subsolver, ...
                     'log_file', log_file, ...
                     'log_level', 1, ...
                     'delta', delta, ...
                     'delta_min', delta_min, ...
                     'sigma', sigma, ...
                     'sigma_min', sigma_min, ...
                     'sigma_max', sigma_max, ...
                     'optTolAbs', optTolAbs, ...
                     'feaTolAbs', feaTolAbs, ...
                     'optTolRel', optTolRel, ...
                     'feaTolRel', feaTolRel, ...
                     'max_iterations', max_iterations, ...
                     'subsolver_options', subsolver_options, ...
                     'delta_dec', delta_dec, ...
                     'lsq_method', lsq_method, ...
                     'sigma_strategy', sigma_strategy, ...
                     'merit_options', struct('lin_explicit', lin_explicit, ...
                                             'rho', rho, ...
                                             'hprod', hprod), ...
                     'lsq_options', struct('delta_max', delta_max, ...
                                           'regularized', regularized));

% Reset random seed
rng(s);

% Check solution ==========================================================

% Recover solution
x = fs.sol.x;
y = fs.sol.y;

x = x(~p.islack);

c = nlp.fcon(x);
g = nlp.gobj(x);
J = nlp.gcon(x);

pr_feas = nlp.prResidual(x,c);
du_feas = nlp.duResidual(x,c,g,J,y);
fprintf('\n%15s : its = %4i : fea = %12.5e : opt = %12.5e\n', ...
    pname, fs.iteration, pr_feas, du_feas);
