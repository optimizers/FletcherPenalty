Input
=====

The basic solver is called via

::

	% nlp is instance of nlpmodel class
	fs = fletcher_solver(nlp, varargin)

**Note**: the original problem (nlp) should be expressed in slack form (equality constraints only).

fletcher_solver inputs
----------------------

================= ======== ================ ===========
Name              Required Default          Description
================= ======== ================ ===========
nlp               required                  Optimization problem to be solved (nlpmodel)
fid               optional 1                Save iteration log in subsolver.log?
logLevel          optional 1                0 (don't print log), 1 (print log)
nlnFlag           optional 1                Keep linear constraints explicit?
sigma             optional 1                Penalty parameter
sigma_min         optional :math:`10^{-6}`  Minimum penalty parameter
sigma_max         optional :math:`10^6`     Maximum penalty parameter
sigma_strategy    optional SIGMA_ADAPTIVE   SIGMA_ADAPTIVE or SIGMA_FIXED
delta             optional :math:`10^{-8}`  Regularization parameter (ignored if lsq.regularized = false)
delta_min         optional :math:`10^{-8}`  Minimum regularization parameter
delta_dec         optional ``@(d) d/10``    Function to update delta
lsq_options       optional ``struct()``     See below
lsq_method        optional LSQ_SNE          Method for solving augmented system
subsolver         optional BCFLASH          Solver for minimizing penalty function
optTolAbs         optional :math:`10^{-6}`  Absolute tolerance for Lagrangian gradient
optTolRel         optional :math:`10^{-6}`  Relative tolerance for Lagrangian gradient
feaTolAbs         optional :math:`10^{-6}`  Absolute tolerance for primal feasibility
feaTolRel         optional :math:`10^{-6}`  Relative tolerance for primal feasibility
max_iterations    optional 500              Maximum number of outer iterations
check_grad        optional false            Run gradient checker every iteration?
subsolver_options optional ``struct()``     See below
merit_options     optional ``struct()``     See below
x0                optional nlp.x0           Initial point
================= ======== ================ ===========

lsq_options
-----------

Depends on ``lsq_method``:

* LSQ_SNE
	================= ================ ===========
	Name              Default          Description
	================= ================ ===========
	regularized       true             Regularize penalty function?
	delta_max         :math:`10^{-4}`  Maximum regularization parameter
	cond_max          :math:`10^8`     Maximum allowable condition number of :math:`A`
	================= ================ ===========

* LSQ_LDL
	================= ================ ===========
	Name              Default          Description
	================= ================ ===========
	regularized       true             Regularize penalty function?
	delta_max         :math:`10^{-4}`  Maximum regularization parameter
	cond_max          :math:`10^8`     Maximum allowable condition number of :math:`A`
	static_p          false            Use static reordering? If true, LDL is faster but possibly unstable.
	================= ================ ===========

* LSQ_MINRES
	================= ================ ===========
	Name              Default          Description
	================= ================ ===========
	regularized       true             Regularize penalty function?
	delta_max         :math:`10^{-4}`  Maximum regularization parameter
	tol               :math:`10^{-8}`  Residual norm at tolerance
	================= ================ ===========

* LSQ_LNLQ
	===================== ================ ===========
	Name                  Default          Description
	===================== ================ ===========
	regularized           true             Regularize penalty function?
	delta_max             :math:`10^{-4}`  Maximum regularization parameter
	termination_condition TOL_RESIDUAL     TOL_RESIDUAL or TOL_ERROR
	tol                   :math:`10^{-8}`  Residual or error norm at tolerance
	===================== ================ ===========

subsolver_options
-----------------

Depends on ``subsolver``:

* BCFLASH
	See BCFLASH for optional arguments.

* KNITRO
	===================== ================ ===========
	Name                  Default          Description
	===================== ================ ===========
	x0                    nlp.x0           Initial point
	opt_file              []               Options file for KNITRO
	===================== ================ ===========

* TRPCG
	None for now.

* IPOPT
	===================== ================ ===========
	Name                  Default          Description
	===================== ================ ===========
	x0                    nlp.x0           Initial point
	mu_strategy           adaptive         Strategy for updating barrier parameter (see IPOPT docs)
	start_feasible        true             Move to feasible point first?
	fid                   ``''``           Output file for IPOPT log
	===================== ================ ===========

* SNOPT
	===================== ================ ===========
	Name                  Default          Description
	===================== ================ ===========
	x0                    nlp.x0           Initial point
	snfid                 ``''``           Output file for SNOPT log
	===================== ================ ===========


merit_options
-------------
	===================== ================ ===========
	Name                  Default          Description
	===================== ================ ===========
	rho                   0                Penalty parameter for :math:`\ell_2`-penalty
	hprod                 1                Hessian product approximation (1-4)
	lin_explicit          false            Keep linear constraints explicit?
	===================== ================ ===========


Output
======

``fs = fletcher_solver(...)`` returns an instance of the fletcher_solver object, which will contain the following additional fields:

===================== ===========
Name                  Description
===================== ===========
sol.x                 Primal solution
sol.y                 Dual solution
sol.f                 Final objective value
sol.info              Additional subsolver-dependent information
exit                  Exit flag
exit_msg              Exit message
===================== ===========
