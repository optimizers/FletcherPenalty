Code Overview
=============

The code is organized into several packages that handle various aspects. 

+solvers
--------

This package contains interfaces to several solvers that can be used to minimize ``\phi_{\sigma}``. Each interface must implement 2 methods (other than a constructor):

* ``[fletcher, info] = solve()``

	This calls the subsolver on the penalty function, and returns the ``fletcher_solver`` object, and an ``info`` struct:

		* ``info.sol.x``: primal solution
		* ``info.sol.y``: dual solution
		* ``info.sol.f``: final objective value

* ``[...] = post_iteration(...)``

	This function should get called at the end of every iteration of the subsolver. It must at some point call ``fletcher_solver.post_iteration()``.

Currently supported solvers:

	* BCFLASH (equality constraints only)
	* TRPCG (equality constraints and explicit linear constraints)
	* KNITRO
	* IPOPT (not recommended)
	* SNOPT (not recommended)

+least_squares
--------------

This package contains interfaces for solving linear systems of the form

.. math ::

	\begin{bmatrix} I & A \\ A^T & -\delta^2 I \end{bmatrix} \begin{bmatrix} p \\ q \end{bmatrix} = \begin{bmatrix} u \\ v \end{bmatrix}.

Each interface must implement 2 methods (other than a constructor):

* ``self = preprocess_local(self, A, Aprod, Atprod, Q, options)``
* ``[q,p] = lsq_local(u, v)``  

Currently supported linear system solvers:
	
	* Semi-normal equations (``lssne``)

		* Requires that Jacobian is explicitly available

	* LDL (``lsldl``)

		* Requires that Jacobian is explicitly available

	* LNLQ (``lslnlq``)
	* MINRES (``lsminres``) (not recommended)

+utils and +helpers
-------------------

Contains various utilities and helper functions used throughout.