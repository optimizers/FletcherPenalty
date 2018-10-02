Introduction
============

This library is for solving problems of the form

.. math ::
	
	\min_{x \in \mathbb{R}^n} \enspace f(x) \enspace\mbox{subject to}\enspace c(x)=0, \enspace \ell \le x \le u,

where :math:`f : \mathbb{R}^n \rightarrow \mathbb{R}` and :math:`c : \mathbb{R}^n \rightarrow \mathbb{R}^m` have two continuous derivatives. The optimization problem is minimized by instead minimizing a smooth exact penalty function introduced by Fletcher in the 1970's:

.. math::

	\min_{x \in \mathbb{R}^n} \enspace \phi_{\sigma}(x) \enspace\mbox{subject to}\enspace \ell \le x \le u.

The penalty function is

.. math ::
	\phi_{\sigma}(x) &:= f(x) - c(x)^T y_{\sigma}(x) + \tfrac{1}{2} \rho \|c(x)\|_2^2 \\
	y_{\sigma}(x) &:= \arg\min_y \tfrac{1}{2} \|g(x) - A(x) y\|^2_{Q(x)} + \sigma c(x)^T y

where :math:`g(x) = \nabla f(x)`, :math:`A(x) = \nabla c^T (x)`, and :math:`Q(x)` is a diagonal matrix function with entries approximating :math:`\min \{ x - l, u-x \}`.

Problem Formulation and Notation
--------------------------------

Problems must be passed with equality constraints and slacks, in the form described above. This library uses the optimizers/model_ library, which contains models which will automatically construct the slack formulation of an arbitrary nonlinear program (``nlpmodel``).

Additional notation to be aware of when implementing your own problem as an nlpmodel:

* The Lagrangian is defined with a **minus** sign:

.. math::

	L(x,y) = f(x) - c(x)^T y

Therefore, ``hlag``, ``hlagprod``, and other methods involving the Lagrangian should be implemented accordingly.

* The Jacobian of :math:`c` is an :math:`m \times n` matrix, so ``gcon``, ``gconprod``, and other methods involving the Jacobian should be implemented accordingly.

.. _model: https://github.com/optimizers/model