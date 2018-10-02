Introduction
============

This library is for solving problems of the form

.. math ::
	
	\min_{x \in \mathbb{R}^n} \enspace f(x) \enspace\mbox{subject to}\enspace c(x)=0, \enspace \ell \le x \le u,

where :math:`f(x)` and :math:`c(x)` have two continuous derivatives. The optimization problem is minimized by instead minimizing a smooth exact penalty function introduced by Fletcher in the 1970's:

.. math::

	\min_{x \in \mathbb{R}^n} \enspace \phi_{\sigma}(x) \enspace\mbox{subject to}\enspace \ell \le x \le u.

The penalty function is

.. math ::
	\phi_{\sigma}(x) &:= f(x) - c(x)^T y_{\sigma}(x) + \tfrac{1}{2} \rho \|c(x)\|_2^2 \\
	y_{\sigma}(x) &:= \arg\min_y \tfrac{1}{2} \|g(x) - A(x) y\|^2_{Q(x)} + \sigma c(x)^T y

where :math:`g(x) = \nabla f(x)`, :math:`A(x) = \nabla c^T (x)`, and :math:`Q(x)` is a diagonal matrix function with entries approximating :math:`\min \{ x - l, u-x \}`.