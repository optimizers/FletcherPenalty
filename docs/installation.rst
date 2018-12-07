Installation
============

To install, clone the repository,

::

	git clone https://github.com/restrin/FletcherPenalty.git

and add it to the matlab path.

Dependencies
------------

**Required**

* Model_ library
* SuiteSparse_ (specifically, the Matlab interface for ``spqr``)
* BCFLASH_

**Optional**

* LinearSystemSolvers_
	* A collection of iterative methods for linear systems
	* LNLQ is required in particular
* CUTEst_
	* Interface for using CUTEst test problems
* Intrelab_ package from Trilinos_
	* Useful for some example problems

.. _Model: https://github.com/optimizers/model
.. _SuiteSparse: http://faculty.cse.tamu.edu/davis/suitesparse.html
.. _BCFLASH: https://github.com/restrin/bcflash
.. _CUTEst: https://github.com/ralna/CUTEst
.. _LinearSystemSolvers: https://github.com/restrin/LinearSystemSolvers
.. _Intrelab: https://github.com/trilinos/Trilinos/tree/master/packages/intrepid/matlab/intrelab
.. _Trilinos: https://trilinos.org/
