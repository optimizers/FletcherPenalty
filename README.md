# FletcherPenalty

## Overview

Fletcher's penalty function is a smooth exact penalty method for solving
problems of the form

`` min_x f(x) subject to c(x) = 0, l <= x <= u``

It is based on the smooth penalty function originally introduced by Fletcher in
the 1970's [1].

## Installation

### Dependencies

FletcherPenalty requires:
  - <a href="https://github.com/optimizers/model">optimizers/model</a> to define the optimization problems,
  - <a href="http://faculty.cse.tamu.edu/davis/suitesparse.html">SuiteSparse</a> for sparse factorizations, and
  - <a href="https://github.com/restrin/bcflash">BCFLASH</a> as one possible subsolver.

### Installation

Once the dependencies are installed, then installation involves cloning this repository and adding it (and the dependencies) to Matlab's path.

## References
[1] R. Fletcher. A class of methods for nonlinear programming with termination and convergence properties. In J. Abadie, editor, *Integer and nonlinear programming*, pages 157-175. North-Holland, Amsterdam, 1970.
