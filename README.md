usdot is a C++/python library for unidimensional semi-discrete (partial) optimal transport.

Installation
===========

After the `git clone`, one has to run `make get_deps` to get the dependencies (that will appear in the `ext/` directory).

Python examples
---------------

The python examples need the library `cpp_import` (`pip import cppimport`) with a `PYTHONPATH` or `sys.path` pointing to the sources. For example:

```export PYTHONPATH=$PYTHONPATH:$PWD/src/python/; python tests/python/test_ot.py```

The solver procedures return an instance of `OtResult` which contains
    * `norm_2_residual_history (array)`: ...
    * `error_message (str)`: empty if ok
    * `barycenters (array)`: ...
    * `boundaries (array)`: left and right cell boudaries
    * `weights (array)`: the Kantorovitch potentials
    * `masses (array)`: the final cell masses (which of course must be the prescribed ones if the solver did not fail)

As can be seen in `src/python/usdot/__init__.py`, they have in common a parameter `ot_parms` which can be `None` or an instance of `OtParms`. `OtParms` contains the attributes
    * max_mass_ratio_error_target (1e-6 by default): target for `max( ( mass - prescribed_mass ) / prescribed_mass )`
    * epsilon: the regularization parameter as explained in the upcoming publication...

They also have parameters 
    * mass_ratio, an optional scalar between 0 and 1 to defined the mass of the diracs divided by the mass of the density
    * dirac_mass_ratios, an optional vector, that can be used to set individual dirac mass values (else, dirac masses are considered as uniform)

Currently, there are 2 solver procedures:
    * `usdot.d2cap( diracs, dens_p, dens_v, ... )` is made to compute ot plan between a set of diracs and continous piecewise affine density. `dens_p` are the "stops", `dens_v` are the values at the stops.
    * `usdot.d2p( diracs, dens_p, dens_v, ... )` is made to compute ot plan between a set of diracs and piecewise polynomial density. Polynomial order is determined by `dens_v.shape[ 1 ] - 1`. Each row gives the coefficients of a polynomial where x = 0 at the beginning of the interval and x = 1 at the end of the interval (intervals are determined by `dens_p`).


C++ examples
------------

The C++ examples can be run using `vfs_build` (available for instance after `pip import vfs`). For example:

```vfs_build run tests/cpp/test_Solver.cpp```

(assuming there is vfs_config.py with the proper compilation flags, as in the tests of the repository)
