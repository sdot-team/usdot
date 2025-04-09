Usdot is a C++/python library for unidimensional semi-discrete partial optimal transport.

Python examples
===============

The python examples need the library `cpp_import` (`pip import cppimport`) with a `PYTHONPATH` or `sys.path` pointing to the sources. For example:

```export PYTHONPATH=$PYTHONPATH:.../cloned_dir/src/python/; python test_ot.py```

The solver procedures return an instance of `OtResult` which contains
    * `error_message (str)`: a string describing the problem if an error occured (empty if no error)
    * `barycenters (array)`: the cell barycenter coordinates
    * `boundaries (array)`: left and right cell boudaries (shape = [nb_cells, 2])
    * `weights (array)`: the Kantorovitch potentials
    * `masses (array)`: the final cell masses (which must correspond to the dirac masses if no error occured)

As can be seen in `src/python/usdot/__init__.py`, they have in common a parameter `ot_parms` which can be `None` or an instance of `OtParms`. `OtParms` contains the attributes
    * `max_mass_ratio_error_target` (1e-6 by default): target for `max( ( mass - prescribed_mass ) / prescribed_mass )`
    * `epsilon`: the regularization parameter (to obtain differentiability).

The other common parameters are:
    * `global_mass_ratio`, an optional scalar between 0 and 1 to define the mass of the diracs divided by the mass of the density
    * `relative_mass_ratios`, an optional 1D array, that can be used to set individual dirac mass values (else, dirac masses are considered as uniform). This vector is automatically normalized (it does not change the global_mass_ratio).

Currently, there are 1 solver procedure:
    * `usdot.from_p1_grid( dirac_positions, density_values, density_beg, density_end, ... )`: compute an ot plan between a set of diracs and a continous piecewise affine density defined on a regular grid. By default, `density_end` is equal to `len( density_values ) - 1`
    <!-- * `usdot.d2p( diracs, dens_p, dens_v, ... )` is made to compute ot plan between a set of diracs and piecewise polynomial density. Polynomial order is determined by `dens_v.shape[ 1 ] - 1`. Each row gives the coefficients of a polynomial where x = 0 at the beginning of the interval and x = 1 at the end of the interval (intervals are determined by `dens_p`). -->


C++ examples
------------

The C++ examples can be run using `vfs_build` (available for instance after `pip import vfs`). For example:

```vfs_build run tests/cpp/test_Solver.cpp```

