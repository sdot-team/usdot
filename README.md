usdot is a C++/python library for unidimensional semi-discrete optimal transport.

Compilation
===========

After the git clone, one has to run 

```make get_deps```

to get the dependencies.

The C++ examples can be run using `vfs_build` (pip import vfs). For example:

```vfs_build run tests/cpp/test_Solver.cpp```

The python examples need `cpp_import` (`pip import cppimport`) with a `PYTHONPATH` or `sys.path` pointing to the sources. For example:

```export PYTHONPATH=$PYTHONPATH:DIR_OF_USDOT/src/python/; python tests/python/test_ot.py```

