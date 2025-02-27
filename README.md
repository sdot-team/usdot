usdot is a C++/python library for unidimensional semi-discrete (partial) optimal transport.

Installation
===========

After the `git clone`, one has to run `make get_deps` to get the dependencies (that will appear in the `ext/` directory).

Python examples
---------------

The python examples need the library `cpp_import` (`pip import cppimport`) with a `PYTHONPATH` or `sys.path` pointing to the sources. For example:

```export PYTHONPATH=$PYTHONPATH:$PWD/src/python/; python tests/python/test_ot.py```

C++ examples
------------

The C++ examples can be run using `vfs_build` (available for instance after `pip import vfs`). For example:

```vfs_build run tests/cpp/test_Solver.cpp```

(assuming there is vfs_config.py with the proper compilation flags, as in the tests of the repository)