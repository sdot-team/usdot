# mamba run -n vfs vfs_build run tests/cpp/test_Density.cpp
# mamba run -n vfs vfs_build run tests/cpp/test_Solver.cpp

all:
	export PYTHONPATH=$$PYTHONPATH:$$PWD/src/python/; python tests/python/test_ot.py
