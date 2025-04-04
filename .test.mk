# mamba run -n vfs vfs_build run tests/cpp/test_Density.cpp
# mamba run -n vfs vfs_build run tests/cpp/test_Solver.cpp
all:
	mamba run -n vfs vfs_build run tests/cpp/test_System.cpp

mu:
	mamba run -n vfs vfs_build run tests/applications/mumble_d2d.cpp

py:
	export PYTHONPATH=$$PYTHONPATH:$$PWD/src/python/; python tests/python/test_ot.py
