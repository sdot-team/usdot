# mamba run -n vfs vfs_build run tests/cpp/test_Density.cpp
# mamba run -n vfs vfs_build run tests/cpp/test_Solver.cpp
all:
	micromamba run -n vfs vfs_build run tests/cpp/test_UsdotSystem.cpp

mu:
	micromamba run -n vfs vfs_build run tests/applications/mumble_d2d.cpp

py:
	# export PYTHONPATH=$$PYTHONPATH:$$PWD/src/python/; python tests/python/test_ot.py
	export PYTHONPATH=$$PYTHONPATH:$$PWD/src/python/; python test.py
