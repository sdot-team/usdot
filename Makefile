eigen_version = 3.4.0

get_deps:
	test -d ext || mkdir ext
	test -d ext/tl || ( cd ext && git clone https://github.com/hleclerc/tl20.git && ln -s tl20/src/cpp/tl tl )
	test -d ext/eigen3 || ( cd ext && curl https://gitlab.com/libeigen/eigen/-/archive/${eigen_version}/eigen-${eigen_version}.zip --output eigen-${eigen_version}.zip && unzip eigen-${eigen_version}.zip && ln -s eigen-${eigen_version} eigen3 )
	test -d ext/pybind11 || ( cd ext && curl -LO https://github.com/pybind/pybind11/archive/refs/tags/v2.13.6.zip && unzip v2.13.6.zip && rm v2.13.6.zip && ln -s pybind11-2.13.6/include/pybind11 pybind11 )
