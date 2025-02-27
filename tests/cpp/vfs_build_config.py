import vfs, sysconfig, numpy

def config( options ):
    options.load_lib( 'https://github.com/catchorg/Catch2.git', lib_names = [ "Catch2Main", "Catch2" ] )

    # vfs.vfs_build_config( options )

    options.add_cpp_flag( '-Wno-deprecated-declarations' )
    options.add_cpp_flag( '-std=c++20' )
    options.add_cpp_flag( '-fPIC' )
    options.add_cpp_flag( '-g3' )
    options.add_cpp_flag( '-O3' )

    options.add_inc_path( '/opt/homebrew/Cellar/eigen/3.4.0_1/include/' )
    options.add_inc_path( '/opt/homebrew/Cellar/boost/1.87.0/include' )
    options.add_inc_path( '../../ext/matplotlib-cpp/' )
    options.add_inc_path( '../../ext/tl20/src/cpp' )
    options.add_inc_path( '../../src/cpp' )
    
    # FLAGS := --cpp-flag=-I/home/leclerc/.vfs_build/ext/Catch2/install/include
    # FLAGS += --cpp-flag=-I/opt/homebrew/Cellar/eigen/3.4.0_1/include/
    # FLAGS += --cpp-flag=-I/opt/homebrew/Cellar/boost/1.87.0/include
    # FLAGS += --cpp-flag=-I/opt/homebrew/Cellar/gmp/6.3.0/include/
    # FLAGS += --cpp-flag=-Iext/tl20/src/cpp/
    # FLAGS += --cpp-flag=-Isrc/cpp/

    # FLAGS += --cpp-flag=-std=c++20
    # FLAGS += --cpp-flag=-O2
    # FLAGS += --cpp-flag=-g3

    # FLAGS += --cpp-flag=-I/usr/X11R6/include
    # FLAGS += --lib-flag=-L/usr/X11R6/lib
    # FLAGS += --lib-flag=-lpthread
    # FLAGS += --lib-flag=-lX11
    # FLAGS += --lib-flag=-lm

    # options.add_lib_path( '/opt/homebrew/opt/gmp/lib' )
    # options.add_lib_name( 'gmp' )

    # options.add_lib_name( sysconfig.get_config_var( 'LIBRARY' ).replace( '.a', '.dylib' ).replace( 'lib', '' ) )

    # options.add_lib_flag( '-rpath' )
    # options.add_lib_flag( sysconfig.get_config_var( 'LIBDIR' ) )

    # options.add_lib_name( 'python3.13' )
    # options.add_lib_path( sysconfig.get_config_var( 'LIBDIR' ) )
    # options.add_inc_path( sysconfig.get_path( "include" ) )
    # options.add_inc_path( numpy.get_include() )
