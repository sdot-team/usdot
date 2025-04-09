#include <partial1D/InitialSolutionDD.h>
#include "catch_main.h"

using TF = InitialSolutionDD::TF;

TEST_CASE( "Init sol DD", "" ) {
    Vec<TF> d0;
    for( TF i = 0; i < 10; ++i )
        d0 << TF( rand() ) / RAND_MAX;

    Vec<TF> d1;
    for( TF i = 0; i < 5; ++i )
        d1 << TF( rand() ) / RAND_MAX;

    InitialSolutionDD inso = make_initial_solution_dd( d0, d1 );
}
