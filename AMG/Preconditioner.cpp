#include "Preconditioner.h"

namespace AlgebraicMultigrid
{

    /*
        Empty preconditioner
    */
    Preconditioner::Preconditioner(BSR* A)
    {

    }

    void inverse()
    {

    }

    void solve(reaVec& x, reaVec& y)
    {
#pragma omp single
        y.assign(x.begin(), x.end());
    }
}