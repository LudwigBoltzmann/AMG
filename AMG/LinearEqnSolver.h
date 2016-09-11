#ifndef ___CG_H___
#define ___CG_H___

#include "Common.h"
#include "MatrixOperation.h"

namespace AlgebraicMultigrid
{
    template<class preconditioner>
    pair_i_r CG(integer n, preconditioner& M,   // dim, preconditioner
        intVec& ptr, intVec& ind, reaVec& val,  // matrix A
        reaVec& x, reaVec& b,                   // rhs and x
        real eps);
}

#endif