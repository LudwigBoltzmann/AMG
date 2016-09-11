#ifndef ___CG_H___
#define ___CG_H___

#include "Common.h"
#include "MatrixOperation.h"

namespace AlgebraicMultigrid
{
    template<class preconditioner>
    pair_i_r CG(integer n, preconditioner& M, intVec& ptr, intVec& ind, reaVec& val, reaVec& x, reaVec& b)
    {
        pair_i_r ret;

        integer one = 1;
        integer zero = 0;


        return ret;
    }
}

#endif