#ifndef ___PRECONDITIONER_H___
#define ___PRECONDITIONER_H___

#include "Common.h"
#include "Matrix.h"

namespace AlgebraicMultigrid
{
    class Preconditioner
    {
    protected:

    public:
        Preconditioner(BSR* A);

        void inverse();
        void solve(reaVec& x, reaVec& y);
    };


}

#endif