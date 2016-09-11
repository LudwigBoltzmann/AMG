#include "LinearEqnSolver.h"


namespace AlgebraicMultigrid
{
    template <class preconditioner>
    pair_i_r CG(integer n, preconditioner& M,   // dim, preconditioner
        intVec& ptr, intVec& ind, reaVec& val,  // matrix A
        reaVec& x, reaVec& b,                   // rhs and x
        real tolerance, integer maxiter)
    {
        pair_i_r ret;

        integer one = 1;
        integer zero = 0;

        reaVec r(n, 0.0);
        reaVec s(n, 0.0);
        reaVec p(n, 0.0);
        reaVec q(n, 0.0);

        residual(n, ptr, ind, val, x, r);
        real normb = norm(b);

        if (normb < eps) {
            ret.first = 0;
            ret.second = normb;
            return ret;
        }
        
        real eps = normb * tolerance;
        real eps2 = tolerance * tolerance;

        real rho1 = 2.0 * eps2;
        real rho2 = 0.0;

        real normr = norm(r);

        integer iter = 0;

        for (; iter < maxiter && norm < abs(normr); iter++) {
            M.solve(r, s);
            rho2 = rho1;
            rho1 = innerProduct(r, s);
        }


    }
}