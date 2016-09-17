#ifndef ___CG_H___
#define ___CG_H___

#include "Common.h"
#include "MatrixOperation.h"

namespace AlgebraicMultigrid
{
    template <class preconditioner>
    void CG(integer n, preconditioner& M,   // dim, preconditioner
        intVec& ptr, intVec& ind, reaVec& val,  // matrix A
        reaVec& x, reaVec& b,                   // rhs and x
        real& tolerance, integer& maxiter)
    {
        reaVec R(n, 0.0);   double* r = R.data();
        reaVec S(n, 0.0);   double* s = S.data();
        reaVec P(n, 0.0);   double* p = P.data();
        reaVec Q(n, 0.0);   double* q = Q.data();

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

            for (int i = 0; i < n; i++) p[i] = rho1 / rho2*p[i] + s[i];

        }


    }
}

#endif