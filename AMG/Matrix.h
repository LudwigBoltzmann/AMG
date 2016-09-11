#ifndef ___MATRIX_H___
#define ___MATRIX_H___

#include "Common.h"

namespace AlgebraicMultigrid
{
    /// Block Sparse Row matrix(row wise)
    class BSR
    {
    private:
        integer     m_nrow;
        integer     m_ncol;
        reaVec      m_val;
        intVec      m_ind;
        intVec      m_ptr;

    public:
        BSR()
            : m_nrow(0), m_ncol(0), m_val(), m_ind(), m_ptr()
        {

        }
    };
}


#endif