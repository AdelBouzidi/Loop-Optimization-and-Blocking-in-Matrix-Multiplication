#include <mcblas3.hpp>
#include <iostream>
/*
 * Routines of cblas level 3
 */

#define A(r,c) A[r*k+c]
#define B(r,c) B[r*n+c]
#define C(r,c) C[r*n+c]

void ssgemm(
    int m, int n, int k,
    float alpha,
    float *A,
    float *B,
    float beta,
    float *C)
{
    int r, c, d;
    int rr, cc, dd;

/*
    const int Tr = 16;
    const int Td = 16;
    const int Tc = 16;

    const int Tr = 32;
    const int Td = 32;
    const int Tc = 32;
    
    const int Tr = 64;
    const int Td = 64;
    const int Tc = 64;  
*/

    const int Tr = 128;
    const int Td = 128;
    const int Tc = 128; 
    
    for (rr = 0; rr < m; rr += Tr)
    {
        for (dd = 0; dd < k; dd += Td)
        {
            for (cc = 0; cc < n; cc += Tc)
            {
                for (r = rr; r < std::min(rr + Tr, m); ++r)
                {
                    for (d = dd; d < std::min(dd + Td, k); ++d)
                    {
                        for (c = cc; c < std::min(cc + Tc, n); ++c)
                        {
                            C(r,c) = beta * C(r,c) + alpha * A(r,d) * B(d,c);
                        }
                    }
                }
            }
        }
    }
}

