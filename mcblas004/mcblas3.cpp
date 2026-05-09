#include <mcblas3.hpp>
#include <iostream>
/*
 * Routines of cblas level 3
 */

void ssgemm(
    int m, int n, int k,
    float alpha,
    float ** A,
    float ** B,
    float beta,
    float ** C)
{
    int r, d, c;
    int rr, dd, cc;

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
    
    for (dd = 0; dd < k; dd += Td)
    {
        for (cc = 0; cc < n; cc += Tc)
        {
            for (rr = 0; rr < m; rr += Tr)
            {
                for (r = rr; r < std::min(rr + Tr, m); ++r)
                {
                    for (d = dd; d < std::min(dd + Td, k); ++d)
                    {
                        for (c = cc; c < std::min(cc + Tc, n); ++c)
                        {
                            C[r][c] = C[r][c] * beta + alpha * A[r][d] * B[d][c];
                        }
                    }
                }
            }
        }
    }
}





/*
#include <mcblas3.hpp>
#include <iostream>


void ssgemm(
    int m, int n, int k,
    float alpha,
    float ** A,
    float ** B,
    float beta,
    float ** C)
{
	// std::cout << "***************************************" << std::endl;
	// std::cout << "   WARNING: Not Yet Implemented..." << std::endl;
	// std::cout << "***************************************" << std::endl;

	int r, c, d;
	int rr, cc, dd;
	const int Tr = 16;
	const int Tc = 16;
	const int Td = 16;

		for (dd = 0;  dd < m; dd += Td)
		{
			for (cc = 0;  cc < m; cc += Tc)
			{
				for (rr = 0;  rr < m; rr += Tr)
				{
				for (c = cc;  c < std::min(cc+Tc, n); ++c)
				{
					for (r = rr;  r < std::min(rr+Tr ,m); ++r)
					{
						float t=0.0;
						for (d = dd;  d < std::min(dd+Td, m); ++d)
						{
							t = t + alpha * A[r][d] * B[d][c];
						}
						C[r][c] = C[r][c] * beta + t;
 					}
				}
			}
		}
	}

}

*/
