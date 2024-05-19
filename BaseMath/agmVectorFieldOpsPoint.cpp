#include "stdDefinitions.h"
#include <math.h>

#include "agmVectorFieldOps.h"
#include "agmScalarFieldOps.h"
#include "DiffCoefs.h"

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::cross_plane(CagmVectorFieldOps *a, const CagmVectorFieldOps *b, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    double tx, ty, tz;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            double ax = a->fieldX[fidx(kx, ky, kz)];
            double ay = a->fieldY[fidx(kx, ky, kz)];
            double az = a->fieldZ[fidx(kx, ky, kz)];
            double bx = b->fieldX[fidx(kx, ky, kz)];
            double by = b->fieldY[fidx(kx, ky, kz)];
            double bz = b->fieldZ[fidx(kx, ky, kz)];

            tx = a->fieldY[fidx(kx, ky, kz)] * b->fieldZ[fidx(kx, ky, kz)] - a->fieldZ[fidx(kx, ky, kz)] * b->fieldY[fidx(kx, ky, kz)];
            ty = a->fieldZ[fidx(kx, ky, kz)] * b->fieldX[fidx(kx, ky, kz)] - a->fieldX[fidx(kx, ky, kz)] * b->fieldZ[fidx(kx, ky, kz)];
            tz = a->fieldX[fidx(kx, ky, kz)] * b->fieldY[fidx(kx, ky, kz)] - a->fieldY[fidx(kx, ky, kz)] * b->fieldX[fidx(kx, ky, kz)];
            fieldX[fidx(kx, ky, kz)] = tx;
            fieldY[fidx(kx, ky, kz)] = ty;
            fieldZ[fidx(kx, ky, kz)] = tz;
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::rot_plane(CagmVectorFieldOps *a, int kz, int scheme)
{
    double zy, yz, xz, zx, yx, xy;
    for (int ky = 0; ky < N[1]; ky++)
        for (int kx = 0; kx < N[0]; kx++)
            if (scheme == 3)
            {
                if (kx == 0)
                {
                    zx = -3 * a->fieldZ[fidx(0, ky, kz)] + 4 * a->fieldZ[fidx(1, ky, kz)] - a->fieldZ[fidx(2, ky, kz)];
                    yx = -3 * a->fieldY[fidx(0, ky, kz)] + 4 * a->fieldY[fidx(1, ky, kz)] - a->fieldY[fidx(2, ky, kz)];
                }
                else if (kx == N[0] - 1)
                {
                    zx = a->fieldZ[fidx(N[0] - 3, ky, kz)] - 4 * a->fieldZ[fidx(N[0] - 2, ky, kz)] + 3 * a->fieldZ[fidx(N[0] - 1, ky, kz)];
                    yx = a->fieldY[fidx(N[0] - 3, ky, kz)] - 4 * a->fieldY[fidx(N[0] - 2, ky, kz)] + 3 * a->fieldY[fidx(N[0] - 1, ky, kz)];
                }
                else
                {
                    zx = a->fieldZ[fidx(kx + 1, ky, kz)] - a->fieldZ[fidx(kx - 1, ky, kz)];
                    yx = a->fieldY[fidx(kx + 1, ky, kz)] - a->fieldY[fidx(kx - 1, ky, kz)];
                }

                if (ky == 0)
                {
                    zy = -3 * a->fieldZ[fidx(kx, 0, kz)] + 4 * a->fieldZ[fidx(kx, 1, kz)] - a->fieldZ[fidx(kx, 2, kz)];
                    xy = -3 * a->fieldX[fidx(kx, 0, kz)] + 4 * a->fieldX[fidx(kx, 1, kz)] - a->fieldX[fidx(kx, 2, kz)];
                }
                else if (ky == N[1] - 1)
                {
                    zy = a->fieldZ[fidx(kx, N[1] - 3, kz)] - 4 * a->fieldZ[fidx(kx, N[1] - 2, kz)] + 3 * a->fieldZ[fidx(kx, N[1] - 1, kz)];
                    xy = a->fieldX[fidx(kx, N[1] - 3, kz)] - 4 * a->fieldX[fidx(kx, N[1] - 2, kz)] + 3 * a->fieldX[fidx(kx, N[1] - 1, kz)];
                }
                else
                {
                    zy = a->fieldZ[fidx(kx, ky + 1, kz)] - a->fieldZ[fidx(kx, ky - 1, kz)];
                    xy = a->fieldX[fidx(kx, ky + 1, kz)] - a->fieldX[fidx(kx, ky - 1, kz)];
                }

                if (kz == 0)
                {
                    xz = -3 * a->fieldX[fidx(kx, ky, 0)] + 4 * a->fieldX[fidx(kx, ky, 1)] - a->fieldX[fidx(kx, ky, 2)];
                    yz = -3 * a->fieldY[fidx(kx, ky, 0)] + 4 * a->fieldY[fidx(kx, ky, 1)] - a->fieldY[fidx(kx, ky, 2)];
                }
                else if (kz == N[2] - 1)
                {
                    xz = a->fieldX[fidx(kx, ky, N[2] - 3)] - 4 * a->fieldX[fidx(kx, ky, N[2] - 2)] + 3 * a->fieldX[fidx(kx, ky, N[2] - 1)];
                    yz = a->fieldY[fidx(kx, ky, N[2] - 3)] - 4 * a->fieldY[fidx(kx, ky, N[2] - 2)] + 3 * a->fieldY[fidx(kx, ky, N[2] - 1)];
                }
                else
                {
                    xz = a->fieldX[fidx(kx, ky, kz + 1)] - a->fieldX[fidx(kx, ky, kz - 1)];
                    yz = a->fieldY[fidx(kx, ky, kz + 1)] - a->fieldY[fidx(kx, ky, kz - 1)];
                }

                fieldX[fidx(kx, ky, kz)] = (zy*step[1] - yz*step[2])*0.5;
                fieldY[fidx(kx, ky, kz)] = (xz*step[2] - zx*step[0])*0.5;
                fieldZ[fidx(kx, ky, kz)] = (yx*step[0] - xy*step[1])*0.5;
            }
            else if (scheme == 5)
		    {
                if (kx == 0)
                {
                    zx = - 25*a->fieldZ[fidx(0, ky, kz)] + 48*a->fieldZ[fidx(1, ky, kz)] - 36*a->fieldZ[fidx(2, ky, kz)] + 16*a->fieldZ[fidx(3, ky, kz)] - 3*a->fieldZ[fidx(4, ky, kz)];
                    yx = - 25*a->fieldY[fidx(0, ky, kz)] + 48*a->fieldY[fidx(1, ky, kz)] - 36*a->fieldY[fidx(2, ky, kz)] + 16*a->fieldY[fidx(3, ky, kz)] - 3*a->fieldY[fidx(4, ky, kz)];
                }
                else if (kx == 1)
                {
                    zx = -  3*a->fieldZ[fidx(0, ky, kz)] - 10*a->fieldZ[fidx(1, ky, kz)] + 18*a->fieldZ[fidx(2, ky, kz)] -  6*a->fieldZ[fidx(3, ky, kz)] +   a->fieldZ[fidx(4, ky, kz)];
                    yx = -  3*a->fieldY[fidx(0, ky, kz)] - 10*a->fieldY[fidx(1, ky, kz)] + 18*a->fieldY[fidx(2, ky, kz)] -  6*a->fieldY[fidx(3, ky, kz)] +   a->fieldY[fidx(4, ky, kz)];
                }
                else if (kx == N[0]-2)
                {
                    zx = -  a->fieldZ[fidx(N[0]-5, ky, kz)] +  6*a->fieldZ[fidx(N[0]-4, ky, kz)] - 18*a->fieldZ[fidx(N[0]-3, ky, kz)] + 10*a->fieldZ[fidx(N[0]-2, ky, kz)] +  3*a->fieldZ[fidx(N[0]-1, ky, kz)];
                    yx = -  a->fieldY[fidx(N[0]-5, ky, kz)] +  6*a->fieldY[fidx(N[0]-4, ky, kz)] - 18*a->fieldY[fidx(N[0]-3, ky, kz)] + 10*a->fieldY[fidx(N[0]-2, ky, kz)] +  3*a->fieldY[fidx(N[0]-1, ky, kz)];
                }
                else if (kx == N[0]-1)
                {
                    zx =  3*a->fieldZ[fidx(N[0]-5, ky, kz)] - 16*a->fieldZ[fidx(N[0]-4, ky, kz)] + 36*a->fieldZ[fidx(N[0]-3, ky, kz)] - 48*a->fieldZ[fidx(N[0]-2, ky, kz)] + 25*a->fieldZ[fidx(N[0]-1, ky, kz)];
                    yx =  3*a->fieldY[fidx(N[0]-5, ky, kz)] - 16*a->fieldY[fidx(N[0]-4, ky, kz)] + 36*a->fieldY[fidx(N[0]-3, ky, kz)] - 48*a->fieldY[fidx(N[0]-2, ky, kz)] + 25*a->fieldY[fidx(N[0]-1, ky, kz)];
                }
                else
                {
                    zx = - a->fieldZ[fidx(kx+2, ky, kz)] + 8*a->fieldZ[fidx(kx+1, ky, kz)] - 8*a->fieldZ[fidx(kx-1, ky, kz)] + a->fieldZ[fidx(kx-2, ky, kz)];
                    yx = - a->fieldY[fidx(kx+2, ky, kz)] + 8*a->fieldY[fidx(kx+1, ky, kz)] - 8*a->fieldY[fidx(kx-1, ky, kz)] + a->fieldY[fidx(kx-2, ky, kz)];
                }

                if (ky == 0)
                {
                    zy = - 25*a->fieldZ[fidx(kx, 0, kz)] + 48*a->fieldZ[fidx(kx, 1, kz)] - 36*a->fieldZ[fidx(kx, 2, kz)] + 16*a->fieldZ[fidx(kx, 3, kz)] - 3*a->fieldZ[fidx(kx, 4, kz)];
                    xy = - 25*a->fieldX[fidx(kx, 0, kz)] + 48*a->fieldX[fidx(kx, 1, kz)] - 36*a->fieldX[fidx(kx, 2, kz)] + 16*a->fieldX[fidx(kx, 3, kz)] - 3*a->fieldX[fidx(kx, 4, kz)];
                }
                else if (ky == 1)
                {
                    zy = -  3*a->fieldZ[fidx(kx, 0, kz)] - 10*a->fieldZ[fidx(kx, 1, kz)] + 18*a->fieldZ[fidx(kx, 2, kz)] -  6*a->fieldZ[fidx(kx, 3, kz)] +   a->fieldZ[fidx(kx, 4, kz)];
                    xy = -  3*a->fieldX[fidx(kx, 0, kz)] - 10*a->fieldX[fidx(kx, 1, kz)] + 18*a->fieldX[fidx(kx, 2, kz)] -  6*a->fieldX[fidx(kx, 3, kz)] +   a->fieldX[fidx(kx, 4, kz)];
                }
                else if (ky == N[1]-2)
                {
                    zy = -  a->fieldZ[fidx(kx, N[1]-5, kz)] +  6*a->fieldZ[fidx(kx, N[1]-4, kz)] - 18*a->fieldZ[fidx(kx, N[1]-3, kz)] + 10*a->fieldZ[fidx(kx, N[1]-2, kz)] +  3*a->fieldZ[fidx(kx, N[1]-1, kz)];
                    xy = -  a->fieldX[fidx(kx, N[1]-5, kz)] +  6*a->fieldX[fidx(kx, N[1]-4, kz)] - 18*a->fieldX[fidx(kx, N[1]-3, kz)] + 10*a->fieldX[fidx(kx, N[1]-2, kz)] +  3*a->fieldX[fidx(kx, N[1]-1, kz)];
                }
                else if (ky == N[1]-1)
                {
                    zy =  3*a->fieldZ[fidx(kx, N[1]-5, kz)] - 16*a->fieldZ[fidx(kx, N[1]-4, kz)] + 36*a->fieldZ[fidx(kx, N[1]-3, kz)] - 48*a->fieldZ[fidx(kx, N[1]-2, kz)] + 25*a->fieldZ[fidx(kx, N[1]-1, kz)];
                    xy =  3*a->fieldX[fidx(kx, N[1]-5, kz)] - 16*a->fieldX[fidx(kx, N[1]-4, kz)] + 36*a->fieldX[fidx(kx, N[1]-3, kz)] - 48*a->fieldX[fidx(kx, N[1]-2, kz)] + 25*a->fieldX[fidx(kx, N[1]-1, kz)];
                }
                else
                {
                    zy = - a->fieldZ[fidx(kx, ky+2, kz)] + 8*a->fieldZ[fidx(kx, ky+1, kz)] - 8*a->fieldZ[fidx(kx, ky-1, kz)] + a->fieldZ[fidx(kx, ky-2, kz)];
                    xy = - a->fieldX[fidx(kx, ky+2, kz)] + 8*a->fieldX[fidx(kx, ky+1, kz)] - 8*a->fieldX[fidx(kx, ky-1, kz)] + a->fieldX[fidx(kx, ky-2, kz)];
                }

                if (kz == 0)
                {
                    xz = - 25*a->fieldX[fidx(kx, ky, 0)] + 48*a->fieldX[fidx(kx, ky, 1)] - 36*a->fieldX[fidx(kx, ky, 2)] + 16*a->fieldX[fidx(kx, ky, 3)] - 3*a->fieldX[fidx(kx, ky, 4)];
                    yz = - 25*a->fieldY[fidx(kx, ky, 0)] + 48*a->fieldY[fidx(kx, ky, 1)] - 36*a->fieldY[fidx(kx, ky, 2)] + 16*a->fieldY[fidx(kx, ky, 3)] - 3*a->fieldY[fidx(kx, ky, 4)];
                }
                else if (kz == 1)
                {
                    xz = -  3*a->fieldX[fidx(kx, ky, 0)] - 10*a->fieldX[fidx(kx, ky, 1)] + 18*a->fieldX[fidx(kx, ky, 2)] -  6*a->fieldX[fidx(kx, ky, 3)] +   a->fieldX[fidx(kx, ky, 4)];
                    yz = -  3*a->fieldY[fidx(kx, ky, 0)] - 10*a->fieldY[fidx(kx, ky, 1)] + 18*a->fieldY[fidx(kx, ky, 2)] -  6*a->fieldY[fidx(kx, ky, 3)] +   a->fieldY[fidx(kx, ky, 4)];
                }
                else if (kz == N[2]-2)
                {
                    xz = -  a->fieldX[fidx(kx, ky, N[2]-5)] +  6*a->fieldX[fidx(kx, ky, N[2]-4)] - 18*a->fieldX[fidx(kx, ky, N[2]-3)] + 10*a->fieldX[fidx(kx, ky, N[2]-2)] +  3*a->fieldX[fidx(kx, ky, N[2]-1)];
                    yz = -  a->fieldY[fidx(kx, ky, N[2]-5)] +  6*a->fieldY[fidx(kx, ky, N[2]-4)] - 18*a->fieldY[fidx(kx, ky, N[2]-3)] + 10*a->fieldY[fidx(kx, ky, N[2]-2)] +  3*a->fieldY[fidx(kx, ky, N[2]-1)];
                }
                else if (kz == N[2]-1)
                {
                    xz =  3*a->fieldX[fidx(kx, ky, N[2]-5)] - 16*a->fieldX[fidx(kx, ky, N[2]-4)] + 36*a->fieldX[fidx(kx, ky, N[2]-3)] - 48*a->fieldX[fidx(kx, ky, N[2]-2)] + 25*a->fieldX[fidx(kx, ky, N[2]-1)];
                    yz =  3*a->fieldY[fidx(kx, ky, N[2]-5)] - 16*a->fieldY[fidx(kx, ky, N[2]-4)] + 36*a->fieldY[fidx(kx, ky, N[2]-3)] - 48*a->fieldY[fidx(kx, ky, N[2]-2)] + 25*a->fieldY[fidx(kx, ky, N[2]-1)];
                }
                else
                {
                    xz = - a->fieldX[fidx(kx, ky, kz+2)] + 8*a->fieldX[fidx(kx, ky, kz+1)] - 8*a->fieldX[fidx(kx, ky, kz-1)] + a->fieldX[fidx(kx, ky, kz-2)];
                    yz = - a->fieldY[fidx(kx, ky, kz+2)] + 8*a->fieldY[fidx(kx, ky, kz+1)] - 8*a->fieldY[fidx(kx, ky, kz-1)] + a->fieldY[fidx(kx, ky, kz-2)];
                }

			    fieldX[fidx(kx, ky, kz)] = (zy*step[1] - yz*step[2]) / 12.0;
			    fieldY[fidx(kx, ky, kz)] = (xz*step[2] - zx*step[0]) / 12.0;
			    fieldZ[fidx(kx, ky, kz)] = (yx*step[0] - xy*step[1]) / 12.0;
		    }
            else // scheme == 7
            {
                if (kx == 0)
                {
                    zx =    s7m3_0 * a->fieldZ[fidx(0, ky, kz)] + s7m3_1 * a->fieldZ[fidx(1, ky, kz)] + s7m3_2 * a->fieldZ[fidx(2, ky, kz)] + s7m3_3 * a->fieldZ[fidx(3, ky, kz)] 
                          + s7m3_4 * a->fieldZ[fidx(4, ky, kz)] + s7m3_5 * a->fieldZ[fidx(5, ky, kz)] + s7m3_6 * a->fieldZ[fidx(6, ky, kz)];
                    yx =    s7m3_0 * a->fieldY[fidx(0, ky, kz)] + s7m3_1 * a->fieldY[fidx(1, ky, kz)] + s7m3_2 * a->fieldY[fidx(2, ky, kz)] + s7m3_3 * a->fieldY[fidx(3, ky, kz)]
                          + s7m3_4 * a->fieldY[fidx(4, ky, kz)] + s7m3_5 * a->fieldY[fidx(5, ky, kz)] + s7m3_6 * a->fieldY[fidx(6, ky, kz)];
                }
                else if (kx == 1)
                {
                    zx =    s7m2_0 * a->fieldZ[fidx(0, ky, kz)] + s7m2_1 * a->fieldZ[fidx(1, ky, kz)] + s7m2_2 * a->fieldZ[fidx(2, ky, kz)] + s7m2_3 * a->fieldZ[fidx(3, ky, kz)]
                          + s7m2_4 * a->fieldZ[fidx(4, ky, kz)] + s7m2_5 * a->fieldZ[fidx(5, ky, kz)] + s7m2_6 * a->fieldZ[fidx(6, ky, kz)];
                    yx =    s7m2_0 * a->fieldY[fidx(0, ky, kz)] + s7m2_1 * a->fieldY[fidx(1, ky, kz)] + s7m2_2 * a->fieldY[fidx(2, ky, kz)] + s7m2_3 * a->fieldY[fidx(3, ky, kz)]
                          + s7m2_4 * a->fieldY[fidx(4, ky, kz)] + s7m2_5 * a->fieldY[fidx(5, ky, kz)] + s7m2_6 * a->fieldY[fidx(6, ky, kz)];
                }
                else if (kx == 2)
                {
                    zx =    s7m1_0 * a->fieldZ[fidx(0, ky, kz)] + s7m1_1 * a->fieldZ[fidx(1, ky, kz)] + s7m1_2 * a->fieldZ[fidx(2, ky, kz)] + s7m1_3 * a->fieldZ[fidx(3, ky, kz)]
                          + s7m1_4 * a->fieldZ[fidx(4, ky, kz)] + s7m1_5 * a->fieldZ[fidx(5, ky, kz)] + s7m1_6 * a->fieldZ[fidx(6, ky, kz)];
                    yx =    s7m1_0 * a->fieldY[fidx(0, ky, kz)] + s7m1_1 * a->fieldY[fidx(1, ky, kz)] + s7m1_2 * a->fieldY[fidx(2, ky, kz)] + s7m1_3 * a->fieldY[fidx(3, ky, kz)]
                          + s7m1_4 * a->fieldY[fidx(4, ky, kz)] + s7m1_5 * a->fieldY[fidx(5, ky, kz)] + s7m1_6 * a->fieldY[fidx(6, ky, kz)];
                }
                else if (kx == N[0]-3)
                {
                    zx = - (s7m1_6 * a->fieldZ[fidx(N[0]-7, ky, kz)] + s7m1_5 * a->fieldZ[fidx(N[0]-6, ky, kz)] + s7m1_4 * a->fieldZ[fidx(N[0]-5, ky, kz)] + s7m1_3 * a->fieldZ[fidx(N[0]-4, ky, kz)]
                          + s7m1_2 * a->fieldZ[fidx(N[0]-3, ky, kz)] + s7m1_1 * a->fieldZ[fidx(N[0]-2, ky, kz)] + s7m1_0 * a->fieldZ[fidx(N[0]-1, ky, kz)]);
                    yx = - (s7m1_6 * a->fieldY[fidx(N[0]-7, ky, kz)] + s7m1_5 * a->fieldY[fidx(N[0]-6, ky, kz)] + s7m1_4 * a->fieldY[fidx(N[0]-5, ky, kz)] + s7m1_3 * a->fieldY[fidx(N[0]-4, ky, kz)]
                          + s7m1_2 * a->fieldY[fidx(N[0]-3, ky, kz)] + s7m1_1 * a->fieldY[fidx(N[0]-2, ky, kz)] + s7m1_0 * a->fieldY[fidx(N[0]-1, ky, kz)]);
                }
                else if (kx == N[0]-2)
                {
                    zx = - (s7m2_6 * a->fieldZ[fidx(N[0]-7, ky, kz)] + s7m2_5 * a->fieldZ[fidx(N[0]-6, ky, kz)] + s7m2_4 * a->fieldZ[fidx(N[0]-5, ky, kz)] + s7m2_3 * a->fieldZ[fidx(N[0]-4, ky, kz)]
                          + s7m2_2 * a->fieldZ[fidx(N[0]-3, ky, kz)] + s7m2_1 * a->fieldZ[fidx(N[0]-2, ky, kz)] + s7m2_0 * a->fieldZ[fidx(N[0]-1, ky, kz)]);
                    yx = - (s7m2_6 * a->fieldY[fidx(N[0]-7, ky, kz)] + s7m2_5 * a->fieldY[fidx(N[0]-6, ky, kz)] + s7m2_4 * a->fieldY[fidx(N[0]-5, ky, kz)] + s7m2_3 * a->fieldY[fidx(N[0]-4, ky, kz)]
                          + s7m2_2 * a->fieldY[fidx(N[0]-3, ky, kz)] + s7m2_1 * a->fieldY[fidx(N[0]-2, ky, kz)] + s7m2_0 * a->fieldY[fidx(N[0]-1, ky, kz)]);
                }
                else if (kx == N[0]-1)
                {
                    zx = - (s7m3_6 * a->fieldZ[fidx(N[0]-7, ky, kz)] + s7m3_5 * a->fieldZ[fidx(N[0]-6, ky, kz)] + s7m3_4 * a->fieldZ[fidx(N[0]-5, ky, kz)] + s7m3_3 * a->fieldZ[fidx(N[0]-4, ky, kz)]
                          + s7m3_2 * a->fieldZ[fidx(N[0]-3, ky, kz)] + s7m3_1 * a->fieldZ[fidx(N[0]-2, ky, kz)] + s7m3_0 * a->fieldZ[fidx(N[0]-1, ky, kz)]);
                    yx = - (s7m3_6 * a->fieldY[fidx(N[0]-7, ky, kz)] + s7m3_5 * a->fieldY[fidx(N[0]-6, ky, kz)] + s7m3_4 * a->fieldY[fidx(N[0]-5, ky, kz)] + s7m3_3 * a->fieldY[fidx(N[0]-4, ky, kz)]
                          + s7m3_2 * a->fieldY[fidx(N[0]-3, ky, kz)] + s7m3_1 * a->fieldY[fidx(N[0]-2, ky, kz)] + s7m3_0 * a->fieldY[fidx(N[0]-1, ky, kz)]);
                }
                else
                {
                    zx =  s7m0_6 * a->fieldZ[fidx(kx+3, ky, kz)] + s7m0_5 * a->fieldZ[fidx(kx+2, ky, kz)] + s7m0_4 * a->fieldZ[fidx(kx+1, ky, kz)]
                        + s7m0_2 * a->fieldZ[fidx(kx-1, ky, kz)] + s7m0_1 * a->fieldZ[fidx(kx-2, ky, kz)] + s7m0_0 * a->fieldZ[fidx(kx-3, ky, kz)];
                    yx =  s7m0_6 * a->fieldY[fidx(kx+3, ky, kz)] + s7m0_5 * a->fieldY[fidx(kx+2, ky, kz)] + s7m0_4 * a->fieldY[fidx(kx+1, ky, kz)]
                        + s7m0_2 * a->fieldY[fidx(kx-1, ky, kz)] + s7m0_1 * a->fieldY[fidx(kx-2, ky, kz)] + s7m0_0 * a->fieldY[fidx(kx-3, ky, kz)];
                }

                // y
                if (ky == 0)
                {
                    zy =    s7m3_0 * a->fieldZ[fidx(kx, 0, kz)] + s7m3_1 * a->fieldZ[fidx(kx, 1, kz)] + s7m3_2 * a->fieldZ[fidx(kx, 2, kz)] + s7m3_3 * a->fieldZ[fidx(kx, 3, kz)] 
                          + s7m3_4 * a->fieldZ[fidx(kx, 4, kz)] + s7m3_5 * a->fieldZ[fidx(kx, 5, kz)] + s7m3_6 * a->fieldZ[fidx(kx, 6, kz)];
                    xy =    s7m3_0 * a->fieldX[fidx(kx, 0, kz)] + s7m3_1 * a->fieldX[fidx(kx, 1, kz)] + s7m3_2 * a->fieldX[fidx(kx, 2, kz)] + s7m3_3 * a->fieldX[fidx(kx, 3, kz)]
                          + s7m3_4 * a->fieldX[fidx(kx, 4, kz)] + s7m3_5 * a->fieldX[fidx(kx, 5, kz)] + s7m3_6 * a->fieldX[fidx(kx, 6, kz)];
                }
                else if (ky == 1)
                {
                    zy =    s7m2_0 * a->fieldZ[fidx(kx, 0, kz)] + s7m2_1 * a->fieldZ[fidx(kx, 1, kz)] + s7m2_2 * a->fieldZ[fidx(kx, 2, kz)] + s7m2_3 * a->fieldZ[fidx(kx, 3, kz)]
                          + s7m2_4 * a->fieldZ[fidx(kx, 4, kz)] + s7m2_5 * a->fieldZ[fidx(kx, 5, kz)] + s7m2_6 * a->fieldZ[fidx(kx, 6, kz)];
                    xy =    s7m2_0 * a->fieldX[fidx(kx, 0, kz)] + s7m2_1 * a->fieldX[fidx(kx, 1, kz)] + s7m2_2 * a->fieldX[fidx(kx, 2, kz)] + s7m2_3 * a->fieldX[fidx(kx, 3, kz)]
                          + s7m2_4 * a->fieldX[fidx(kx, 4, kz)] + s7m2_5 * a->fieldX[fidx(kx, 5, kz)] + s7m2_6 * a->fieldX[fidx(kx, 6, kz)];
                }
                else if (ky == 2)
                {
                    zy =    s7m1_0 * a->fieldZ[fidx(kx, 0, kz)] + s7m1_1 * a->fieldZ[fidx(kx, 1, kz)] + s7m1_2 * a->fieldZ[fidx(kx, 2, kz)] + s7m1_3 * a->fieldZ[fidx(kx, 3, kz)]
                          + s7m1_4 * a->fieldZ[fidx(kx, 4, kz)] + s7m1_5 * a->fieldZ[fidx(kx, 5, kz)] + s7m1_6 * a->fieldZ[fidx(kx, 6, kz)];
                    xy =    s7m1_0 * a->fieldX[fidx(kx, 0, kz)] + s7m1_1 * a->fieldX[fidx(kx, 1, kz)] + s7m1_2 * a->fieldX[fidx(kx, 2, kz)] + s7m1_3 * a->fieldX[fidx(kx, 3, kz)]
                          + s7m1_4 * a->fieldX[fidx(kx, 4, kz)] + s7m1_5 * a->fieldX[fidx(kx, 5, kz)] + s7m1_6 * a->fieldX[fidx(kx, 6, kz)];
                }
                else if (ky == N[1]-3)
                {
                    zy = - (s7m1_6 * a->fieldZ[fidx(kx, N[1]-7, kz)] + s7m1_5 * a->fieldZ[fidx(kx, N[1]-6, kz)] + s7m1_4 * a->fieldZ[fidx(kx, N[1]-5, kz)] + s7m1_3 * a->fieldZ[fidx(kx, N[1]-4, kz)]
                          + s7m1_2 * a->fieldZ[fidx(kx, N[1]-3, kz)] + s7m1_1 * a->fieldZ[fidx(kx, N[1]-2, kz)] + s7m1_0 * a->fieldZ[fidx(kx, N[1]-1, kz)]);
                    xy = - (s7m1_6 * a->fieldX[fidx(kx, N[1]-7, kz)] + s7m1_5 * a->fieldX[fidx(kx, N[1]-6, kz)] + s7m1_4 * a->fieldX[fidx(kx, N[1]-5, kz)] + s7m1_3 * a->fieldX[fidx(kx, N[1]-4, kz)]
                          + s7m1_2 * a->fieldX[fidx(kx, N[1]-3, kz)] + s7m1_1 * a->fieldX[fidx(kx, N[1]-2, kz)] + s7m1_0 * a->fieldX[fidx(kx, N[1]-1, kz)]);
                }
                else if (ky == N[1]-2)
                {
                    zy = - (s7m2_6 * a->fieldZ[fidx(kx, N[1]-7, kz)] + s7m2_5 * a->fieldZ[fidx(kx, N[1]-6, kz)] + s7m2_4 * a->fieldZ[fidx(kx, N[1]-5, kz)] + s7m2_3 * a->fieldZ[fidx(kx, N[1]-4, kz)]
                          + s7m2_2 * a->fieldZ[fidx(kx, N[1]-3, kz)] + s7m2_1 * a->fieldZ[fidx(kx, N[1]-2, kz)] + s7m2_0 * a->fieldZ[fidx(kx, N[1]-1, kz)]);
                    xy = - (s7m2_6 * a->fieldX[fidx(kx, N[1]-7, kz)] + s7m2_5 * a->fieldX[fidx(kx, N[1]-6, kz)] + s7m2_4 * a->fieldX[fidx(kx, N[1]-5, kz)] + s7m2_3 * a->fieldX[fidx(kx, N[1]-4, kz)]
                          + s7m2_2 * a->fieldX[fidx(kx, N[1]-3, kz)] + s7m2_1 * a->fieldX[fidx(kx, N[1]-2, kz)] + s7m2_0 * a->fieldX[fidx(kx, N[1]-1, kz)]);
                }
                else if (ky == N[1]-1)
                {
                    zy = - (s7m3_6 * a->fieldZ[fidx(kx, N[1]-7, kz)] + s7m3_5 * a->fieldZ[fidx(kx, N[1]-6, kz)] + s7m3_4 * a->fieldZ[fidx(kx, N[1]-5, kz)] + s7m3_3 * a->fieldZ[fidx(kx, N[1]-4, kz)]
                          + s7m3_2 * a->fieldZ[fidx(kx, N[1]-3, kz)] + s7m3_1 * a->fieldZ[fidx(kx, N[1]-2, kz)] + s7m3_0 * a->fieldZ[fidx(kx, N[1]-1, kz)]);
                    xy = - (s7m3_6 * a->fieldX[fidx(kx, N[1]-7, kz)] + s7m3_5 * a->fieldX[fidx(kx, N[1]-6, kz)] + s7m3_4 * a->fieldX[fidx(kx, N[1]-5, kz)] + s7m3_3 * a->fieldX[fidx(kx, N[1]-4, kz)]
                          + s7m3_2 * a->fieldX[fidx(kx, N[1]-3, kz)] + s7m3_1 * a->fieldX[fidx(kx, N[1]-2, kz)] + s7m3_0 * a->fieldX[fidx(kx, N[1]-1, kz)]);
                }
                else
                {
                    zy =  s7m0_6 * a->fieldZ[fidx(kx, ky+3, kz)] + s7m0_5 * a->fieldZ[fidx(kx, ky+2, kz)] + s7m0_4 * a->fieldZ[fidx(kx, ky+1, kz)]
                        + s7m0_2 * a->fieldZ[fidx(kx, ky-1, kz)] + s7m0_1 * a->fieldZ[fidx(kx, ky-2, kz)] + s7m0_0 * a->fieldZ[fidx(kx, ky-3, kz)];
                    xy =  s7m0_6 * a->fieldX[fidx(kx, ky+3, kz)] + s7m0_5 * a->fieldX[fidx(kx, ky+2, kz)] + s7m0_4 * a->fieldX[fidx(kx, ky+1, kz)]
                        + s7m0_2 * a->fieldX[fidx(kx, ky-1, kz)] + s7m0_1 * a->fieldX[fidx(kx, ky-2, kz)] + s7m0_0 * a->fieldX[fidx(kx, ky-3, kz)];
                }

                // z
                if (kz == 0)
                {
                    yz =    s7m3_0 * a->fieldY[fidx(kx, ky, 0)] + s7m3_1 * a->fieldY[fidx(kx, ky, 1)] + s7m3_2 * a->fieldY[fidx(kx, ky, 2)] + s7m3_3 * a->fieldY[fidx(kx, ky, 3)] 
                          + s7m3_4 * a->fieldY[fidx(kx, ky, 4)] + s7m3_5 * a->fieldY[fidx(kx, ky, 5)] + s7m3_6 * a->fieldY[fidx(kx, ky, 6)];
                    xz =    s7m3_0 * a->fieldX[fidx(kx, ky, 0)] + s7m3_1 * a->fieldX[fidx(kx, ky, 1)] + s7m3_2 * a->fieldX[fidx(kx, ky, 2)] + s7m3_3 * a->fieldX[fidx(kx, ky, 3)]
                          + s7m3_4 * a->fieldX[fidx(kx, ky, 4)] + s7m3_5 * a->fieldX[fidx(kx, ky, 5)] + s7m3_6 * a->fieldX[fidx(kx, ky, 6)];
                }
                else if (kz == 1)
                {
                    yz =    s7m2_0 * a->fieldY[fidx(kx, ky, 0)] + s7m2_1 * a->fieldY[fidx(kx, ky, 1)] + s7m2_2 * a->fieldY[fidx(kx, ky, 2)] + s7m2_3 * a->fieldY[fidx(kx, ky, 3)]
                          + s7m2_4 * a->fieldY[fidx(kx, ky, 4)] + s7m2_5 * a->fieldY[fidx(kx, ky, 5)] + s7m2_6 * a->fieldY[fidx(kx, ky, 6)];
                    xz =    s7m2_0 * a->fieldX[fidx(kx, ky, 0)] + s7m2_1 * a->fieldX[fidx(kx, ky, 1)] + s7m2_2 * a->fieldX[fidx(kx, ky, 2)] + s7m2_3 * a->fieldX[fidx(kx, ky, 3)]
                          + s7m2_4 * a->fieldX[fidx(kx, ky, 4)] + s7m2_5 * a->fieldX[fidx(kx, ky, 5)] + s7m2_6 * a->fieldX[fidx(kx, ky, 6)];
                }
                else if (kz == 2)
                {
                    yz =    s7m1_0 * a->fieldY[fidx(kx, ky, 0)] + s7m1_1 * a->fieldY[fidx(kx, ky, 1)] + s7m1_2 * a->fieldY[fidx(kx, ky, 2)] + s7m1_3 * a->fieldY[fidx(kx, ky, 3)]
                          + s7m1_4 * a->fieldY[fidx(kx, ky, 4)] + s7m1_5 * a->fieldY[fidx(kx, ky, 5)] + s7m1_6 * a->fieldY[fidx(kx, ky, 6)];
                    xz =    s7m1_0 * a->fieldX[fidx(kx, ky, 0)] + s7m1_1 * a->fieldX[fidx(kx, ky, 1)] + s7m1_2 * a->fieldX[fidx(kx, ky, 2)] + s7m1_3 * a->fieldX[fidx(kx, ky, 3)]
                          + s7m1_4 * a->fieldX[fidx(kx, ky, 4)] + s7m1_5 * a->fieldX[fidx(kx, ky, 5)] + s7m1_6 * a->fieldX[fidx(kx, ky, 6)];
                }
                else if (kz == N[2]-3)
                {
                    yz = - (s7m1_6 * a->fieldY[fidx(kx, ky, N[2]-7)] + s7m1_5 * a->fieldY[fidx(kx, ky, N[2]-6)] + s7m1_4 * a->fieldY[fidx(kx, ky, N[2]-5)] + s7m1_3 * a->fieldY[fidx(kx, ky, N[2]-4)]
                          + s7m1_2 * a->fieldY[fidx(kx, ky, N[2]-3)] + s7m1_1 * a->fieldY[fidx(kx, ky, N[2]-2)] + s7m1_0 * a->fieldY[fidx(kx, ky, N[2]-1)]);
                    xz = - (s7m1_6 * a->fieldX[fidx(kx, ky, N[2]-7)] + s7m1_5 * a->fieldX[fidx(kx, ky, N[2]-6)] + s7m1_4 * a->fieldX[fidx(kx, ky, N[2]-5)] + s7m1_3 * a->fieldX[fidx(kx, ky, N[2]-4)]
                          + s7m1_2 * a->fieldX[fidx(kx, ky, N[2]-3)] + s7m1_1 * a->fieldX[fidx(kx, ky, N[2]-2)] + s7m1_0 * a->fieldX[fidx(kx, ky, N[2]-1)]);
                }
                else if (kz == N[2]-2)
                {
                    yz = - (s7m2_6 * a->fieldY[fidx(kx, ky, N[2]-7)] + s7m2_5 * a->fieldY[fidx(kx, ky, N[2]-6)] + s7m2_4 * a->fieldY[fidx(kx, ky, N[2]-5)] + s7m2_3 * a->fieldY[fidx(kx, ky, N[2]-4)]
                          + s7m2_2 * a->fieldY[fidx(kx, ky, N[2]-3)] + s7m2_1 * a->fieldY[fidx(kx, ky, N[2]-2)] + s7m2_0 * a->fieldY[fidx(kx, ky, N[2]-1)]);
                    xz = - (s7m2_6 * a->fieldX[fidx(kx, ky, N[2]-7)] + s7m2_5 * a->fieldX[fidx(kx, ky, N[2]-6)] + s7m2_4 * a->fieldX[fidx(kx, ky, N[2]-5)] + s7m2_3 * a->fieldX[fidx(kx, ky, N[2]-4)]
                          + s7m2_2 * a->fieldX[fidx(kx, ky, N[2]-3)] + s7m2_1 * a->fieldX[fidx(kx, ky, N[2]-2)] + s7m2_0 * a->fieldX[fidx(kx, ky, N[2]-1)]);
                }
                else if (kz == N[2]-1)
                {
                    yz = - (s7m3_6 * a->fieldY[fidx(kx, ky, N[2]-7)] + s7m3_5 * a->fieldY[fidx(kx, ky, N[2]-6)] + s7m3_4 * a->fieldY[fidx(kx, ky, N[2]-5)] + s7m3_3 * a->fieldY[fidx(kx, ky, N[2]-4)]
                          + s7m3_2 * a->fieldY[fidx(kx, ky, N[2]-3)] + s7m3_1 * a->fieldY[fidx(kx, ky, N[2]-2)] + s7m3_0 * a->fieldY[fidx(kx, ky, N[2]-1)]);
                    xz = - (s7m3_6 * a->fieldX[fidx(kx, ky, N[2]-7)] + s7m3_5 * a->fieldX[fidx(kx, ky, N[2]-6)] + s7m3_4 * a->fieldX[fidx(kx, ky, N[2]-5)] + s7m3_3 * a->fieldX[fidx(kx, ky, N[2]-4)]
                          + s7m3_2 * a->fieldX[fidx(kx, ky, N[2]-3)] + s7m3_1 * a->fieldX[fidx(kx, ky, N[2]-2)] + s7m3_0 * a->fieldX[fidx(kx, ky, N[2]-1)]);
                }
                else
                {
                    yz =  s7m0_6 * a->fieldY[fidx(kx, ky, kz+3)] + s7m0_5 * a->fieldY[fidx(kx, ky, kz+2)] + s7m0_4 * a->fieldY[fidx(kx, ky, kz+1)]
                        + s7m0_2 * a->fieldY[fidx(kx, ky, kz-1)] + s7m0_1 * a->fieldY[fidx(kx, ky, kz-2)] + s7m0_0 * a->fieldY[fidx(kx, ky, kz-3)];
                    xz =  s7m0_6 * a->fieldX[fidx(kx, ky, kz+3)] + s7m0_5 * a->fieldX[fidx(kx, ky, kz+2)] + s7m0_4 * a->fieldX[fidx(kx, ky, kz+1)]
                        + s7m0_2 * a->fieldX[fidx(kx, ky, kz-1)] + s7m0_1 * a->fieldX[fidx(kx, ky, kz-2)] + s7m0_0 * a->fieldX[fidx(kx, ky, kz-3)];
                }

                fieldX[fidx(kx, ky, kz)] = (zy*step[1] - yz*step[2]) / 60.0;
			    fieldY[fidx(kx, ky, kz)] = (xz*step[2] - zx*step[0]) / 60.0;
			    fieldZ[fidx(kx, ky, kz)] = (yx*step[0] - xy*step[1]) / 60.0;
            }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::grad_plane(CagmScalarFieldOps *a, int kz, int scheme)
{
    // check equiv. sizes!
    int kx, ky;
    double dx, dy, dz;
    for (int ky = 0; ky < N[1]; ky++)
        for (int kx = 0; kx < N[0]; kx++)
            if (scheme == 3)
            {
                if (kx == 0)
                    dx = -3 * a->field[fidx(0, ky, kz)] + 4 * a->field[fidx(1, ky, kz)] - a->field[fidx(2, ky, kz)];
                else if (kx == N[0] - 1)
                    dx = a->field[fidx(N[0] - 3, ky, kz)] - 4 * a->field[fidx(N[0] - 2, ky, kz)] + 3 * a->field[fidx(N[0] - 1, ky, kz)];
                else
                    dx = a->field[fidx(kx + 1, ky, kz)] - a->field[fidx(kx - 1, ky, kz)];

                if (ky == 0)
                    dy = -3 * a->field[fidx(kx, 0, kz)] + 4 * a->field[fidx(kx, 1, kz)] - a->field[fidx(kx, 2, kz)];
                else if (ky == N[1] - 1)
                    dy = a->field[fidx(kx, N[1] - 3, kz)] - 4 * a->field[fidx(kx, N[1] - 2, kz)] + 3 * a->field[fidx(kx, N[1] - 1, kz)];
                else
                    dy = a->field[fidx(kx, ky + 1, kz)] - a->field[fidx(kx, ky - 1, kz)];

                if (kz == 0)
                    dz = -3 * a->field[fidx(kx, ky, 0)] + 4 * a->field[fidx(kx, ky, 1)] - a->field[fidx(kx, ky, 2)];
                else if (kz == N[2] - 1)
                    dz = a->field[fidx(kx, ky, N[2] - 3)] - 4 * a->field[fidx(kx, ky, N[2] - 2)] + 3 * a->field[fidx(kx, ky, N[2] - 1)];
                else
                    dz = a->field[fidx(kx, ky, kz + 1)] - a->field[fidx(kx, ky, kz - 1)];

                fieldX[fidx(kx, ky, kz)] = dx*0.5*step[0];
                fieldY[fidx(kx, ky, kz)] = dy*0.5*step[1];
                fieldZ[fidx(kx, ky, kz)] = dz*0.5*step[2];
            }
            else if (scheme == 5)
		    {
                if (kx == 0)
                    dx = - 25*a->field[fidx(0, ky, kz)] + 48*a->field[fidx(1, ky, kz)] - 36*a->field[fidx(2, ky, kz)] + 16*a->field[fidx(3, ky, kz)] - 3*a->field[fidx(4, ky, kz)];
                else if (kx == 1)
                    dx = -  3*a->field[fidx(0, ky, kz)] - 10*a->field[fidx(1, ky, kz)] + 18*a->field[fidx(2, ky, kz)] -  6*a->field[fidx(3, ky, kz)] +   a->field[fidx(4, ky, kz)];
                else if (kx == N[0]-2)
                    dx = -  a->field[fidx(N[0]-5, ky, kz)] +  6*a->field[fidx(N[0]-4, ky, kz)] - 18*a->field[fidx(N[0]-3, ky, kz)] + 10*a->field[fidx(N[0]-2, ky, kz)] +  3*a->field[fidx(N[0]-1, ky, kz)];
                else if (kx == N[0]-1)
                    dx =  3*a->field[fidx(N[0]-5, ky, kz)] - 16*a->field[fidx(N[0]-4, ky, kz)] + 36*a->field[fidx(N[0]-3, ky, kz)] - 48*a->field[fidx(N[0]-2, ky, kz)] + 25*a->field[fidx(N[0]-1, ky, kz)];
                else
                    dx = - a->field[fidx(kx+2, ky, kz)] + 8*a->field[fidx(kx+1, ky, kz)] - 8*a->field[fidx(kx-1, ky, kz)] + a->field[fidx(kx-2, ky, kz)];

                if (ky == 0)
                    dy = - 25*a->field[fidx(kx, 0, kz)] + 48*a->field[fidx(kx, 1, kz)] - 36*a->field[fidx(kx, 2, kz)] + 16*a->field[fidx(kx, 3, kz)] - 3*a->field[fidx(kx, 4, kz)];
                else if (ky == 1)
                    dy = -  3*a->field[fidx(kx, 0, kz)] - 10*a->field[fidx(kx, 1, kz)] + 18*a->field[fidx(kx, 2, kz)] -  6*a->field[fidx(kx, 3, kz)] +   a->field[fidx(kx, 4, kz)];
                else if (ky == N[1]-2)
                    dy = -  a->field[fidx(kx, N[1]-5, kz)] +  6*a->field[fidx(kx, N[1]-4, kz)] - 18*a->field[fidx(kx, N[1]-3, kz)] + 10*a->field[fidx(kx, N[1]-2, kz)] +  3*a->field[fidx(kx, N[1]-1, kz)];
                else if (ky == N[1]-1)
                    dy =  3*a->field[fidx(kx, N[1]-5, kz)] - 16*a->field[fidx(kx, N[1]-4, kz)] + 36*a->field[fidx(kx, N[1]-3, kz)] - 48*a->field[fidx(kx, N[1]-2, kz)] + 25*a->field[fidx(kx, N[1]-1, kz)];
                else
                    dy = - a->field[fidx(kx, ky+2, kz)] + 8*a->field[fidx(kx, ky+1, kz)] - 8*a->field[fidx(kx, ky-1, kz)] + a->field[fidx(kx, ky-2, kz)];

                if (kz == 0)
                    dz = - 25*a->field[fidx(kx, ky, 0)] + 48*a->field[fidx(kx, ky, 1)] - 36*a->field[fidx(kx, ky, 2)] + 16*a->field[fidx(kx, ky, 3)] - 3*a->field[fidx(kx, ky, 4)];
                else if (kz == 1)
                    dz = -  3*a->field[fidx(kx, ky, 0)] - 10*a->field[fidx(kx, ky, 1)] + 18*a->field[fidx(kx, ky, 2)] -  6*a->field[fidx(kx, ky, 3)] +   a->field[fidx(kx, ky, 4)];
                else if (kz == N[2]-2)
                    dz = -  a->field[fidx(kx, ky, N[2]-5)] +  6*a->field[fidx(kx, ky, N[2]-4)] - 18*a->field[fidx(kx, ky, N[2]-3)] + 10*a->field[fidx(kx, ky, N[2]-2)] +  3*a->field[fidx(kx, ky, N[2]-1)];
                else if (kz == N[2]-1)
                    dz =  3*a->field[fidx(kx, ky, N[2]-5)] - 16*a->field[fidx(kx, ky, N[2]-4)] + 36*a->field[fidx(kx, ky, N[2]-3)] - 48*a->field[fidx(kx, ky, N[2]-2)] + 25*a->field[fidx(kx, ky, N[2]-1)];
                else
                    dz = - a->field[fidx(kx, ky, kz+2)] + 8*a->field[fidx(kx, ky, kz+1)] - 8*a->field[fidx(kx, ky, kz-1)] + a->field[fidx(kx, ky, kz-2)];

			    fieldX[fidx(kx, ky, kz)] = dx*step[0] / 12.0;
			    fieldY[fidx(kx, ky, kz)] = dy*step[1] / 12.0;
			    fieldZ[fidx(kx, ky, kz)] = dz*step[2] / 12.0;
		    }
            else // scheme == 7
            {
                // x
                if (kx == 0)
                    dx = s7m3_0 * a->field[fidx(0, ky, kz)] + s7m3_1 * a->field[fidx(1, ky, kz)] + s7m3_2 * a->field[fidx(2, ky, kz)] + s7m3_3 * a->field[fidx(3, ky, kz)]
                       + s7m3_4 * a->field[fidx(4, ky, kz)] + s7m3_5 * a->field[fidx(5, ky, kz)] + s7m3_6 * a->field[fidx(6, ky, kz)];
                else if (kx == 1)
                    dx = s7m2_0 * a->field[fidx(0, ky, kz)] + s7m2_1 * a->field[fidx(1, ky, kz)] + s7m2_2 * a->field[fidx(2, ky, kz)] + s7m2_3 * a->field[fidx(3, ky, kz)]
                       + s7m2_4 * a->field[fidx(4, ky, kz)] + s7m2_5 * a->field[fidx(5, ky, kz)] + s7m2_6 * a->field[fidx(6, ky, kz)];
                else if (kx == 2)
                    dx = s7m1_0 * a->field[fidx(0, ky, kz)] + s7m1_1 * a->field[fidx(1, ky, kz)] + s7m1_2 * a->field[fidx(2, ky, kz)] + s7m1_3 * a->field[fidx(3, ky, kz)]
                       + s7m1_4 * a->field[fidx(4, ky, kz)] + s7m1_5 * a->field[fidx(5, ky, kz)] + s7m1_6 * a->field[fidx(6, ky, kz)];
                else if (kx == N[0]-3)
                    dx = - (s7m1_6 * a->field[fidx(N[0]-7, ky, kz)] + s7m1_5 * a->field[fidx(N[0]-6, ky, kz)] + s7m1_4 * a->field[fidx(N[0]-5, ky, kz)] + s7m1_3 * a->field[fidx(N[0]-4, ky, kz)]
                          + s7m1_2 * a->field[fidx(N[0]-3, ky, kz)] + s7m1_1 * a->field[fidx(N[0]-2, ky, kz)] + s7m1_0 * a->field[fidx(N[0]-1, ky, kz)]);
                else if (kx == N[0]-2)
                    dx = - (s7m2_6 * a->field[fidx(N[0]-7, ky, kz)] + s7m2_5 * a->field[fidx(N[0]-6, ky, kz)] + s7m2_4 * a->field[fidx(N[0]-5, ky, kz)] + s7m2_3 * a->field[fidx(N[0]-4, ky, kz)]
                          + s7m2_2 * a->field[fidx(N[0]-3, ky, kz)] + s7m2_1 * a->field[fidx(N[0]-2, ky, kz)] + s7m2_0 * a->field[fidx(N[0]-1, ky, kz)]);
                else if (kx == N[0]-1)
                    dx = - (s7m3_6 * a->field[fidx(N[0]-7, ky, kz)] + s7m3_5 * a->field[fidx(N[0]-6, ky, kz)] + s7m3_4 * a->field[fidx(N[0]-5, ky, kz)] + s7m3_3 * a->field[fidx(N[0]-4, ky, kz)]
                          + s7m3_2 * a->field[fidx(N[0]-3, ky, kz)] + s7m3_1 * a->field[fidx(N[0]-2, ky, kz)] + s7m3_0 * a->field[fidx(N[0]-1, ky, kz)]);
                else
                    dx = s7m0_6 * a->field[fidx(kx+3, ky, kz)] + s7m0_5 * a->field[fidx(kx+2, ky, kz)] + s7m0_4 * a->field[fidx(kx+1, ky, kz)]
                       + s7m0_2 * a->field[fidx(kx-1, ky, kz)] + s7m0_1 * a->field[fidx(kx-2, ky, kz)] + s7m0_0 * a->field[fidx(kx-3, ky, kz)];

                // y
                if (ky == 0)
                    dy = s7m3_0 * a->field[fidx(kx, 0, kz)] + s7m3_1 * a->field[fidx(kx, 1, kz)] + s7m3_2 * a->field[fidx(kx, 2, kz)] + s7m3_3 * a->field[fidx(kx, 3, kz)]
                       + s7m3_4 * a->field[fidx(kx, 4, kz)] + s7m3_5 * a->field[fidx(kx, 5, kz)] + s7m3_6 * a->field[fidx(kx, 6, kz)];
                else if (ky == 1)
                    dy = s7m2_0 * a->field[fidx(kx, 0, kz)] + s7m2_1 * a->field[fidx(kx, 1, kz)] + s7m2_2 * a->field[fidx(kx, 2, kz)] + s7m2_3 * a->field[fidx(kx, 3, kz)]
                       + s7m2_4 * a->field[fidx(kx, 4, kz)] + s7m2_5 * a->field[fidx(kx, 5, kz)] + s7m2_6 * a->field[fidx(kx, 6, kz)];
                else if (ky == 2)
                    dy = s7m1_0 * a->field[fidx(kx, 0, kz)] + s7m1_1 * a->field[fidx(kx, 1, kz)] + s7m1_2 * a->field[fidx(kx, 2, kz)] + s7m1_3 * a->field[fidx(kx, 3, kz)]
                       + s7m1_4 * a->field[fidx(kx, 4, kz)] + s7m1_5 * a->field[fidx(kx, 5, kz)] + s7m1_6 * a->field[fidx(kx, 6, kz)];
                else if (ky == N[1]-3)
                    dy = - (s7m1_6 * a->field[fidx(kx, N[1]-7, kz)] + s7m1_5 * a->field[fidx(kx, N[1]-6, kz)] + s7m1_4 * a->field[fidx(kx, N[1]-5, kz)] + s7m1_3 * a->field[fidx(kx, N[1]-4, kz)]
                          + s7m1_2 * a->field[fidx(kx, N[1]-3, kz)] + s7m1_1 * a->field[fidx(kx, N[1]-2, kz)] + s7m1_0 * a->field[fidx(kx, N[1]-1, kz)]);
                else if (ky == N[1]-2)
                    dy = - (s7m2_6 * a->field[fidx(kx, N[1]-7, kz)] + s7m2_5 * a->field[fidx(kx, N[1]-6, kz)] + s7m2_4 * a->field[fidx(kx, N[1]-5, kz)] + s7m2_3 * a->field[fidx(kx, N[1]-4, kz)]
                          + s7m2_2 * a->field[fidx(kx, N[1]-3, kz)] + s7m2_1 * a->field[fidx(kx, N[1]-2, kz)] + s7m2_0 * a->field[fidx(kx, N[1]-1, kz)]);
                else if (ky == N[1]-1)
                    dy = - (s7m3_6 * a->field[fidx(kx, N[1]-7, kz)] + s7m3_5 * a->field[fidx(kx, N[1]-6, kz)] + s7m3_4 * a->field[fidx(kx, N[1]-5, kz)] + s7m3_3 * a->field[fidx(kx, N[1]-4, kz)]
                          + s7m3_2 * a->field[fidx(kx, N[1]-3, kz)] + s7m3_1 * a->field[fidx(kx, N[1]-2, kz)] + s7m3_0 * a->field[fidx(kx, N[1]-1, kz)]);
                else
                    dy = s7m0_6 * a->field[fidx(kx, ky+3, kz)] + s7m0_5 * a->field[fidx(kx, ky+2, kz)] + s7m0_4 * a->field[fidx(kx, ky+1, kz)]
                       + s7m0_2 * a->field[fidx(kx, ky-1, kz)] + s7m0_1 * a->field[fidx(kx, ky-2, kz)] + s7m0_0 * a->field[fidx(kx, ky-3, kz)];

                // z
                if (kz == 0)
                    dz = s7m3_0 * a->field[fidx(kx, ky, 0)] + s7m3_1 * a->field[fidx(kx, ky, 1)] + s7m3_2 * a->field[fidx(kx, ky, 2)] + s7m3_3 * a->field[fidx(kx, ky, 3)]
                       + s7m3_4 * a->field[fidx(kx, ky, 4)] + s7m3_5 * a->field[fidx(kx, ky, 5)] + s7m3_6 * a->field[fidx(kx, ky, 6)];
                else if (kz == 1)
                    dz = s7m2_0 * a->field[fidx(kx, ky, 0)] + s7m2_1 * a->field[fidx(kx, ky, 1)] + s7m2_2 * a->field[fidx(kx, ky, 2)] + s7m2_3 * a->field[fidx(kx, ky, 3)]
                       + s7m2_4 * a->field[fidx(kx, ky, 4)] + s7m2_5 * a->field[fidx(kx, ky, 5)] + s7m2_6 * a->field[fidx(kx, ky, 6)];
                else if (kz == 2)
                    dz = s7m1_0 * a->field[fidx(kx, ky, 0)] + s7m1_1 * a->field[fidx(kx, ky, 1)] + s7m1_2 * a->field[fidx(kx, ky, 2)] + s7m1_3 * a->field[fidx(kx, ky, 3)]
                       + s7m1_4 * a->field[fidx(kx, ky, 4)] + s7m1_5 * a->field[fidx(kx, ky, 5)] + s7m1_6 * a->field[fidx(kx, ky, 6)];
                else if (kz == N[2]-3)
                    dz = - (s7m1_6 * a->field[fidx(kx, ky, N[2]-7)] + s7m1_5 * a->field[fidx(kx, ky, N[2]-6)] + s7m1_4 * a->field[fidx(kx, ky, N[2]-5)] + s7m1_3 * a->field[fidx(kx, ky, N[2]-4)]
                          + s7m1_2 * a->field[fidx(kx, ky, N[2]-3)] + s7m1_1 * a->field[fidx(kx, ky, N[2]-2)] + s7m1_0 * a->field[fidx(kx, ky, N[2]-1)]);
                else if (kz == N[2]-2)
                    dz = - (s7m2_6 * a->field[fidx(kx, ky, N[2]-7)] + s7m2_5 * a->field[fidx(kx, ky, N[2]-6)] + s7m2_4 * a->field[fidx(kx, ky, N[2]-5)] + s7m2_3 * a->field[fidx(kx, ky, N[2]-4)]
                          + s7m2_2 * a->field[fidx(kx, ky, N[2]-3)] + s7m2_1 * a->field[fidx(kx, ky, N[2]-2)] + s7m2_0 * a->field[fidx(kx, ky, N[2]-1)]);
                else if (kz == N[2]-1)
                    dz = - (s7m3_6 * a->field[fidx(kx, ky, N[2]-7)] + s7m3_5 * a->field[fidx(kx, ky, N[2]-6)] + s7m3_4 * a->field[fidx(kx, ky, N[2]-5)] + s7m3_3 * a->field[fidx(kx, ky, N[2]-4)]
                          + s7m3_2 * a->field[fidx(kx, ky, N[2]-3)] + s7m3_1 * a->field[fidx(kx, ky, N[2]-2)] + s7m3_0 * a->field[fidx(kx, ky, N[2]-1)]);
                else
                    dz = s7m0_6 * a->field[fidx(kx, ky, kz+3)] + s7m0_5 * a->field[fidx(kx, ky, kz+2)] + s7m0_4 * a->field[fidx(kx, ky, kz+1)]
                       + s7m0_2 * a->field[fidx(kx, ky, kz-1)] + s7m0_1 * a->field[fidx(kx, ky, kz-2)] + s7m0_0 * a->field[fidx(kx, ky, kz-3)];

			    fieldX[fidx(kx, ky, kz)] = dx*step[0] / 60.0;
			    fieldY[fidx(kx, ky, kz)] = dy*step[1] / 60.0;
			    fieldZ[fidx(kx, ky, kz)] = dz*step[2] / 60.0;
            }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::mult_plane(double c, CagmVectorFieldOps *a, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            fieldX[fidx(kx, ky, kz)] = a->fieldX[fidx(kx, ky, kz)] * c;
            fieldY[fidx(kx, ky, kz)] = a->fieldY[fidx(kx, ky, kz)] * c;
            fieldZ[fidx(kx, ky, kz)] = a->fieldZ[fidx(kx, ky, kz)] * c;
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::mult_plane(CagmScalarFieldOps *c, CagmVectorFieldOps *a, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            fieldX[fidx(kx, ky, kz)] = (a->fieldX[fidx(kx, ky, kz)]) * (c->field[fidx(kx, ky, kz)]);
            fieldY[fidx(kx, ky, kz)] = (a->fieldY[fidx(kx, ky, kz)]) * (c->field[fidx(kx, ky, kz)]);
            fieldZ[fidx(kx, ky, kz)] = (a->fieldZ[fidx(kx, ky, kz)]) * (c->field[fidx(kx, ky, kz)]);
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::add_plane(CagmVectorFieldOps *a, CagmVectorFieldOps *b, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            fieldX[fidx(kx, ky, kz)] = a->fieldX[fidx(kx, ky, kz)] + b->fieldX[fidx(kx, ky, kz)];
            fieldY[fidx(kx, ky, kz)] = a->fieldY[fidx(kx, ky, kz)] + b->fieldY[fidx(kx, ky, kz)];
            fieldZ[fidx(kx, ky, kz)] = a->fieldZ[fidx(kx, ky, kz)] + b->fieldZ[fidx(kx, ky, kz)];
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::sub_plane(CagmVectorFieldOps *a, CagmVectorFieldOps *b, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            fieldX[fidx(kx, ky, kz)] = a->fieldX[fidx(kx, ky, kz)] - b->fieldX[fidx(kx, ky, kz)];
            fieldY[fidx(kx, ky, kz)] = a->fieldY[fidx(kx, ky, kz)] - b->fieldY[fidx(kx, ky, kz)];
            fieldZ[fidx(kx, ky, kz)] = a->fieldZ[fidx(kx, ky, kz)] - b->fieldZ[fidx(kx, ky, kz)];
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::neg_plane(CagmVectorFieldOps *a, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            fieldX[fidx(kx, ky, kz)] = -a->fieldX[fidx(kx, ky, kz)];
            fieldY[fidx(kx, ky, kz)] = -a->fieldY[fidx(kx, ky, kz)];
            fieldZ[fidx(kx, ky, kz)] = -a->fieldZ[fidx(kx, ky, kz)];
        }

    return 0;
}

//-----------------------------------------------------------------------
double CagmVectorFieldOps::max2_plane(int kz)
{
    double tmax = 0;
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            double t = fieldX[fidx(kx, ky, kz)] * fieldX[fidx(kx, ky, kz)] + fieldY[fidx(kx, ky, kz)] * fieldY[fidx(kx, ky, kz)] + fieldZ[fidx(kx, ky, kz)] * fieldZ[fidx(kx, ky, kz)];
            if (t > tmax)
                tmax = t;
        }

    return tmax;
}
