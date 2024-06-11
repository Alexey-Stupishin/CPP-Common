#include "stdDefinitions.h"
#include <math.h>
#include <cfloat>

#include "agsFieldsCommon.h"
#include "agmScalarFieldOps.h"
#include "agmVectorFieldOps.h"
#include "DiffCoefs.h"

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::div_plane(CagmVectorFieldOps *a, int kz, int scheme)
{
    double dx, dy, dz;
    for (int ky = NL[1]; ky < NH[1]; ky++)
        for (int kx = NL[0]; kx < NH[0]; kx++)
            if (scheme == 3)
            {
                if (kx == 0)
                    dx = -3 * a->fieldX[fidx(0, ky, kz)] + 4 * a->fieldX[fidx(1, ky, kz)] - a->fieldX[fidx(2, ky, kz)];
                else if (kx == N[0] - 1)
                    dx = a->fieldX[fidx(N[0] - 3, ky, kz)] - 4 * a->fieldX[fidx(N[0] - 2, ky, kz)] + 3 * a->fieldX[fidx(N[0] - 1, ky, kz)];
                else
                    dx = a->fieldX[fidx(kx + 1, ky, kz)] - a->fieldX[fidx(kx - 1, ky, kz)];

                if (ky == 0)
                    dy = -3 * a->fieldY[fidx(kx, 0, kz)] + 4 * a->fieldY[fidx(kx, 1, kz)] - a->fieldY[fidx(kx, 2, kz)];
                else if (ky == N[1] - 1)
                    dy = a->fieldY[fidx(kx, N[1] - 3, kz)] - 4 * a->fieldY[fidx(kx, N[1] - 2, kz)] + 3 * a->fieldY[fidx(kx, N[1] - 1, kz)];
                else
                    dy = a->fieldY[fidx(kx, ky + 1, kz)] - a->fieldY[fidx(kx, ky - 1, kz)];

                if (kz == 0)
                    dz = -3 * a->fieldZ[fidx(kx, ky, 0)] + 4 * a->fieldZ[fidx(kx, ky, 1)] - a->fieldZ[fidx(kx, ky, 2)];
                else if (kz == N[2] - 1)
                    dz = a->fieldZ[fidx(kx, ky, N[2] - 3)] - 4 * a->fieldZ[fidx(kx, ky, N[2] - 2)] + 3 * a->fieldZ[fidx(kx, ky, N[2] - 1)];
                else
                    dz = a->fieldZ[fidx(kx, ky, kz + 1)] - a->fieldZ[fidx(kx, ky, kz - 1)];

                field[fidx(kx, ky, kz)] = (dx*step[0] + dy*step[1] + dz*step[2]) * 0.5;
            }
            else if (scheme == 5)
            {
                if (kx == 0)
                    dx = - 25*a->fieldX[fidx(0, ky, kz)] + 48*a->fieldX[fidx(1, ky, kz)] - 36*a->fieldX[fidx(2, ky, kz)] + 16*a->fieldX[fidx(3, ky, kz)] - 3*a->fieldX[fidx(4, ky, kz)];
                else if (kx == 1)
                    dx = -  3*a->fieldX[fidx(0, ky, kz)] - 10*a->fieldX[fidx(1, ky, kz)] + 18*a->fieldX[fidx(2, ky, kz)] -  6*a->fieldX[fidx(3, ky, kz)] +   a->fieldX[fidx(4, ky, kz)];
                else if (kx == N[0]-2)
                    dx = -  a->fieldX[fidx(N[0]-5, ky, kz)] +  6*a->fieldX[fidx(N[0]-4, ky, kz)] - 18*a->fieldX[fidx(N[0]-3, ky, kz)] + 10*a->fieldX[fidx(N[0]-2, ky, kz)] +  3*a->fieldX[fidx(N[0]-1, ky, kz)];
                else if (kx == N[0]-1)
                    dx =  3*a->fieldX[fidx(N[0]-5, ky, kz)] - 16*a->fieldX[fidx(N[0]-4, ky, kz)] + 36*a->fieldX[fidx(N[0]-3, ky, kz)] - 48*a->fieldX[fidx(N[0]-2, ky, kz)] + 25*a->fieldX[fidx(N[0]-1, ky, kz)];
                else
                    dx = - a->fieldX[fidx(kx+2, ky, kz)] + 8*a->fieldX[fidx(kx+1, ky, kz)] - 8*a->fieldX[fidx(kx-1, ky, kz)] + a->fieldX[fidx(kx-2, ky, kz)];

                if (ky == 0)
                    dy = - 25*a->fieldY[fidx(kx, 0, kz)] + 48*a->fieldY[fidx(kx, 1, kz)] - 36*a->fieldY[fidx(kx, 2, kz)] + 16*a->fieldY[fidx(kx, 3, kz)] - 3*a->fieldY[fidx(kx, 4, kz)];
                else if (ky == 1)
                    dy = -  3*a->fieldY[fidx(kx, 0, kz)] - 10*a->fieldY[fidx(kx, 1, kz)] + 18*a->fieldY[fidx(kx, 2, kz)] -  6*a->fieldY[fidx(kx, 3, kz)] +   a->fieldY[fidx(kx, 4, kz)];
                else if (ky == N[1]-2)
                    dy = -  a->fieldY[fidx(kx, N[1]-5, kz)] +  6*a->fieldY[fidx(kx, N[1]-4, kz)] - 18*a->fieldY[fidx(kx, N[1]-3, kz)] + 10*a->fieldY[fidx(kx, N[1]-2, kz)] +  3*a->fieldY[fidx(kx, N[1]-1, kz)];
                else if (ky == N[1]-1)
                    dy =  3*a->fieldY[fidx(kx, N[1]-5, kz)] - 16*a->fieldY[fidx(kx, N[1]-4, kz)] + 36*a->fieldY[fidx(kx, N[1]-3, kz)] - 48*a->fieldY[fidx(kx, N[1]-2, kz)] + 25*a->fieldY[fidx(kx, N[1]-1, kz)];
                else
                    dy = - a->fieldY[fidx(kx, ky+2, kz)] + 8*a->fieldY[fidx(kx, ky+1, kz)] - 8*a->fieldY[fidx(kx, ky-1, kz)] + a->fieldY[fidx(kx, ky-2, kz)];

                if (kz == 0)
                    dz = - 25*a->fieldZ[fidx(kx, ky, 0)] + 48*a->fieldZ[fidx(kx, ky, 1)] - 36*a->fieldZ[fidx(kx, ky, 2)] + 16*a->fieldZ[fidx(kx, ky, 3)] - 3*a->fieldZ[fidx(kx, ky, 4)];
                else if (kz == 1)
                    dz = -  3*a->fieldZ[fidx(kx, ky, 0)] - 10*a->fieldZ[fidx(kx, ky, 1)] + 18*a->fieldZ[fidx(kx, ky, 2)] -  6*a->fieldZ[fidx(kx, ky, 3)] +   a->fieldZ[fidx(kx, ky, 4)];
                else if (kz == N[2]-2)
                    dz = -  a->fieldZ[fidx(kx, ky, N[2]-5)] +  6*a->fieldZ[fidx(kx, ky, N[2]-4)] - 18*a->fieldZ[fidx(kx, ky, N[2]-3)] + 10*a->fieldZ[fidx(kx, ky, N[2]-2)] +  3*a->fieldZ[fidx(kx, ky, N[2]-1)];
                else if (kz == N[2]-1)
                    dz =  3*a->fieldZ[fidx(kx, ky, N[2]-5)] - 16*a->fieldZ[fidx(kx, ky, N[2]-4)] + 36*a->fieldZ[fidx(kx, ky, N[2]-3)] - 48*a->fieldZ[fidx(kx, ky, N[2]-2)] + 25*a->fieldZ[fidx(kx, ky, N[2]-1)];
                else
                    dz = - a->fieldZ[fidx(kx, ky, kz+2)] + 8*a->fieldZ[fidx(kx, ky, kz+1)] - 8*a->fieldZ[fidx(kx, ky, kz-1)] + a->fieldZ[fidx(kx, ky, kz-2)];

			    field[fidx(kx, ky, kz)] = (dx*step[0] + dy*step[1] + dz*step[2]) / 12.0;
            }
            else // scheme == 7
            {
                if (kx == 0)
                    dx =    s7m3_0 * a->fieldX[fidx(0, ky, kz)] + s7m3_1 * a->fieldX[fidx(1, ky, kz)] + s7m3_2 * a->fieldX[fidx(2, ky, kz)] + s7m3_3 * a->fieldX[fidx(3, ky, kz)]
                          + s7m3_4 * a->fieldX[fidx(4, ky, kz)] + s7m3_5 * a->fieldX[fidx(5, ky, kz)] + s7m3_6 * a->fieldX[fidx(6, ky, kz)];
                else if (kx == 1)
                    dx =    s7m2_0 * a->fieldX[fidx(0, ky, kz)] + s7m2_1 * a->fieldX[fidx(1, ky, kz)] + s7m2_2 * a->fieldX[fidx(2, ky, kz)] + s7m2_3 * a->fieldX[fidx(3, ky, kz)]
                          + s7m2_4 * a->fieldX[fidx(4, ky, kz)] + s7m2_5 * a->fieldX[fidx(5, ky, kz)] + s7m2_6 * a->fieldX[fidx(6, ky, kz)];
                else if (kx == 2)
                    dx =    s7m1_0 * a->fieldX[fidx(0, ky, kz)] + s7m1_1 * a->fieldX[fidx(1, ky, kz)] + s7m1_2 * a->fieldX[fidx(2, ky, kz)] + s7m1_3 * a->fieldX[fidx(3, ky, kz)]
                          + s7m1_4 * a->fieldX[fidx(4, ky, kz)] + s7m1_5 * a->fieldX[fidx(5, ky, kz)] + s7m1_6 * a->fieldX[fidx(6, ky, kz)];
                else if (kx == N[0]-3)
                    dx = - (s7m1_6 * a->fieldX[fidx(N[0]-7, ky, kz)] + s7m1_5 * a->fieldX[fidx(N[0]-6, ky, kz)] + s7m1_4 * a->fieldX[fidx(N[0]-5, ky, kz)] + s7m1_3 * a->fieldX[fidx(N[0]-4, ky, kz)]
                          + s7m1_2 * a->fieldX[fidx(N[0]-3, ky, kz)] + s7m1_1 * a->fieldX[fidx(N[0]-2, ky, kz)] + s7m1_0 * a->fieldX[fidx(N[0]-1, ky, kz)]);
                else if (kx == N[0]-2)
                    dx = - (s7m2_6 * a->fieldX[fidx(N[0]-7, ky, kz)] + s7m2_5 * a->fieldX[fidx(N[0]-6, ky, kz)] + s7m2_4 * a->fieldX[fidx(N[0]-5, ky, kz)] + s7m2_3 * a->fieldX[fidx(N[0]-4, ky, kz)]
                          + s7m2_2 * a->fieldX[fidx(N[0]-3, ky, kz)] + s7m2_1 * a->fieldX[fidx(N[0]-2, ky, kz)] + s7m2_0 * a->fieldX[fidx(N[0]-1, ky, kz)]);
                else if (kx == N[0]-1)
                    dx = - (s7m3_6 * a->fieldX[fidx(N[0]-7, ky, kz)] + s7m3_5 * a->fieldX[fidx(N[0]-6, ky, kz)] + s7m3_4 * a->fieldX[fidx(N[0]-5, ky, kz)] + s7m3_3 * a->fieldX[fidx(N[0]-4, ky, kz)]
                          + s7m3_2 * a->fieldX[fidx(N[0]-3, ky, kz)] + s7m3_1 * a->fieldX[fidx(N[0]-2, ky, kz)] + s7m3_0 * a->fieldX[fidx(N[0]-1, ky, kz)]);
                else
                    dx = s7m0_6 * a->fieldX[fidx(kx+3, ky, kz)] + s7m0_5 * a->fieldX[fidx(kx+2, ky, kz)] + s7m0_4 * a->fieldX[fidx(kx+1, ky, kz)] 
                       + s7m0_2 * a->fieldX[fidx(kx-1, ky, kz)] + s7m0_1 * a->fieldX[fidx(kx-2, ky, kz)] + s7m0_0 * a->fieldX[fidx(kx-3, ky, kz)];

                if (ky == 0)
                    dy =    s7m3_0 * a->fieldY[fidx(kx, 0, kz)] + s7m3_1 * a->fieldY[fidx(kx, 1, kz)] + s7m3_2 * a->fieldY[fidx(kx, 2, kz)] + s7m3_3 * a->fieldY[fidx(kx, 3, kz)]
                          + s7m3_4 * a->fieldY[fidx(kx, 4, kz)] + s7m3_5 * a->fieldY[fidx(kx, 5, kz)] + s7m3_6 * a->fieldY[fidx(kx, 6, kz)];
                else if (ky == 1)
                    dy =    s7m2_0 * a->fieldY[fidx(kx, 0, kz)] + s7m2_1 * a->fieldY[fidx(kx, 1, kz)] + s7m2_2 * a->fieldY[fidx(kx, 2, kz)] + s7m2_3 * a->fieldY[fidx(kx, 3, kz)]
                          + s7m2_4 * a->fieldY[fidx(kx, 4, kz)] + s7m2_5 * a->fieldY[fidx(kx, 5, kz)] + s7m2_6 * a->fieldY[fidx(kx, 6, kz)];
                else if (ky == 2)
                    dy =    s7m1_0 * a->fieldY[fidx(kx, 0, kz)] + s7m1_1 * a->fieldY[fidx(kx, 1, kz)] + s7m1_2 * a->fieldY[fidx(kx, 2, kz)] + s7m1_3 * a->fieldY[fidx(kx, 3, kz)]
                          + s7m1_4 * a->fieldY[fidx(kx, 4, kz)] + s7m1_5 * a->fieldY[fidx(kx, 5, kz)] + s7m1_6 * a->fieldY[fidx(kx, 6, kz)];
                else if (ky == N[1]-3)
                    dy = - (s7m1_6 * a->fieldY[fidx(kx, N[1]-7, kz)] + s7m1_5 * a->fieldY[fidx(kx, N[1]-6, kz)] + s7m1_4 * a->fieldY[fidx(kx, N[1]-5, kz)] + s7m1_3 * a->fieldY[fidx(kx, N[1]-4, kz)]
                          + s7m1_2 * a->fieldY[fidx(kx, N[1]-3, kz)] + s7m1_1 * a->fieldY[fidx(kx, N[1]-2, kz)] + s7m1_0 * a->fieldY[fidx(kx, N[1]-1, kz)]);
                else if (ky == N[1]-2)
                    dy = - (s7m2_6 * a->fieldY[fidx(kx, N[1]-7, kz)] + s7m2_5 * a->fieldY[fidx(kx, N[1]-6, kz)] + s7m2_4 * a->fieldY[fidx(kx, N[1]-5, kz)] + s7m2_3 * a->fieldY[fidx(kx, N[1]-4, kz)]
                          + s7m2_2 * a->fieldY[fidx(kx, N[1]-3, kz)] + s7m2_1 * a->fieldY[fidx(kx, N[1]-2, kz)] + s7m2_0 * a->fieldY[fidx(kx, N[1]-1, kz)]);
                else if (ky == N[1]-1)
                    dy = - (s7m3_6 * a->fieldY[fidx(kx, N[1]-7, kz)] + s7m3_5 * a->fieldY[fidx(kx, N[1]-6, kz)] + s7m3_4 * a->fieldY[fidx(kx, N[1]-5, kz)] + s7m3_3 * a->fieldY[fidx(kx, N[1]-4, kz)]
                          + s7m3_2 * a->fieldY[fidx(kx, N[1]-3, kz)] + s7m3_1 * a->fieldY[fidx(kx, N[1]-2, kz)] + s7m3_0 * a->fieldY[fidx(kx, N[1]-1, kz)]);
                else
                    dy = s7m0_6 * a->fieldY[fidx(kx, ky+3, kz)] + s7m0_5 * a->fieldY[fidx(kx, ky+2, kz)] + s7m0_4 * a->fieldY[fidx(kx, ky+1, kz)] 
                       + s7m0_2 * a->fieldY[fidx(kx, ky-1, kz)] + s7m0_1 * a->fieldY[fidx(kx, ky-2, kz)] + s7m0_0 * a->fieldY[fidx(kx, ky-3, kz)];

                if (kz == 0)
                    dz =    s7m3_0 * a->fieldZ[fidx(kx, ky, 0)] + s7m3_1 * a->fieldZ[fidx(kx, ky, 1)] + s7m3_2 * a->fieldZ[fidx(kx, ky, 2)] + s7m3_3 * a->fieldZ[fidx(kx, ky, 3)]
                          + s7m3_4 * a->fieldZ[fidx(kx, ky, 4)] + s7m3_5 * a->fieldZ[fidx(kx, ky, 5)] + s7m3_6 * a->fieldZ[fidx(kx, ky, 6)];
                else if (kz == 1)
                    dz =    s7m2_0 * a->fieldZ[fidx(kx, ky, 0)] + s7m2_1 * a->fieldZ[fidx(kx, ky, 1)] + s7m2_2 * a->fieldZ[fidx(kx, ky, 2)] + s7m2_3 * a->fieldZ[fidx(kx, ky, 3)]
                          + s7m2_4 * a->fieldZ[fidx(kx, ky, 4)] + s7m2_5 * a->fieldZ[fidx(kx, ky, 5)] + s7m2_6 * a->fieldZ[fidx(kx, ky, 6)];
                else if (kz == 2)
                    dz =    s7m1_0 * a->fieldZ[fidx(kx, ky, 0)] + s7m1_1 * a->fieldZ[fidx(kx, ky, 1)] + s7m1_2 * a->fieldZ[fidx(kx, ky, 2)] + s7m1_3 * a->fieldZ[fidx(kx, ky, 3)]
                          + s7m1_4 * a->fieldZ[fidx(kx, ky, 4)] + s7m1_5 * a->fieldZ[fidx(kx, ky, 5)] + s7m1_6 * a->fieldZ[fidx(kx, ky, 6)];
                else if (kz == N[2]-3)
                    dz = - (s7m1_6 * a->fieldZ[fidx(kx, ky, N[2]-7)] + s7m1_5 * a->fieldZ[fidx(kx, ky, N[2]-6)] + s7m1_4 * a->fieldZ[fidx(kx, ky, N[2]-5)] + s7m1_3 * a->fieldZ[fidx(kx, ky, N[2]-4)]
                          + s7m1_2 * a->fieldZ[fidx(kx, ky, N[2]-3)] + s7m1_1 * a->fieldZ[fidx(kx, ky, N[2]-2)] + s7m1_0 * a->fieldZ[fidx(kx, ky, N[2]-1)]);
                else if (kz == N[2]-2)
                    dz = - (s7m2_6 * a->fieldZ[fidx(kx, ky, N[2]-7)] + s7m2_5 * a->fieldZ[fidx(kx, ky, N[2]-6)] + s7m2_4 * a->fieldZ[fidx(kx, ky, N[2]-5)] + s7m2_3 * a->fieldZ[fidx(kx, ky, N[2]-4)]
                          + s7m2_2 * a->fieldZ[fidx(kx, ky, N[2]-3)] + s7m2_1 * a->fieldZ[fidx(kx, ky, N[2]-2)] + s7m2_0 * a->fieldZ[fidx(kx, ky, N[2]-1)]);
                else if (kz == N[2]-1)
                    dz = - (s7m3_6 * a->fieldZ[fidx(kx, ky, N[2]-7)] + s7m3_5 * a->fieldZ[fidx(kx, ky, N[2]-6)] + s7m3_4 * a->fieldZ[fidx(kx, ky, N[2]-5)] + s7m3_3 * a->fieldZ[fidx(kx, ky, N[2]-4)]
                          + s7m3_2 * a->fieldZ[fidx(kx, ky, N[2]-3)] + s7m3_1 * a->fieldZ[fidx(kx, ky, N[2]-2)] + s7m3_0 * a->fieldZ[fidx(kx, ky, N[2]-1)]);
                else
                    dz = s7m0_6 * a->fieldZ[fidx(kx, ky, kz+3)] + s7m0_5 * a->fieldZ[fidx(kx, ky, kz+2)] + s7m0_4 * a->fieldZ[fidx(kx, ky, kz+1)] 
                       + s7m0_2 * a->fieldZ[fidx(kx, ky, kz-1)] + s7m0_1 * a->fieldZ[fidx(kx, ky, kz-2)] + s7m0_0 * a->fieldZ[fidx(kx, ky, kz-3)];

                field[fidx(kx, ky, kz)] = (dx*step[0] + dy*step[1] + dz*step[2]) / 60.0;
            }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::dot_plane(CagmVectorFieldOps *a, CagmVectorFieldOps *b, int kz, CagmScalarFieldOps *Weight)
{
    double w = 1;
    for (int ky = NL[1]; ky < NH[1]; ky++)
        for (int kx = NL[0]; kx < NH[0]; kx++)
        {
            if (Weight)
                w = Weight->field[fidx(kx, ky, kz)];
            field[fidx(kx, ky, kz)] = (a->fieldX[fidx(kx, ky, kz)] * b->fieldX[fidx(kx, ky, kz)]
                                     + a->fieldY[fidx(kx, ky, kz)] * b->fieldY[fidx(kx, ky, kz)]
                                     + a->fieldZ[fidx(kx, ky, kz)] * b->fieldZ[fidx(kx, ky, kz)]
                ) * w;
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::abs2_plane(CagmVectorFieldOps *a, int kz, CagmScalarFieldOps *Weight)
{
    return dot_plane(a, a, kz, Weight);
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::abs_plane(CagmVectorFieldOps *a, int kz, CagmScalarFieldOps *Weight)
{
    double w = 1;
    for (int ky = NL[1]; ky < NH[1]; ky++)
        for (int kx = NL[0]; kx < NH[0]; kx++)
        {
            if (Weight)
                w = Weight->field[fidx(kx, ky, kz)];
            field[fidx(kx, ky, kz)] = sqrt((a->fieldX[fidx(kx, ky, kz)] * a->fieldX[fidx(kx, ky, kz)]
                                          + a->fieldY[fidx(kx, ky, kz)] * a->fieldY[fidx(kx, ky, kz)]
                                          + a->fieldZ[fidx(kx, ky, kz)] * a->fieldZ[fidx(kx, ky, kz)]
                                           ) * w
                                          );
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::sqrt_plane(CagmScalarFieldOps *a, int kz)
{
    for (int ky = NL[1]; ky < NH[1]; ky++)
        for (int kx = NL[0]; kx < NH[0]; kx++)
            field[fidx(kx, ky, kz)] = sqrt(a->field[fidx(kx, ky, kz)]);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::inv_plane(CagmScalarFieldOps *a, int kz)
{
    for (int ky = NL[1]; ky < NH[1]; ky++)
        for (int kx = NL[0]; kx < NH[0]; kx++)
            field[fidx(kx, ky, kz)] = (a->field[fidx(kx, ky, kz)] < tolerance_zero ? 0.0 : 1.0 / (a->field[fidx(kx, ky, kz)]+tolerance_denom));

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::invabs_plane(CagmVectorFieldOps *a, int kz, CagmScalarFieldOps *Weight)
{
    double w = 1;
    for (int ky = NL[1]; ky < NH[1]; ky++)
        for (int kx = NL[0]; kx < NH[0]; kx++)
        {
            if (Weight)
                w = Weight->field[fidx(kx, ky, kz)];
            double a2 = sqrt((a->fieldX[fidx(kx, ky, kz)] * a->fieldX[fidx(kx, ky, kz)]
                            + a->fieldY[fidx(kx, ky, kz)] * a->fieldY[fidx(kx, ky, kz)]
                            + a->fieldZ[fidx(kx, ky, kz)] * a->fieldZ[fidx(kx, ky, kz)]
                             ) * w
                            );

           field[fidx(kx, ky, kz)] = (a2 < tolerance_zero ? 0.0 : 1.0 / (a2+tolerance_denom));
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::mult_plane(double c, CagmScalarFieldOps *a, int kz)
{
    for (int ky = NL[1]; ky < NH[1]; ky++)
        for (int kx = NL[0]; kx < NH[0]; kx++)
            field[fidx(kx, ky, kz)] = a->field[fidx(kx, ky, kz)] * c;

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::mult_plane(CagmScalarFieldOps *b, CagmScalarFieldOps *a, int kz)
{
    for (int ky = NL[1]; ky < NH[1]; ky++)
        for (int kx = NL[0]; kx < NH[0]; kx++)
            if (a->field[fidx(kx, ky, kz)] == 0)
                field[fidx(kx, ky, kz)] = 0;
            else
                field[fidx(kx, ky, kz)] = a->field[fidx(kx, ky, kz)] * b->field[fidx(kx, ky, kz)];

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::add_plane(CagmScalarFieldOps *a, CagmScalarFieldOps *b, int kz)
{
    for (int ky = NL[1]; ky < NH[1]; ky++)
        for (int kx = NL[0]; kx < NH[0]; kx++)
            field[fidx(kx, ky, kz)] = a->field[fidx(kx, ky, kz)] + b->field[fidx(kx, ky, kz)];

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::sub_plane(CagmScalarFieldOps *a, CagmScalarFieldOps *b, int kz)
{
    for (int ky = NL[1]; ky < NH[1]; ky++)
        for (int kx = NL[0]; kx < NH[0]; kx++)
            field[fidx(kx, ky, kz)] = a->field[fidx(kx, ky, kz)] - b->field[fidx(kx, ky, kz)];

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::neg_plane(CagmScalarFieldOps *a, int kz)
{
    for (int ky = NL[1]; ky < NH[1]; ky++)
        for (int kx = NL[0]; kx < NH[0]; kx++)
            field[fidx(kx, ky, kz)] = -a->field[fidx(kx, ky, kz)];

    return 0;
}

//-----------------------------------------------------------------------
double CagmScalarFieldOps::sum_plane(int kz, CagmScalarFieldOps *weight)
{
    double wsum = 0;
    for (int ky = NL[1]; ky < NH[1]; ky++)
        for (int kx = NL[0]; kx < NH[0]; kx++)
        {
            if (weight)
    			wsum += field[fidx(kx, ky, kz)] * weight->field[fidx(kx, ky, kz)];
            else
    			wsum += field[fidx(kx, ky, kz)];
        }

    return wsum;
}

//-----------------------------------------------------------------------
double CagmScalarFieldOps::max_plane(int kz)
{
    double wmax = - DBL_MAX;
    for (int ky = NL[1]; ky < NH[1]; ky++)
        for (int kx = NL[0]; kx < NH[0]; kx++)
            if (field[fidx(kx, ky, kz)] > wmax)
                wmax = field[fidx(kx, ky, kz)];

    return wmax;
}
