#include "stdDefinitions.h"
#include <math.h>
#include <cfloat>

#include "mfoGlobals.h"

#include "agsFieldsCommon.h"
#include "agmScalarFieldOps.h"
#include "agmVectorFieldOps.h"

//-----------------------------------------------------------------------
CagmScalarFieldOps::CagmScalarFieldOps(int *_N, double *_step, int *_NL, int *_NH)
    : CubeXD(_N, 1, _step, _NL, _NH)
{
    Initialize();
}

//-----------------------------------------------------------------------
CagmScalarFieldOps::CagmScalarFieldOps(CubeXD *from)
    : CubeXD(from)
{
    Initialize();
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::Initialize()
{
    field = new double * [N[1]*N[2]];

    tolerance_zero = WiegelmannInversionTolerance;
    tolerance_denom = WiegelmannInversionDenom;

    return 0;
}

//-----------------------------------------------------------------------
CagmScalarFieldOps::~CagmScalarFieldOps()
{
    Delete();
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::Delete()
{
    delete [] field;

    return 0;
}

//-----------------------------------------------------------------------
double *CagmScalarFieldOps::getAddress(int kx, int ky, int kz)
{
    return &field[fidx(kx, ky, kz)];
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::stretch(CagmScalarFieldOps*src)
{
    double cx = (double)(src->N[0] - 1)/(double)(N[0] - 1);
    double cy = (double)(src->N[1] - 1)/(double)(N[1] - 1);
    double cz = (double)(src->N[2] - 1)/(double)(N[2] - 1);

	int x1, y1, z1;
    double tx, ty, tz;
    for (int kz = NL[2]; kz < NH[2]; kz++)
    {
        l_div_position((kz*cz), (src->N[2]), z1, tz);
        for (int ky = NL[1]; ky < NH[1]; ky++)
        {
            l_div_position(ky*cy, src->N[1], y1, ty);
            for (int kx = NL[0]; kx < NH[0]; kx++)
			{
                l_div_position(kx*cx, src->N[0], x1, tx);
                field[fidx(kx, ky, kz)] = l_get_point(src->N[1], src->field, x1, y1, z1, tx, ty, tz);
            }
        }
    }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::div(CagmVectorFieldOps *a)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        div_plane(a, kz);

    return 0;
}

    
//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::dot(CagmVectorFieldOps *a, CagmVectorFieldOps *b, CagmScalarFieldOps *Weight)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        dot_plane(a, b, kz, Weight);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::abs2(CagmVectorFieldOps *a, CagmScalarFieldOps *Weight)
{
    return dot(a, a, Weight);
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::abs(CagmVectorFieldOps *a)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        for (int ky = NL[1]; ky < NH[1]; ky++)
            for (int kx = NL[0]; kx < NH[0]; kx++)
				field[fidx(kx, ky, kz)] = sqrt(a->fieldX[fidx(kx, ky, kz)]*a->fieldX[fidx(kx, ky, kz)] +
                                               a->fieldY[fidx(kx, ky, kz)]*a->fieldY[fidx(kx, ky, kz)] +
                                               a->fieldZ[fidx(kx, ky, kz)]*a->fieldZ[fidx(kx, ky, kz)]);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::projection(CagmVectorFieldOps *a, double *d)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        for (int ky = NL[1]; ky < NH[1]; ky++)
            for (int kx = NL[0]; kx < NH[0]; kx++)
				field[fidx(kx, ky, kz)] = a->fieldX[fidx(kx, ky, kz)] * d[0] +
                                          a->fieldY[fidx(kx, ky, kz)] * d[1] +
                                          a->fieldZ[fidx(kx, ky, kz)] * d[2];

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::inv(CagmScalarFieldOps *a)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        inv_plane(a, kz);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::inv(void)
{
    return inv(this);
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::mult(double c, CagmScalarFieldOps *a)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        mult_plane(c, a, kz);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::mult(double c)
{
    return mult(c, this);
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::mult(CagmScalarFieldOps *c, CagmScalarFieldOps *a)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        mult_plane(c, a, kz);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::mult(CagmScalarFieldOps *c)
{
    return mult(c, this);
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::add(CagmScalarFieldOps *a, CagmScalarFieldOps *b)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        add_plane(a, b, kz);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::add(CagmScalarFieldOps *a)
{
    return add(this, a);
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::sub(CagmScalarFieldOps *a, CagmScalarFieldOps *b)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        sub_plane(a, b, kz);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::sub(CagmScalarFieldOps *a)
{
    return sub(this, a);
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::neg(CagmScalarFieldOps *a)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        neg_plane(a, kz);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::neg()
{
    return neg(this);
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::acos()
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        for (int ky = NL[1]; ky < NH[1]; ky++)
            for (int kx = NL[0]; kx < NH[0]; kx++)
				field[fidx(kx, ky, kz)] = ::acos(field[fidx(kx, ky, kz)]);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::zero()
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        for (int ky = NL[1]; ky < NH[1]; ky++)
            for (int kx = NL[0]; kx < NH[0]; kx++)
                field[fidx(kx, ky, kz)] = 0;

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::zeroZ0()
{
    for (int ky = NL[1]; ky < NH[1]; ky++)
        for (int kx = NL[0]; kx < NH[0]; kx++)
            field[fidx(kx, ky, 0)] = 0;

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::setZlevel(int level, double w)
{
    for (int ky = NL[1]; ky < NH[1]; ky++)
        for (int kx = NL[0]; kx < NH[0]; kx++)
            field[ky+level*N[1]][kx] = w;

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::setPlane(CagmScalarFieldOps *plane, int wplane, int from, int to)
{
	int kx, ky, kz;
    if (wplane & PLANE_Z)
    {
        for (ky = NL[1]; ky < NH[1]; ky++)
            for (kx = NL[0]; kx < NH[0]; kx++)
			    field[fidx(kx, ky, to)] = plane->field[ky+from*plane->N[1]][kx];
    }
    if (wplane & PLANE_Y)
    {
        for (kz = NL[2]; kz < NH[2]; kz++)
            for (kx = NL[0]; kx < NH[0]; kx++)
			    field[fidx(kx, to, kz)] = plane->field[from+kz*plane->N[1]][kx];
    }
    if (wplane & PLANE_X)
    {
        for (kz = NL[2]; kz < NH[2]; kz++)
            for (ky = NL[1]; ky < NH[1]; ky++)
			    field[fidx(to, ky, kz)] = plane->field[ky+kz*plane->N[1]][from];
    }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::LOS(CagmVectorFieldOps *a, double *dircos)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        for (int ky = NL[1]; ky < NH[1]; ky++)
            for (int kx = NL[0]; kx < NH[0]; kx++)
				field[fidx(kx, ky, kz)] =   a->fieldX[fidx(kx, ky, kz)]*dircos[0]
                                          + a->fieldY[fidx(kx, ky, kz)]*dircos[1]
                                          + a->fieldZ[fidx(kx, ky, kz)]*dircos[2]
                ;

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::rotate2D(CagmScalarFieldOps *ac, CagmScalarFieldOps *as, CagmScalarFieldOps *acn, CagmScalarFieldOps *asn, double cosz)
{
    double sinz = sqrt(1-cosz*cosz);
    int N[3];
    ac->dimensions(N);
    for (int kz = ac->NL[2]; kz < ac->NH[2]; kz++)
        for (int ky = ac->NL[1]; ky < ac->NH[1]; ky++)
            for (int kx = ac->NL[0]; kx < ac->NH[0]; kx++)
            {
				acn->field[fidx(kx, ky, kz)] =  ac->field[fidx(kx, ky, kz)]*cosz + as->field[fidx(kx, ky, kz)]*sinz;
				asn->field[fidx(kx, ky, kz)] = -ac->field[fidx(kx, ky, kz)]*sinz + as->field[fidx(kx, ky, kz)]*cosz;
            }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::sqDiff(CagmScalarFieldOps *a1, CagmScalarFieldOps *a2)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        for (int ky = NL[1]; ky < NH[1]; ky++)
            for (int kx = NL[0]; kx < NH[0]; kx++)
            {
                if ( ::fabs(a1->field[fidx(kx, ky, kz)]) <= ::fabs(a2->field[fidx(kx, ky, kz)]) )
				    field[fidx(kx, ky, kz)] = 0;
                else
				    field[fidx(kx, ky, kz)] = sqrt(a1->field[fidx(kx, ky, kz)]*a1->field[fidx(kx, ky, kz)] - a2->field[fidx(kx, ky, kz)]*a2->field[fidx(kx, ky, kz)]);
            }

    return 0;
}
//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::relax(CagmScalarFieldOps *cond, CagmScalarFieldOps *weight)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        for (int ky = NL[1]; ky < NH[1]; ky++)
            for (int kx = NL[0]; kx < NH[0]; kx++)
				field[fidx(kx, ky, kz)] = weight->field[fidx(kx, ky, kz)]*(cond->field[fidx(kx, ky, kz)] - field[fidx(kx, ky, kz)]) + field[fidx(kx, ky, kz)];

    return 0;
}

//-----------------------------------------------------------------------
double CagmScalarFieldOps::sum(CagmScalarFieldOps *weight)
{
    double wsum = 0;
    for (int kz = NL[2]; kz < NH[2]; kz++)
        wsum += sum_plane(kz, weight);

    return wsum;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::power(double pw)
{
    return power(this, pw);
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::power(CagmScalarFieldOps *a, double pw)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        for (int ky = NL[1]; ky < NH[1]; ky++)
            for (int kx = NL[0]; kx < NH[0]; kx++)
                field[fidx(kx, ky, kz)] = pow(a->field[fidx(kx, ky, kz)], pw);

    return 0;
}

//-----------------------------------------------------------------------
double CagmScalarFieldOps::maxval(void)
{
    double wmax = - DBL_MAX, wmax_p;
	int kz;
    for (kz = NL[2]; kz < NH[2]; kz++)
    {
        wmax_p = max_plane(kz);
        if (wmax_p > wmax)
            wmax = wmax_p;
    }

    return wmax;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::limWeight(int limType, CagmScalarFieldOps *calc, CagmScalarFieldOps *cond)
{
    if (limType != 0)
    {
        for (int kz = NL[2]; kz < NH[2]; kz++)
            for (int ky = NL[1]; ky < NH[1]; ky++)
                for (int kx = NL[0]; kx < NH[0]; kx++)
                    if (limType < 0 && calc->field[fidx(kx, ky, kz)] < cond->field[fidx(kx, ky, kz)])
                        field[fidx(kx, ky, kz)] = 0;
                    else if (limType > 0 && calc->field[fidx(kx, ky, kz)] > cond->field[fidx(kx, ky, kz)])
                        field[fidx(kx, ky, kz)] = 0;
    }

    return 0;
}
