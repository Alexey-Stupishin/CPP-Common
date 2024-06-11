#include "stdDefinitions.h"
#include <math.h>

#include "agsFieldsCommon.h"
#include "agmVectorFieldOps.h"
#include "agmScalarFieldOps.h"

#include "agmRotate3D.h"

//-----------------------------------------------------------------------
CagmVectorFieldOps::CagmVectorFieldOps(int *_N, double *_step, int *_NL, int *_NH)
    : CubeXD(_N, 3, _step, _NL, _NH)
{
    Initialize();
}

//-----------------------------------------------------------------------
CagmVectorFieldOps::CagmVectorFieldOps(CubeXD *from)
    : CubeXD(from)
{
    Initialize();
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::Initialize()
{
    fieldX = new double * [N[1]*N[2]];
    fieldY = new double * [N[1]*N[2]];
    fieldZ = new double * [N[1]*N[2]];

    return 0;
}

//-----------------------------------------------------------------------
CagmVectorFieldOps::~CagmVectorFieldOps()
{
    Delete();
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::Delete()
{
    delete [] fieldX;
    delete [] fieldY;
    delete [] fieldZ;

    return 0;
}

//-----------------------------------------------------------------------
double *CagmVectorFieldOps::getAddress(int v, int kx, int ky, int kz)
{
    if (v == 0)
        return &fieldX[fidx(kx, ky, kz)];
    else if (v == 1)
        return &fieldY[fidx(kx, ky, kz)];
    else
        return &fieldZ[fidx(kx, ky, kz)];
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::cross(CagmVectorFieldOps *a, const CagmVectorFieldOps *b)
{
	// check equiv. sizes!
	int kz;
    for (kz = 0; kz < N[2]; kz++)
        cross_plane(a, b, kz);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::cross(CagmVectorFieldOps *a)
{
    return cross(this, a);
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::rot(CagmVectorFieldOps *a)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        rot_plane(a, kz);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::grad(CagmScalarFieldOps *a)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        grad_plane(a, kz);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::mult(double c, CagmVectorFieldOps *a)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        mult_plane(c, a, kz);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::mult(double c)
{
    return mult(c, this);
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::mult(CagmScalarFieldOps *c, CagmVectorFieldOps *a)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        mult_plane(c, a, kz);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::mult(CagmVectorFieldOps *a, CagmScalarFieldOps *c)
{
    return mult(c, a);
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::mult(CagmScalarFieldOps *c)
{
    return mult(c, this);
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::mult(CagmVectorFieldOps *a, double *d)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        for (int ky = NL[1]; ky < NH[1]; ky++)
            for (int kx = NL[0]; kx < NH[0]; kx++)
			{
				fieldX[fidx(kx, ky, kz)] = (a->fieldX[fidx(kx, ky, kz)]) * d[0];
				fieldY[fidx(kx, ky, kz)] = (a->fieldY[fidx(kx, ky, kz)]) * d[1];
				fieldZ[fidx(kx, ky, kz)] = (a->fieldZ[fidx(kx, ky, kz)]) * d[2];
			}

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::add(CagmVectorFieldOps *a, CagmVectorFieldOps *b)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        add_plane(a, b, kz);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::add(CagmVectorFieldOps *a)
{
    return add(this, a);
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::sub(CagmVectorFieldOps *a, CagmVectorFieldOps *b)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        sub_plane(a, b, kz);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::sub(CagmVectorFieldOps *a)
{
    return sub(this, a);
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::neg(CagmVectorFieldOps *a)
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        neg_plane(a, kz);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::neg()
{
    return neg(this);
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::zero()
{
    for (int kz = NL[2]; kz < NH[2]; kz++)
        for (int ky = NL[1]; ky < NH[1]; ky++)
            for (int kx = NL[0]; kx < NH[0]; kx++)
			{
				fieldX[fidx(kx, ky, kz)] = 0;
				fieldY[fidx(kx, ky, kz)] = 0;
				fieldZ[fidx(kx, ky, kz)] = 0;
			}

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::shift(int n)
{
    int from0 = 0;
    int from = from0;
    int to0 = N[2];
    int to = to0;
    if (n > 0)
        from += n;
    else
        to += n;

	int kx, ky, kz;
    for (kz = from; kz < to; kz++)
        for (ky = 0; ky < N[1]; ky++)
            for (kx = 0; kx < N[0]; kx++)
			{
				fieldX[fidx(kx, ky, kz)] = fieldX[fidx(kx, ky, kz-n)];
				fieldY[fidx(kx, ky, kz)] = fieldY[fidx(kx, ky, kz-n)];
				fieldZ[fidx(kx, ky, kz)] = fieldY[fidx(kx, ky, kz-n)];
			}
    for (kz = from0; kz < from; kz++)
        for (ky = 0; ky < N[1]; ky++)
            for (kx = 0; kx < N[0]; kx++)
			{
				fieldX[fidx(kx, ky, kz)] = 0;
				fieldY[fidx(kx, ky, kz)] = 0;
				fieldZ[fidx(kx, ky, kz)] = 0;
			}
    for (kz = to+1; kz < to0; kz++)
        for (ky = 0; ky < N[1]; ky++)
            for (kx = 0; kx < N[0]; kx++)
			{
				fieldX[fidx(kx, ky, kz)] = 0;
				fieldY[fidx(kx, ky, kz)] = 0;
				fieldZ[fidx(kx, ky, kz)] = 0;
			}

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::setZlevel(int wplane, int level, double w)
{
	int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
		{
            if (wplane & PLANE_X)
			    fieldX[ky+level*N[1]][kx] = w;
            if (wplane & PLANE_Y)
			    fieldY[ky+level*N[1]][kx] = w;
            if (wplane & PLANE_Z)
			    fieldZ[ky+level*N[1]][kx] = w;
		}

    return 0;
}


//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::setVector(double *d)
{
	int kx, ky, kz;
    for (kz = 0; kz < N[2]; kz++)
        for (ky = 0; ky < N[1]; ky++)
            for (kx = 0; kx < N[0]; kx++)
			{
				fieldX[fidx(kx, ky, kz)] = d[0];
				fieldY[fidx(kx, ky, kz)] = d[1];
				fieldZ[fidx(kx, ky, kz)] = d[2];
			}

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::getPlane(CagmVectorFieldOps *plane, int wplane, int from, int to)
{
	int kx, ky, kz;
    if (wplane & PLANE_Z)
    {
        for (ky = NL[1]; ky < NH[1]; ky++)
            for (kx = NL[0]; kx < NH[0]; kx++)
		    {
			    plane->fieldX[ky+to*plane->N[1]][kx] = fieldX[fidx(kx, ky, from)];
			    plane->fieldY[ky+to*plane->N[1]][kx] = fieldY[fidx(kx, ky, from)];
			    plane->fieldZ[ky+to*plane->N[1]][kx] = fieldZ[fidx(kx, ky, from)];
		    }
    }
    if (wplane & PLANE_Y)
    {
        for (kz = NL[2]; kz < NH[2]; kz++)
            for (kx = NL[0]; kx < NH[0]; kx++)
		    {
			    plane->fieldX[to+kz*plane->N[1]][kx] = fieldX[fidx(kx, from, kz)];
			    plane->fieldY[to+kz*plane->N[1]][kx] = fieldY[fidx(kx, from, kz)];
			    plane->fieldZ[to+kz*plane->N[1]][kx] = fieldZ[fidx(kx, from, kz)];
		    }
    }
    if (wplane & PLANE_X)
    {
        for (kz = NL[2]; kz < NH[2]; kz++)
            for (ky = NL[1]; ky < NH[1]; ky++)
		    {
			    plane->fieldX[ky+kz*plane->N[1]][to] = fieldX[fidx(from, ky, kz)];
			    plane->fieldY[ky+kz*plane->N[1]][to] = fieldY[fidx(from, ky, kz)];
			    plane->fieldZ[ky+kz*plane->N[1]][to] = fieldZ[fidx(from, ky, kz)];
		    }
    }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::setPlaneComp(CagmVectorFieldOps *plane, int wplane, int wcomp, int from, int to)
{
	int kx, ky, kz;
    if (wplane & PLANE_Z)
    {
        for (ky = NL[1]; ky < NH[1]; ky++)
            for (kx = NL[0]; kx < NH[0]; kx++)
		    {
                if (wcomp & PLANE_X)
			        fieldX[fidx(kx, ky, to)] = plane->fieldX[ky+from*plane->N[1]][kx];
                if (wcomp & PLANE_Y)
			        fieldY[fidx(kx, ky, to)] = plane->fieldY[ky+from*plane->N[1]][kx];
                if (wcomp & PLANE_Z)
    			    fieldZ[fidx(kx, ky, to)] = plane->fieldZ[ky+from*plane->N[1]][kx];
		    }
    }
    if (wplane & PLANE_Y)
    {
        for (kz = NL[2]; kz < NH[2]; kz++)
            for (kx = NL[0]; kx < NH[0]; kx++)
		    {
                if (wcomp & PLANE_X)
			        fieldX[fidx(kx, to, kz)] = plane->fieldX[from+kz*plane->N[1]][kx];
                if (wcomp & PLANE_Y)
			        fieldY[fidx(kx, to, kz)] = plane->fieldY[from+kz*plane->N[1]][kx];
                if (wcomp & PLANE_Z)
			        fieldZ[fidx(kx, to, kz)] = plane->fieldZ[from+kz*plane->N[1]][kx];
		    }
    }
    if (wplane & PLANE_X)
    {
        for (kz = NL[2]; kz < NH[2]; kz++)
            for (ky = NL[1]; ky < NH[1]; ky++)
		    {
                if (wcomp & PLANE_X)
			        fieldX[fidx(to, ky, kz)] = plane->fieldX[ky+kz*plane->N[1]][from];
                if (wcomp & PLANE_Y)
			        fieldY[fidx(to, ky, kz)] = plane->fieldY[ky+kz*plane->N[1]][from];
                if (wcomp & PLANE_Z)
			        fieldZ[fidx(to, ky, kz)] = plane->fieldZ[ky+kz*plane->N[1]][from];
		    }
    }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::setPlane(CagmVectorFieldOps *plane, int wplane, int from, int to)
{
    return setPlaneComp(plane, wplane, PLANE_X + PLANE_Y + PLANE_Z, from, to);
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::getComponent(CagmScalarFieldOps *comp, int wcomp)
{
	int kx, ky, kz;
    for (kz = 0; kz < N[2]; kz++)
        for (ky = 0; ky < N[1]; ky++)
            for (kx = 0; kx < N[0]; kx++)
			{
                if (wcomp == PLANE_X)
			        comp->field[fidx(kx, ky, kz)] = fieldX[fidx(kx, ky, kz)];
                else if (wcomp == PLANE_Y)
			        comp->field[fidx(kx, ky, kz)] = fieldY[fidx(kx, ky, kz)];
                else if (wcomp == PLANE_Z)
			        comp->field[fidx(kx, ky, kz)] = fieldZ[fidx(kx, ky, kz)];
			}

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::setComponent(CagmScalarFieldOps *comp, int wcomp)
{
	int kx, ky, kz;
    for (kz = 0; kz < N[2]; kz++)
        for (ky = 0; ky < N[1]; ky++)
            for (kx = 0; kx < N[0]; kx++)
			{
                if (wcomp == PLANE_X)
			        fieldX[fidx(kx, ky, kz)] = comp->field[fidx(kx, ky, kz)];
                else if (wcomp == PLANE_Y)
			        fieldY[fidx(kx, ky, kz)] = comp->field[fidx(kx, ky, kz)];
                else if (wcomp == PLANE_Z)
			        fieldZ[fidx(kx, ky, kz)] = comp->field[fidx(kx, ky, kz)];
			}

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::getTransv(CagmScalarFieldOps *a)
{
	int kx, ky, kz;
    for (kz = 0; kz < N[2]; kz++)
        for (ky = 0; ky < N[1]; ky++)
            for (kx = 0; kx < N[0]; kx++)
				a->field[fidx(kx, ky, kz)] = sqrt(fieldX[fidx(kx, ky, kz)]*fieldX[fidx(kx, ky, kz)] + fieldY[fidx(kx, ky, kz)]*fieldY[fidx(kx, ky, kz)]);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::getBounds(CagmVectorFieldOps *boundsx, CagmVectorFieldOps *boundsy, CagmVectorFieldOps *boundsz)
{
	// check equiv. sizes!
    getPlane(boundsx, PLANE_X, NL[0], 0);
    getPlane(boundsx, PLANE_X, NH[0]-1, 1);
    getPlane(boundsy, PLANE_Y, NL[1], 0);
    getPlane(boundsy, PLANE_Y, NH[1]-1, 1);
    getPlane(boundsz, PLANE_Z, NL[2], 0);
    getPlane(boundsz, PLANE_Z, NH[2]-1, 1);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::setBounds(CagmVectorFieldOps *boundsx, CagmVectorFieldOps *boundsy, CagmVectorFieldOps *boundsz)
{
	// check equiv. sizes!
    setPlane(boundsx, PLANE_X, 0, NL[0]);
    setPlane(boundsx, PLANE_X, 1, NH[0]-1);
    setPlane(boundsy, PLANE_Y, 0, NL[1]);
    setPlane(boundsy, PLANE_Y, 1, NH[1]-1);
    setPlane(boundsz, PLANE_Z, 0, NL[2]);
    setPlane(boundsz, PLANE_Z, 1, NH[2]-1);

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::stretch(CagmVectorFieldOps *src, Interpolator /* inetrp */, double /* p1 */, double /* p2 */, double /* p3 */)
{
    double cx = (double)(src->N[0] - 1)/(double)(N[0] - 1);
    double cy = (double)(src->N[1] - 1)/(double)(N[1] - 1);
    double cz = (double)(src->N[2] - 1)/(double)(N[2] - 1);

    //if (inetrp == Lanczos)
    //{
    //    int win = (int)p1;
    //    int size = (int)p2;
    //}

	int kx, ky, kz;
	int x1, y1, z1;
    double tx, ty, tz;
    for (kz = NL[2]; kz < NH[2]; kz++)
    {
        l_div_position((kz*cz), (src->N[2]), z1, tz);
        for (ky = NL[1]; ky < NH[1]; ky++)
        {
            l_div_position(ky*cy, src->N[1], y1, ty);
            for (kx = NL[0]; kx < NH[0]; kx++)
			{
                l_div_position(kx*cx, src->N[0], x1, tx);
                fieldX[fidx(kx, ky, kz)] = l_get_point(src->N[1], src->fieldX, x1, y1, z1, tx, ty, tz);
                fieldY[fidx(kx, ky, kz)] = l_get_point(src->N[1], src->fieldY, x1, y1, z1, tx, ty, tz);
                fieldZ[fidx(kx, ky, kz)] = l_get_point(src->N[1], src->fieldZ, x1, y1, z1, tx, ty, tz);
            }
        }
    }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::conv(CagmVectorFieldOps *src, CagmScalarFieldOps *win)
{
    int kx, ky, kz;
    int wx, wy, wz;
    int tx;
    int cx = (int)floor((win->NL[0] + win->NH[0]) / 2);
    int cy = (int)floor((win->NL[1] + win->NH[1]) / 2);
    int cz = (int)floor((win->NL[2] + win->NH[2]) / 2);
    //int M = (win->NH[0] - win->NL[0])*(win->NH[1] - win->NL[1])*(win->NH[2] - win->NL[2]);
    for (kz = NL[2]+cz; kz < NH[2]-cz; kz++)
        for (ky = NL[1]+cy; ky < NH[1]-cy; ky++)
            for (kx = NL[0]+cx; kx < NH[0]-cx; kx++)
            {
                double sx = 0, sy = 0, sz = 0;
                for (wz = win->NL[2]; wz < win->NH[2]; wz++)
                {
                    for (wy = win->NL[1]; wy < win->NH[1]; wy++)
                    {
                        for (wx = win->NL[0]; wx < win->NH[0]; wx++)
                        {
                            tx = kx + wx - cx;
                            sx += src->fieldX[fidx(kx, ky, kz)] * win->field[wy + wz*win->N[1]][wx];
                            sy += src->fieldY[fidx(kx, ky, kz)] * win->field[wy + wz*win->N[1]][wx];
                            sz += src->fieldZ[fidx(kx, ky, kz)] * win->field[wy + wz*win->N[1]][wx];
                        }
                    }
                }
                fieldX[fidx(kx, ky, kz)] = sx; // / M;
                fieldY[fidx(kx, ky, kz)] = sy; // / M;
                fieldZ[fidx(kx, ky, kz)] = sz; // / M;
            }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::inCube(const double *coord, const double absBoundAchieve, const double relBoundAchieve)
{
    uint32_t out = 0;
    if (coord[0] < relBoundAchieve*N[0] || coord[0] < absBoundAchieve || 
        coord[0] > N[0]-relBoundAchieve*N[0] - 1 || coord[0] > N[0] - 1 - absBoundAchieve)
        out = out | PLANE_X;
    if (coord[1] < relBoundAchieve*N[1] || coord[1] < absBoundAchieve || 
        coord[1] > N[1]-relBoundAchieve*N[1] - 1 || coord[1] > N[1] - 1 - absBoundAchieve)
        out = out | PLANE_Y;
    if (coord[2] < relBoundAchieve*N[2] || coord[2] < absBoundAchieve || 
        coord[2] > N[2]-relBoundAchieve*N[2] - 1 || coord[2] > N[2] - 1 - absBoundAchieve)
        out = out | PLANE_Z;

    return out;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::getPoint(const double *coord, double *vect)
{
    uint32_t out = 0;
    if (coord[0] < 0 || coord[0] > N[0]-1)
        out = out | PLANE_X;
    if (coord[1] < 0 || coord[1] > N[1]-1)
        out = out | PLANE_Y;
    if (coord[2] < 0 || coord[2] > N[2]-1)
        out = out | PLANE_Z;

    if (!out)
    {
        int x1, y1, z1;
        double tx, ty, tz;
        l_div_position(coord[0], N[0], x1, tx);
        l_div_position(coord[1], N[1], y1, ty);
        l_div_position(coord[2], N[2], z1, tz);
        vect[0] = l_get_point(N[1], fieldX, x1, y1, z1, tx, ty, tz);
        vect[1] = l_get_point(N[1], fieldY, x1, y1, z1, tx, ty, tz);
        vect[2] = l_get_point(N[1], fieldZ, x1, y1, z1, tx, ty, tz);
    }

    return out;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::rotate3D(CagmRotate3D *rotator, bool direction)
{
    int kx, ky, kz;
    double v[3];
    for (kz = 0; kz < N[2]; kz++)
        for (ky = 0; ky < N[1]; ky++)
            for (kx = 0; kx < N[0]; kx++)
            {
                rotator->rotate(fieldX[fidx(kx, ky, kz)], fieldY[fidx(kx, ky, kz)], fieldZ[fidx(kx, ky, kz)], v, v+1, v+2, direction);
                fieldX[fidx(kx, ky, kz)] = v[0];
                fieldY[fidx(kx, ky, kz)] = v[1];
                fieldZ[fidx(kx, ky, kz)] = v[2];
            }

    return 0;
}
