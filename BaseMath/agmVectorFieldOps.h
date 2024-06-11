#pragma once

#include "math.h"
#include "agmRKF45.h"

#include "agm3DgridDefines.h"
#include "CubeXD.h"

class CagmScalarFieldOps;
class CagmRotate3D;

//------------------------------------------------------------------
class CagmVectorFieldOps : public CubeXD
{
friend class CagmScalarFieldOps;

public:
    typedef enum {None = 0, Boundary = 2, OutOfCube = 4, BufferOverload = 8, RKF45Problem = 16
                 }  Status;
    typedef enum {Linear = 0, Lanczos = 1
                 }  Interpolator;

    static double v_norm(double *v1, double *v2)
    {
        return sqrt(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
    }

    static uint32_t GetAllocSize(int *N)
    {
        return 3*sizeof(double)*N[1]*N[2] + sizeof(CagmVectorFieldOps);
    }

protected:
    double **fieldX, **fieldY, **fieldZ;

public:
	CagmVectorFieldOps(int *_N, double *_step = nullptr, int *_NphysL = nullptr, int *_NphysH = nullptr);
	CagmVectorFieldOps(CubeXD *);
	virtual ~CagmVectorFieldOps();

    double *getAddress(int v, int kx, int ky, int kz);

    uint32_t cross_plane(CagmVectorFieldOps *a, const CagmVectorFieldOps *b, int z);
    uint32_t rot_plane(CagmVectorFieldOps *a, int kz, int scheme = 3);
	uint32_t grad_plane(CagmScalarFieldOps *a, int kz, int scheme = 3);
    uint32_t mult_plane(CagmScalarFieldOps *c, CagmVectorFieldOps *a, int kz);
    uint32_t mult_plane(double c, CagmVectorFieldOps *a, int kz);
    uint32_t add_plane(CagmVectorFieldOps *a, CagmVectorFieldOps *b, int kz);
    uint32_t sub_plane(CagmVectorFieldOps *a, CagmVectorFieldOps *b, int kz);
    uint32_t neg_plane(CagmVectorFieldOps *a, int kz);
    double max2_plane(int kz);

	uint32_t cross(CagmVectorFieldOps *a, const CagmVectorFieldOps *b);
    uint32_t cross(CagmVectorFieldOps *a);
    uint32_t rot(CagmVectorFieldOps *a);
    uint32_t grad(CagmScalarFieldOps *a);
	uint32_t mult(double c, CagmVectorFieldOps *a);
    uint32_t mult(double c);
	uint32_t mult(CagmScalarFieldOps *c, CagmVectorFieldOps *a);
	uint32_t mult(CagmVectorFieldOps *a, CagmScalarFieldOps *c);
	uint32_t mult(CagmScalarFieldOps *c);
	uint32_t mult(CagmVectorFieldOps *a, double *d);
    uint32_t add(CagmVectorFieldOps *a, CagmVectorFieldOps *b);
    uint32_t add(CagmVectorFieldOps *a);
	uint32_t sub(CagmVectorFieldOps *a, CagmVectorFieldOps *b);
    uint32_t sub(CagmVectorFieldOps *a);
    uint32_t neg(CagmVectorFieldOps *a);
    uint32_t neg();
    uint32_t zero();

    uint32_t shift(int n);
    uint32_t setZlevel(int wplane, int level, double w = 1.0);
    uint32_t setVector(double *w);

    uint32_t getPlane(CagmVectorFieldOps *plane, int wplane, int from, int to);
    uint32_t setPlaneComp(CagmVectorFieldOps *plane, int wplane, int wcomp, int from, int to);
    uint32_t setPlane(CagmVectorFieldOps *plane, int wplane, int from, int to);

    uint32_t getComponent(CagmScalarFieldOps *comp, int wcomp);
    uint32_t setComponent(CagmScalarFieldOps *comp, int wcomp);
    uint32_t getTransv(CagmScalarFieldOps *a);

    uint32_t getBounds(CagmVectorFieldOps *boundsx, CagmVectorFieldOps *boundsy, CagmVectorFieldOps *boundsz);
    uint32_t setBounds(CagmVectorFieldOps *boundsx, CagmVectorFieldOps *boundsy, CagmVectorFieldOps *boundsz);

	uint32_t stretch(CagmVectorFieldOps *a, Interpolator = Linear, double p1 = 0, double p2 = 0, double p3 = 0);
    uint32_t conv(CagmVectorFieldOps *src, CagmScalarFieldOps *win);
    uint32_t inCube(const double * coord, const double absBoundAchieve = 0, const double relBoundAchieve = 0);
    uint32_t getPoint(const double *coord, double *vect);

    Status getOneFullLine(CagmRKF45 *rkf45, double *start, int direction, double step, double boundAchieve, double boundAchieveBottom,
        int maxResult, int *length, double *coord, int *status);

    uint32_t rotate3D(CagmRotate3D *, bool);

protected:
	uint32_t Initialize();
    Status getOneLine(CagmRKF45 *rkf45, CagmRKF45Vect *rkfv, double step, double *coord, int maxlen, int *length, CagmRKF45::Status *status, bool noDuplicate = false);
    uint32_t Delete();
};
