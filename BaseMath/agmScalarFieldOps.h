#pragma once

#include "agm3DgridDefines.h"
#include "CubeXD.h"

class CagmVectorFieldOps;

//------------------------------------------------------------------
class CagmScalarFieldOps : public CubeXD
{
friend class CagmVectorFieldOps;

public:
    static uint32_t GetAllocSize(int *N)
    {
        return sizeof(double)*N[1]*N[2] + sizeof(CagmScalarFieldOps);
    }

protected:
    double **field;
    double tolerance_zero;
    double tolerance_denom;

public:
    CagmScalarFieldOps(int *_N, double *_step = nullptr, int *_NL = nullptr, int *_NH = nullptr);
	CagmScalarFieldOps(CubeXD *);
	virtual ~CagmScalarFieldOps();

    double *getAddress(int kx, int ky, int kz);

    uint32_t div_plane(CagmVectorFieldOps *a, int kz, int scheme = 3);
    uint32_t dot_plane(CagmVectorFieldOps *a, CagmVectorFieldOps *b, int kz, CagmScalarFieldOps *Weight = nullptr);
    uint32_t abs2_plane(CagmVectorFieldOps *a, int kz, CagmScalarFieldOps *Weight = nullptr);
    uint32_t abs_plane(CagmVectorFieldOps *a, int kz, CagmScalarFieldOps *Weight = nullptr);
    uint32_t sqrt_plane(CagmScalarFieldOps *a, int kz);
    uint32_t inv_plane(CagmScalarFieldOps *a, int kz);
    uint32_t invabs_plane(CagmVectorFieldOps *a, int kz, CagmScalarFieldOps *Weight = nullptr);
    uint32_t mult_plane(CagmScalarFieldOps *b, CagmScalarFieldOps *a, int kz);
    uint32_t mult_plane(double c, CagmScalarFieldOps *a, int kz);
    uint32_t add_plane(CagmScalarFieldOps *a, CagmScalarFieldOps *b, int kz);
    uint32_t sub_plane(CagmScalarFieldOps *a, CagmScalarFieldOps *b, int kz);
    uint32_t neg_plane(CagmScalarFieldOps *a, int kz);
    double sum_plane(int kz, CagmScalarFieldOps *weight = nullptr);
    double max_plane(int kz);

    uint32_t stretch(CagmScalarFieldOps*src);

    uint32_t div(CagmVectorFieldOps *a);
    uint32_t dot(CagmVectorFieldOps *a, CagmVectorFieldOps *b, CagmScalarFieldOps *Weight = nullptr);
    uint32_t abs2(CagmVectorFieldOps *a, CagmScalarFieldOps *Weight = nullptr);
    uint32_t abs(CagmVectorFieldOps *a);
	uint32_t projection(CagmVectorFieldOps *a, double *d);
	uint32_t inv(CagmScalarFieldOps *a);
    uint32_t inv(void);
	uint32_t mult(double c, CagmScalarFieldOps *a);
    uint32_t mult(double c);
	uint32_t mult(CagmScalarFieldOps *c, CagmScalarFieldOps *a);
    uint32_t mult(CagmScalarFieldOps *c);
    uint32_t add(CagmScalarFieldOps *a, CagmScalarFieldOps *b);
    uint32_t add(CagmScalarFieldOps *a);
    uint32_t sub(CagmScalarFieldOps *a, CagmScalarFieldOps *b);
    uint32_t sub(CagmScalarFieldOps *a);
    uint32_t neg(CagmScalarFieldOps *a);
    uint32_t neg();
    uint32_t acos();
    uint32_t power(CagmScalarFieldOps *a, double pw);
    uint32_t power(double pw);
    uint32_t zero();
    uint32_t zeroZ0();
    uint32_t setZlevel(int level, double w);
    uint32_t setPlane(CagmScalarFieldOps *plane, int wplane, int from, int to);

	uint32_t LOS(CagmVectorFieldOps *a, double *dircos);
    static uint32_t rotate2D(CagmScalarFieldOps *ac, CagmScalarFieldOps *as, CagmScalarFieldOps *acn, CagmScalarFieldOps *asn, double cosz);
    uint32_t sqDiff(CagmScalarFieldOps *a1, CagmScalarFieldOps *a2);
    uint32_t relax(CagmScalarFieldOps *mult, CagmScalarFieldOps *weight);

    double sum(CagmScalarFieldOps *weight = nullptr);
    double maxval(void);
    uint32_t limWeight(int limType, CagmScalarFieldOps *calc, CagmScalarFieldOps *cond);

protected:
	uint32_t Initialize();
    uint32_t Delete();
};

// static void indices(int p, int Np, int q, int Nq, int *p1, int *p2, double *pf, int *q1, int *q2, double *qf);

