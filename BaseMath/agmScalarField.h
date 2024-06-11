#pragma once
#include "stdDefinitions.h"

#include "agmScalarFieldOps.h"
#include "agsFieldsCommon.h"
#include "binUtilities.h"
#include "agmBaseMath.h"

class CagmVectorField;

enum {SWF_NONE = 0, SWF_COS = 1, SWF_UDF = 99};

//------------------------------------------------------------------
class CagmScalarField : public CagmScalarFieldOps
{
friend class CagmVectorField;

public:
    static uint32_t GetAllocSize(int *N)
    {
        return sizeof(double)*N[0]*N[1]*N[2] + sizeof(CagmScalarField) + CagmScalarFieldOps::GetAllocSize(N);
    }

public:
    double *allocField;

protected:
	bool isRef;

protected:
    uint32_t Alloc()
    {
        allocField = new double[N[0]*N[1]*N[2]];
        memset(allocField, 0, sizeof(double)*N[0]*N[1]*N[2]);

        int ky, kz;
        for (ky = 0; ky < N[1]; ky++)
            for (kz = 0; kz < N[2]; kz++)
                field[ky + kz*N[1]] = allocField + (ky + kz*N[1])*N[0];

        return 0;
    }

    uint32_t Delete()
    {
		if (!isRef)
            delete [] allocField;

        return 0;
    }

public:
	CagmScalarField(int *_N, double *_step = nullptr, int *_NL = nullptr, int *_NH = nullptr)
		: CagmScalarFieldOps(_N, _step, _NL, _NH)
        , allocField(nullptr)
        , isRef(false)
    {
        Alloc();
	}

	CagmScalarField(CbinDataStruct* data, int idx = 0)
        : CagmScalarFieldOps(data->GetDimensions())
        , allocField(nullptr)  
        , isRef(false)
    {
        Alloc();
        data->Copy(allocField, nullptr, nullptr, idx);
	}

    CagmScalarField(CubeXD *_from, bool copy = false)
        : CagmScalarFieldOps(_from)
        , allocField(nullptr)  
        , isRef(false)
    {
        dim_ = 1;
        Alloc();

        if (!copy)
            return;

        CagmScalarField *from = (CagmScalarField *)_from;

        memcpy(allocField, from->allocField, sizeof(double)*N[0]*N[1]*N[2]);
    }

    uint32_t setField(double *Bc)
    {
        memcpy(allocField, Bc, N[0]*N[1]*N[2]*sizeof(double));

        return 0;
    }

    uint32_t Copy2Bottom(CagmScalarField *source)
    {
        for (int kz = source->NL[2]; kz < source->NH[2]; kz++)
            for (int ky = source->NL[1]; ky < source->NH[1]; ky++)
                for (int kx = source->NL[0]; kx < source->NH[0]; kx++)
                    allocField[(ky + kz*N[1])*N[0] + kx] = source->allocField[(ky + kz*source->N[1])*source->N[0] + kx];
        
        return 0;
    }

	virtual ~CagmScalarField()
    {
        Delete();
    }

    uint32_t GetFieldAddress(double **p)
    {
        *p = allocField;

        return 0;
    }

    uint32_t GetData(CbinDataStruct* data)
    {
        return data->Create(N, allocField);
    }

    uint32_t Weight(int type, double bound, int *xb, int *yb, int *zb, double rel_bound_weight = 1.0)
    {
        int kx, ky, kz;
        double *wx = new double[N[0]];
        double *wy = new double[N[1]];
        double *wz = new double[N[2]];

        double zD = (NH[2]-NL[2])*bound;
        zb[0] = 0; zb[1] = 0;
        for (kz = NL[2]; kz < NH[2]; kz++)
        {
            if (type == SWF_COS && kz > NH[2]-zD-1)
                wz[kz] = cos(v_pi_c*0.5* rel_bound_weight *(kz-NH[2]+zD+1)/zD);
            else
            {
                wz[kz] = v_1;
                if (kz > zb[1])
                    zb[1] = kz;
            }
        }

        double yD = (NH[1]-NL[1])*bound;
        yb[0] = N[1]; yb[1] = 0;
        for (ky = NL[1]; ky < NH[1]; ky++)
        {
            if (type == SWF_COS && ky > NH[1]-yD-1)
                wy[ky] = cos(v_pi_c*0.5* rel_bound_weight *(ky-NH[1]+1+yD)/yD);
            else if (type == SWF_COS && ky < NL[1]+yD)
                wy[ky] = cos(v_pi_c*0.5* rel_bound_weight *(v_1 - (ky+NL[1])/yD));
            else
            {
                wy[ky] = v_1;
                if (ky < yb[0])
                    yb[0] = ky;
                if (ky > yb[1])
                    yb[1] = ky;
            }
        }

        double xD = (NH[0]-NL[0])*bound;
        xb[0] = N[0]; xb[1] = 0;
        for (kx = NL[0]; kx < NH[0]; kx++)
        {
            if (type == SWF_COS && kx > NH[0]-xD-1)
                wx[kx] = cos(v_pi_c*0.5* rel_bound_weight *(kx-NH[0]+1+xD)/xD);
            else if (type == SWF_COS && kx < NL[0]+xD)
                wx[kx] = cos(v_pi_c*0.5* rel_bound_weight *(v_1 - (kx+NL[0])/xD));
            else
            {
                wx[kx] = v_1;
                if (kx < xb[0])
                    xb[0] = kx;
                if (kx > xb[1])
                    xb[1] = kx;
            }
        }

        for (kz = NL[2]; kz < NH[2]; kz++)
            for (ky = NL[1]; ky < NH[1]; ky++)
                for (kx = NL[0]; kx < NH[0]; kx++)
                    allocField[kx + (ky + kz*N[1])*N[0]] = wx[kx]*wy[ky]*wz[kz];

        delete [] wx;
        delete [] wy;
        delete [] wz;

        return 0;
    }
};
