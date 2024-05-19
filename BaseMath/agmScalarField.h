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

    uint32_t Copy(const CagmScalarField& from)
    {
        if (from.isRef)
        {
            int ky, kz;
            for (kz = 0; kz < N[2]; kz++)
                for (ky = 0; ky < N[1]; ky++)
                    memcpy(allocField + (ky + kz*N[1])*N[0], from.allocField + (ky + kz*N[1])*N[0], sizeof(double)*N[0]);
        }
        else
            memcpy(allocField, from.allocField, sizeof(double)*N[0]*N[1]*N[2]);

        SetSteps((double *)from.step);

        return 0;
    }

    uint32_t Copy(CbinDataStruct* data, int idx = 0)
    {
        data->Copy(allocField, nullptr, nullptr, idx);

        double step[] = {1.0, 1.0, 1.0};
        SetSteps(step);

        return 0;
    }

    uint32_t Copy(double *S) // use carefully!
    {
        memcpy(allocField, S, sizeof(double)*N[0]*N[1]*N[2]);

        double step[] = {1.0, 1.0, 1.0};
        SetSteps(step);

        return 0;
    }

    uint32_t Delete()
    {
		if (!isRef)
            delete [] allocField;

        return 0;
    }

public:
	CagmScalarField(int *_N, bool isAlloc = true, int *_DphysL = nullptr, int *_DphysH = nullptr)
		: CagmScalarFieldOps(_N, _DphysL, _DphysH),
          allocField(nullptr)  
    {
		isRef = !isAlloc;
        if (isAlloc)
            Alloc();
	}

	CagmScalarField(CbinDataStruct* data, int idx = 0)
        : CagmScalarFieldOps(data->GetDimensions()),
          allocField(nullptr)  
        {
    		isRef = false;
            Alloc();
            Copy(data, idx);
		}

	CagmScalarField(int *N, double *S)
        : CagmScalarFieldOps(N),
          allocField(nullptr)  
        {
    		isRef = false;
            Alloc();
            Copy(S);
		}

    CagmScalarField(const CagmScalarField& from) // copy constructor, creates only solid copy
        : CagmScalarFieldOps(from),
          allocField(nullptr)  
    {
        isRef = false;
        Alloc();
        Copy(from);
    }

	CagmScalarField(CagmScalarFieldOps *source, int *M, int *Mmin, int *_DphysL = nullptr, int *_DphysH = nullptr)
		: CagmScalarFieldOps(M, _DphysL, _DphysH),
          allocField(nullptr)  
    {
		isRef = true;
        SetMargins(source, Mmin, _DphysL, _DphysH);
        SetSteps(source->step);
	}

    CagmScalarField& operator=(const CagmScalarField& from) // creates only solid copy
    {
        if (this == &from)
            return *this;

        this->CagmScalarFieldOps::operator=(from);

        Delete();
        Alloc();
        Copy(from);

        return *this;
    }

    uint32_t CreateConvWindow();

    uint32_t setField(double *Bc)
    {
        memcpy(allocField, Bc, N[0]*N[1]*N[2]*sizeof(double));

        return 0;
    }

    uint32_t Copy(CagmScalarField *source, int *Mmin)
    {
        // M == N !
	    int kx, ky, kz;
        for (kz = 0; kz < N[2]; kz++)
            for (ky = 0; ky < N[1]; ky++)
                for (kx = 0; kx < N[0]; kx++)
                    allocField[(ky + kz*N[1])*N[0] + kx] = source->allocField[(ky+Mmin[1] + (kz+Mmin[2])*source->N[1])*source->N[0] + kx+Mmin[0]];

        SetSteps(source->step);
        
        return 0;
    }

    uint32_t Copy2Bottom(CagmScalarField *source)
    {
        // N[0] = source->N[0]; N[1] = source->N[1]; 
	    int kx, ky, kz;
        for (kz = 0; kz < source->N[2]; kz++)
            for (ky = 0; ky < source->N[1]; ky++)
                for (kx = 0; kx < source->N[0]; kx++)
                    allocField[(ky + kz*N[1])*N[0] + kx] = source->allocField[(ky + kz*source->N[1])*source->N[0] + kx];

        SetSteps(source->step);
        
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

        double zD = N[2]*bound;
        zb[0] = 0; zb[1] = 0;
        for (kz = 0; kz < N[2]; kz++)
        {
            if (type == SWF_COS && kz > N[2]-zD-1)
                wz[kz] = cos(v_pi_c*0.5* rel_bound_weight *(kz-N[2]+zD+1)/zD);
            else
            {
                wz[kz] = v_1;
                if (kz > zb[1])
                    zb[1] = kz;
            }
        }

        double yD = N[1]*bound;
        yb[0] = N[1]; yb[1] = 0;
        for (ky = 0; ky < N[1]; ky++)
        {
            if (type == SWF_COS && ky > N[1]-yD-1)
                wy[ky] = cos(v_pi_c*0.5* rel_bound_weight *(ky-N[1]+1+yD)/yD);
            else if (type == SWF_COS && ky < yD)
                wy[ky] = cos(v_pi_c*0.5* rel_bound_weight *(v_1 - ky/yD));
            else
            {
                wy[ky] = v_1;
                if (ky < yb[0])
                    yb[0] = ky;
                if (ky > yb[1])
                    yb[1] = ky;
            }
        }

        double xD = N[0]*bound;
        xb[0] = N[0]; xb[1] = 0;
        for (kx = 0; kx < N[0]; kx++)
        {
            if (type == SWF_COS && kx > N[0]-xD-1)
                wx[kx] = cos(v_pi_c*0.5* rel_bound_weight *(kx-N[0]+1+xD)/xD);
            else if (type == SWF_COS && kx < xD)
                wx[kx] = cos(v_pi_c*0.5* rel_bound_weight *(v_1 - kx/xD));
            else
            {
                wx[kx] = v_1;
                if (kx < xb[0])
                    xb[0] = kx;
                if (kx > xb[1])
                    xb[1] = kx;
            }
        }

        for (kz = 0; kz < N[2]; kz++)
            for (ky = 0; ky < N[1]; ky++)
                for (kx = 0; kx < N[0]; kx++)
                    allocField[kx + (ky + kz*N[1])*N[0]] = wx[kx]*wy[ky]*wz[kz];

        delete [] wx;
        delete [] wy;
        delete [] wz;

        return 0;
    }

    uint32_t divD(CagmVectorField *a);
    uint32_t dotD(CagmVectorField *a, CagmVectorField *b);
    uint32_t abs2D(CagmVectorField *a);
    uint32_t absD(CagmVectorField *a);
    uint32_t invD(CagmScalarField *a);
    uint32_t invD(void);
    uint32_t multD(double c, CagmScalarField *a);
    uint32_t multD(double c);
    uint32_t multD(CagmScalarField *c, CagmScalarField *a);
    uint32_t multD(CagmScalarField *c);
    uint32_t addD(CagmScalarField *a, CagmScalarField *b);
    uint32_t addD(CagmScalarField *a);
    uint32_t subD(CagmScalarField *a, CagmScalarField *b);
    uint32_t subD(CagmScalarField *a);
    uint32_t negD(CagmScalarField *a);
    uint32_t negD();
    uint32_t zeroD();
};
