#pragma once
#include "stdDefinitions.h"

#include "agmVectorFieldOps.h"
#include "agsFieldsCommon.h"
#include "binUtilities.h"

class CagmRotate3D;
class CagmMetricsLim;
class CagmMetricsCosLim;

class CagmScalarField;

//------------------------------------------------------------------
class CagmVectorField : public CagmVectorFieldOps
{
friend class CagmScalarField;

public:
    static uint32_t GetAllocSize(int *N)
    {
        return 3*sizeof(double)*N[0]*N[1]*N[2] + sizeof(CagmVectorField) + CagmVectorFieldOps::GetAllocSize(N);
    }

public:
    double *allocFieldX,  *allocFieldY,  *allocFieldZ;

protected:
	bool isRef;

protected:
    uint32_t Alloc()
    {
        allocFieldX = new double[N[0]*N[1]*N[2]];
        allocFieldY = new double[N[0]*N[1]*N[2]];
        allocFieldZ = new double[N[0]*N[1]*N[2]];
        memset(allocFieldX, 0, sizeof(double)*N[0]*N[1]*N[2]);
        memset(allocFieldY, 0, sizeof(double)*N[0]*N[1]*N[2]);
        memset(allocFieldZ, 0, sizeof(double)*N[0]*N[1]*N[2]);

        for (int ky = 0; ky < N[1]; ky++)
            for (int kz = 0; kz < N[2]; kz++)
            {
                fieldX[ky + kz*N[1]] = allocFieldX + (ky + kz*N[1])*N[0];
                fieldY[ky + kz*N[1]] = allocFieldY + (ky + kz*N[1])*N[0];
                fieldZ[ky + kz*N[1]] = allocFieldZ + (ky + kz*N[1])*N[0];
            }

        return 0;
    }

    uint32_t Delete()
    {
		if (!isRef)
        {
            delete [] allocFieldX;
            delete [] allocFieldY;
            delete [] allocFieldZ;
        }

        return 0;
    }

public:
	CagmVectorField(int *_N, double *_step = nullptr, int *_NL = nullptr, int *_NH = nullptr)
        : CagmVectorFieldOps(_N, _step, _NL, _NH)
        , allocFieldX(nullptr)  
        , allocFieldY(nullptr)  
        , allocFieldZ(nullptr)  
        , isRef(false)
    {
        Alloc();
	}

	CagmVectorField(CbinDataStruct* data)
        : CagmVectorFieldOps(data->GetDimensions())
        , allocFieldX(nullptr)  
        , allocFieldY(nullptr)  
        , allocFieldZ(nullptr)  
        , isRef(false)
    {
        Alloc();
        data->Copy(allocFieldX, allocFieldY, allocFieldZ);
	}

    CagmVectorField(CubeXD *_from, bool copy = false)
        : CagmVectorFieldOps(_from)
        , allocFieldX(nullptr)  
        , allocFieldY(nullptr)  
        , allocFieldZ(nullptr)  
    {
        dim_ = 3;
        Alloc();

        if (!copy)
            return;

        CagmVectorField *from = (CagmVectorField *)_from;

        isRef = from->isRef;
        Copy(from);
    }

	CagmVectorField(int *N, double *Bx, double *By, double *Bz, bool _isRef = false)
        : CagmVectorFieldOps(N)
        , allocFieldX(nullptr)  
        , allocFieldY(nullptr)  
        , allocFieldZ(nullptr)  
        , isRef(_isRef)
    {
        if (isRef)
        {
            for (int kz = 0; kz < N[2]; kz++)
                for (int ky = 0; ky < N[1]; ky++)
                {
                    fieldX[ky + kz*N[1]] = Bx + (ky + kz*N[1])*N[0];
                    fieldY[ky + kz*N[1]] = By + (ky + kz*N[1])*N[0];
                    fieldZ[ky + kz*N[1]] = Bz + (ky + kz*N[1])*N[0];
                }
        }
        else
        {
            Alloc();

            CopyComp(Bx, 0);
            CopyComp(By, 1);
            CopyComp(Bz, 2);
        }
	}

    uint32_t Copy(CagmVectorField *from)
    {
        if (from->isRef)
        {
            for (int ky = 0; ky < N[1]; ky++)
                for (int kz = 0; kz < N[2]; kz++)
                {
                    memcpy(allocFieldX + (ky + kz*N[1])*N[0], from->allocFieldX + (ky + kz*N[1])*N[0], sizeof(double)*N[0]);
                    memcpy(allocFieldY + (ky + kz*N[1])*N[0], from->allocFieldY + (ky + kz*N[1])*N[0], sizeof(double)*N[0]);
                    memcpy(allocFieldZ + (ky + kz*N[1])*N[0], from->allocFieldZ + (ky + kz*N[1])*N[0], sizeof(double)*N[0]);
                }
        }
        else
        {
            memcpy(allocFieldX, from->allocFieldX, sizeof(double)*N[0] * N[1] * N[2]);
            memcpy(allocFieldY, from->allocFieldY, sizeof(double)*N[0] * N[1] * N[2]);
            memcpy(allocFieldZ, from->allocFieldZ, sizeof(double)*N[0] * N[1] * N[2]);
        }

        return 0;
	}

    uint32_t CopyComp(double *Bc, int comp)
    {
        double *c = nullptr;
        if (comp == 0)
            c = allocFieldX;
        else if (comp == 1)
            c = allocFieldY;
        else if (comp == 2)
            c = allocFieldZ;

        memcpy(c, Bc, N[0]*N[1]*N[2]*sizeof(double));

        return 0;
    }

    uint32_t GetComp(double *Bc, int comp)
    {
        double *c = nullptr;
        if (comp == 0)
            c = allocFieldX;
        else if (comp == 1)
            c = allocFieldY;
        else if (comp == 2)
            c = allocFieldZ;

        memcpy(Bc, c, N[0]*N[1]*N[2]*sizeof(double));

        return 0;
    }

	virtual ~CagmVectorField()
    {
        Delete();
    }

    uint32_t GetFieldAddress(double **px, double **py, double **pz)
    {
        *px = allocFieldX;
        *py = allocFieldY;
        *pz = allocFieldZ;

        return 0;
    }

	uint32_t GetData(CbinDataStruct* data)
    {
        return data->Create(N, allocFieldX, allocFieldY, allocFieldZ);
    }

    uint32_t condWeight(int WiegelmannProcCondBase, CagmVectorField *baseField, CagmVectorField *baseWeight, int WiegelmannProcCondBase2, CagmVectorField *baseField2, CagmVectorField *baseWeight2,
        int WiegelmannProcCondAbs, CagmScalarField *absField, CagmScalarField *absWeight, int WiegelmannProcCondAbs2, CagmScalarField *absField2, CagmScalarField *absWeight2,
        int WiegelmannProcCondLOS, CagmScalarField *losField, CagmScalarField *losWeight, int WiegelmannProcCondLOS2, CagmScalarField *losField2, CagmScalarField *losWeight2, 
        CagmRotate3D *rotator);
};
