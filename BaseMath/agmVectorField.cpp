#include "stdDefinitions.h"
#include "agsFieldsCommon.h"
#define _USE_MATH_DEFINES
#include "math.h"

#include "agmVectorField.h"
#include "agmScalarField.h"
#include "agmMetrics.h"

#include "agmRotate3D.h"

//-----------------------------------------------------------------------
uint32_t CagmVectorField::condWeight(int WiegelmannProcCondBase, CagmVectorField *baseField, CagmVectorField *baseWeight, int WiegelmannProcCondBase2, CagmVectorField *baseField2, CagmVectorField *baseWeight2,
                                  int WiegelmannProcCondAbs, CagmScalarField *absField, CagmScalarField *absWeight, int WiegelmannProcCondAbs2, CagmScalarField *absField2, CagmScalarField *absWeight2,
                                  int WiegelmannProcCondLOS, CagmScalarField *losField, CagmScalarField *losWeight, int WiegelmannProcCondLOS2, CagmScalarField *losField2, CagmScalarField *losWeight2, 
                                  CagmRotate3D *rotator)
{
    int iAbs = (absField && absWeight ? WiegelmannProcCondAbs : 0);
    int iAbs2 = (absField2 && absWeight2 ? WiegelmannProcCondAbs2 : 0);
    int iLOS = (losField && losWeight ? WiegelmannProcCondLOS : 0);
    int iLOS2 = (losField2 && losWeight2 ? WiegelmannProcCondLOS2 : 0);
    int ixyz = (baseField && baseWeight ? WiegelmannProcCondBase : 0);
    int ixyz2 = (baseField2 && baseWeight2 ? WiegelmannProcCondBase2 : 0);
    if (iLOS != 0 || iLOS2 != 0 || iAbs != 0 || iAbs2 != 0)
    {
        CagmScalarField w(this);
        CagmScalarField Babs(this);
        Babs.abs(this);

        if (iAbs != 0)
        {
            w = *absWeight;
            if (iAbs == 1)
                w.limWeight(1, &Babs, absField);
            Babs.relax(absField, &w);
        }

        if (iAbs2 == 1)
        {
            w = *absWeight2;
            w.limWeight(-1, &Babs, absField2);
            Babs.relax(absField2, &w);
        }

        if (iLOS != 0 || iLOS2 != 0)
        {
            CagmVectorField *Brot = nullptr;
            bool bRotate = false;
            if (rotator->isEye())
                Brot = this;
            else
            {
                bRotate = true;
                Brot = new CagmVectorField(this);
                Brot->rotate3D(rotator, true);
            }

            CagmScalarField *proj = new CagmScalarField(Brot);
            proj->projection(this, rotator->vcos);
            
            CagmScalarField Blos(Brot);
            CagmScalarField Btr(Brot);
            CagmScalarField Bv(Brot);

            Brot->getTransv(&Btr);
            Brot->getComponent(&Blos, PLANE_Z);

            if (iLOS != 0)
            {
                w = *losWeight;
                if (iLOS == 1)
                    w.limWeight(1, &Blos, losField);
                Blos.relax(losField, &w);
            }
            if (iLOS2 == 1)
            {
                w = *losWeight2;
                w.limWeight(-1, &Blos, losField2);
                Blos.relax(losField2, &w);
            }

            CagmScalarField Btrn(Brot);
            Btrn.sqDiff(&Babs, &Blos);

            Brot->setComponent(&Blos, PLANE_Z);

            Btr.inv();
            Btrn.mult(&Btr);
            Brot->getComponent(&Bv, PLANE_X);
            Bv.mult(&Btrn);
            Brot->setComponent(&Bv, PLANE_X);
            Brot->getComponent(&Bv, PLANE_Y);
            Bv.mult(&Btrn);
            Brot->setComponent(&Bv, PLANE_Y);

            Brot->rotate3D(rotator, false);

            if (bRotate)
            {
                *this = *Brot;
                delete Brot;
            }
        }
        else
        {
            CagmScalarField Babsn(this);
            Babsn.abs(this);
            Babsn.inv();
            Babs.mult(&Babsn);

            this->mult(&Babs);
        }

        return 0;
    }

    if (ixyz != 0 || ixyz2 != 0)
    {
        CagmScalarField Bv(this);
        CagmScalarField bv(this);
        CagmScalarField wv(this);
        CagmScalarField bv2(this);
        CagmScalarField wv2(this);

        this->getComponent(&Bv, PLANE_X);
        if (ixyz != 0)
        {
            baseField->getComponent(&bv, PLANE_X);
            baseWeight->getComponent(&wv, PLANE_X);
            if (ixyz == 1)
                wv.limWeight(1, &bv, &wv);
            Bv.relax(&bv, &wv);
        }
        if (ixyz2 == 1)
        {
            baseField2->getComponent(&bv, PLANE_X);
            baseWeight2->getComponent(&wv, PLANE_X);
            wv.limWeight(-1, &bv, &wv);
            Bv.relax(&bv, &wv);
        }
        this->setComponent(&Bv, PLANE_X);

        this->getComponent(&Bv, PLANE_Y);
        if (ixyz != 0)
        {
            baseField->getComponent(&bv, PLANE_Y);
            baseWeight->getComponent(&wv, PLANE_Y);
            if (ixyz == 1)
                wv.limWeight(1, &bv, &wv);
            Bv.relax(&bv, &wv);
        }
        if (ixyz2 == 1)
        {
            baseField2->getComponent(&bv, PLANE_Y);
            baseWeight2->getComponent(&wv, PLANE_Y);
            wv.limWeight(-1, &bv, &wv);
            Bv.relax(&bv, &wv);
        }
        this->setComponent(&Bv, PLANE_Y);

        this->getComponent(&Bv, PLANE_Z);
        if (ixyz != 0)
        {
            baseField->getComponent(&bv, PLANE_Z);
            baseWeight->getComponent(&wv, PLANE_Z);
            if (ixyz == 1)
                wv.limWeight(1, &bv, &wv);
            Bv.relax(&bv, &wv);
        }
        if (ixyz2 == 1)
        {
            baseField2->getComponent(&bv, PLANE_Z);
            baseWeight2->getComponent(&wv, PLANE_Z);
            wv.limWeight(-1, &bv, &wv);
            Bv.relax(&bv, &wv);
        }
        this->setComponent(&Bv, PLANE_Z);
    }

    return 0;
}
