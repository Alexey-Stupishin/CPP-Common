#include "stdDefinitions.h"
#include "agsFieldsCommon.h"
#define _USE_MATH_DEFINES
#include "math.h"

#include "agmVectorField.h"
#include "agmScalarField.h"

//-----------------------------------------------------------------------
uint32_t CagmVectorField::FIA2XYZ(CagmVectorField *FIA)
{
	int kx, ky, kz;
    for (kz = 0; kz < N[2]; kz++)
        for (ky = 0; ky < N[1]; ky++)
            for (kx = 0; kx < N[0]; kx++)
            {
                double cos_incl = cos(Deg2Rad(FIA->fieldY[fidx(kx, ky, kz)]));
                double sin_incl = sin(Deg2Rad(FIA->fieldY[fidx(kx, ky, kz)]));
                double cos_azim = cos(Deg2Rad(FIA->fieldZ[fidx(kx, ky, kz)]));
                double sin_azim = sin(Deg2Rad(FIA->fieldZ[fidx(kx, ky, kz)]));
                fieldZ[fidx(kx, ky, kz)] = FIA->fieldX[fidx(kx, ky, kz)] * cos_incl;
                double transv = FIA->fieldX[fidx(kx, ky, kz)] * sin_incl;
                fieldX[fidx(kx, ky, kz)] = transv * cos_azim;
                fieldY[fidx(kx, ky, kz)] = transv * sin_azim;
            }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorField::InvertAzimuth(CagmVectorField *B)
{
	int kx, ky, kz;
    for (kz = 0; kz < N[2]; kz++)
        for (ky = 0; ky < N[1]; ky++)
            for (kx = 0; kx < N[0]; kx++)
            {
                fieldX[fidx(kx, ky, kz)] = - B->fieldX[fidx(kx, ky, kz)];
                fieldY[fidx(kx, ky, kz)] = - B->fieldY[fidx(kx, ky, kz)];
            }

    return 0;
}

#define lay(k,n) (k & (1<<n) ? 1: 0)

//-----------------------------------------------------------------------
double *CagmVectorField::disambigGetF(CagmVectorField *_B000, CagmVectorField *_B180, CagmVectorField *_Bref, double *_step)
{
    step[0] = _step[0];
    step[1] = _step[1];
    step[2] = 0;

    B000 = _B000;
    B180 = _B180;
    Bref = _Bref;

    int sizeB = N[0]*N[1];
    double *Bexx = new double[sizeB*2];
    double *Bexy = new double[sizeB*2];
    for (int ky = 0; ky < N[1]; ky++)
        for (int kx = 0; kx < N[0]; kx++)
        {
            int pos = kx + ky*N[0];
            Bexx[pos] =       _B000->fieldX[fidx2(kx, ky)];
            Bexx[pos+sizeB] = _B180->fieldX[fidx2(kx, ky)];
            Bexy[pos] =       _B000->fieldY[fidx2(kx, ky)];
            Bexy[pos+sizeB] = _B180->fieldY[fidx2(kx, ky)];
        }

    int sizeF = (N[0]+1)*(N[1]+1);
    delete[] daF;
    daF = new double[sizeF*16];
    for (int k = 0; k < sizeF*16; k++)
        daF[k] = 0;

    for (int ky = 0; ky < N[1]-1; ky++)
        for (int kx = 0; kx < N[0] - 1; kx++)
        {
            double divPhL = (_Bref->fieldX[fidx2(kx+1, ky)] - _Bref->fieldX[fidx2(kx, ky)])/step[0];
            double divPhT = (_Bref->fieldY[fidx2(kx, ky+1)] - _Bref->fieldY[fidx2(kx, ky)])/step[1];
            double divPhR = (_Bref->fieldX[fidx2(kx+1, ky+1)] - _Bref->fieldX[fidx2(kx, ky+1)])/step[0];
            double divPhB = (_Bref->fieldY[fidx2(kx+1, ky+1)] - _Bref->fieldY[fidx2(kx+1, ky)])/step[1];
            int pos = kx   + ky    *N[0];
            int pdx = kx+1 + ky    *N[0];
            int pdy = kx   + (ky+1)*N[0];
            int pdd = kx+1 + (ky+1)*N[0];
            for (int k = 0; k < 16; k++)
            {
                double divBhL = (Bexx[pdx + sizeB*lay(k,0)] - Bexx[pos + sizeB*lay(k,1)])/step[0] - divPhL;
                double divBhT = (Bexy[pdy + sizeB*lay(k,3)] - Bexy[pos + sizeB*lay(k,1)])/step[1] - divPhT;
                double divBhR = (Bexx[pdd + sizeB*lay(k,2)] - Bexx[pdy + sizeB*lay(k,3)])/step[0] - divPhR;
                double divBhB = (Bexy[pdd + sizeB*lay(k,2)] - Bexy[pdx + sizeB*lay(k,0)])/step[1] - divPhB;
                double Fdiv = (fabs(divBhL + divBhT) + fabs(divBhR + divBhB) + fabs(divBhR + divBhT) + fabs(divBhL + divBhB))/4;

                double dYxL = (Bexy[pdx + sizeB*lay(k,0)] - Bexy[pos + sizeB*lay(k,1)])/2/step[0];
                double dXyT = (Bexx[pdy + sizeB*lay(k,3)] - Bexx[pos + sizeB*lay(k,1)])/2/step[1];
                double dYxR = (Bexy[pdd + sizeB*lay(k,2)] - Bexy[pdy + sizeB*lay(k,3)])/2/step[0];
                double dXyB = (Bexx[pdd + sizeB*lay(k,2)] - Bexx[pdx + sizeB*lay(k,0)])/2/step[1];

                double JzLT = dYxL - dXyT;
                double JzRB = dYxR - dXyB;
                double JzLB = dYxL - dXyB;
                double JzRT = dYxR - dXyT;

                double FJ = (fabs(JzLT) + fabs(JzRB) + fabs(JzLB) + fabs(JzRT))/4;

                daF[kx+1 + (ky+1)*(N[0]+1) + sizeF*k] = Fdiv + FJ;
            }
        }

    delete [] Bexx;
    delete [] Bexy;

    return daF;
}

//-----------------------------------------------------------------------
double *CagmVectorField::disambigGetT(double ktfactor_vp, double ktfactor_v0, double ktfactor_M, double ktfactor_p, double ktfactor_init)
{
    int sizeB = N[0] * N[1];
    int sizeF = (N[0]+1) * (N[1]+1);
    int N01 = N[0] + 1;
    //double *rule = new double[sizeB];

    //CagmScalarField s(N);
    //s.abs(B000);
    //double BSmax =  s.max_plane(0);
    //for (int ky = 0; ky < N[1]; ky++)
    //    for (int kx = 0; kx < N[0]; kx++)
    //        rule[kx + ky*N[0]] = s.field[fidx(kx, ky, 0)]/BSmax;

    double * kTex = new double[sizeF];
    for (int ky = 0; ky < N[1]+1; ky++)
        for (int kx = 0; kx < N[0]+1; kx++)
        {
            int pos = kx + ky*(N[0] + 1);
            double Fmax = 0;
            for (int k = 0; k < 16; k++)
                Fmax = daF[pos + k*sizeF] > Fmax ? daF[pos + k*sizeF] : Fmax;
            kTex[pos] = Fmax;
        }

//kT = kT*(1+param.kTFactor_M);
//kTmax = xmax2(kT);
//kk = kT/kTmax;
//kT = kTmax .* kk.^param.kTFactor_p;

    delete[] daT;
    daT = new double[sizeB];
    double kTmax = 0;
    for (int ky = 0; ky < N[1]; ky++)
        for (int kx = 0; kx < N[0]; kx++)
        {
            int pos = kx   + ky    *N[0];
            int pd0 = kx   + ky    *N01;
            int pdx = kx+1 + ky    *N01;
            int pdy = kx   + (ky+1)*N01;
            int pdd = kx+1 + (ky+1)*N01;
            double KTmax = kTex[pd0];
            KTmax = kTex[pdx] > KTmax ? kTex[pdx] : KTmax;
            KTmax = kTex[pdy] > KTmax ? kTex[pdy] : KTmax;
            daT[pos] = kTex[pdd] > KTmax ? kTex[pdd] : KTmax;
            kTmax = daT[pos] > kTmax ? daT[pos] : kTmax;
        }

    //double alpha = log((ktfactor_vp-ktfactor_M)/(ktfactor_v0-ktfactor_M)) / log(ktfactor_p);
    //if (isnan(alpha))
    //    alpha = 0;
    for (int ky = 0; ky < N[1]; ky++)
        for (int kx = 0; kx < N[0]; kx++)
        {
            int pos = kx + ky*N[0];
            daT[pos] = kTmax * pow(daT[pos]/kTmax, ktfactor_p);
            //double exTFactor = pow(rule[pos], alpha);
            //daT[pos] *= ktfactor_init * (ktfactor_M + (ktfactor_v0-ktfactor_M)*exTFactor);
        }

    //delete [] rule;
    delete [] kTex;

    return daT;
}
