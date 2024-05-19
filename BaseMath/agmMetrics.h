#pragma once

class CagmMetricsJB
{
public:
    double *m;
    double *mj;
    double *mB;

    CagmMetricsJB(int length)
    {
        m = new double[length];
        mj = new double[length];
        mB = new double[length];
        for (int k = 0; k < length; k++)
        {
            m[k] = 0;
            mj[k] = 0;
            mB[k] = 0;
        }
    }
    ~CagmMetricsJB()
    {
        delete [] m;
        delete [] mj;
        delete [] mB;
    }

    void set(CagmMetricsJB *_m, int pos)
    {
        m[pos] = _m->m[0];
        mj[pos] = _m->mj[0];
        mB[pos] = _m->mB[0];
    }
};

class CagmMetricsLim
{
public:
    CagmMetricsJB *mW;
    CagmMetricsJB *mL;

    CagmMetricsLim(int length)
    {
        mW = new CagmMetricsJB(length);
        mL = new CagmMetricsJB(length);
    }
    ~CagmMetricsLim()
    {
        delete mW;
        delete mL;
    }

    void set(CagmMetricsLim *_m, int pos)
    {
        mW->set(_m->mW, pos);
        mL->set(_m->mL, pos);
    }
};

class CagmMetricsCos
{
public:
    double *c;
    double *B4c;

    CagmMetricsCos(int length)
    {
        c = new double[length];
        B4c = new double[length];
        for (int k = 0; k < length; k++)
        {
            c[k] = 0;
            B4c[k] = 0;
        }
    }
    ~CagmMetricsCos()
    {
        delete [] c;
        delete [] B4c;
    }

    void set(CagmMetricsCos *_m, int pos)
    {
        c[pos] = _m->c[0];
        B4c[pos] = _m->B4c[0];
    }
};

class CagmMetricsCosLim
{
public:
    CagmMetricsCos *mW;
    CagmMetricsCos *mL;

    CagmMetricsCosLim(int length)
    {
        mW = new CagmMetricsCos(length);
        mL = new CagmMetricsCos(length);
    }
    ~CagmMetricsCosLim()
    {
        delete mW;
        delete mL;
    }

    void set(CagmMetricsCosLim *_m, int pos)
    {
        mW->set(_m->mW, pos);
        mL->set(_m->mL, pos);
    }
};

class CagmMetrics
{
public:
    int *depth;
    int *iterN;
    double *steps;
    int *stepNs;
    double *dLs;
    int *z_size;
    int *voxels;
    
    double *dtime;
    double *LBottom;

    double *residualT, *dthetaT, *dthetaBT;

    double *residual, *dtheta, *dthetaB;
    
    int szN2;
    double *B2sum;
    double *F2max;
    double *sS;
    double *sSW;
    double *sSJ;
    double *sJ;
    double *sSB;
    double *sB;

    double *L, *LF, *LD;

    // just for info, not array
    double max_dL_incr; 
    int counter;

    CagmMetrics(int length, int _szN2)
        : szN2(_szN2)
        , counter(0)
    {
        depth = new int[length];
        iterN = new int[length];
        steps = new double[length];
        stepNs = new int[length];
        dLs = new double[length];
        z_size = new int[length];
        voxels = new int[length];

        LBottom = new double[length];
        dtime = new double[length];

        dthetaT = new double[length];
        dthetaBT = new double[length];
        residualT = new double[length];

        for (int k = 0; k < length; k++)
        {
            depth[k] = 0;
            iterN[k] = 0;
            steps[k] = 0;
            stepNs[k] = 0;
            dLs[k] = 0;
            z_size[k] = 0;
            voxels[k] = 0;

            LBottom[k] = 0;

            dthetaT[k] = 0;
            dthetaBT[k] = 0;
            residualT[k] = 0;
        }

        dtheta = new double[length*_szN2];
        dthetaB = new double[length*_szN2];
        residual = new double[length*_szN2];

        B2sum = new double[length*_szN2];
        F2max = new double[length*_szN2];
        sS = new double[length*_szN2];
        sSW = new double[length*_szN2];
        sSJ = new double[length*_szN2];
        sJ = new double[length*_szN2];
        sSB = new double[length*_szN2];
        sB = new double[length*_szN2];

        L = new double[length*_szN2];
        LF = new double[length*_szN2];
        LD = new double[length*_szN2];
    }
    ~CagmMetrics()
    {
        delete [] depth;
        delete [] iterN;
        delete[] steps;
        delete[] stepNs;
        delete [] dLs;
        delete[] z_size;
        delete[] voxels;

        delete [] dtime;
        delete[] LBottom;

        delete[] dthetaT;
        delete[] dthetaBT;
        delete[] residualT;

        delete [] B2sum;
        delete [] F2max;
        delete [] sS;
        delete [] sSW;
        delete [] sSJ;
        delete [] sJ;
        delete [] sSB;
        delete [] sB;

        delete [] L;
        delete [] LF;
        delete [] LD;

        delete[] dtheta;
        delete[] dthetaB;
        delete[] residual;
    }

    void setBase(double _dLs, double _step, int _stepN, int _depth, int _iterN, double _dtime, int _z_size, int _voxels)
    {
        depth[0] = _depth;
        iterN[0] = _iterN;
        steps[0] = _step;
        stepNs[0] = _stepN;
        dLs[0] = _dLs;
        z_size[0] = _z_size;
        voxels[0] = _voxels;

        dtime[0] = _dtime;
    }

    // call from callback
    void set(CagmMetrics *_m, int pos)
    {
        counter = _m->counter;
        depth[pos] = _m->depth[0];
        iterN[pos] = _m->iterN[0];
        steps[pos] = _m->steps[0];
        stepNs[pos] = _m->stepNs[0];
        dLs[pos] = _m->dLs[0];
        z_size[pos] = _m->z_size[0];
        voxels[pos] = _m->voxels[0];

        dtime[pos] = _m->dtime[0];
        LBottom[pos] = _m->LBottom[0];

        memcpy(L + pos*szN2, _m->L, sizeof(double)*_m->szN2);

        memcpy(dtheta + pos*szN2, _m->dtheta, sizeof(double)*_m->szN2);
        memcpy(dthetaB + pos*szN2, _m->dthetaB, sizeof(double)*_m->szN2);
        memcpy(residual + pos*szN2, _m->residual, sizeof(double)*_m->szN2);
        memcpy(B2sum + pos*szN2, _m->B2sum, sizeof(double)*_m->szN2);
        memcpy(F2max + pos*szN2, _m->F2max, sizeof(double)*_m->szN2);
        memcpy(sS + pos*szN2, _m->sS, sizeof(double)*_m->szN2);
        memcpy(sSW + pos*szN2, _m->sSW, sizeof(double)*_m->szN2);
        memcpy(sSJ + pos*szN2, _m->sSJ, sizeof(double)*_m->szN2);
        memcpy(sJ + pos*szN2, _m->sJ, sizeof(double)*_m->szN2);
        memcpy(sSB + pos*szN2, _m->sSB, sizeof(double)*_m->szN2);
        memcpy(sB + pos*szN2, _m->sB, sizeof(double)*_m->szN2);
        memcpy(LF + pos*szN2, _m->LF, sizeof(double)*_m->szN2);
        memcpy(LD + pos*szN2, _m->LD, sizeof(double)*_m->szN2);

        dthetaT[pos] = _m->dthetaT[0];
        dthetaBT[pos] = _m->dthetaBT[0];
        residualT[pos] = _m->residualT[0];
    }
};
