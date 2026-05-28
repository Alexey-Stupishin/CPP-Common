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
    double *B2sumW;
    
    double *dtime;
    double *LBottom;

    double *residualT, *dthetaT, *dthetaBT;

    int *step_failed;

    double *residual, *dtheta, *dthetaB;
    
    double *sS;
    double *sSW;
    double *sSJ;
    double *sJ;
    double *sSB;
    double *sB;

    double *LF, *LD;

    // just for info, not arrays
    int counter;
    int szN2;
    double max_dL_incr; 

    CagmMetrics(int length, int _szN2, bool calc_metrics = false)
        : szN2(_szN2)
        , counter(0)
        , dtheta(nullptr)
        , dthetaB(nullptr)
        , residual(nullptr)
        , sS(nullptr)
        , sSW(nullptr)
        , sSJ(nullptr)
        , sJ(nullptr)
        , sSB(nullptr)
        , sB(nullptr)
        , LF(nullptr) 
        , LD(nullptr)
    {
        depth = new int[length];
        iterN = new int[length];
        steps = new double[length];
        stepNs = new int[length];
        dLs = new double[length];
        z_size = new int[length];
        voxels = new int[length];
        B2sumW = new double[length];

        LBottom = new double[length];
        dtime = new double[length];

        dthetaT = new double[length];
        dthetaBT = new double[length];
        residualT = new double[length];

        step_failed = new int[length];

        for (int k = 0; k < length; k++)
        {
            depth[k] = 0;
            iterN[k] = 0;
            steps[k] = 0;
            stepNs[k] = 0;
            dLs[k] = 0;
            z_size[k] = 0;
            voxels[k] = 0;
            B2sumW[k] = 0;

            LBottom[k] = 0;
            dtime[k] = 0;

            dthetaT[k] = 0;
            dthetaBT[k] = 0;
            residualT[k] = 0;

            step_failed[k] = 0;
        }

        if (calc_metrics)
        {
            dtheta = new double[length*_szN2];
            dthetaB = new double[length*_szN2];
            residual = new double[length*_szN2];

            sS = new double[length*_szN2];
            sSW = new double[length*_szN2];
            sSJ = new double[length*_szN2];
            sJ = new double[length*_szN2];
            sSB = new double[length*_szN2];
            sB = new double[length*_szN2];

            LF = new double[length*_szN2];
            LD = new double[length*_szN2];
        }
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
        delete [] B2sumW;

        delete [] dtime;
        delete[] LBottom;

        delete[] dthetaT;
        delete[] dthetaBT;
        delete[] residualT;

        delete[] step_failed;

        delete [] sS;
        delete [] sSW;
        delete [] sSJ;
        delete [] sJ;
        delete [] sSB;
        delete [] sB;

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
        step_failed[pos] = _m->step_failed[0];

        z_size[pos] = _m->z_size[0];
        voxels[pos] = _m->voxels[0];
        B2sumW[pos] = _m->B2sumW[0];

        dtime[pos] = _m->dtime[0];
        LBottom[pos] = _m->LBottom[0];

        dthetaT[pos] = _m->dthetaT[0];
        dthetaBT[pos] = _m->dthetaBT[0];
        residualT[pos] = _m->residualT[0];

        if (dtheta) memcpy(dtheta + pos*szN2, _m->dtheta, sizeof(double)*_m->szN2);
        if (dthetaB) memcpy(dthetaB + pos*szN2, _m->dthetaB, sizeof(double)*_m->szN2);
        if (residual) memcpy(residual + pos*szN2, _m->residual, sizeof(double)*_m->szN2);
        if (sS) memcpy(sS + pos*szN2, _m->sS, sizeof(double)*_m->szN2);
        if (sSW) memcpy(sSW + pos*szN2, _m->sSW, sizeof(double)*_m->szN2);
        if (sSJ) memcpy(sSJ + pos*szN2, _m->sSJ, sizeof(double)*_m->szN2);
        if (sJ) memcpy(sJ + pos*szN2, _m->sJ, sizeof(double)*_m->szN2);
        if (sSB) memcpy(sSB + pos*szN2, _m->sSB, sizeof(double)*_m->szN2);
        if (sB) (sB + pos*szN2, _m->sB, sizeof(double)*_m->szN2);
        if (LF) memcpy(LF + pos*szN2, _m->LF, sizeof(double)*_m->szN2);
        if (LD) memcpy(LD + pos*szN2, _m->LD, sizeof(double)*_m->szN2);
    }
};
