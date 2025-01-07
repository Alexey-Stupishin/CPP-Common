#pragma once

#include <math.h>

class CubeXD
{
protected:
    int dim_;
    int N[3];
    double step[3];
    int NL[3], NH[3];

public:
    CubeXD(int *_N, int _dim_, double *_step = nullptr, int *_NL = nullptr, int *_NH = nullptr)
        : dim_(_dim_)
    {
        N[0] = _N[0]; N[1] = _N[1]; N[2] = _N[2];

        setConfig(_step, _NL, _NH);
    }

    CubeXD(CubeXD *sample)
    {
        init(sample);
    }

    int dim()
    {
        return dim_;
    }

    void dimensions(int *_N)
    {
        _N[0] = N[0]; _N[1] = N[1]; _N[2] = N[2];
    }

    int size()
    {
        return N[0]*N[1]*N[2];
    }

    int getMatryoshkaDepth(int minChunk, double factor)
    {
        int minN = N[0];
        if (minN > N[1])
            minN = N[1];

        double d = (double)minN / minChunk;
        return (int)ceil(log(d) / log(factor));
    }

protected:
    void init(CubeXD *sample)
    {
        dim_ = sample->dim_;
        N[0] = sample->N[0]; N[1] = sample->N[1]; N[2] = sample->N[2];
        setConfig(sample);
    }

    void setConfig(double *_step = nullptr, int *_NL = nullptr, int *_NH = nullptr)
    {
        setSteps(_step);
        setNPhys(_NL, _NH);
    }

    void setConfig(CubeXD *sample)
    {
        setConfig(sample->step, sample->NL, sample->NH);
    }

    void setSteps(double *_step)
    {
        if (_step)
            { step[0] = _step[0]; step[1] = _step[1]; step[2] = _step[2]; }
        else
            { step[0] = 1.0; step[1] = 1.0; step[2] = 1.0; }
    }

    void setNPhys(int *_NL, int *_NH)
    {
        if (_NL)
            { NL[0] = _NL[0]; NL[1] = _NL[1]; NL[2] = _NL[2]; }
        else
            { NL[0] =    0; NL[1] =    0; NL[2] =    0; }

        if (_NH)
            { NH[0] = _NH[0]; NH[1] = _NH[1]; NH[2] = _NH[2]; }
        else
            { NH[0] = N[0]; NH[1] = N[1]; NH[2] = N[2]; }
    }

    void setPhysLims(int minChunk, double factor)
    {
        int depth = getMatryoshkaDepth(minChunk, factor);

        for (int k = 0; k < 3; k++)
        {
            int coef = 1 << (depth-1);
            int points = (int)floor((N[k]-1)/coef);
            int new_lng = (points+1)*coef + 1;
            int diff = N[k] - new_lng;
            NL[k] = k == 2 ? 0 : (int)floor(diff / 2);
            NH[k] = NL[k] + new_lng;
        }
    }
};