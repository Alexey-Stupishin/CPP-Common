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
    CubeXD(int *_N, int _dim_, double *_step = nullptr, int *_NL = nullptr, int *_NH = nullptr);
    CubeXD(CubeXD *sample);

    int dim();
    void dimensions(int *_N);
    int height();
    int size();
    int getMatryoshkaDepth(int minChunk = 25, double factor = 2);

    int getGlobalID(int kx, int ky, int kz);
    int getGlobalID(double x, double y, double z);
    int getGlobalID(double *coord);
    void parseGlobalID(int idx, int *kx, int *ky, int *kz);

    int reindexXY(int index);
    void transposeXY(int *data);
    void transposeXY(double *data);
    void transposeIDXY(int *data);

protected:
    void init(CubeXD *sample);

    void setConfig(double *_step = nullptr, int *_NL = nullptr, int *_NH = nullptr);
    void setConfig(CubeXD *sample);

    void setSteps(double *_step);
    void setNPhys(int *_NL, int *_NH);
    void setPhysLims(int minChunk, double factor);
};