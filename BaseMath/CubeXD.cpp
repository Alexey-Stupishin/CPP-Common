#include "stdDefinitions.h"

#include "CubeXD.h"

//-------------------------------------------------------------------------------
CubeXD::CubeXD(int *_N, int _dim_, double *_step, int *_NL, int *_NH)
    : dim_(_dim_)
{
    N[0] = _N[0]; N[1] = _N[1]; N[2] = _N[2];

    setConfig(_step, _NL, _NH);
}

//-------------------------------------------------------------------------------
CubeXD::CubeXD(CubeXD *sample)
{
    init(sample);
}

//-------------------------------------------------------------------------------
int CubeXD::dim()
{
    return dim_;
}

//-------------------------------------------------------------------------------
void CubeXD::dimensions(int *_N)
{
    _N[0] = N[0]; _N[1] = N[1]; _N[2] = N[2];
}

//-------------------------------------------------------------------------------
int CubeXD::height()
{
    return N[2];
}

//-------------------------------------------------------------------------------
int CubeXD::size()
{
    return N[0] * N[1] * N[2];
}

//-------------------------------------------------------------------------------
int CubeXD::getMatryoshkaDepth(int minChunk, double factor)
{
    int minN = N[0];
    if (minN > N[1])
        minN = N[1];

    double d = (double)minN / minChunk;
    return (int)ceil(log(d) / log(factor));
}

//-------------------------------------------------------------------------------
int CubeXD::getGlobalID(int kx, int ky, int kz)
{
    return (kz*N[1] + ky)*N[0] + kx;
}

//-------------------------------------------------------------------------------
int CubeXD::getGlobalID(double x, double y, double z)
{
    int kx = (int)floor(x);
    int ky = (int)floor(y);
    int kz = (int)floor(z);

    return getGlobalID(kx, ky, kz);
}

//-------------------------------------------------------------------------------
int CubeXD::getGlobalID(double *coord)
{
    return getGlobalID(coord[0], coord[1], coord[2]);
}

//-------------------------------------------------------------------------------
void CubeXD::parseGlobalID(int idx, int *kx, int *ky, int *kz)
{
    *kx = idx % N[0];
    int kyz = (idx - *kx) / N[0];
    *ky = kyz % N[1];
    *kz = (kyz - *ky) / N[1];
}

//-------------------------------------------------------------------------------
int CubeXD::reindexXY(int index)
{
    int kx = index % N[0];
    int kyz = index / N[0];
    int ky = kyz % N[1];
    int kz = kyz / N[1];
    return (kz*N[0] + kx)*N[1] + ky;
}

//-------------------------------------------------------------------------------
void CubeXD::transposeXY(int *data)
{
    int V = N[0] * N[1] * N[2];
    int *reform = new int[V];

    for (int k = 0; k < V; k++)
        reform[reindexXY(k)] = data[k];

    memcpy(data, reform, V*sizeof(int));

    delete[] reform;
}

//-------------------------------------------------------------------------------
void CubeXD::transposeXY(double *data)
{
    int V = N[0] * N[1] * N[2];
    double *reform = new double[V];

    for (int k = 0; k < V; k++)
        reform[reindexXY(k)] = data[k];

    memcpy(data, reform, V*sizeof(double));

    delete[] reform;
}

//-------------------------------------------------------------------------------
void CubeXD::transposeIDXY(int *data)
{
    int V = N[0] * N[1] * N[2];
    int *reform_index = new int[V];

    for (int k = 0; k < V; k++)
        reform_index[k] = reindexXY(data[k]);
    transposeXY(reform_index);
    memcpy(data, reform_index, V*sizeof(int));

    delete[] reform_index;
}

//-------------------------------------------------------------------------------
void CubeXD::init(CubeXD *sample)
{
    dim_ = sample->dim_;
    N[0] = sample->N[0]; N[1] = sample->N[1]; N[2] = sample->N[2];
    setConfig(sample);
}

//-------------------------------------------------------------------------------
void CubeXD::setConfig(double *_step, int *_NL, int *_NH)
{
    setSteps(_step);
    setNPhys(_NL, _NH);
}

//-------------------------------------------------------------------------------
void CubeXD::setConfig(CubeXD *sample)
{
    setConfig(sample->step, sample->NL, sample->NH);
}

//-------------------------------------------------------------------------------
void CubeXD::setSteps(double *_step)
{
    if (_step)
    {
        step[0] = _step[0]; step[1] = _step[1]; step[2] = _step[2];
    }
    else
    {
        step[0] = 1.0; step[1] = 1.0; step[2] = 1.0;
    }
}

//-------------------------------------------------------------------------------
void CubeXD::setNPhys(int *_NL, int *_NH)
{
    if (_NL)
    {
        NL[0] = _NL[0]; NL[1] = _NL[1]; NL[2] = _NL[2];
    }
    else
    {
        NL[0] = 0; NL[1] = 0; NL[2] = 0;
    }

    if (_NH)
    {
        NH[0] = _NH[0]; NH[1] = _NH[1]; NH[2] = _NH[2];
    }
    else
    {
        NH[0] = N[0]; NH[1] = N[1]; NH[2] = N[2];
    }
}

//-------------------------------------------------------------------------------
void CubeXD::setPhysLims(int minChunk, double factor)
{
    int depth = getMatryoshkaDepth(minChunk, factor);

    for (int k = 0; k < 3; k++)
    {
        int coef = 1 << (depth - 1);
        int points = (int)floor((N[k] - 1) / coef);
        int new_lng = (points + 1)*coef + 1;
        int diff = N[k] - new_lng;
        NL[k] = k == 2 ? 0 : (int)floor(diff / 2);
        NH[k] = NL[k] + new_lng;
    }
}
