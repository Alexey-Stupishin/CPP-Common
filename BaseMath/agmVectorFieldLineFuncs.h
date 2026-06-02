#pragma once
#include "stdDefinitions.h"

class T_Lines
{
protected:
    double dir;
    double cubeLimLow[3], cubeLimHigh[3];
    CagmVectorFieldOps *vfield;
    int N[3];

public:
    T_Lines(double _dir, CagmVectorFieldOps *_vfield, double _absBoundAchieve, double _relBoundAchieve)
        : dir(_dir)
        , vfield(_vfield)
    {
        vfield->dimensions(N);
        cubeLimLow[0]  = _absBoundAchieve + _relBoundAchieve*N[0];
        cubeLimHigh[0] = N[0] - 1 - cubeLimLow[0];
        cubeLimLow[1]  = _absBoundAchieve + _relBoundAchieve*N[1];
        cubeLimHigh[1] = N[1] - 1 - cubeLimLow[1];
        cubeLimLow[2]  = _absBoundAchieve;
        cubeLimHigh[2] = N[2] - 1 - (_absBoundAchieve + _relBoundAchieve*N[2]);
    }

    double inverse()
    {
        dir = -dir;
        return dir;
    }

    bool inBoundCube(const CagmRKF45Vect& v)
    {
        return v.e[0] < cubeLimLow[0] || v.e[0] > cubeLimHigh[0]
            || v.e[1] < cubeLimLow[1] || v.e[1] > cubeLimHigh[1]
            || v.e[2] < cubeLimLow[2] || v.e[2] > cubeLimHigh[2]
            ;
    }

    uint32_t derivatives(const CagmRKF45Vect& v, CagmRKF45Vect& vp)
    {
        vp.e[0] = 0; vp.e[1] = 0; vp.e[2] = 0;

        if ((v.e[0] >= 0 && v.e[0]  <= N[0] - 1
          && v.e[1] >= 0 && v.e[1]  <= N[1] - 1
          && v.e[2] >= 0 && v.e[2]  <= N[2] - 1))
        {
            double field[3];
            vfield->getPoint(((CagmRKF45Vect&)v).v(), field);
            vp.e[0] = field[0] * dir;
            vp.e[1] = field[1] * dir;
            vp.e[2] = field[2] * dir;
            double n = 1.0/sqrt(vp.e[0]*vp.e[0] + vp.e[1]*vp.e[1] + vp.e[2]*vp.e[2]);
            vp.e[0] *= n;
            vp.e[1] *= n;
            vp.e[2] *= n;
            return 0;
        }

        return 1;
    }
};
