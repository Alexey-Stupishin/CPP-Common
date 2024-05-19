#pragma once

#include "agmBaseMath.h"

class CagmRotate3D
{
public:
    double vcos[3];
    double lat, lon;

protected:
    double VCS2VAR[3][3], VAR2VCS[3][3];
    double sinlat, sinlon, coslat, coslon;
    double eyeTol;

public:
	CagmRotate3D(double *_vcos, double _eyeTol = 1e-6)
    {
        init(_vcos, _eyeTol);
    }
	CagmRotate3D(double _lon, double _lat, double _eyeTol = 1e-6)
    {
        init(_lon, _lat, _eyeTol);
    }
	CagmRotate3D()
    {
    }
	virtual ~CagmRotate3D()
    {}

    void copyRotator(CagmRotate3D *rotator)
    {
        memcpy(vcos, rotator->vcos, 3*sizeof(double));
        lat = rotator->lat;
        lon = rotator->lon;

        memcpy(VCS2VAR, rotator->VCS2VAR, 9 * sizeof(double));
        memcpy(VAR2VCS, rotator->VAR2VCS, 9 * sizeof(double));
        sinlat = rotator->sinlat;
        sinlon = rotator->sinlon;
        coslat = rotator->coslat;
        coslon = rotator->coslon;
        eyeTol = rotator->eyeTol;
    }

	virtual void init(double *_vcos, double _eyeTol = 1e-6)
    {
        eyeTol = _eyeTol;

        for (int d= 0; d < 3; d++)
            vcos[d] = _vcos[d];

        sinlon = -vcos[1];
        coslon = sqrt(1 - sinlon*sinlon);
        coslat = vcos[2]/coslon;
        sinlat = -vcos[0]/coslon;

        lat = asin(sinlat)*v_180/v_pi_c;
        lon = asin(sinlon)*v_180/v_pi_c;

        createMatrices();
    }
	virtual void init(double _lon, double _lat, double _eyeTol = 1e-6)
    {
        eyeTol = _eyeTol;

        lon = _lon;
        lat = _lat;

        getDirCos(this->lon, this->lat, this->vcos, &(this->sinlon), &(this->coslon), &(this->sinlat), &(this->coslat));

        createMatrices();
    }

    static void getDirCos(double _lon, double _lat, double *_vcos
                        , double *_sinlon = nullptr, double *_coslon = nullptr, double *_sinlat = nullptr, double *_coslat = nullptr)
    {
        double sinlon = sin(_lon*v_pi_c / v_180);
        double coslon = cos(_lon*v_pi_c / v_180);
        double sinlat = sin(_lat*v_pi_c / v_180);
        double coslat = cos(_lat*v_pi_c / v_180);
        if (_sinlon)
            *_sinlon = sinlon;
        if (_coslon)
            *_coslon = coslon;
        if (_sinlat)
            *_sinlat = sinlat;
        if (_coslat)
            *_coslat = coslat;

        _vcos[0] = -sinlat*coslon;
        _vcos[1] = -sinlon;
        _vcos[2] = coslat*coslon;
    }

    bool isEye()
    {
        return fabs(vcos[0]) < eyeTol && fabs(vcos[1]) < eyeTol && fabs(1 - vcos[2]) < eyeTol;
    }

    void rotate(double *v, double *vr, bool direction)
    {
        if (direction)
            rotateP2V(v, vr);
        else
            rotateV2P(v, vr);
    }

    void rotate(double x, double y, double z, double *xr, double *yr, double *zr, bool direction)
    {
        if (direction)
            rotateP2V(x, y, z, xr, yr, zr);
        else
            rotateV2P(x, y, z, xr, yr, zr);
    }

    void rotateP2V(double x, double y, double z, double *xr, double *yr, double *zr)
    {
        double xt = VAR2VCS[0][0] * x + VAR2VCS[0][1] * y + VAR2VCS[0][2] * z;
        double yt = VAR2VCS[1][0] * x + VAR2VCS[1][1] * y + VAR2VCS[1][2] * z;
        double zt = VAR2VCS[2][0] * x + VAR2VCS[2][1] * y + VAR2VCS[2][2] * z;
        *xr = xt;
        *yr = yt;
        *zr = zt;
    }
    void rotateP2V(double *v, double *vr)
    {
        rotateP2V(v[0], v[1], v[2], vr, vr + 1, vr + 2);
    }
    void rotateP2V(double *v)
    {
        rotateP2V(v[0], v[1], v[2], v, v + 1, v + 2);
    }

    void rotateV2P(double x, double y, double z, double *xr, double *yr, double *zr)
    {
        *xr = VCS2VAR[0][0] * x + VCS2VAR[0][1] * y + VCS2VAR[0][2] * z;
        *yr = VCS2VAR[1][0] * x + VCS2VAR[1][1] * y + VCS2VAR[1][2] * z;
        *zr = VCS2VAR[2][0] * x + VCS2VAR[2][1] * y + VCS2VAR[2][2] * z;
    }
    void rotateV2P(double *v, double *vr)
    {
        rotateV2P(v[0], v[1], v[2], vr, vr + 1, vr + 2);
    }
protected:
    void createMatrices()
    {
        VAR2VCS[0][0] = VCS2VAR[0][0] = coslat;
        VAR2VCS[1][0] = VCS2VAR[0][1] = -sinlat*sinlon;
        VAR2VCS[2][0] = VCS2VAR[0][2] = -sinlat*coslon;
        VAR2VCS[0][1] = VCS2VAR[1][0] = 0;
        VAR2VCS[1][1] = VCS2VAR[1][1] = coslon;
        VAR2VCS[2][1] = VCS2VAR[1][2] = -sinlon;
        VAR2VCS[0][2] = VCS2VAR[2][0] = sinlat;
        VAR2VCS[1][2] = VCS2VAR[2][1] = coslat*sinlon;
        VAR2VCS[2][2] = VCS2VAR[2][2] = coslat*coslon;
    }
};
