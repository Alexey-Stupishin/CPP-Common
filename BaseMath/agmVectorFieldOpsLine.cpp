#include "stdDefinitions.h"
#include <math.h>

#include "agsFieldsCommon.h"
#include "agmVectorFieldOps.h"
#include "agmRKF45.h"
#include "agmVectorFieldLineFuncs.h"

static void v_copyCoord(double *coord, CagmRKF45Vect *rkfv)
{
    coord[0] = (*rkfv)[0];
    coord[1] = (*rkfv)[1];
    coord[2] = (*rkfv)[2];
}

//-----------------------------------------------------------------------
CagmVectorFieldOps::Status CagmVectorFieldOps::getOneLine(CagmRKF45 *rkf45, CagmRKF45Vect *rkfv, double _step, double *coord, int maxlen, int *length, CagmRKF45::Status *status, bool noDuplicate)
{
    *length = 0;
    CagmVectorFieldOps::Status outstatus = CagmVectorFieldOps::Status::None;
    double t = 0;
    double s = _step;
    if (!noDuplicate)
    {
        v_copyCoord(coord, rkfv);
        (*length)++;
    }
    for (int i = 1; i < maxlen; i++)
    {
        if (i >= maxlen)
        {
            outstatus = CagmVectorFieldOps::Status::BufferOverload;
            break;
        }

        *status = rkf45->calculate(t, *rkfv, t + s, false);
        if (!CagmRKF45::isStatusCritical(*status))
        {
            if (*status != CagmRKF45::Status::EndNoMove)
                v_copyCoord(coord+3*((*length)++), rkfv);

            if (CagmRKF45::isFinished(*status))
            {
                outstatus = (*status == CagmRKF45::Status::EndByCond ?
                    CagmVectorFieldOps::Status::Boundary : CagmVectorFieldOps::Status::Stable);
                break;
            }
        }
        else
        {
            outstatus = CagmVectorFieldOps::Status::RKF45Problem;
            break;
        }
    }

    return outstatus;
}

//-----------------------------------------------------------------------
CagmVectorFieldOps::Status CagmVectorFieldOps::getOneFullLine(CagmRKF45 *rkf45, double *start, int direction, double _step, double absBoundAchieve, double relBoundAchieve, 
    int maxLength, int *length, double *coord, int *code)
{
    _step = fabs(_step);
    int dirc = (direction >= 0 ? 1 : -1);
    CagmRKF45Vect rkfv(3, start);

    CagmVectorFieldOps::Status status;

    *length = 0;
    *code = 0;

    double B0[3];
    uint32_t rc = getPoint(start, B0);
    if (rc != 0)
    {
        *code = (int)CagmVectorFieldOps::Status::OutOfCube;
        return CagmVectorFieldOps::Status::OutOfCube;
    }

    CagmRKF45::Status rkfstatus;

    T_Lines data(dirc, this, absBoundAchieve, relBoundAchieve);

    rkf45->reinit(&data);
    int currlen;
    status = getOneLine(rkf45, &rkfv, _step, coord, maxLength, &currlen, &rkfstatus);

    bool needBack = (direction == 0);
    bool noDuplicate = false;
    if (currlen > 1)
    {
        if (needBack)
        {
            for (int krev = 0; krev < currlen/2; krev++)
            {
                double tmp = coord[krev*3];   coord[krev*3]   = coord[(currlen-1-krev)*3];   coord[(currlen-1-krev)*3]   = tmp;
                       tmp = coord[krev*3+1]; coord[krev*3+1] = coord[(currlen-1-krev)*3+1]; coord[(currlen-1-krev)*3+1] = tmp;
                       tmp = coord[krev*3+2]; coord[krev*3+2] = coord[(currlen-1-krev)*3+2]; coord[(currlen-1-krev)*3+2] = tmp;
            }
            noDuplicate = true;
        }
    }
    else
        currlen = 0;

    int rest = maxLength - currlen;
    double *posc = coord + currlen*3;
    *length += currlen;
    *code = (int)status + (((int)rkfstatus)<<16);

    if (status & CagmVectorFieldOps::Status::BufferOverload)
        return status;

    if (needBack)
    {
        data.inverse();
        rkfv = start;
        rkf45->reinit(&data);
        status = getOneLine(rkf45, &rkfv, _step, posc, rest, &currlen, &rkfstatus, noDuplicate);
        if (currlen <= 1)
            currlen = 0;

        *length += currlen;
        *code += (((int)status)<<8) + (((int)rkfstatus)<<24);
    }

    return status;
}
