#pragma once
#include "stdDefinitions.h"

class T_Lines
{
public:
    double dir;
    double absBoundAchieve;
    double absBoundAchieveBottom;
    double relBoundAchieve;
    CagmVectorFieldOps *vfield;
    int *N;
    double **fieldX, **fieldY, **fieldZ;
};
