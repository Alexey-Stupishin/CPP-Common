#pragma once

#include <time.h>

class aguTimeTicToc
{
protected:
    clock_t start;

public:
    aguTimeTicToc()
    {
        tic();
    }

    clock_t tic()
    {
        start = clock();
        return start;
    }

    double toc()
    {
        return (double)(clock() - start) / CLOCKS_PER_SEC;
    }
};
