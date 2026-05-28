#pragma once
#include <chrono>
#include <ctime>

class aguTimeTicToc
{
protected:
    std::chrono::time_point<std::chrono::steady_clock> timestamp;
    clock_t clocktimestamp;

public:
    aguTimeTicToc()
    {
        tic();
        tic_clock();
    }

    std::chrono::time_point<std::chrono::steady_clock> tic()
    {
        return timestamp = std::chrono::steady_clock::now(); 
    }

    clock_t tic_clock()
    {
        return clocktimestamp = clock(); 
    }

    double toc()
    {
        auto now = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed = now - timestamp;
        return elapsed.count();
    }

    double toc_clock()
    {
        return (double)(clock() - clocktimestamp) / CLOCKS_PER_SEC;
    }
};

//#include <time.h>
//clock_t start = clock();
//return (double)(clock() - start) / CLOCKS_PER_SEC;

//FILETIME FileTime;
//ULARGE_INTEGER time;
//GetSystemTimeAsFileTime(&FileTime);
//time.HighPart = FileTime.dwHighDateTime;
//time.LowPart = FileTime.dwLowDateTime;

//ULONGLONG udt = time.QuadPart-time0.QuadPart;
//double seconds = (double)udt*1e-7;
