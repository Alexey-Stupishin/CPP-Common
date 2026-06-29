#pragma once

#include <random>

inline void Durstenfeld(int n, int *x)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    for (int k = n - 1; k > 0; k--)
    {
        std::uniform_int_distribution<int> distrib(0, k);
        int idx = distrib(gen);
        int t = x[k]; x[k] = x[idx]; x[idx] = t;
    }
}

inline void Durstenfeld(int n, int *x, int *y)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    for (int k = n - 1; k > 0; k--)
    {
        std::uniform_int_distribution<int> distrib(0, k);
        int idx = distrib(gen);
        int t = x[k]; x[k] = x[idx]; x[idx] = t;
        t = y[k]; y[k] = y[idx]; y[idx] = t;
    }
}
