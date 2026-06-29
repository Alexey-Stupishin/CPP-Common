#pragma once
#include <math.h>

class CagmRKF45Vect
{
public:
    int n;
    double *e;

public:
    CagmRKF45Vect()
        : n(0),
        e(nullptr)
    {
    }

    CagmRKF45Vect(int _n, double _v = 0)
        : n(_n),
          e(nullptr) 
    {
        init(n, _v);
    }

    CagmRKF45Vect(int _n, double *v)
        : n(_n),
        e(nullptr)
    {
        if (n > 0)
        {
            e = new double[n];
            for (int i = 0; i < n; i++)
                e[i] = v[i];
        }
    }

    CagmRKF45Vect(const CagmRKF45Vect& v)
        : e(nullptr) 
    {
        n = v.GetLength();
        if (n > 0)
            e = new double[n];
        for (int i = 0; i < n; i++)
            e[i] = v[i];
    }

    void init(int _n, double _v = 0)
    {
        n = _n;
        if (n > 0)
        {
            e = new double[n];
            for (int i = 0; i < n; i++)
                e[i] = _v;
        }
    }

    void clear()
    {
        delete [] e;
    }

    virtual ~CagmRKF45Vect()
    {
        clear();
    }

    CagmRKF45Vect& vmin(const CagmRKF45Vect& v)
    {
        for (int i = 0; i < n; i++)
            e[i] = e[i] < v[i] ? e[i] : v[i];

        return *this;
    }

    CagmRKF45Vect& vmax(const CagmRKF45Vect& v)
    {
        for (int i = 0; i < n; i++)
            e[i] = e[i] > v[i] ? e[i] : v[i];

        return *this;
    }

    static bool inRange(const CagmRKF45Vect& min, const CagmRKF45Vect& max, double lim)
    {
        for (int i = 0; i < min.n; i++)
        {
            if (fabs(min[i] - max[i]) > lim)
                return false;
        }

        return true;
    }

    CagmRKF45Vect& operator=(const CagmRKF45Vect& v)
    {
        delete [] e;
        n = v.GetLength();
        if (n > 0)
            e = new double[n];
        for (int i = 0; i < n; i++)
            e[i] = v[i];

        return *this;
    }

    CagmRKF45Vect& operator=(double *v)
    {
        for (int i = 0; i < n; i++)
            e[i] = v[i];

        return *this;
    }

    CagmRKF45Vect& operator+=(const CagmRKF45Vect& v)
    {
        for (int i = 0; i < n; i++)
            e[i] += v[i];

        return *this;
    }

    CagmRKF45Vect& operator*=(const double d)
    {
        for (int i = 0; i < n; i++)
            e[i] *= d;

        return *this;
    }

    int GetLength() const { return n; }

    double& operator[](const int i) const
    {
        // exceptions?
        return e[i];
    }

    const double* v() // not safe!
    {
        return e;
    }
};
