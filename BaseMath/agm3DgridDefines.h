#pragma once

#define fidx(kx, ky, kz) (ky)+(kz)*N[1]][(kx)

#define fidx3(kx, ky, kz) (((ky) + (kz)*N[1])*N[0] + (kx))

#define l_div_position(p, N, k1, tk) \
{ \
    if (p >= double((N)-1) || fabs(p-(double((N)-1))) < 1e-5) \
    { \
        k1 = (N)-2; \
        tk = 1; \
    } \
    else \
    { \
        k1 = (int)floor(p); \
        if (k1 < 0) \
        { \
            k1 = 0;  \
            tk = 0; \
        } \
        else \
            tk = (p) - k1; \
    } \
}

#define l_get_point(Ny, tfield, x1, y1, z1, tx, ty, tz) \
    ((1-(tz))* ((1-(ty))* ((1-(tx))*tfield[(y1)     +  (z1)   *Ny][x1] + (tx)*tfield[(y1)     +  (z1)   *Ny][x1+1]) + \
                    (ty)* ((1-(tx))*tfield[((y1)+1) +  (z1)   *Ny][x1] + (tx)*tfield[((y1)+1) +  (z1)   *Ny][x1+1]))  \
   +     (tz)* ((1-(ty))* ((1-(tx))*tfield[(y1)     + ((z1)+1)*Ny][x1] + (tx)*tfield[(y1)     + ((z1)+1)*Ny][x1+1]) + \
                    (ty)* ((1-(tx))*tfield[((y1)+1) + ((z1)+1)*Ny][x1] + (tx)*tfield[((y1)+1) + ((z1)+1)*Ny][x1+1]))  \
    )
