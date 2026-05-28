#pragma once

#ifdef _WINDOWS
#pragma warning (disable : 4996)
#endif

size_t formatThouthands(uint64_t n, char *buffer)
{
    buffer[0] = 0;
    char buf_t[32];
    buf_t[0] = 0;
    do
    {
        int tri = n % 1000;
        if (n < 1000)
            sprintf(buffer, "%3d", tri);
        else
            sprintf(buffer, "%03d", tri);
        strcat(buffer, " ");
        strcat(buffer, buf_t);
        strcpy(buf_t, buffer);

        n /= 1000;
    } while (n > 0);

    return strlen(buffer);
}
