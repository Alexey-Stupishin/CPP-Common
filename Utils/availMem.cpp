#include "stdDefinitions.h"
#include "availMem.h"

#ifdef _WINDOWS
#include <windows.h>
#include "psapi.h"
#pragma warning (disable : 4996)
#endif

uint64_t processMemory(char *buffer)
{
#ifdef _WINDOWS
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
    SIZE_T virtualMemUsedByMe = pmc.PrivateUsage;

    if (buffer)
        sprintf(buffer, "%7.3f", (double)virtualMemUsedByMe / 1024 / 1024 / 1024);

    return virtualMemUsedByMe;
#else
    return 0;
#endif
}

uint64_t availableMemory(char *buffer)
{
#ifdef _WINDOWS
    MEMORYSTATUSEX memstatus;
    memstatus.dwLength = sizeof(MEMORYSTATUSEX);
    GlobalMemoryStatusEx(&memstatus);

    if (buffer)
    {
        double avail = (double)memstatus.ullAvailPhys / 1024 / 1024 / 1024;
        sprintf(buffer, "%7.3f", avail);
    }

    return memstatus.ullAvailPhys;
#else
    return 0;
#endif
}

    //Linux:
    //https://forums.raspberrypi.com/viewtopic.php?t=235578
    //
    //const std::string info[] = {"Cached:", "Buffers:", "MemFree:", "MemTotal:"};
    //int intInfo[4];
    // 
    //std::string token;
    //std::ifstream file("/proc/meminfo");
 
    //while(file >> token){

    //    for(short i = 4; i != -1; --i)
    //    {   
    //        if(token == info[i])
    //        {   
    //            file >> intInfo[i];
    //            
    //            if(i == 0){
    //                file.close(); 
    //                return std::to_string((intInfo[3] - intInfo[2] - (intInfo[1] + intInfo[0])) / 1024) + "m ";
    //            }

    //        }
    //    }

    //   // ignore rest of the line
    //    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');   
    //}
