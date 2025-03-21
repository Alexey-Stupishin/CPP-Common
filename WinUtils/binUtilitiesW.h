#pragma once

#include "binUtilities.h"
#include "aslMap.h"

enum {MFO_DATA_FILE_OK = 0, MFO_DATA_FILE_NOT_EXIST = -99, MFO_DATA_FILE_BIN_WRONG_FORMAT = -100, MFO_DATA_FILE_TOO_SHORT = -1, MFO_DATA_FILE_BIN_READ_ERROR = -2, MFO_DATA_FILE_HDF_READ_ERROR = -3};
#define HDR_BUFSIZE 5

class CbinDataStructW : public CbinDataStruct
{
public:
    CbinDataStructW() {}
    virtual ~CbinDataStructW() {}

    int ReadFile(char *filename, aslMap<int, 64> *map);
    int Read(FILE *fid, aslMap<int, 64> *map);
    int ReadHDF(char *filename, aslMap<int, 64> *map);
};
