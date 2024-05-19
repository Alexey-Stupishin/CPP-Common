#pragma once

#ifdef _WINDOWS

#include "binUtilities.h"
#include "agmScalarField.h"
#include "agmVectorField.h"

class debug_data_trace_win : public CbinDataStruct
{
protected:
    FILE * fid;

public:
    debug_data_trace_win(char * filename)
        : CbinDataStruct()
    {
        fid = fopen(filename, "wb");
        CbinDataStruct::WriteHeader(fid);
    }

    virtual ~debug_data_trace_win()
    {
        CbinDataStruct::WriteFooter(fid);
        fclose(fid);
    }

    void write_vector(CagmVectorField *v)
    {
        CbinDataStruct data;
        v->GetData(&data);
        data.Write(fid);
    }

    void write_scalar(CagmScalarField *v)
    {
        CbinDataStruct data;
        v->GetData(&data);
        data.Write(fid);
    }
};

class debug_state_trace_win : public CbinDataStruct
{
protected:
    FILE * fid;

public:
    debug_state_trace_win(char * filename)
        : CbinDataStruct()
    {
        fid = fopen(filename, "wb");
        CbinDataStruct::WriteHeader(fid);
    }

    virtual ~debug_state_trace_win()
    {
        CbinDataStruct::WriteFooter(fid);
        fclose(fid);
    }
};

#endif
