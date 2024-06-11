#pragma once

#ifdef _WINDOWS

#include "binUtilities.h"
#include "agmScalarField.h"
#include "agmVectorField.h"

#ifdef _WINDOWS
#pragma warning(disable:4996)
#endif

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
        CbinDataStruct _data;
        v->GetData(&_data);
        _data.Write(fid);
    }

    void write_scalar(CagmScalarField *v)
    {
        CbinDataStruct _data;
        v->GetData(&_data);
        _data.Write(fid);
    }

    void write(CubeXD *v)
    {
        if (v->dim() == 1)
            write_scalar((CagmScalarField *)v);
        else
            write_vector((CagmVectorField *)v);
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
