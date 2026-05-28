#pragma once

#ifdef _WINDOWS

#include <Windows.h>
#include "baseCritSect.h"
#include <stdio.h>
#include "TimeTicToc.h"

class BWndDebug
{
protected:
    FILE *m_file;
    char m_sFileProp[8];
    BOOL m_bFlush;
    HWND m_hWndListBox;
    BOOL m_bLastString;
    BOOL m_bTimeStamp;

    aguTimeTicToc tic;

    char s[1024];

    CbaseCritSect cs;

    LARGE_INTEGER frequency_time, timeBegin;

public:
    BWndDebug(LPCSTR _lpsFileName = "", LPCSTR _sFileProp = "wb", BOOL _bFlush = TRUE, 
                         HWND _hWnd = NULL, BOOL _bLastString = TRUE, BOOL _bTimeStamp = FALSE);
    virtual ~BWndDebug();
    void Init(LPCSTR _lpsFileName = "", LPCSTR _sFileProp = "wb", BOOL _bFlush = TRUE, 
                         HWND _hWnd = NULL, BOOL _bLastString = TRUE, BOOL _bTimeStamp = FALSE);
    void Close();
    void Say(char *format, ...);
    void memory(char *where);
};

#endif

