
#ifdef IPHREEQC_NO_FORTRAN_MODULE
#include <ctype.h>   /* isgraph */
#include <stdlib.h>  /* malloc */
#include <memory.h>  /* memcpy */
#include <assert.h>  /* assert */
#include <stdio.h>   /* sprintf */
#include "phrqtype.h"
#include "IPhreeqc.h"

#include "fwrap.h"

char *
f2cstring(char* fstring, size_t len)
{
    char *cstr, *str;
    long i;

    str = fstring;
    for (i = (long) len - 1; i >= 0 && !isgraph((int) str[i]); i--);
    cstr = (char *) malloc((size_t) (i + 2));
    if (!cstr) return 0;

    cstr[i + 1] = '\0';
    if ((i + 1) > 0) memcpy(cstr, str, (size_t) (i + 1));
    return cstr;
}

void
padfstring(char *dest, const char *src, unsigned int len)
{
    size_t sofar;

    for (sofar = 0; (sofar < len) && (*src != '\0'); ++sofar)
        *dest++ = *src++;

    while (sofar++ < len)
        *dest++ = ' ';
}

IPQ_RESULT
AccumulateLineF(int *id, char *line, size_t line_length)
{
	IPQ_RESULT n;
	char* cline;

	cline = f2cstring(line, line_length);
	if (!cline)
	{
		::AddError(*id, "AccumulateLine: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	n = ::AccumulateLine(*id, cline);
	free(cline);
	return n;
}

int
AddErrorF(int *id, char *error_msg, size_t len)
{
	int n;
	char* cmsg;

	cmsg = f2cstring(error_msg, len);
	if (!cmsg)
	{
		::AddError(*id, "AddError: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	n = ::AddError(*id, cmsg);
	free(cmsg);
	return n;
}

int
AddWarningF(int *id, char *warn_msg, size_t len)
{
	int n;
	char* cmsg;

	cmsg = f2cstring(warn_msg, len);
	if (!cmsg)
	{
		::AddError(*id, "AddWarning: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	n = ::AddWarning(*id, cmsg);
	free(cmsg);
	return n;
}

IPQ_RESULT
ClearAccumulatedLinesF(int *id)
{
	return ::ClearAccumulatedLines(*id);
}

int
CreateIPhreeqcF(void)
{
	return ::CreateIPhreeqc();
}

int
DestroyIPhreeqcF(int *id)
{
	return ::DestroyIPhreeqc(*id);
}

int
GetComponentCountF(int *id)
{
	return ::GetComponentCount(*id);
}

void
GetComponentF(int *id, int *n, char* comp, size_t line_length)
{
	padfstring(comp, ::GetComponent(*id, (*n) - 1), (unsigned int) line_length);
}

int
GetCurrentSelectedOutputUserNumberF(int *id)
{
	return ::GetCurrentSelectedOutputUserNumber(*id);
}

void
GetDumpFileNameF(int *id, char* fname, size_t fname_length)
{
	padfstring(fname, ::GetDumpFileName(*id), (unsigned int) fname_length);
}

int
GetDumpFileOnF(int *id)
{
	return ::GetDumpFileOn(*id);
}

/*
GetDumpStringF
*/

int
GetDumpStringLineCountF(int *id)
{
	return ::GetDumpStringLineCount(*id);
}

void
GetDumpStringLineF(int *id, int* n, char* line, size_t line_length)
{
	padfstring(line, ::GetDumpStringLine(*id, (*n) - 1), (unsigned int) line_length);
}

int
GetDumpStringOnF(int *id)
{
	return ::GetDumpStringOn(*id);
}

void
GetErrorFileNameF(int *id, char* fname, size_t fname_length)
{
	padfstring(fname, ::GetErrorFileName(*id), (unsigned int) fname_length);
}

int
GetErrorFileOnF(int *id)
{
	return ::GetErrorFileOn(*id);
}

int
GetErrorOnF(int *id)
{
	return ::GetErrorOn(*id);
}

/*
GetErrorStringF
*/

int
GetErrorStringLineCountF(int *id)
{
	return ::GetErrorStringLineCount(*id);
}

void
GetErrorStringLineF(int *id, int* n, char* line, size_t line_length)
{
	padfstring(line, ::GetErrorStringLine(*id, (*n) - 1), (unsigned int) line_length);
}

int
GetErrorStringOnF(int *id)
{
	return ::GetErrorStringOn(*id);
}

void 
GetLogFileNameF(int *id, char* fname, size_t fname_length)
{
	padfstring(fname, ::GetLogFileName(*id), (unsigned int) fname_length);
}

int
GetLogFileOnF(int *id)
{
	return ::GetLogFileOn(*id);
}

int
GetLogStringLineCountF(int *id)
{
	return ::GetLogStringLineCount(*id);
}

void
GetLogStringLineF(int *id, int* n, char* line, size_t line_length)
{
	padfstring(line, ::GetLogStringLine(*id, (*n) - 1), (unsigned int) line_length);
}

int
GetLogStringOnF(int *id)
{
	return ::GetLogStringOn(*id);
}

int
GetNthSelectedOutputUserNumberF(int *id, int* n)
{
	return ::GetNthSelectedOutputUserNumber(*id, (*n) - 1);
}

void
GetOutputFileNameF(int *id, char* fname, size_t fname_length)
{
	padfstring(fname, ::GetOutputFileName(*id), (unsigned int) fname_length);
}

int
GetOutputStringLineCountF(int *id)
{
	return ::GetOutputStringLineCount(*id);
}

void
GetOutputStringLineF(int *id, int* n, char* line, size_t line_length)
{
	padfstring(line, ::GetOutputStringLine(*id, (*n) - 1), (unsigned int) line_length);
}

int
GetOutputStringOnF(int *id)
{
	return ::GetOutputStringOn(*id);
}

int
GetOutputFileOnF(int *id)
{
	return ::GetOutputFileOn(*id);
}

int
GetSelectedOutputColumnCountF(int *id)
{
	return ::GetSelectedOutputColumnCount(*id);
}

int
GetSelectedOutputCountF(int *id)
{
	return ::GetSelectedOutputCount(*id);
}

void
GetSelectedOutputFileNameF(int *id, char* fname, size_t fname_length)
{
	padfstring(fname, ::GetSelectedOutputFileName(*id), (unsigned int) fname_length);
}

int
GetSelectedOutputFileOnF(int *id)
{
	return ::GetSelectedOutputFileOn(*id);
}

/*
GetSelectedOutputStringF
*/

int
GetSelectedOutputStringLineCountF(int *id)
{
	return ::GetSelectedOutputStringLineCount(*id);
}

void
GetSelectedOutputStringLineF(int *id, int* n, char* line, size_t line_length)
{
	padfstring(line, ::GetSelectedOutputStringLine(*id, (*n) - 1), (unsigned int) line_length);
}

int
GetSelectedOutputStringOnF(int *id)
{
	return ::GetSelectedOutputStringOn(*id);
}

int
GetSelectedOutputRowCountF(int *id)
{
	int rows = ::GetSelectedOutputRowCount(*id);
	if (rows > 0)
	{
		rows -= 1;
	}
	return rows;
}

IPQ_RESULT
GetSelectedOutputValueF(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, size_t svalue_length)
{
	IPQ_RESULT result;
	VAR v;
	VarInit(&v);
	char buffer[100];

	int adjcol = *col - 1;
	result = ::GetSelectedOutputValue(*id, *row, adjcol, &v);

	switch (v.type)
	{
	case TT_EMPTY:
		*vtype = v.type;
		break;
	case TT_ERROR:
		*vtype = v.type;
		break;
	case TT_LONG:
		*vtype = TT_DOUBLE;
		*dvalue = (double)v.lVal;
		::snprintf(buffer, sizeof(buffer), "%ld", v.lVal);
		padfstring(svalue, buffer, (unsigned int) svalue_length);
		break;
	case TT_DOUBLE:
		*vtype = v.type;
		*dvalue = v.dVal;
		::snprintf(buffer, sizeof(buffer), "%23.15e", v.dVal);
		padfstring(svalue, buffer, (unsigned int) svalue_length);
		break;
	case TT_STRING:
		*vtype = v.type;
		padfstring(svalue, v.sVal, (unsigned int) svalue_length);
		break;
	default:
		assert(0);
	}
	::VarClear(&v);
	return result;
}

void
GetVersionStringF(char* version, size_t version_length)
{
	padfstring(version, ::GetVersionString(), (unsigned int) version_length);
}

/*
GetWarningStringF
*/

int
GetWarningStringLineCountF(int *id)
{
	return ::GetWarningStringLineCount(*id);
}

void
GetWarningStringLineF(int *id, int* n, char* line, size_t line_length)
{
	padfstring(line, ::GetWarningStringLine(*id, (*n) - 1), (unsigned int) line_length);
}

int
LoadDatabaseF(int *id, char* filename, size_t filename_length)
{
	char* cfilename;

	cfilename = f2cstring(filename, filename_length);
	if (!cfilename)
	{
		::AddError(*id, "LoadDatabase: Out of memory.\n");
		return VR_OUTOFMEMORY;
	}

	int n = ::LoadDatabase(*id, cfilename);
	free(cfilename);
	return n;
}

int
LoadDatabaseStringF(int *id, char* input, size_t input_length)
{
	char* cinput;

	cinput = f2cstring(input, input_length);
	if (!cinput)
	{
		::AddError(*id, "LoadDatabaseString: Out of memory.\n");
		return VR_OUTOFMEMORY;
	}

	int n = ::LoadDatabaseString(*id, cinput);
	free(cinput);
	return n;
}

void
OutputAccumulatedLinesF(int *id)
{
	::OutputAccumulatedLines(*id);
}

void
OutputErrorStringF(int *id)
{
	::OutputErrorString(*id);
}

void
OutputWarningStringF(int *id)
{
	::OutputWarningString(*id);
}

int
RunAccumulatedF(int *id)
{
	return ::RunAccumulated(*id);
}

int
RunFileF(int *id, char* filename, size_t filename_length)
{
	char* cfilename;

	cfilename = f2cstring(filename, filename_length);
	if (!cfilename)
	{
		::AddError(*id, "RunFile: Out of memory.\n");
		return (int)VR_OUTOFMEMORY;
	}

	int n = ::RunFile(*id, cfilename);
	free(cfilename);
	return n;
}

int
RunStringF(int *id, char* input, size_t input_length)
{
	char* cinput;

	cinput = f2cstring(input, input_length);
	if (!cinput)
	{
		::AddError(*id, "RunString: Out of memory.\n");
		return (int)VR_OUTOFMEMORY;
	}

	int n = ::RunString(*id, cinput);
	free(cinput);
	return n;
}

IPQ_RESULT
SetBasicFortranCallbackF(int *id, double (*fcn)(double *x1, double *x2, const char *str, size_t l))
{
	return ::SetBasicFortranCallback(*id, fcn);
}

IPQ_RESULT
SetCurrentSelectedOutputUserNumberF(int *id, int *n)
{
	return ::SetCurrentSelectedOutputUserNumber(*id, *n);
}

IPQ_RESULT
SetDumpFileNameF(int *id, char* fname, size_t fname_length)
{
	char* cinput;

	cinput = f2cstring(fname, fname_length);
	if (!cinput)
	{
		::AddError(*id, "SetDumpFileName: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	IPQ_RESULT n = ::SetDumpFileName(*id, cinput);
	free(cinput);
	return n;
}

IPQ_RESULT
SetDumpFileOnF(int *id, int* dump_on)
{
	return ::SetDumpFileOn(*id, *dump_on);
}

IPQ_RESULT
SetDumpStringOnF(int *id, int* dump_string_on)
{
	return ::SetDumpStringOn(*id, *dump_string_on);
}

IPQ_RESULT
SetErrorFileNameF(int *id, char* fname, size_t fname_length)
{
	char* cinput;

	cinput = f2cstring(fname, fname_length);
	if (!cinput)
	{
		::AddError(*id, "SetErrorFileName: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	IPQ_RESULT n = ::SetErrorFileName(*id, cinput);
	free(cinput);
	return n;
}

IPQ_RESULT
SetErrorFileOnF(int *id, int* error_file_on)
{
	return ::SetErrorFileOn(*id, *error_file_on);
}

IPQ_RESULT
SetErrorOnF(int *id, int* error_on)
{
	return ::SetErrorOn(*id, *error_on);
}

IPQ_RESULT
SetErrorStringOnF(int *id, int* error_string_on)
{
	return ::SetErrorStringOn(*id, *error_string_on);
}

IPQ_RESULT
SetLogFileNameF(int *id, char* fname, size_t fname_length)
{
	char* cinput;

	cinput = f2cstring(fname, fname_length);
	if (!cinput)
	{
		::AddError(*id, "SetLogFileName: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	IPQ_RESULT n = ::SetLogFileName(*id, cinput);
	free(cinput);
	return n;
}

IPQ_RESULT
SetLogFileOnF(int *id, int* log_on)
{
	return ::SetLogFileOn(*id, *log_on);
}

IPQ_RESULT
SetLogStringOnF(int *id, int* log_string_on)
{
	return ::SetLogStringOn(*id, *log_string_on);
}

IPQ_RESULT
SetOutputFileNameF(int *id, char* fname, size_t fname_length)
{
	char* cinput;

	cinput = f2cstring(fname, fname_length);
	if (!cinput)
	{
		::AddError(*id, "SetOutputFileName: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	IPQ_RESULT n = ::SetOutputFileName(*id, cinput);
	free(cinput);
	return n;
}

IPQ_RESULT
SetOutputFileOnF(int *id, int* output_on)
{
	return ::SetOutputFileOn(*id, *output_on);
}

IPQ_RESULT
SetOutputStringOnF(int *id, int* output_string_on)
{
	return ::SetOutputStringOn(*id, *output_string_on);
}

IPQ_RESULT
SetSelectedOutputFileNameF(int *id, char* fname, size_t fname_length)
{
	char* cinput;

	cinput = f2cstring(fname, fname_length);
	if (!cinput)
	{
		::AddError(*id, "SetSelectedOutputFileName: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	IPQ_RESULT n = ::SetSelectedOutputFileName(*id, cinput);
	free(cinput);
	return n;
}

IPQ_RESULT
SetSelectedOutputFileOnF(int *id, int* sel_on)
{
	return ::SetSelectedOutputFileOn(*id, *sel_on);
}

IPQ_RESULT
SetSelectedOutputStringOnF(int *id, int* selected_output_string_on)
{
	return ::SetSelectedOutputStringOn(*id, *selected_output_string_on);
}
#endif