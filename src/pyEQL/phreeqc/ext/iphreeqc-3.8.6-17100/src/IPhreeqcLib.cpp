#include <cassert>
#include <iostream>
#include <map>

#include "IPhreeqc.h"
#include "IPhreeqc.hpp"
#include "thread.h"

class IPhreeqcLib
{
public:
	//static void CleanupIPhreeqcInstances(void);
	static int CreateIPhreeqc(void);
	static IPQ_RESULT DestroyIPhreeqc(int n);
	static IPhreeqc* GetInstance(int n);
};

IPQ_RESULT
AccumulateLine(int id, const char *line)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		switch (IPhreeqcPtr->AccumulateLine(line))
		{
		case VR_OK:
			return IPQ_OK;
		case VR_OUTOFMEMORY:
			return IPQ_OUTOFMEMORY;
		default:
			assert(false);
		}
	}
	return IPQ_BADINSTANCE;
}

int
AddError(int id, const char* error_msg)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return (int)IPhreeqcPtr->AddError(error_msg);
	}
	return IPQ_BADINSTANCE;
}

int
AddWarning(int id, const char* warn_msg)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return (int)IPhreeqcPtr->AddWarning(warn_msg);
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
ClearAccumulatedLines(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->ClearAccumulatedLines();
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

int
CreateIPhreeqc(void)
{
	return IPhreeqcLib::CreateIPhreeqc();
}

IPQ_RESULT
DestroyIPhreeqc(int id)
{
	return IPhreeqcLib::DestroyIPhreeqc(id);
}

// TODO Maybe GetAccumulatedLines

const char*
GetComponent(int id, int n)
{
	static const char err_msg[] = "GetComponent: Invalid instance id.\n";

	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetComponent(n);
	}
	return err_msg;
}

int
GetComponentCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return (int)IPhreeqcPtr->GetComponentCount();
	}
	return IPQ_BADINSTANCE;
}

int
GetCurrentSelectedOutputUserNumber(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		int n = IPhreeqcPtr->GetCurrentSelectedOutputUserNumber();
		switch (n)
		{
		case VR_INVALIDARG:
			return IPQ_INVALIDARG;
		}
		return n;
	}
	return IPQ_BADINSTANCE;
}

const char*
GetDumpFileName(int id)
{
	static const char empty[] = "";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetDumpFileName();
	}
	return empty;
}

int
GetDumpFileOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetDumpFileOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

const char*
GetDumpString(int id)
{
	static const char empty[] = "";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetDumpString();
	}
	return empty;
}

const char*
GetDumpStringLine(int id, int n)
{
	static const char err_msg[] = "GetDumpStringLine: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetDumpStringLine(n);
	}
	return err_msg;
}

int
GetDumpStringLineCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetDumpStringLineCount();
	}
	return 0;
}

int
GetDumpStringOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetDumpStringOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

const char*
GetErrorFileName(int id)
{
	static const char empty[] = "";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetErrorFileName();
	}
	return empty;
}

int
GetErrorFileOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetErrorFileOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

int
GetErrorOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetErrorOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

const char*
GetErrorString(int id)
{
	static const char err_msg[] = "GetErrorString: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetErrorString();
	}
	return err_msg;
}

const char*
GetErrorStringLine(int id, int n)
{
	static const char err_msg[] = "GetErrorStringLine: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetErrorStringLine(n);
	}
	return err_msg;
}

int
GetErrorStringLineCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return (int)IPhreeqcPtr->GetErrorStringLineCount();
	}
	return IPQ_BADINSTANCE;
}

int
GetErrorStringOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetErrorStringOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

const char*
GetLogFileName(int id)
{
	static const char empty[] = "";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetLogFileName();
	}
	return empty;
}

int
GetLogFileOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetLogFileOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

const char*
GetLogString(int id)
{
	static const char empty[] = "";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetLogString();
	}
	return empty;
}

const char*
GetLogStringLine(int id, int n)
{
	static const char err_msg[] = "GetLogStringLine: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetLogStringLine(n);
	}
	return err_msg;
}

int
GetLogStringLineCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetLogStringLineCount();
	}
	return 0;
}

int
GetLogStringOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetLogStringOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

int
GetNthSelectedOutputUserNumber(int id, int n)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		int nth = IPhreeqcPtr->GetNthSelectedOutputUserNumber(n);
		switch (nth)
		{
		case VR_INVALIDARG:  return IPQ_INVALIDARG;
		}
		return nth;
	}
	return IPQ_BADINSTANCE;
}

const char*
GetOutputFileName(int id)
{
	static const char empty[] = "";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetOutputFileName();
	}
	return empty;
}

int
GetOutputFileOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetOutputFileOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

const char*
GetOutputString(int id)
{
	static const char empty[] = "";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetOutputString();
	}
	return empty;
}

const char*
GetOutputStringLine(int id, int n)
{
	static const char err_msg[] = "GetOutputStringLine: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetOutputStringLine(n);
	}
	return err_msg;
}

int
GetOutputStringLineCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetOutputStringLineCount();
	}
	return 0;
}

int
GetOutputStringOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetOutputStringOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

int
GetSelectedOutputColumnCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetSelectedOutputColumnCount();
	}
	return IPQ_BADINSTANCE;
}

int
GetSelectedOutputCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetSelectedOutputCount();
	}
	return IPQ_BADINSTANCE;
}

const char*
GetSelectedOutputFileName(int id)
{
	static const char empty[] = "";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetSelectedOutputFileName();
	}
	return empty;
}

int
GetSelectedOutputFileOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetSelectedOutputFileOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

int
GetSelectedOutputRowCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetSelectedOutputRowCount();
	}
	return IPQ_BADINSTANCE;
}

const char*
GetSelectedOutputString(int id)
{
	static const char empty[] = "";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetSelectedOutputString();
	}
	return empty;
}

const char*
GetSelectedOutputStringLine(int id, int n)
{
	static const char err_msg[] = "GetSelectedOutputStringLine: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetSelectedOutputStringLine(n);
	}
	return err_msg;
}

int
GetSelectedOutputStringLineCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetSelectedOutputStringLineCount();
	}
	return 0;
}

int
GetSelectedOutputStringOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetSelectedOutputStringOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
GetSelectedOutputValue(int id, int row, int col, VAR* pVAR)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		switch (IPhreeqcPtr->GetSelectedOutputValue(row, col, pVAR))
		{
		case VR_OK:          return IPQ_OK;
		case VR_OUTOFMEMORY: return IPQ_OUTOFMEMORY;
		case VR_BADVARTYPE:  return IPQ_BADVARTYPE;
		case VR_INVALIDARG:  return IPQ_INVALIDARG;
		case VR_INVALIDROW:  return IPQ_INVALIDROW;
		case VR_INVALIDCOL:  return IPQ_INVALIDCOL;
		default:
			assert(false);
		}
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
GetSelectedOutputValue2(int id, int row, int col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		switch (IPhreeqcPtr->GetSelectedOutputValue2(row, col, vtype, dvalue, svalue, svalue_length))
		{
		case VR_OK:          return IPQ_OK;
		case VR_OUTOFMEMORY: return IPQ_OUTOFMEMORY;
		case VR_BADVARTYPE:  return IPQ_BADVARTYPE;
		case VR_INVALIDARG:  return IPQ_INVALIDARG;
		case VR_INVALIDROW:  return IPQ_INVALIDROW;
		case VR_INVALIDCOL:  return IPQ_INVALIDCOL;
		default:
			assert(false);
		}
	}
	return IPQ_BADINSTANCE;
}

const char*
GetVersionString(void)
{
	return IPhreeqc::GetVersionString();
}

const char*
GetWarningString(int id)
{
	static const char err_msg[] = "GetWarningString: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetWarningString();
	}
	return err_msg;
}

const char*
GetWarningStringLine(int id, int n)
{
	static const char err_msg[] = "GetWarningStringLine: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetWarningStringLine(n);
	}
	return err_msg;
}

int
GetWarningStringLineCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return (int)IPhreeqcPtr->GetWarningStringLineCount();
	}
	return IPQ_BADINSTANCE;
}

int
LoadDatabase(int id, const char* filename)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->LoadDatabase(filename);
	}
	return IPQ_BADINSTANCE;
}

int
LoadDatabaseString(int id, const char* input)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->LoadDatabaseString(input);
	}
	return IPQ_BADINSTANCE;
}

void
OutputAccumulatedLines(int id)
{
	static const char err_msg[] = "OutputAccumulatedLines: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->OutputAccumulatedLines();
		return;
	}
#if !defined(R_SO)
	std::cout << err_msg << std::endl;
#endif
}

void
OutputErrorString(int id)
{
	static const char err_msg[] = "OutputErrorString: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->OutputErrorString();
		return;
	}
#if !defined(R_SO)
	std::cout << err_msg << std::endl;
#endif
}

void
OutputWarningString(int id)
{
	static const char err_msg[] = "OutputWarningString: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->OutputWarningString();
		return;
	}
#if !defined(R_SO)
	std::cout << err_msg << std::endl;
#endif
}

int
RunAccumulated(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->RunAccumulated();
	}
	return IPQ_BADINSTANCE;
}

int
RunFile(int id, const char* filename)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->RunFile(filename);
	}
	return IPQ_BADINSTANCE;
}

int
RunString(int id, const char* input)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->RunString(input);
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetBasicCallback(int id, double (*fcn)(double x1, double x2, const char *str, void *cookie), void *cookie1)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetBasicCallback(fcn, cookie1);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}
#if !defined(R_SO)
#ifdef IPHREEQC_NO_FORTRAN_MODULE
IPQ_RESULT
SetBasicFortranCallback(int id, double (*fcn)(double *x1, double *x2, const char *str, size_t l))
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetBasicFortranCallback(fcn);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}
#else
IPQ_RESULT
SetBasicFortranCallback(int id, double (*fcn)(double *x1, double *x2, const char *str, int l))
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetBasicFortranCallback(fcn);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}
#endif /* IPHREEQC_NO_FORTRAN_MODULE */
#endif /* !defined(R_SO) */
IPQ_RESULT
SetCurrentSelectedOutputUserNumber(int id, int n)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		switch (IPhreeqcPtr->SetCurrentSelectedOutputUserNumber(n))
		{
		case VR_INVALIDARG: return IPQ_INVALIDARG;
		case VR_OK:         return IPQ_OK;
		default:            assert(false);
		}
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetDumpFileName(int id, const char* filename)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetDumpFileName(filename);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetDumpFileOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetDumpFileOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetDumpStringOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetDumpStringOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetErrorFileName(int id, const char* filename)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetErrorFileName(filename);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetErrorFileOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetErrorFileOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetErrorOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetErrorOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetErrorStringOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetErrorStringOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetLogFileName(int id, const char* filename)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetLogFileName(filename);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetLogFileOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetLogFileOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetLogStringOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetLogStringOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetOutputFileName(int id, const char* filename)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetOutputFileName(filename);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetOutputFileOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetOutputFileOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetOutputStringOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetOutputStringOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetSelectedOutputFileName(int id, const char* filename)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetSelectedOutputFileName(filename);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetSelectedOutputFileOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetSelectedOutputFileOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetSelectedOutputStringOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetSelectedOutputStringOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

// helper functions
//

int
IPhreeqcLib::CreateIPhreeqc(void)
{
	int n = IPQ_OUTOFMEMORY;
	IPhreeqc* IPhreeqcPtr;
	try
	{
		IPhreeqcPtr = new IPhreeqc;
		n = (int) IPhreeqcPtr->Index;
	}
	catch (const std::bad_alloc&)
	{
		return IPQ_OUTOFMEMORY;
	}
	return n;
}

IPQ_RESULT
IPhreeqcLib::DestroyIPhreeqc(int id)
{
	IPQ_RESULT retval = IPQ_BADINSTANCE;
	if (id >= 0)
	{
		if (IPhreeqc *ptr = IPhreeqcLib::GetInstance(id))
		{
			delete ptr;
			retval = IPQ_OK;
		}
	}
	return retval;
}

IPhreeqc*
IPhreeqcLib::GetInstance(int id)
{
	IPhreeqc* instance = 0;
	mutex_lock(&map_lock);
	std::map<size_t, IPhreeqc*>::iterator it = IPhreeqc::Instances.find(size_t(id));
	if (it != IPhreeqc::Instances.end())
	{
		instance = (*it).second;
	}
	mutex_unlock(&map_lock);
	return instance;
}
//// static method
//void IPhreeqcLib::CleanupIPhreeqcInstances(void)
//{
//	std::map<size_t, IPhreeqc*>::iterator it = IPhreeqc::Instances.begin();
//	std::vector<IPhreeqc*> ip_list;
//	for ( ; it != IPhreeqc::Instances.end(); it++)
//	{
//		ip_list.push_back(it->second);
//	}
//	for (size_t i = 0; i < ip_list.size(); i++)
//	{
//		delete ip_list[i];
//	}
//}
