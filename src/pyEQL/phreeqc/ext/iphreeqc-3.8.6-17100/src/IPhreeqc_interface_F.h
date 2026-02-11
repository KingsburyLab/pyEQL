#ifndef __IPHREEQC_INTERFACE__H
#define __IPHREEQC_INTERFACE__H

#include "PHRQ_exports.h"

#ifdef SKIP
#if defined(FC_FUNC)
#define AccumulateLineF                     FC_FUNC (accumulatelinef,                     ACCUMULATELINEF)
#define AddErrorF                           FC_FUNC (adderrorf,                           ADDERRORF)
#define AddWarningF                         FC_FUNC (addwarningf,                         ADDWARNINGF)
#define ClearAccumulatedLinesF              FC_FUNC (clearaccumulatedlinesf,              CLEARACCUMULATEDLINESF)
#define CreateIPhreeqcF                     FC_FUNC (createiphreeqcf,                     CREATEIPHREEQCF)
#define DestroyIPhreeqcF                    FC_FUNC (destroyiphreeqcf,                    DESTROYIPHREEQCF)
#define GetComponentF                       FC_FUNC (getcomponentf,                       GETCOMPONENTF)
#define GetComponentCountF                  FC_FUNC (getcomponentcountf,                  GETCOMPONENTCOUNTF)
#define GetCurrentSelectedOutputUserNumberF FC_FUNC (getcurrentselectedoutputusernumberf, GETCURRENTSELECTEDOUTPUTUSERNUMBERF)
#define GetDumpFileNameF                    FC_FUNC (getdumpfilenamef,                    GETDUMPFILENAMEF)
#define GetDumpFileOnF                      FC_FUNC (getdumpfileonf,                      GETDUMPFILEONF)
#define GetDumpStringLineF                  FC_FUNC (getdumpstringlinef,                  GETDUMPSTRINGLINEF)
#define GetDumpStringLineCountF             FC_FUNC (getdumpstringlinecountf,             GETDUMPSTRINGLINECOUNTF)
#define GetDumpStringOnF                    FC_FUNC (getdumpstringonf,                    GETDUMPSTRINGONF)
#define GetErrorFileNameF                   FC_FUNC (geterrorfilenamef,                   GETERRORFILENAMEF)
#define GetErrorFileOnF                     FC_FUNC (geterrorfileonf,                     GETERRORFILEONF)
#define GetErrorOnF                         FC_FUNC (geterroronf,                         GETERRORONF)
#define GetErrorStringLineF                 FC_FUNC (geterrorstringlinef,                 GETERRORSTRINGLINEF)
#define GetErrorStringLineCountF            FC_FUNC (geterrorstringlinecountf,            GETERRORSTRINGLINECOUNTF)
#define GetErrorStringOnF                   FC_FUNC (geterrorstringonf,                   GETERRORSTRINGONF)
#define GetLogFileNameF                     FC_FUNC (getlogfilenamef,                     GETLOGFILENAMEF)
#define GetLogFileOnF                       FC_FUNC (getlogfileonf,                       GETLOGFILEONF)
#define GetLogStringLineF                   FC_FUNC (getlogstringlinef,                   GETLOGSTRINGLINEF)
#define GetLogStringLineCountF              FC_FUNC (getlogstringlinecountf,              GETLOGSTRINGLINECOUNTF)
#define GetLogStringOnF                     FC_FUNC (getlogstringonf,                     GETLOGSTRINGONF)
#define GetNthSelectedOutputUserNumberF     FC_FUNC (getnthselectedoutputusernumberf,     GETNTHSELECTEDOUTPUTUSERNUMBERF)
#define GetOutputFileNameF                  FC_FUNC (getoutputfilenamef,                  GETOUTPUTFILENAMEF)
#define GetOutputFileOnF                    FC_FUNC (getoutputfileonf,                    GETOUTPUTFILEONF)
#define GetOutputStringLineF                FC_FUNC (getoutputstringlinef,                GETOUTPUTSTRINGLINEF)
#define GetOutputStringLineCountF           FC_FUNC (getoutputstringlinecountf,           GETOUTPUTSTRINGLINECOUNTF)
#define GetOutputStringOnF                  FC_FUNC (getoutputstringonf,                  GETOUTPUTSTRINGONF)
#define GetSelectedOutputColumnCountF       FC_FUNC (getselectedoutputcolumncountf,       GETSELECTEDOUTPUTCOLUMNCOUNTF)
#define GetSelectedOutputCountF             FC_FUNC (getselectedoutputcountf,             GETSELECTEDOUTPUTCOUNTF)
#define GetSelectedOutputFileNameF          FC_FUNC (getselectedoutputfilenamef,          GETSELECTEDOUTPUTFILENAMEF)
#define GetSelectedOutputFileOnF            FC_FUNC (getselectedoutputfileonf,            GETSELECTEDOUTPUTFILEONF)
#define GetSelectedOutputRowCountF          FC_FUNC (getselectedoutputrowcountf,          GETSELECTEDOUTPUTROWCOUNTF)
#define GetSelectedOutputStringLineF        FC_FUNC (getselectedoutputstringlinef,        GETSELECTEDOUTPUTSTRINGLINEF)
#define GetSelectedOutputStringLineCountF   FC_FUNC (getselectedoutputstringlinecountf,   GETSELECTEDOUTPUTSTRINGLINECOUNTF)
#define GetSelectedOutputStringOnF          FC_FUNC (getselectedoutputstringonf,          GETSELECTEDOUTPUTSTRINGONF)
#define GetSelectedOutputValueF             FC_FUNC (getselectedoutputvaluef,             GETSELECTEDOUTPUTVALUEF)
#define GetVersionStringF                   FC_FUNC (getversionstringf,                   GETVERSIONSTRINGF)
#define GetWarningStringLineF               FC_FUNC (getwarningstringlinef,               GETWARNINGSTRINGLINEF)
#define GetWarningStringLineCountF          FC_FUNC (getwarningstringlinecountf,          GETWARNINGSTRINGLINECOUNTF)
#define LoadDatabaseF                       FC_FUNC (loaddatabasef,                       LOADDATABASEF)
#define LoadDatabaseStringF                 FC_FUNC (loaddatabasestringf,                 LOADDATABASESTRINGF)
#define OutputAccumulatedLinesF             FC_FUNC (outputaccumulatedlinesf,             OUTPUTACCUMULATEDLINESF)
#define OutputErrorStringF                  FC_FUNC (outputerrorstringf,                  OUTPUTERRORSTRINGF)
#define OutputWarningStringF                FC_FUNC (outputwarningstringf,                OUTPUTWARNINGSTRINGF)
#define RunAccumulatedF                     FC_FUNC (runaccumulatedf,                     RUNACCUMULATEDF)
#define RunFileF                            FC_FUNC (runfilef,                            RUNFILEF)
#define RunStringF                          FC_FUNC (runstringf,                          RUNSTRINGF)
#define SetBasicFortranCallbackF            FC_FUNC (setbasicfortrancallbackf,            SETFOTRANBASICCALLBACKF)
#define SetCurrentSelectedOutputUserNumberF FC_FUNC (setcurrentselectedoutputusernumberf, SETCURRENTSELECTEDOUTPUTUSERNUMBERF)
#define SetDumpFileNameF                    FC_FUNC (setdumpfilenamef,                    SETDUMPFILENAMEF)
#define SetDumpFileOnF                      FC_FUNC (setdumpfileonf,                      SETDUMPFILEONF)
#define SetDumpStringOnF                    FC_FUNC (setdumpstringonf,                    SETDUMPSTRINGONF)
#define SetErrorFileNameF                   FC_FUNC (seterrorfilenamef,                   SETERRORFILENAMEF)
#define SetErrorFileOnF                     FC_FUNC (seterrorfileonf,                     SETERRORFILEONF)
#define SetErrorOnF                         FC_FUNC (seterroronf,                         SETERRORONF)
#define SetErrorStringOnF                   FC_FUNC (seterrorstringonf,                   SETERRORSTRINGONF)
#define SetLogFileNameF                     FC_FUNC (setlogfilenamef,                     SETLOGFILENAMEF)
#define SetLogFileOnF                       FC_FUNC (setlogfileonf,                       SETLOGFILEONF)
#define SetLogStringOnF                     FC_FUNC (setlogstringonf,                     SETLOGSTRINGONF)
#define SetOutputFileNameF                  FC_FUNC (setoutputfilenamef,                  SETOUTPUTFILENAMEF)
#define SetOutputFileOnF                    FC_FUNC (setoutputfileonf,                    SETOUTPUTFILEONF)
#define SetOutputStringOnF                  FC_FUNC (setoutputstringonf,                  SETOUTPUTSTRINGONF)
#define SetSelectedOutputFileNameF          FC_FUNC (setselectedoutputfilenamef,          SETSELECTEDOUTPUTFILENAMEF)
#define SetSelectedOutputFileOnF            FC_FUNC (setselectedoutputfileonf,            SETSELECTEDOUTPUTFILEONF)
#define SetSelectedOutputStringOnF          FC_FUNC (setselectedoutputstringonf,          SETSELECTEDOUTPUTSTRINGONF)
#endif /* FC_FUNC */
#endif

#if defined(__cplusplus)
extern "C" {
#endif

  IPQ_DLL_EXPORT IPQ_RESULT AccumulateLineF(int *id, char *line);
  IPQ_DLL_EXPORT int        AddErrorF(int *id, char *error_msg);
  IPQ_DLL_EXPORT int        AddWarningF(int *id, char *warn_msg);
  IPQ_DLL_EXPORT IPQ_RESULT ClearAccumulatedLinesF(int *id);
  IPQ_DLL_EXPORT int        CreateIPhreeqcF(void);
  IPQ_DLL_EXPORT int        DestroyIPhreeqcF(int *id);
  IPQ_DLL_EXPORT void       GetComponentF(int *id, int* n, char* line, int* line_length);
  IPQ_DLL_EXPORT int        GetComponentCountF(int *id);
  IPQ_DLL_EXPORT int        GetCurrentSelectedOutputUserNumberF(int *id);
  IPQ_DLL_EXPORT void       GetDumpFileNameF(int *id, char* filename, int* filename_length);
  IPQ_DLL_EXPORT int        GetDumpFileOnF(int *id);
  IPQ_DLL_EXPORT void       GetDumpStringLineF(int *id, int* n, char* line, int* line_length);
  IPQ_DLL_EXPORT int        GetDumpStringLineCountF(int *id);
  IPQ_DLL_EXPORT int        GetDumpStringOnF(int *id);
  IPQ_DLL_EXPORT void       GetErrorFileNameF(int *id, char* filename, int* filename_length);
  IPQ_DLL_EXPORT int        GetErrorFileOnF(int *id);
  IPQ_DLL_EXPORT int        GetErrorOnF(int *id);
  IPQ_DLL_EXPORT void       GetErrorStringLineF(int *id, int* n, char* line, int* line_length);
  IPQ_DLL_EXPORT int        GetErrorStringLineCountF(int *id);
  IPQ_DLL_EXPORT int        GetErrorStringOnF(int *id);
  IPQ_DLL_EXPORT void       GetLogFileNameF(int *id, char* filename, int* filename_length);
  IPQ_DLL_EXPORT int        GetLogFileOnF(int *id);
  IPQ_DLL_EXPORT void       GetLogStringLineF(int *id, int* n, char* line, int* line_length);
  IPQ_DLL_EXPORT int        GetLogStringLineCountF(int *id);
  IPQ_DLL_EXPORT int        GetLogStringOnF(int *id);
  IPQ_DLL_EXPORT int        GetNthSelectedOutputUserNumberF(int *id, int* n);
  IPQ_DLL_EXPORT void       GetOutputFileNameF(int *id, char* filename, int* filename_length);
  IPQ_DLL_EXPORT int        GetOutputFileOnF(int *id);
  IPQ_DLL_EXPORT void       GetOutputStringLineF(int *id, int* n, char* line, int* line_length);
  IPQ_DLL_EXPORT int        GetOutputStringLineCountF(int *id);
  IPQ_DLL_EXPORT int        GetOutputStringOnF(int *id);
  IPQ_DLL_EXPORT int        GetSelectedOutputColumnCountF(int *id);
  IPQ_DLL_EXPORT int        GetSelectedOutputCountF(int *id);
  IPQ_DLL_EXPORT void       GetSelectedOutputFileNameF(int *id, char* filename, int* filename_length);
  IPQ_DLL_EXPORT int        GetSelectedOutputFileOnF(int *id);
  IPQ_DLL_EXPORT int        GetSelectedOutputRowCountF(int *id);
  IPQ_DLL_EXPORT void       GetSelectedOutputStringLineF(int *id, int* n, char* line, int* line_length);
  IPQ_DLL_EXPORT int        GetSelectedOutputStringLineCountF(int *id);
  IPQ_DLL_EXPORT int        GetSelectedOutputStringOnF(int *id);
  IPQ_DLL_EXPORT IPQ_RESULT GetSelectedOutputValueF(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, int* svalue_length);
  IPQ_DLL_EXPORT void       GetVersionStringF(char* version, int* version_length);
  IPQ_DLL_EXPORT void       GetWarningStringLineF(int *id, int* n, char* line, int* line_length);
  IPQ_DLL_EXPORT int        GetWarningStringLineCountF(int *id);
  IPQ_DLL_EXPORT int        LoadDatabaseF(int *id, char* filename);
  IPQ_DLL_EXPORT int        LoadDatabaseStringF(int *id, char* input);
  IPQ_DLL_EXPORT void       OutputAccumulatedLinesF(int *id);
  IPQ_DLL_EXPORT void       OutputErrorStringF(int *id);
  IPQ_DLL_EXPORT void       OutputWarningStringF(int *id);
  IPQ_DLL_EXPORT int        RunAccumulatedF(int *id);
  IPQ_DLL_EXPORT int        RunFileF(int *id, char* filename);
  IPQ_DLL_EXPORT int        RunStringF(int *id, char* input);
#ifdef IPHREEQC_NO_FORTRAN_MODULE
  IPQ_DLL_EXPORT IPQ_RESULT SetBasicFortranCallbackF(int *id, double (*fcn)(double *x1, double *x2, const char *str, size_t l));
#else
  IPQ_DLL_EXPORT IPQ_RESULT SetBasicFortranCallbackF(int *id, double (*fcn)(double *x1, double *x2, const char *str, int l));
#endif
  IPQ_DLL_EXPORT IPQ_RESULT SetCurrentSelectedOutputUserNumberF(int *id, int *n);
  IPQ_DLL_EXPORT IPQ_RESULT SetDumpFileNameF(int *id, char* fname);
  IPQ_DLL_EXPORT IPQ_RESULT SetDumpFileOnF(int *id, int* dump_on);
  IPQ_DLL_EXPORT IPQ_RESULT SetDumpStringOnF(int *id, int* dump_string_on);
  IPQ_DLL_EXPORT IPQ_RESULT SetErrorFileNameF(int *id, char* fname);
  IPQ_DLL_EXPORT IPQ_RESULT SetErrorFileOnF(int *id, int* error_file_on);
  IPQ_DLL_EXPORT IPQ_RESULT SetErrorOnF(int *id, int* error_on);
  IPQ_DLL_EXPORT IPQ_RESULT SetErrorStringOnF(int *id, int* error_string_on);
  IPQ_DLL_EXPORT IPQ_RESULT SetLogFileNameF(int *id, char* fname);
  IPQ_DLL_EXPORT IPQ_RESULT SetLogFileOnF(int *id, int* log_file_on);
  IPQ_DLL_EXPORT IPQ_RESULT SetLogStringOnF(int *id, int* log_string_on);
  IPQ_DLL_EXPORT IPQ_RESULT SetOutputFileNameF(int *id, char* fname);
  IPQ_DLL_EXPORT IPQ_RESULT SetOutputFileOnF(int *id, int* output_on);
  IPQ_DLL_EXPORT IPQ_RESULT SetOutputStringOnF(int *id, int* output_string_on);
  IPQ_DLL_EXPORT IPQ_RESULT SetSelectedOutputFileNameF(int *id, char* fname);
  IPQ_DLL_EXPORT IPQ_RESULT SetSelectedOutputFileOnF(int *id, int* selected_output_file_on);
  IPQ_DLL_EXPORT IPQ_RESULT SetSelectedOutputStringOnF(int *id, int* selected_output_string_on);

#if defined(__cplusplus)
}
#endif
void padfstring(char *dest, const char *src, unsigned int len);
#endif  /* __IPHREEQC_INTERFACE__H */
