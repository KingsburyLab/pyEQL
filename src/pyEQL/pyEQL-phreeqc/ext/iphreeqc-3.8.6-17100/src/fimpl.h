IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(accumulateline, ACCUMULATELINE, accumulateline_, ACCUMULATELINE_)(int *id, char *line, size_t len)
{
	return AccumulateLineF(id, line, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(adderror, ADDERROR, adderror_, ADDERROR_)(int *id, char *error_msg, size_t len)
{
	return AddErrorF(id, error_msg, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(addwarning, ADDWARNING, addwarning_, ADDWARNING_)(int *id, char *warn_msg, size_t len)
{
	return AddWarningF(id, warn_msg, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(clearaccumulatedlines, CLEARACCUMULATEDLINES, clearaccumulatedlines_, CLEARACCUMULATEDLINES_)(int *id)
{
	return ClearAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(createiphreeqc, CREATEIPHREEQC, createiphreeqc_, CREATEIPHREEQC_)(void)
{
	return CreateIPhreeqcF();
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(destroyiphreeqc, DESTROYIPHREEQC, destroyiphreeqc_, DESTROYIPHREEQC_)(int *id)
{
	return DestroyIPhreeqcF(id);
}
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(getcomponent, GETCOMPONENT, getcomponent_, GETCOMPONENT_)(int *id, int *n, char* line, size_t line_length)
{
	GetComponentF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getcomponentcount, GETCOMPONENTCOUNT, getcomponentcount_, GETCOMPONENTCOUNT_)(int *id)
{
	return GetComponentCountF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getcurrentselectedoutputusernumber, GETCURRENTSELECTEDOUTPUTUSERNUMBER, getcurrentselectedoutputusernumber_, GETCURRENTSELECTEDOUTPUTUSERNUMBER_)(int *id)
{
	return GetCurrentSelectedOutputUserNumberF(id);
}
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(getdumpfilename, GETDUMPFILENAME, getdumpfilename_, GETDUMPFILENAME_)(int *id, char *filename, size_t len)
{
	GetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getdumpfileon, GETDUMPFILEON, getdumpfileon_, GETDUMPFILEON_)(int *id)
{
	return GetDumpFileOnF(id);
}
// GetDumpString
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(getdumpstringline, GETDUMPSTRINGLINE, getdumpstringline_, GETDUMPSTRINGLINE_)(int *id, int *n, char* line, size_t line_length)
{
	GetDumpStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getdumpstringlinecount, GETDUMPSTRINGLINECOUNT, getdumpstringlinecount_, GETDUMPSTRINGLINECOUNT_)(int *id)
{
	return GetDumpStringLineCountF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getdumpstringon, GETDUMPSTRINGON, getdumpstringon_, GETDUMPSTRINGON_)(int *id)
{
	return GetDumpStringOnF(id);
}
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(geterrorfilename, GETERRORFILENAME, geterrorfilename_, GETERRORFILENAME_)(int *id, char *filename, size_t len)
{
	GetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(geterrorfileon, GETERRORFILEON, geterrorfileon_, GETERRORFILEON_)(int *id)
{
	return GetErrorFileOnF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(geterroron, GETERRORON, geterroron_, GETERRORON_)(int *id)
{
	return GetErrorOnF(id);
}
// GetErrorString
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(geterrorstringline, GETERRORSTRINGLINE, geterrorstringline_, GETERRORSTRINGLINE_)(int *id, int *n, char* line, size_t line_length)
{
	GetErrorStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(geterrorstringlinecount, GETERRORSTRINGLINECOUNT, geterrorstringlinecount_, GETERRORSTRINGLINECOUNT_)(int *id)
{
	return GetErrorStringLineCountF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(geterrorstringon, GETERRORSTRINGON, geterrorstringon_, GETERRORSTRINGON_)(int *id)
{
	return GetErrorStringOnF(id);
}
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(getlogfilename, GETLOGFILENAME, getlogfilename_, GETLOGFILENAME_)(int *id, char *filename, size_t len)
{
	GetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getlogfileon, GETLOGFILEON, getlogfileon_, GETLOGFILEON_)(int *id)
{
	return GetLogFileOnF(id);
}
// GetLogString
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(getlogstringline, GETLOGSTRINGLINE, getlogstringline_, GETLOGSTRINGLINE_)(int *id, int *n, char* line, size_t line_length)
{
	GetLogStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getlogstringlinecount, GETLOGSTRINGLINECOUNT, getlogstringlinecount_, GETLOGSTRINGLINECOUNT_)(int *id)
{
	return GetLogStringLineCountF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getlogstringon, GETLOGSTRINGON, getlogstringon_, GETLOGSTRINGON_)(int *id)
{
	return GetLogStringOnF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getnthselectedoutputusernumber, GETNTHSELECTEDOUTPUTUSERNUMBER, getnthselectedoutputusernumber_, GETNTHSELECTEDOUTPUTUSERNUMBER_)(int *id, int *n)
{
	return GetNthSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(getoutputfilename, GETOUTPUTFILENAME, getoutputfilename_, GETOUTPUTFILENAME_)(int *id, char *filename, size_t len)
{
	GetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getoutputfileon, GETOUTPUTFILEON, getoutputfileon_, GETOUTPUTFILEON_)(int *id)
{
	return GetOutputFileOnF(id);
}
// GetOutputString
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(getoutputstringline, GETOUTPUTSTRINGLINE, getoutputstringline_, GETOUTPUTSTRINGLINE_)(int *id, int *n, char* line, size_t line_length)
{
	GetOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getoutputstringlinecount, GETOUTPUTSTRINGLINECOUNT, getoutputstringlinecount_, GETOUTPUTSTRINGLINECOUNT_)(int *id)
{
	return GetOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getoutputstringon, GETOUTPUTSTRINGON, getoutputstringon_, GETOUTPUTSTRINGON_)(int *id)
{
	return GetOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getselectedoutputcolumncount, GETSELECTEDOUTPUTCOLUMNCOUNT, getselectedoutputcolumncount_, GETSELECTEDOUTPUTCOLUMNCOUNT_)(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getselectedoutputcount, GETSELECTEDOUTPUTCOUNT, getselectedoutputcount_, GETSELECTEDOUTPUTCOUNT_)(int *id)
{
	return GetSelectedOutputCountF(id);
}
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(getselectedoutputfilename, GETSELECTEDOUTPUTFILENAME, getselectedoutputfilename_, GETSELECTEDOUTPUTFILENAME_)(int *id, char *filename, size_t len)
{
	GetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getselectedoutputfileon, GETSELECTEDOUTPUTFILEON, getselectedoutputfileon_, GETSELECTEDOUTPUTFILEON_)(int *id)
{
	return GetSelectedOutputFileOnF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getselectedoutputrowcount, GETSELECTEDOUTPUTROWCOUNT, getselectedoutputrowcount_, GETSELECTEDOUTPUTROWCOUNT_)(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
// GetSelectedOutputString
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(getselectedoutputstringline, GETSELECTEDOUTPUTSTRINGLINE, getselectedoutputstringline_, GETSELECTEDOUTPUTSTRINGLINE_)(int *id, int *n, char* line, size_t line_length)
{
	GetSelectedOutputStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getselectedoutputstringlinecount, GETSELECTEDOUTPUTSTRINGLINECOUNT, getselectedoutputstringlinecount_, GETSELECTEDOUTPUTSTRINGLINECOUNT_)(int *id)
{
	return GetSelectedOutputStringLineCountF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getselectedoutputstringon, GETSELECTEDOUTPUTSTRINGON, getselectedoutputstringon_, GETSELECTEDOUTPUTSTRINGON_)(int *id)
{
	return GetSelectedOutputStringOnF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getselectedoutputvalue, GETSELECTEDOUTPUTVALUE, getselectedoutputvalue_, GETSELECTEDOUTPUTVALUE_)(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, size_t svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(getversionstring, GETVERSIONSTRING, getversionstring_, GETVERSIONSTRING_)(char* version, size_t version_length)
{
	GetVersionStringF(version, version_length);
}
// GetWarningString
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(getwarningstringline, GETWARNINGSTRINGLINE, getwarningstringline_, GETWARNINGSTRINGLINE_)(int *id, int *n, char* line, size_t line_length)
{
	GetWarningStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(getwarningstringlinecount, GETWARNINGSTRINGLINECOUNT, getwarningstringlinecount_, GETWARNINGSTRINGLINECOUNT_)(int *id)
{
	return GetWarningStringLineCountF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(loaddatabase, LOADDATABASE, loaddatabase_, LOADDATABASE_)(int *id, char *filename, size_t len)
{
	return LoadDatabaseF(id, filename, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(loaddatabasestring, LOADDATABASESTRING, loaddatabasestring_, LOADDATABASESTRING_)(int *id, char *input, size_t len)
{
	return LoadDatabaseStringF(id, input, len);
}
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(outputaccumulatedlines, OUTPUTACCUMULATEDLINES, outputaccumulatedlines_, OUTPUTACCUMULATEDLINES_)(int *id)
{
	OutputAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(outputerrorstring, OUTPUTERRORSTRING, outputerrorstring_, OUTPUTERRORSTRING_)(int *id)
{
	OutputErrorStringF(id);
}
IPQ_DLL_EXPORT void IPQ_DECL IPQ_CASE_UND(outputwarningstring, OUTPUTWARNINGSTRING, outputwarningstring_, OUTPUTWARNINGSTRING_)(int *id)
{
	OutputWarningStringF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(runaccumulated, RUNACCUMULATED, runaccumulated_, RUNACCUMULATED_)(int *id)
{
	return RunAccumulatedF(id);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(runfile, RUNFILE, runfile_, RUNFILE_)(int *id, char *filename, size_t len)
{
	return RunFileF(id, filename, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(runstring, RUNSTRING, runstring_, RUNSTRING_)(int *id, char *input, size_t len)
{
	return RunStringF(id, input, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(setbasicfortrancallback, SETBASICFORTRANCALLBACK, setbasicfortrancallback_, SETBASICFORTRANCALLBACK_)(int *id, double (*fcn)(double *x1, double *x2, const char *str, size_t l))
{
	return SetBasicFortranCallbackF(id, fcn);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(setcurrentselectedoutputusernumber, SETCURRENTSELECTEDOUTPUTUSERNUMBER, setcurrentselectedoutputusernumber_, SETCURRENTSELECTEDOUTPUTUSERNUMBER_)(int *id, int *n)
{
	return SetCurrentSelectedOutputUserNumberF(id, n);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(setdumpfilename, SETDUMPFILENAME, setdumpfilename_, SETDUMPFILENAME_)(int *id, char *filename, size_t len)
{
	return SetDumpFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(setdumpfileon, SETDUMPFILEON, setdumpfileon_, SETDUMPFILEON_)(int *id, int *dump_on)
{
	return SetDumpFileOnF(id, dump_on);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(setdumpstringon, SETDUMPSTRINGON, setdumpstringon_, SETDUMPSTRINGON_)(int *id, int *dump_string_on)
{
	return SetDumpStringOnF(id, dump_string_on);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(seterrorfilename, SETERRORFILENAME, seterrorfilename_, SETERRORFILENAME_)(int *id, char *filename, size_t len)
{
	return SetErrorFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(seterrorfileon, SETERRORFILEON, seterrorfileon_, SETERRORFILEON_)(int *id, int *error_on)
{
	return SetErrorFileOnF(id, error_on);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(seterroron, SETERRORON, seterroron_, SETERRORON_)(int *id, int *error_on)
{
	return SetErrorOnF(id, error_on);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(seterrorstringon, SETERRORSTRINGON, seterrorstringon_, SETERRORSTRINGON_)(int *id, int *error_string_on)
{
	return SetErrorStringOnF(id, error_string_on);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(setlogfilename, SETLOGFILENAME, setlogfilename_, SETLOGFILENAME_)(int *id, char *filename, size_t len)
{
	return SetLogFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(setlogfileon, SETLOGFILEON, setlogfileon_, SETLOGFILEON_)(int *id, int *log_on)
{
	return SetLogFileOnF(id, log_on);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(setlogstringon, SETLOGSTRINGON, setlogstringon_, SETLOGSTRINGON_)(int *id, int *log_on)
{
	return SetLogStringOnF(id, log_on);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(setoutputfilename, SETOUTPUTFILENAME, setoutputfilename_, SETOUTPUTFILENAME_)(int *id, char *filename, size_t len)
{
	return SetOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(setoutputfileon, SETOUTPUTFILEON, setoutputfileon_, SETOUTPUTFILEON_)(int *id, int *output_on)
{
	return SetOutputFileOnF(id, output_on);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(setoutputstringon, SETOUTPUTSTRINGON, setoutputstringon_, SETOUTPUTSTRINGON_)(int *id, int *output_on)
{
	return SetOutputStringOnF(id, output_on);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(setselectedoutputfilename, SETSELECTEDOUTPUTFILENAME, setselectedoutputfilename_, SETSELECTEDOUTPUTFILENAME_)(int *id, char *filename, size_t len)
{
	return SetSelectedOutputFileNameF(id, filename, len);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(setselectedoutputfileon, SETSELECTEDOUTPUTFILEON, setselectedoutputfileon_, SETSELECTEDOUTPUTFILEON_)(int *id, int *selout_file_on)
{
	return SetSelectedOutputFileOnF(id, selout_file_on);
}
IPQ_DLL_EXPORT int  IPQ_DECL IPQ_CASE_UND(setselectedoutputstringon, SETSELECTEDOUTPUTSTRINGON, setselectedoutputstringon_, SETSELECTEDOUTPUTSTRINGON_)(int *id, int *selout_string_on)
{
	return SetSelectedOutputStringOnF(id, selout_string_on);
}
