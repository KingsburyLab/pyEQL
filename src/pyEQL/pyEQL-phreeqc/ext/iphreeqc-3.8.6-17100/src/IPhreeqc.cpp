#include <memory>                       // auto_ptr
#include <map>
#include <string.h>
#include "IPhreeqc.hpp"                 // IPhreeqc
#include "Phreeqc.h"                    // Phreeqc
#include "thread.h"
#include "Version.h"

#include "Debug.h"                      // ASSERT
#include "ErrorReporter.hxx"            // CErrorReporter
#include "CSelectedOutput.hxx"          // CSelectedOutput
#include "SelectedOutput.h"             // SelectedOutput
#include "dumper.h"                     // dumper

// statics
std::map<size_t, IPhreeqc*> IPhreeqc::Instances;
size_t IPhreeqc::InstancesIndex = 0;

std::string IPhreeqc::Version(VERSION_STRING);

static const char empty[] = "";


IPhreeqc::IPhreeqc(void)
: DatabaseLoaded(false)
, ClearAccumulated(false)
, UpdateComponents(true)
, OutputFileOn(false)
, LogFileOn(false)
, ErrorFileOn(false)
, DumpOn(false)
, DumpStringOn(false)
, OutputStringOn(false)
, LogStringOn(false)
, ErrorStringOn(true)
, ErrorReporter(0)
, WarningStringOn(true)
, WarningReporter(0)
, CurrentSelectedOutputUserNumber(1)
, PhreeqcPtr(0)
, input_file(0)
, database_file(0)
{
	this->ErrorReporter   = new CErrorReporter<std::ostringstream>;
	this->WarningReporter = new CErrorReporter<std::ostringstream>;
	this->PhreeqcPtr      = new Phreeqc(this);

	ASSERT(this->PhreeqcPtr->phast == 0);
	this->UnLoadDatabase();

	mutex_lock(&map_lock);
	this->Index = IPhreeqc::InstancesIndex++;
	std::map<size_t, IPhreeqc*>::value_type instance(this->Index, this);
	/*std::pair<std::map<size_t, IPhreeqc*>::iterator, bool> pr = */IPhreeqc::Instances.insert(instance);
	mutex_unlock(&map_lock);

	this->SelectedOutputStringOn[1] = false;

	this->SelectedOutputFileOnMap[1] = false;
	this->SelectedOutputFileNameMap[1] = this->sel_file_name(1);
	
	this->OutputFileName = create_file_name("phreeqc", "out");
	this->ErrorFileName  = create_file_name("phreeqc", "err");
	this->LogFileName    = create_file_name("phreeqc", "log");

	this->DumpFileName   = create_file_name("dump", "out");
	this->PhreeqcPtr->dump_info.Set_file_name(this->DumpFileName);
}

IPhreeqc::~IPhreeqc(void)
{
#if !defined(NDEBUG)
	this->OutputFileOn = false;
#endif
	delete this->PhreeqcPtr;
	delete this->WarningReporter;
	delete this->ErrorReporter;

	std::map< int, CSelectedOutput* >::iterator sit = this->SelectedOutputMap.begin();
	for (; sit != this->SelectedOutputMap.end(); ++sit)
	{
		delete (*sit).second;
	}
	this->SelectedOutputMap.clear();

	mutex_lock(&map_lock);
	std::map<size_t, IPhreeqc*>::iterator it = IPhreeqc::Instances.find(this->Index);
	if (it != IPhreeqc::Instances.end())
	{
		IPhreeqc::Instances.erase(it);
	}
	mutex_unlock(&map_lock);
}

VRESULT IPhreeqc::AccumulateLine(const char *line)
{
	try
	{
		if (this->ClearAccumulated)
		{
			this->ClearAccumulatedLines();
			this->ClearAccumulated = false;
		}

		this->ErrorReporter->Clear();
		this->WarningReporter->Clear();
		this->StringInput.append(line);
		this->StringInput.append("\n");
		return VR_OK;
	}
	catch (...)
	{
		this->AddError("AccumulateLine: An unhandled exception occurred.\n");
		throw;
	}
	return VR_OUTOFMEMORY;
}

size_t IPhreeqc::AddError(const char* str)
{
	return this->ErrorReporter->AddError(str);
}

size_t IPhreeqc::AddWarning(const char* str)
{
	return this->WarningReporter->AddError(str);
}

void IPhreeqc::ClearAccumulatedLines(void)
{
	this->StringInput.erase();
}

const std::string& IPhreeqc::GetAccumulatedLines(void)
{
	return this->StringInput;
}

const char* IPhreeqc::GetComponent(int n)
{
	static const char empty[] = "";
	this->ListComponents();
	if (n < 0 || n >= (int)this->Components.size())
	{
		return empty;
	}
	std::list< std::string >::iterator it = this->Components.begin();
	for(int i = 0; i < n; ++i)
	{
		++it;
	}
	return (*it).c_str();
}

size_t IPhreeqc::GetComponentCount(void)
{
	this->ListComponents();
	return this->Components.size();
}

int IPhreeqc::GetCurrentSelectedOutputUserNumber(void)const
{
	return this->CurrentSelectedOutputUserNumber;
}

const char* IPhreeqc::GetDumpFileName(void)const
{
	return this->DumpFileName.c_str();
}

bool IPhreeqc::GetDumpFileOn(void)const
{
	return this->DumpOn;
}

const char* IPhreeqc::GetDumpString(void)const
{
	static const char err_msg[] = "GetDumpString: DumpStringOn not set.\n";
	if (!this->DumpStringOn)
	{
		return err_msg;
	}
	return this->DumpString.c_str();
}

const char* IPhreeqc::GetDumpStringLine(int n)
{
	static const char empty[] = "";
	if (n < 0 || n >= this->GetDumpStringLineCount())
	{
		return empty;
	}
	return this->DumpLines[n].c_str();
}

int IPhreeqc::GetDumpStringLineCount(void)const
{
	return (int)this->DumpLines.size();
}

bool IPhreeqc::GetDumpStringOn(void)const
{
	return this->DumpStringOn;
}

const char* IPhreeqc::GetErrorFileName(void)const
{
	return this->ErrorFileName.c_str();
}

bool IPhreeqc::GetErrorFileOn(void)const
{
	return this->ErrorFileOn;
}

bool IPhreeqc::GetErrorOn(void)const
{
	return this->Get_error_on();
}

const char* IPhreeqc::GetErrorString(void)
{
	static const char err_msg[] = "GetErrorString: ErrorStringOn not set.\n";
	static const char err_msg2[] = "GetErrorString: ErrorOn not set.\n";
	if (!this->error_on)
	{
		return err_msg2;
	}
	if (!this->ErrorStringOn)
	{
		return err_msg;
	}
	this->ErrorString = ((CErrorReporter<std::ostringstream>*)this->ErrorReporter)->GetOS()->str();
	return this->ErrorString.c_str();
}

const char* IPhreeqc::GetErrorStringLine(int n)
{
	static const char empty[] = "";
	if (n < 0 || n >= this->GetErrorStringLineCount())
	{
		return empty;
	}
	return this->ErrorLines[n].c_str();
}

int IPhreeqc::GetErrorStringLineCount(void)const
{
	return (int)this->ErrorLines.size();
}

bool IPhreeqc::GetErrorStringOn(void)const
{
	return this->ErrorStringOn;
}

int IPhreeqc::GetId(void)const
{
	return (int)this->Index;
}

const char* IPhreeqc::GetLogFileName(void)const
{
	return this->LogFileName.c_str();
}

bool IPhreeqc::GetLogFileOn(void)const
{
	return this->LogFileOn;
}

const char* IPhreeqc::GetLogString(void)const
{
	static const char err_msg[] = "GetLogString: LogStringOn not set.\n";
	if (!this->LogStringOn)
	{
		return err_msg;
	}
	return this->LogString.c_str();
}

const char* IPhreeqc::GetLogStringLine(int n)const
{
	static const char empty[] = "";
	if (n < 0 || n >= this->GetLogStringLineCount())
	{
		return empty;
	}
	return this->LogLines[n].c_str();
}

int IPhreeqc::GetLogStringLineCount(void)const
{
	return (int)this->LogLines.size();
}

bool IPhreeqc::GetLogStringOn(void)const
{
	return this->LogStringOn;
}

int IPhreeqc::GetNthSelectedOutputUserNumber(int n)const
{
	int nth = VR_INVALIDARG;
	std::map< int, SelectedOutput >::const_iterator ci = this->PhreeqcPtr->SelectedOutput_map.begin();
	for (int i = 0; ci != this->PhreeqcPtr->SelectedOutput_map.end(); ++ci, ++i)
	{
		if (i == n)
		{
			nth = (*ci).first;
			break;
		}
	}
	return nth;
}

const char* IPhreeqc::GetOutputFileName(void)const
{
	return this->OutputFileName.c_str();
}

bool IPhreeqc::GetOutputFileOn(void)const
{
	return this->OutputFileOn;
}

const char* IPhreeqc::GetOutputString(void)const
{
	static const char err_msg[] = "GetOutputString: OutputStringOn not set.\n";
	if (!this->OutputStringOn)
	{
		return err_msg;
	}
	return this->OutputString.c_str();
}

const char* IPhreeqc::GetOutputStringLine(int n)const
{
	static const char empty[] = "";
	if (n < 0 || n >= this->GetOutputStringLineCount())
	{
		return empty;
	}
	return this->OutputLines[n].c_str();
}

int IPhreeqc::GetOutputStringLineCount(void)const
{
	return (int)this->OutputLines.size();
}

bool IPhreeqc::GetOutputStringOn(void)const
{
	return this->OutputStringOn;
}

int IPhreeqc::GetSelectedOutputColumnCount(void)const
{
	std::map< int, CSelectedOutput* >::const_iterator ci = this->SelectedOutputMap.find(this->CurrentSelectedOutputUserNumber);
	if (ci != this->SelectedOutputMap.end())
	{
		return (int)(*ci).second->GetColCount();
	}
	return 0;
}

int IPhreeqc::GetSelectedOutputCount(void)const
{
	ASSERT(this->PhreeqcPtr->SelectedOutput_map.size() == this->SelectedOutputMap.size());
	return (int) this->PhreeqcPtr->SelectedOutput_map.size();
}

const char* IPhreeqc::GetSelectedOutputFileName(void)const
{
	static const char empty[] = "";
	std::map< int, std::string >::const_iterator ci = this->SelectedOutputFileNameMap.find(this->CurrentSelectedOutputUserNumber);
	if (ci != this->SelectedOutputFileNameMap.end())
	{
		return (*ci).second.c_str();
	}
	return empty;
}

bool IPhreeqc::GetSelectedOutputFileOn(void)const
{
	return this->get_sel_out_file_on(this->CurrentSelectedOutputUserNumber);
}

int IPhreeqc::GetSelectedOutputRowCount(void)const
{
	std::map< int, CSelectedOutput* >::const_iterator ci = this->SelectedOutputMap.find(this->CurrentSelectedOutputUserNumber);
	if (ci != this->SelectedOutputMap.end())
	{
		return (int)(*ci).second->GetRowCount();
	}
	return 0;
}

const char* IPhreeqc::GetSelectedOutputString(void)const
{
	static const char err_msg[] = "GetSelectedOutputString: SelectedOutputStringOn not set.\n";
	std::map< int, bool >::const_iterator ci = this->SelectedOutputStringOn.find(this->CurrentSelectedOutputUserNumber);
	if (ci == this->SelectedOutputStringOn.end())
	{
		return err_msg;
	}
	std::map< int, std::string >::const_iterator cit = this->SelectedOutputStringMap.find(this->CurrentSelectedOutputUserNumber);
	if (cit != this->SelectedOutputStringMap.end())
	{
		return (*cit).second.c_str();
	}
	return empty;
}

const char* IPhreeqc::GetSelectedOutputStringLine(int n)
{
	static const char empty[] = "";
	if (n < 0 || n >= this->GetSelectedOutputStringLineCount())
	{
		return empty;
	}
	// don't need to check CurrentSelectedOutputUserNumber since
	// GetSelectedOutputStringLineCount will catch
	return this->SelectedOutputLinesMap[this->CurrentSelectedOutputUserNumber][n].c_str();
}

int IPhreeqc::GetSelectedOutputStringLineCount(void)const
{
	std::map< int, std::vector < std::string > >::const_iterator cit = this->SelectedOutputLinesMap.find(this->CurrentSelectedOutputUserNumber);
	if (cit != this->SelectedOutputLinesMap.end())
	{
		return (int) (*cit).second.size();
	}
	return 0;
}

bool IPhreeqc::GetSelectedOutputStringOn(void)const
{
	return this->get_sel_out_string_on(this->CurrentSelectedOutputUserNumber);
}

VRESULT IPhreeqc::GetSelectedOutputValue(int row, int col, VAR* pVAR)
{
	VRESULT v;
	this->ErrorReporter->Clear();
	if (!pVAR)
	{
		this->AddError("GetSelectedOutputValue: VR_INVALIDARG pVAR is NULL.\n");
		this->update_errors();
		return VR_INVALIDARG;
	}

	std::map< int, CSelectedOutput* >::const_iterator ci = this->SelectedOutputMap.find(this->CurrentSelectedOutputUserNumber);
	if (ci != this->SelectedOutputMap.end())
	{
		v = (*ci).second->Get(row, col, pVAR);
		switch (v)
		{
		case VR_OK:
			break;
		case VR_INVALIDROW:
			this->AddError("GetSelectedOutputValue: VR_INVALIDROW Row index out of range.\n");
			this->update_errors();
			break;
		case VR_INVALIDCOL:
			this->AddError("GetSelectedOutputValue: VR_INVALIDCOL Column index out of range.\n");
			this->update_errors();
			break;
		case VR_OUTOFMEMORY:
			this->AddError("GetSelectedOutputValue: VR_OUTOFMEMORY Out of memory.\n");
			this->update_errors();
			break;
		case VR_BADVARTYPE:
			this->AddError("GetSelectedOutputValue: VR_BADVARTYPE pVar must be initialized(VarInit) and/or cleared(VarClear).\n");
			this->update_errors();
			break;
		case VR_INVALIDARG:
			// not possible
			break;
		default:
			assert(0);
		}
	}
	else
	{
		char buffer[120];
		v = VR_INVALIDARG;
		::snprintf(buffer, sizeof(buffer), "GetSelectedOutputValue: VR_INVALIDARG Invalid selected-output user number %d.\n", this->CurrentSelectedOutputUserNumber);
		this->AddError(buffer);
		this->update_errors();
	}
	return v;
}

VRESULT IPhreeqc::GetSelectedOutputValue2(int row, int col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	VRESULT result;
	VAR v;
	VarInit(&v);
	char buffer[100];

	result = this->GetSelectedOutputValue(row, col, &v);

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
		::strncpy(svalue, buffer, svalue_length);
		break;
	case TT_DOUBLE:
		*vtype = v.type;
		*dvalue = v.dVal;
		::snprintf(buffer, sizeof(buffer), "%23.15e", v.dVal);
		::strncpy(svalue, buffer, svalue_length);
		break;
	case TT_STRING:
		*vtype = v.type;
		::strncpy(svalue, v.sVal, svalue_length);
		break;
	default:
		assert(0);
	}
	::VarClear(&v);
	return result;
}

const char* IPhreeqc::GetVersionString(void)
{
	return IPhreeqc::Version.c_str();
}

const char* IPhreeqc::GetWarningString(void)
{
	this->WarningString = ((CErrorReporter<std::ostringstream>*)this->WarningReporter)->GetOS()->str();
	return this->WarningString.c_str();
}

const char* IPhreeqc::GetWarningStringLine(int n)
{
	static const char empty[] = "";
	if (n < 0 || n >= this->GetWarningStringLineCount())
	{
		return empty;
	}
	return this->WarningLines[n].c_str();
}

int IPhreeqc::GetWarningStringLineCount(void)const
{
	return (int)this->WarningLines.size();
}

std::list< std::string > IPhreeqc::ListComponents(void)
{
	if (this->UpdateComponents)
	{
		this->Components.clear();
		this->PhreeqcPtr->list_components(this->Components);
		this->PhreeqcPtr->list_EquilibriumPhases(this->EquilibriumPhasesList);
		this->PhreeqcPtr->list_GasComponents(this->GasComponentsList);
		this->PhreeqcPtr->list_KineticReactions(this->KineticReactionsList);
		this->PhreeqcPtr->list_SolidSolutions(this->SolidSolutionComponentsList,this->SolidSolutionNamesList);
		this->PhreeqcPtr->list_Surfaces(this->SurfaceTypeList, this->SurfaceNamesList);
		this->PhreeqcPtr->list_Exchangers(this->ExchangeNamesList);
		this->UpdateComponents = false;
	}
	return this->Components;
}

int IPhreeqc::load_db(const char* filename)
{
	try
	{
		// cleanup
		//
		this->UnLoadDatabase();

		// open file
		//
		std::ifstream ifs;
		ifs.open(filename);

		if (!ifs.is_open())
		{
			std::ostringstream oss;
			oss << "LoadDatabase: Unable to open:" << "\"" << filename << "\".";
			this->PhreeqcPtr->error_msg(oss.str().c_str(), STOP); // throws
		}

		// read input
		//
		ASSERT(this->PhreeqcPtr->phrq_io->get_istream() == NULL);
		this->PhreeqcPtr->phrq_io->push_istream(&ifs, false);
		this->PhreeqcPtr->read_database();
	}
	catch (const IPhreeqcStop&)
	{
		this->close_input_files();
	}
	catch (...)
	{
		const char *errmsg = "LoadDatabase: An unhandled exception occurred.\n";
		try
		{
			this->PhreeqcPtr->error_msg(errmsg, STOP); // throws IPhreeqcStop
		}
		catch (const IPhreeqcStop&)
		{
			// do nothing
		}
		throw;
	}
	this->PhreeqcPtr->phrq_io->clear_istream();
	this->DatabaseLoaded = (this->PhreeqcPtr->get_input_errors() == 0);
	return this->PhreeqcPtr->get_input_errors();
}

int IPhreeqc::LoadDatabase(const char* filename)
{
	// save I/O state
	bool bSaveErrorFileOn  = this->ErrorFileOn;
	bool bSaveOutputOn     = this->OutputFileOn;
	bool bSaveLogFileOn    = this->LogFileOn;
	this->ErrorFileOn      = false;
	this->OutputFileOn     = false;
	this->LogFileOn        = false;

	int n = this->load_db(filename);
	if (n == 0)
	{
		n = this->test_db();
	}

	// restore I/O state
	this->ErrorFileOn  = bSaveErrorFileOn;
	this->OutputFileOn = bSaveOutputOn;
	this->LogFileOn    = bSaveLogFileOn;

	return n;
}

int IPhreeqc::load_db_str(const char* input)
{
	try
	{
		// cleanup
		//
		this->UnLoadDatabase();

		std::string s(input);
		std::istringstream iss(s);

		// read input
		//
		ASSERT(this->PhreeqcPtr->phrq_io->get_istream() == NULL);
		this->PhreeqcPtr->phrq_io->push_istream(&iss, false);
		this->PhreeqcPtr->read_database();
	}
	catch (const IPhreeqcStop&)
	{
		this->close_input_files();
	}
	catch(...)
	{
		const char *errmsg = "LoadDatabaseString: An unhandled exception occurred.\n";
		try
		{
			this->PhreeqcPtr->error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (const IPhreeqcStop&)
		{
			// do nothing
		}
		throw;
	}

	this->PhreeqcPtr->phrq_io->clear_istream();
	this->DatabaseLoaded = (this->PhreeqcPtr->get_input_errors() == 0);
	return this->PhreeqcPtr->get_input_errors();
}

int IPhreeqc::LoadDatabaseString(const char* input)
{
	// save I/O state
	bool bSaveErrorFileOn  = this->ErrorFileOn;
	bool bSaveOutputOn     = this->OutputFileOn;
	bool bSaveLogFileOn    = this->LogFileOn;
	this->ErrorFileOn      = false;
	this->OutputFileOn     = false;
	this->LogFileOn        = false;

	int n = this->load_db_str(input);
	if (n == 0)
	{
		n = this->test_db();
	}

	// restore I/O state
	this->ErrorFileOn  = bSaveErrorFileOn;
	this->OutputFileOn = bSaveOutputOn;
	this->LogFileOn    = bSaveLogFileOn;

	return n;
}

void IPhreeqc::OutputAccumulatedLines(void)
{
#if !defined(R_SO)
	std::cout << this->StringInput.c_str() << std::endl;
#endif
}

void IPhreeqc::OutputErrorString(void)
{
#if !defined(R_SO)
	std::cout << this->GetErrorString() << std::endl;
#endif
}

void IPhreeqc::OutputWarningString(void)
{
#if !defined(R_SO)
	std::cout << this->GetWarningString() << std::endl;
#endif
}

int IPhreeqc::RunAccumulated(void)
{
	static const char *sz_routine = "RunAccumulated";
	try
	{
		// DON'T clear accumulated (just set
		// flag to on next AccumulateLine)

		// these may throw
		this->open_output_files(sz_routine);
		this->check_database(sz_routine);

		this->PhreeqcPtr->input_error = 0;
		this->io_error_count = 0;

		// create input stream
		std::istringstream iss(this->GetAccumulatedLines());

		// this may throw
		this->do_run(sz_routine, &iss, NULL, NULL, NULL);
	}
	catch (const IPhreeqcStop&)
	{
		// do nothing
	}
	catch(std::exception &e)
	{
		std::string errmsg("RunAccumulated: ");
		errmsg += e.what();
		try
		{
			this->PhreeqcPtr->error_msg(errmsg.c_str(), STOP); // throws PhreeqcStop
		}
		catch (const IPhreeqcStop&)
		{
			// do nothing
		}
		throw;
	}
	catch(...)
	{
		const char *errmsg = "RunAccumulated: An unhandled exception occurred.\n";
		try
		{
			this->PhreeqcPtr->error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (const IPhreeqcStop&)
		{
			// do nothing
		}
		throw;
	}

	this->ClearAccumulated = true;
	this->close_output_files();
	this->update_errors();
	this->PhreeqcPtr->phrq_io->clear_istream();

	return this->PhreeqcPtr->get_input_errors();
}

int IPhreeqc::RunFile(const char* filename)
{
	static const char *sz_routine = "RunFile";
	try
	{
		// clear accumulated
		//
		this->ClearAccumulatedLines();
		this->ClearAccumulated = false;

		// these may throw
		this->open_output_files(sz_routine);
		this->check_database(sz_routine);

		this->PhreeqcPtr->input_error = 0;
		this->io_error_count = 0;

		// open file
		//
		std::ifstream ifs;
		ifs.open(filename);
		if (!ifs.is_open())
		{
			std::ostringstream oss;
			oss << "RunFile: Unable to open:" << "\"" << filename << "\".";
			this->PhreeqcPtr->error_msg(oss.str().c_str(), STOP); // throws
		}

		// this may throw
		this->do_run(sz_routine, &ifs, NULL, NULL, NULL);
	}
	catch (const IPhreeqcStop&)
	{
		this->close_input_files();
	}
	catch(std::exception &e)
	{
		std::string errmsg("RunFile: ");
		errmsg += e.what();
		try
		{
			this->PhreeqcPtr->error_msg(errmsg.c_str(), STOP); // throws PhreeqcStop
		}
		catch (const IPhreeqcStop&)
		{
			// do nothing
		}
		throw;
	}
	catch(...)
	{
		const char *errmsg = "RunFile: An unhandled exception occurred.\n";
		try
		{
			this->PhreeqcPtr->error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (const IPhreeqcStop&)
		{
			// do nothing
		}
		throw;
	}

	this->close_output_files();
	this->update_errors();
	this->PhreeqcPtr->phrq_io->clear_istream();

	return this->PhreeqcPtr->get_input_errors();
}

int IPhreeqc::RunString(const char* input)
{
	static const char *sz_routine = "RunString";
	try
	{
		// clear accumulated
		//
		this->ClearAccumulatedLines();
		this->ClearAccumulated = false;

		// these may throw
		this->open_output_files(sz_routine);
		this->check_database(sz_routine);

		this->PhreeqcPtr->input_error = 0;
		this->io_error_count = 0;

		// create input stream
		std::string s(input);
		std::istringstream iss(s);

		// this may throw
		this->do_run(sz_routine, &iss, NULL, NULL, NULL);
	}
	catch (const IPhreeqcStop&)
	{
		this->close_input_files();
	}
	catch(std::exception &e)
	{
		std::string errmsg("RunString: ");
		errmsg += e.what();
		try
		{
			this->PhreeqcPtr->error_msg(errmsg.c_str(), STOP); // throws PhreeqcStop
		}
		catch (const IPhreeqcStop&)
		{
			// do nothing
		}
		throw;
	}
	catch(...)
	{
		const char *errmsg = "RunString: An unhandled exception occurred.\n";
		try
		{
			this->PhreeqcPtr->error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (const IPhreeqcStop&)
		{
			// do nothing
		}
		throw;
	}

	this->close_output_files();
	this->update_errors();
	this->PhreeqcPtr->phrq_io->clear_istream();

	return this->PhreeqcPtr->get_input_errors();
}

void IPhreeqc::SetBasicCallback(double (*fcn)(double x1, double x2, const char *str, void *cookie), void *cookie1)
{
	this->PhreeqcPtr->register_basic_callback(fcn, cookie1);
}
#ifdef IPHREEQC_NO_FORTRAN_MODULE
void IPhreeqc::SetBasicFortranCallback(double (*fcn)(double *x1, double *x2, const char *str, size_t l))
{
	this->PhreeqcPtr->register_fortran_basic_callback(fcn);
}
#else
void IPhreeqc::SetBasicFortranCallback(double (*fcn)(double *x1, double *x2, const char *str, int l))
{
	this->PhreeqcPtr->register_fortran_basic_callback(fcn);
}
#endif
VRESULT IPhreeqc::SetCurrentSelectedOutputUserNumber(int n)
{
	if (0 <= n)
	{
		this->CurrentSelectedOutputUserNumber = n;
		return VR_OK;
	}
	return VR_INVALIDARG;
}

void IPhreeqc::SetDumpFileName(const char *filename)
{
	if (filename && ::strlen(filename))
	{
		this->DumpFileName = filename;
		this->PhreeqcPtr->dump_info.Set_file_name(this->DumpFileName);
	}
}

void IPhreeqc::SetDumpFileOn(bool bValue)
{
	this->DumpOn = bValue;
}

void IPhreeqc::SetDumpStringOn(bool bValue)
{
	this->DumpStringOn = bValue;
}

void IPhreeqc::SetErrorFileName(const char *filename)
{
	if (filename && ::strlen(filename))
	{
		this->ErrorFileName = filename;
	}
}

void IPhreeqc::SetErrorFileOn(bool bValue)
{
	this->ErrorFileOn = bValue;
}

void IPhreeqc::SetErrorOn(bool bValue)
{
	this->Set_error_on(bValue);
}

void IPhreeqc::SetErrorStringOn(bool bValue)
{
	this->ErrorStringOn = bValue;
}

void IPhreeqc::SetLogFileName(const char *filename)
{
	if (filename && ::strlen(filename))
	{
		this->LogFileName = filename;
	}
}

void IPhreeqc::SetLogFileOn(bool bValue)
{
	this->LogFileOn = bValue;
}

void IPhreeqc::SetLogStringOn(bool bValue)
{
	this->LogStringOn = bValue;
}

void IPhreeqc::SetOutputFileName(const char *filename)
{
	if (filename && ::strlen(filename))
	{
		this->OutputFileName = filename;
	}
}

void IPhreeqc::SetOutputStringOn(bool bValue)
{
	this->OutputStringOn = bValue;
}

void IPhreeqc::SetOutputFileOn(bool bValue)
{
	this->OutputFileOn = bValue;
}

void IPhreeqc::SetSelectedOutputFileName(const char *filename)
{
	if (filename && ::strlen(filename))
	{
		// Can't use this->PhreeqcPtr->SelectedOutput_map since it's necessary
		// to override the default filename "selected_output_%d.%d.sel" 
		this->SelectedOutputFileNameMap[this->CurrentSelectedOutputUserNumber] = std::string(filename);
	}
}

void IPhreeqc::SetSelectedOutputFileOn(bool bValue)
{
	if (0 <= this->CurrentSelectedOutputUserNumber)
	{
		this->SelectedOutputFileOnMap[this->CurrentSelectedOutputUserNumber] = bValue;
	}
}

void IPhreeqc::SetSelectedOutputStringOn(bool bValue)
{
	this->SelectedOutputStringOn[this->CurrentSelectedOutputUserNumber] = bValue;
}

int IPhreeqc::test_db(void)
{
	std::ostringstream oss;
	int sn = this->PhreeqcPtr->next_user_number(Keywords::KEY_SOLUTION);
	//oss << "SOLUTION " << sn <<"; DELETE; -solution " << sn;
	oss << "SOLUTION " << sn << ";";
	if (this->PhreeqcPtr->llnl_temp.size() > 0)
	{
		oss << "-temp " << this->PhreeqcPtr->llnl_temp[0] << ";";
	}
	oss << "DELETE; -solution " << sn;
	this->PhreeqcPtr->set_reading_database(TRUE);
	int n = this->RunString(oss.str().c_str());
	this->PhreeqcPtr->set_reading_database(FALSE);
	return n;
}

void IPhreeqc::UnLoadDatabase(void)
{
	// init IPhreeqc
	//
	this->DatabaseLoaded   = false;
	this->UpdateComponents = true;
	this->Components.clear();

	// clear accumulated
	//
	this->ClearAccumulatedLines();
	this->ClearAccumulated = false;

	// clear error state
	//
	ASSERT(this->ErrorReporter);
	this->ErrorReporter->Clear();
	this->ErrorString.clear();

	// clear warning state
	//
	ASSERT(this->WarningReporter);
	this->WarningReporter->Clear();
	this->WarningString.clear();

	// clear selectedoutput
	//
	this->CurrentSelectedOutputUserNumber = 1;
	this->SelectedOutputFileOnMap.clear();
	this->SelectedOutputFileOnMap[1] = false;
	this->SelectedOutputStringOn.clear();
	this->SelectedOutputStringOn[1] = false;

	std::map< int, CSelectedOutput* >::iterator itt = this->SelectedOutputMap.begin();
	for (; itt != this->SelectedOutputMap.end(); ++itt)
	{
		delete (*itt).second;
	}
	this->SelectedOutputMap.clear();
	this->SelectedOutputStringMap.clear();
	this->SelectedOutputLinesMap.clear();


	// clear dump string
	//
	this->DumpString.clear();
	this->DumpLines.clear();

	// initialize phreeqc
	//
	this->PhreeqcPtr->clean_up();
	this->PhreeqcPtr->init();
	this->PhreeqcPtr->do_initialize();
	this->PhreeqcPtr->input_error = 0;
	this->io_error_count = 0;
}

int IPhreeqc::EndRow(void)
{
	if (this->PhreeqcPtr->current_selected_output)
	{
		std::map< int, CSelectedOutput* >::iterator it =
			this->SelectedOutputMap.find(this->PhreeqcPtr->current_selected_output->Get_n_user());

		if (it != this->SelectedOutputMap.end())
		{
			// ensure all user_punch headings are included
			ASSERT(this->PhreeqcPtr->n_user_punch_index >= 0);
			if (this->PhreeqcPtr->current_user_punch)
			{
				for (size_t i = this->PhreeqcPtr->n_user_punch_index; i < this->PhreeqcPtr->current_user_punch->Get_headings().size(); ++i)
				{
					(*it).second->PushBackEmpty(this->PhreeqcPtr->current_user_punch->Get_headings()[i].c_str());
				}
			}
			return (*it).second->EndRow();
		}
	}
	return 0;
}

void IPhreeqc::check_database(const char* sz_routine)
{
 	this->ErrorReporter->Clear();
	this->WarningReporter->Clear();
	std::map< int, CSelectedOutput* >::iterator it = this->SelectedOutputMap.begin();
	for (; it != this->SelectedOutputMap.end(); ++it)
	{
		delete (*it).second;
	}
	this->SelectedOutputMap.clear();
	this->SelectedOutputStringMap.clear();
	this->SelectedOutputLinesMap.clear();

	// release
	this->LogString.clear();
	this->LogLines.clear();
	this->OutputString.clear();
	this->OutputLines.clear();


	if (!this->DatabaseLoaded)
	{
		std::ostringstream oss;
		oss << sz_routine << ": No database is loaded";
		this->PhreeqcPtr->input_error = 1;
		this->PhreeqcPtr->error_msg(oss.str().c_str(), STOP); // throws
	}
}

void IPhreeqc::do_run(const char* sz_routine, std::istream* pis, PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie)
{
	char token[MAX_LENGTH];

/*
 *   Maybe should be in read_input
 */
	this->PhreeqcPtr->first_read_input = TRUE;

/*
 *   call pre-run callback
 */
	if (pfn_pre)
	{
		pfn_pre(cookie);
	}

/*
 *   set read callback
 */
#if (__GNUC__ && (__cplusplus >= 201103L)) || (_MSC_VER >= 1600)
	std::unique_ptr<std::istringstream> auto_iss=NULL;
#else
	std::auto_ptr<std::istringstream> auto_iss(0);
#endif
	if (!pis)
	{
		auto_iss.reset(new std::istringstream(this->GetAccumulatedLines()));
		this->PhreeqcPtr->phrq_io->push_istream(auto_iss.get(), false);
	}
	else
	{
		ASSERT(this->PhreeqcPtr->phrq_io->get_istream() == NULL);
		this->PhreeqcPtr->phrq_io->push_istream(pis, false);
	}


/*
 *   Read input data for simulation
 */
	for (this->PhreeqcPtr->simulation = 1; ; this->PhreeqcPtr->simulation++)
	{

#ifdef PHREEQ98
   		AddSeries = !connect_simulations;
#endif
		::snprintf(token, sizeof(token), "Reading input data for simulation %d.", this->PhreeqcPtr->simulation);

		// bool save_punch_in = this->PhreeqcPtr->SelectedOutput_map.size() > 0;

		this->PhreeqcPtr->dup_print(token, TRUE);
		if (this->PhreeqcPtr->read_input() == EOF)
			break;
		
		if (this->PhreeqcPtr->simulation == 1)
		{
			// force headings for selected output (on every call to do_run)
			// might want to split tidy_punch to avoid duplicate searches like master_bsearch
			std::map< int, SelectedOutput >::iterator pit = this->PhreeqcPtr->SelectedOutput_map.begin();
			for (; pit != this->PhreeqcPtr->SelectedOutput_map.end(); ++pit)
			{
				(*pit).second.Set_new_def(true);
				this->PhreeqcPtr->keycount[Keywords::KEY_SELECTED_OUTPUT] = 1;
			}
		}

		// bool bWarning = false;
		std::map< int, SelectedOutput >::iterator mit = this->PhreeqcPtr->SelectedOutput_map.begin();
		for (; mit != this->PhreeqcPtr->SelectedOutput_map.end(); ++mit)
		{
			if (this->SelectedOutputMap.find(mit->first) == this->SelectedOutputMap.end())
			{
				// int -> CSelectedOutput*
				std::map< int, CSelectedOutput* >::value_type item((*mit).first, new CSelectedOutput());
				this->SelectedOutputMap.insert(item);

				// int -> std::string
				this->SelectedOutputStringMap.insert(
					std::map< int, std::string >::value_type((*mit).first, std::string()));
			}
			else
			{
				ASSERT(this->SelectedOutputMap.find((*mit).first) != this->SelectedOutputMap.end());
				ASSERT(this->SelectedOutputStringMap.find((*mit).first) != this->SelectedOutputStringMap.end());
			}
		}
		ASSERT(this->PhreeqcPtr->SelectedOutput_map.size() == this->SelectedOutputMap.size());
		ASSERT(this->PhreeqcPtr->SelectedOutput_map.size() == this->SelectedOutputStringMap.size());
		if (!this->PhreeqcPtr->title_x.empty())
		{
			::snprintf(token, sizeof(token), "TITLE");
			this->PhreeqcPtr->dup_print(token, TRUE);
			if (this->PhreeqcPtr->pr.headings == TRUE)
			{
				this->PhreeqcPtr->output_msg(this->PhreeqcPtr->title_x.c_str());
				this->PhreeqcPtr->output_msg("\n\n");
			}
		}

#ifdef SWIG_SHARED_OBJ
		if (this->PhreeqcPtr->SelectedOutput_map.size() > 0)
		{
			//
			// (punch.in == TRUE) when any "RUN" has contained
			// a SELECTED_OUTPUT block since the last LoadDatabase call.
			//
			// Since LoadDatabase initializes punch.in to FALSE
			// (via UnLoadDatabase...do_initialize)
			// and punch.in is set to TRUE in read_selected_output
			//
			// This causes the SELECTED_OUTPUT to contain the same headings
			// until another SELECTED_OUTPUT is defined which sets the variable
			// punch.new_def to TRUE
			//
			// WHAT IF A USER_PUNCH IS DEFINED?? IS punch.new_def SET TO
			// TRUE ???
			//
			//
			std::map< int, SelectedOutput >::iterator ai = this->PhreeqcPtr->SelectedOutput_map.begin();
			for (; ai != this->PhreeqcPtr->SelectedOutput_map.end(); ++ai)
			{
				if (!this->SelectedOutputFileOnMap[(*ai).first])
				{
					ASSERT((*ai).second.Get_punch_ostream() == 0);
				}
			}

			if (this->PhreeqcPtr->pr.punch == FALSE)
			{
				// No selected_output for this simulation
				// this happens when
				//    PRINT;  -selected_output false
				// is given as input
				// Note: this also disables the CSelectedOutput object
				ASSERT(TRUE);
			}
			else
			{
				std::map< int, SelectedOutput >::iterator it = this->PhreeqcPtr->SelectedOutput_map.begin();
				for (; it != this->PhreeqcPtr->SelectedOutput_map.end(); ++it)
				{
					if (this->SelectedOutputFileOnMap[(*it).first] && !(*it).second.Get_punch_ostream())
					{
						//
						// LoadDatabase
						// do_run -- containing SELECTED_OUTPUT ****TODO**** check -file option
						// another do_run without SELECTED_OUTPUT
						//
						ASSERT(!this->SelectedOutputFileNameMap[(*it).first].empty());
						std::string filename = this->SelectedOutputFileNameMap[(*it).first];
						if (!punch_open(filename.c_str(), std::ios_base::out, (*it).first))
						{
							std::ostringstream oss;
							oss << sz_routine << ": Unable to open:" << "\"" << filename << "\".\n";
							this->PhreeqcPtr->warning_msg(oss.str().c_str());
						}
						else
						{
							ASSERT(this->Get_punch_ostream() != NULL);
							ASSERT((*it).second.Get_punch_ostream() == NULL);

							int n_user = (*it).first;
							this->PhreeqcPtr->SelectedOutput_map[n_user].Set_punch_ostream(this->Get_punch_ostream());
							this->Set_punch_ostream(NULL);
							
							// output selected_output headings
							(*it).second.Set_new_def(TRUE);
							this->PhreeqcPtr->tidy_punch();
						}
					}
				}
			}
		}
		else
		{
			ASSERT(TRUE);
		}
		

		std::map< int, SelectedOutput >::iterator it = this->PhreeqcPtr->SelectedOutput_map.begin();
		for (; it != this->PhreeqcPtr->SelectedOutput_map.end(); ++it)
		{
			if (this->SelectedOutputFileOnMap[(*it).first])
			{
				ASSERT((*it).second.Get_punch_ostream());
			}
			else
			{
				ASSERT(!(*it).second.Get_punch_ostream());
			}
		}

		// Consider this addition
		{
			this->PhreeqcPtr->pr.all = (this->OutputFileOn || this->OutputStringOn) ? TRUE : FALSE;
			//this->PhreeqcPtr->pr.punch = (this->SelectedOutputFileOn || this->SelectedOutputStringOn) ? TRUE : FALSE;
		}
		/* the converse is not necessarily true */

		this->PhreeqcPtr->n_user_punch_index = -1;
#endif // SWIG_SHARED_OBJ
		{
			this->PhreeqcPtr->pr.all = (this->OutputFileOn || this->OutputStringOn) ? TRUE : FALSE;
		}

		this->PhreeqcPtr->tidy_model();
#ifdef PHREEQ98
                if (!phreeq98_debug)
				{
#endif

/*
 *   Calculate distribution of species for initial solutions
 */
		if (this->PhreeqcPtr->new_solution)
			this->PhreeqcPtr->initial_solutions(TRUE);
/*
 *   Calculate distribution for exchangers
 */
		if (this->PhreeqcPtr->new_exchange)
			this->PhreeqcPtr->initial_exchangers(TRUE);
/*
 *   Calculate distribution for surfaces
 */
		if (this->PhreeqcPtr->new_surface)
			this->PhreeqcPtr->initial_surfaces(TRUE);
/*
 *   Calculate initial gas composition
 */
		if (this->PhreeqcPtr->new_gas_phase)
			this->PhreeqcPtr->initial_gas_phases(TRUE);
/*
 *   Calculate reactions
 */
		this->PhreeqcPtr->reactions();
/*
 *   Calculate inverse models
 */
		this->PhreeqcPtr->inverse_models();
/*
 *   Calculate advection
 */
		if (this->PhreeqcPtr->use.Get_advect_in())
		{
			this->PhreeqcPtr->dup_print("Beginning of advection calculations.", TRUE);
			this->PhreeqcPtr->advection();
		}
/*
 *   Calculate transport
 */
		if (this->PhreeqcPtr->use.Get_trans_in())
		{
			this->PhreeqcPtr->dup_print("Beginning of transport calculations.", TRUE);
			this->PhreeqcPtr->transport();
		}
/*
 *   run
 */
		this->PhreeqcPtr->run_as_cells();
/*
 *   Calculate mixes
 */
		this->PhreeqcPtr->do_mixes();
/*
 *   Copy
 */
		if (this->PhreeqcPtr->new_copy) this->PhreeqcPtr->copy_entities();
/*
 *   dump
 */
		dumper dump_info_save(this->PhreeqcPtr->dump_info);
		if (this->DumpOn)
		{
			this->PhreeqcPtr->dump_entities();
			this->DumpFileName = this->PhreeqcPtr->dump_info.Get_file_name();
		}
		if (this->DumpStringOn)
		{
			this->PhreeqcPtr->dump_info = dump_info_save;
			if (this->PhreeqcPtr->dump_info.Get_bool_any())
			{
				std::ostringstream oss;
				this->PhreeqcPtr->dump_ostream(oss);
				if (this->PhreeqcPtr->dump_info.Get_append())
				{
					this->DumpString += oss.str();
				}
				else
				{
					this->DumpString = oss.str();
				}

				/* Fill dump lines */
				this->DumpLines.clear();
				std::istringstream iss(this->DumpString);
				std::string line;
				while (std::getline(iss, line))
				{
					this->DumpLines.push_back(line);
				}
			}
		}
/*
 *   delete
 */
		this->PhreeqcPtr->delete_entities();

/*
 *   End of simulation
 */
		this->PhreeqcPtr->dup_print( "End of simulation.", TRUE);
#ifdef PHREEQ98
                } /* if (!phreeq98_debug) */
#endif
	}

/*
 *   Display successful status
 */
	this->PhreeqcPtr->do_status();

/*
 *   call post-run callback
 */
	if (pfn_post)
	{
		pfn_post(cookie);
	}

	this->UpdateComponents = true;
	this->update_errors();

	// update lines
	//

	if (this->LogStringOn)
	{
		// output lines
		std::istringstream iss(this->LogString);
		std::string line;
		while (std::getline(iss, line))
		{
			this->LogLines.push_back(line);
		}
	}

	if (this->OutputStringOn)
	{
		// output lines
		std::istringstream iss(this->OutputString);
		std::string line;
		while (std::getline(iss, line))
		{
			this->OutputLines.push_back(line);
		}
	}

	for (int i = 0; i < this->GetSelectedOutputCount(); ++i)
	{
		int n = this->GetNthSelectedOutputUserNumber(i);
		if (this->get_sel_out_string_on(n))
		{
			std::map< int, std::string >::iterator mit = this->SelectedOutputStringMap.find(n);
			ASSERT(mit != this->SelectedOutputStringMap.end());
			if (mit != this->SelectedOutputStringMap.end())
			{
				// output lines
				std::istringstream iss((*mit).second);
				std::string line;
				while (std::getline(iss, line))
				{
					this->SelectedOutputLinesMap[(*mit).first].push_back(line);
				}
			}
		}
	}
}

void IPhreeqc::update_errors(void)
{
	this->ErrorLines.clear();
	this->ErrorString = ((CErrorReporter<std::ostringstream>*)this->ErrorReporter)->GetOS()->str();
	if (this->ErrorString.size())
	{
		std::istringstream iss(this->ErrorString);
		std::string line;
		while (std::getline(iss, line))
		{
			this->ErrorLines.push_back(line);
		}
	}

	this->WarningLines.clear();
	this->WarningString = ((CErrorReporter<std::ostringstream>*)this->WarningReporter)->GetOS()->str();
	if (this->WarningString.size())
	{
		std::istringstream iss(this->WarningString);
		std::string line;
		while (std::getline(iss, line))
		{
			this->WarningLines.push_back(line);
		}
	}
}

void IPhreeqc::log_msg(const char * str)
{
	if (this->LogStringOn && this->log_on)
	{
		this->LogString += str;
	}
	ASSERT(!(this->LogFileOn != (this->log_ostream != 0)));
	this->PHRQ_io::log_msg(str);
}

void IPhreeqc::error_msg(const char *str, bool stop)
{
	ASSERT(this->error_ostream == &std::cerr || (!(this->ErrorFileOn != (this->error_ostream != 0))));

	if (this->error_ostream != NULL && this->error_on)
	{
		(*this->error_ostream) << str;
	}
	bool save_error_on = this->error_on;
	this->error_on = false;
	this->PHRQ_io::error_msg(str);
	this->error_on = save_error_on;

	if (this->ErrorStringOn && this->error_on)
	{
		this->AddError(str);
	}
	if (stop)
	{
		if (this->error_ostream && this->error_on)
		{
			(*this->error_ostream) << "Stopping.\n";
			this->error_ostream->flush();
		}
		throw IPhreeqcStop();
	}
}

void IPhreeqc::warning_msg(const char *str)
{
	ASSERT(this->error_ostream == &std::cerr || (!(this->ErrorFileOn != (this->error_ostream != 0))));

	if (this->error_ostream != NULL && this->error_on)
	{
		(*this->error_ostream) << str << "\n";
	}
	bool save_error_on = this->error_on;
	this->error_on = false;
	this->PHRQ_io::warning_msg(str);
	this->error_on = save_error_on;

	std::ostringstream oss;
	oss << str << std::endl;
	if (this->WarningStringOn)
	{
		this->AddWarning(oss.str().c_str());
	}
}

void IPhreeqc::output_msg(const char * str)
{
	if (this->OutputStringOn && this->output_on)
	{
		this->OutputString += str;
	}
	ASSERT(!(this->OutputFileOn != (this->output_ostream != 0)));
	this->PHRQ_io::output_msg(str);
}

void IPhreeqc::screen_msg(const char *err_str)
{
	// no-op
}

void IPhreeqc::punch_msg(const char *str)
{
	if (this->get_sel_out_string_on(this->PhreeqcPtr->current_selected_output->Get_n_user()) && this->punch_on)
	{
		ASSERT(this->PhreeqcPtr->current_selected_output);
		if (this->PhreeqcPtr->current_selected_output)
		{
			this->SelectedOutputStringMap[this->PhreeqcPtr->current_selected_output->Get_n_user()] += str;
		}
	}
	ASSERT(!(this->SelectedOutputFileOnMap[this->PhreeqcPtr->current_selected_output->Get_n_user()] != (this->PhreeqcPtr->current_selected_output->Get_punch_ostream() != 0)));
	this->PHRQ_io::punch_msg(str);
}

void IPhreeqc::open_output_files(const char* sz_routine)
{
	if (this->OutputFileOn)
	{
		if (this->output_ostream != NULL)
		{
			safe_close(&this->output_ostream);
		}
		if (!this->output_ostream)
		{
			if ( (this->output_ostream = new std::ofstream(this->OutputFileName.c_str())) == NULL)
			{
				std::ostringstream oss;
				oss << sz_routine << ": Unable to open:" << "\"" << this->OutputFileName << "\".\n";
				this->warning_msg(oss.str().c_str());
			}
		}
	}
	if (this->ErrorFileOn)
	{
		if (this->error_ostream != NULL)
		{
			safe_close(&this->error_ostream);
		}
		if (!this->error_ostream)
		{
			if ( (this->error_ostream = new std::ofstream(this->ErrorFileName.c_str())) == NULL)
			{
				std::ostringstream oss;
				oss << sz_routine << ": Unable to open:" << "\"" << this->ErrorFileName << "\".\n";
				this->warning_msg(oss.str().c_str());
			}
		}
	}
	if (this->LogFileOn)
	{
		if (this->log_ostream != NULL)
		{
			safe_close(&this->log_ostream);
		}
		if (!this->log_ostream)
		{
			if ( (this->log_ostream = new std::ofstream(this->LogFileName.c_str())) == NULL)
			{
				std::ostringstream oss;
				oss << sz_routine << ": Unable to open:" << "\"" << this->LogFileName << "\".\n";
				this->warning_msg(oss.str().c_str());
			}
		}
	}
}

int IPhreeqc::close_input_files(void)
{
	int i = 0;
	if (this->database_file)
	{
		i |= fclose(this->database_file);
	}
	if (this->input_file)
	{
		i |= fclose(this->input_file);
	}
	this->input_file = this->database_file = NULL;
	return (i);
}

int IPhreeqc::close_output_files(void)
{
	int ret = 0;

	safe_close(&this->output_ostream);
	safe_close(&this->log_ostream);
	safe_close(&this->dump_ostream);
	safe_close(&this->error_ostream);

	std::map< int, SelectedOutput >::iterator it = this->PhreeqcPtr->SelectedOutput_map.begin();
	for (; it != this->PhreeqcPtr->SelectedOutput_map.end(); ++it)
	{
		std::ostream *ptr = (*it).second.Get_punch_ostream();
		safe_close(&ptr);
		(*it).second.Set_punch_ostream(NULL);
	}

	this->punch_ostream = 0;

	return ret;
}

void IPhreeqc::fpunchf(const char *name, const char *format, double d)
{
	try
	{
		this->PHRQ_io::fpunchf(name, format, d);
		if (this->get_sel_out_string_on(this->PhreeqcPtr->current_selected_output->Get_n_user()) && this->punch_on)
		{
			ASSERT(this->SelectedOutputStringMap.find(this->PhreeqcPtr->current_selected_output->Get_n_user()) != this->SelectedOutputStringMap.end());
			PHRQ_io::fpunchf_helper(&(this->SelectedOutputStringMap[this->PhreeqcPtr->current_selected_output->Get_n_user()]), format, d);
		}
		ASSERT(this->SelectedOutputMap.find(this->PhreeqcPtr->current_selected_output->Get_n_user()) != this->SelectedOutputMap.end());
		this->SelectedOutputMap[this->PhreeqcPtr->current_selected_output->Get_n_user()]->PushBackDouble(name, d);
	}
	catch (const std::bad_alloc&)
	{
		this->PhreeqcPtr->malloc_error();
	}
}

void IPhreeqc::fpunchf(const char *name, const char *format, char *s)
{
	try
	{
		this->PHRQ_io::fpunchf(name, format, s);
		if (this->get_sel_out_string_on(this->PhreeqcPtr->current_selected_output->Get_n_user()) && this->punch_on)
		{
			ASSERT(this->SelectedOutputStringMap.find(this->PhreeqcPtr->current_selected_output->Get_n_user()) != this->SelectedOutputStringMap.end());
			PHRQ_io::fpunchf_helper(&(this->SelectedOutputStringMap[this->PhreeqcPtr->current_selected_output->Get_n_user()]), format, s);
		}
		ASSERT(this->SelectedOutputMap.find(this->PhreeqcPtr->current_selected_output->Get_n_user()) != this->SelectedOutputMap.end());
		this->SelectedOutputMap[this->PhreeqcPtr->current_selected_output->Get_n_user()]->PushBackString(name, s);
	}
	catch (const std::bad_alloc&)
	{
		this->PhreeqcPtr->malloc_error();
	}
}

void IPhreeqc::fpunchf(const char *name, const char *format, int i)
{
	try
	{
		this->PHRQ_io::fpunchf(name, format, i);
		if (this->get_sel_out_string_on(this->PhreeqcPtr->current_selected_output->Get_n_user()) && this->punch_on)
		{
			ASSERT(this->SelectedOutputStringMap.find(this->PhreeqcPtr->current_selected_output->Get_n_user()) != this->SelectedOutputStringMap.end());
			PHRQ_io::fpunchf_helper(&(this->SelectedOutputStringMap[this->PhreeqcPtr->current_selected_output->Get_n_user()]), format, i);
		}
		ASSERT(this->SelectedOutputMap.find(this->PhreeqcPtr->current_selected_output->Get_n_user()) != this->SelectedOutputMap.end());
		this->SelectedOutputMap[this->PhreeqcPtr->current_selected_output->Get_n_user()]->PushBackLong(name, (long)i);
	}
	catch (const std::bad_alloc&)
	{
		this->PhreeqcPtr->malloc_error();
	}
}

void IPhreeqc::fpunchf_end_row(const char *format)
{
	this->EndRow();
}

bool IPhreeqc::get_sel_out_file_on(int n)const
{
	// if not found in list SelectedOutputFileOn is false
	std::map< int, bool >::const_iterator ci = this->SelectedOutputFileOnMap.find(n);
	if (ci != this->SelectedOutputFileOnMap.end())
	{
		return (*ci).second;
	}
	return false;
}

bool IPhreeqc::get_sel_out_string_on(int n)const
{
	// if not found in list SelectedOutputStringOn is false
	std::map< int, bool >::const_iterator ci = this->SelectedOutputStringOn.find(this->CurrentSelectedOutputUserNumber);
	if (ci != this->SelectedOutputStringOn.end())
	{
		return (*ci).second;
	}
	return false;
}

bool IPhreeqc::punch_open(const char *file_name, std::ios_base::openmode mode, int n_user)
{
	if (this->PhreeqcPtr->SelectedOutput_map[n_user].Get_have_punch_name() &&
		!this->PhreeqcPtr->SelectedOutput_map[n_user].Get_file_name().empty())
	{
		this->SelectedOutputFileNameMap[n_user] = this->PhreeqcPtr->SelectedOutput_map[n_user].Get_file_name();
	}
	else if (this->SelectedOutputFileNameMap[n_user].empty())
	{
		this->SelectedOutputFileNameMap[n_user] = this->sel_file_name(n_user);
	}
	if (this->get_sel_out_file_on(n_user))
	{
		ASSERT(!this->SelectedOutputFileNameMap[n_user].empty());
		this->PhreeqcPtr->SelectedOutput_map[n_user].Set_file_name(this->SelectedOutputFileNameMap[n_user]);
		return this->ofstream_open(&punch_ostream, this->SelectedOutputFileNameMap[n_user].c_str(), mode);
	}
	return true;
}

bool IPhreeqc::output_open(const char *file_name, std::ios_base::openmode mode)
{
	if (this->OutputFileOn)
	{
		return this->PHRQ_io::output_open(file_name, mode);
	}
	return true;
}

std::string IPhreeqc::sel_file_name(int n_user)
{
	std::ostringstream oss;
	oss << "selected_" << n_user << "." << this->Index << ".out";
	return oss.str();
}

std::string IPhreeqc::create_file_name(const char *prefix, const char *suffix)
{
	std::ostringstream oss;
	oss << prefix << "." << this->Index << "." << suffix;
	return oss.str();
}
