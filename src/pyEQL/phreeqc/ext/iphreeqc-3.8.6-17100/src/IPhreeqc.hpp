/*! @file IPhreeqc.hpp
	@brief C++ Documentation
*/

#ifndef INC_IPHREEQC_HPP
#define INC_IPHREEQC_HPP

#include <exception>
#include <list>
#include <vector>
#include <map>
#include <cstdarg>
#include "IPhreeqcCallbacks.h"      /* PFN_PRERUN_CALLBACK, PFN_POSTRUN_CALLBACK, PFN_CATCH_CALLBACK */
#include "Var.h"                    /* VRESULT */
#include "PHRQ_io.h"

#include "PHRQ_exports.h"

class Phreeqc;
class IErrorReporter;
class CSelectedOutput;
class SelectedOutput;

/**
 * @class IPhreeqcStop
 *
 * @brief This class is derived from std::exception and is thrown
 * when an unrecoverable error has occurred.
 */
class IPQ_DLL_EXPORT IPhreeqcStop : public std::exception 
{
public:
  virtual const char *what() const throw () {return "Failure in IPhreeqc\n";}
};

/**
 * @class IPhreeqc
 *
 * @brief Provides an interface to PHREEQC (Version 3)--A Computer
 * Program for Speciation, Batch-Reaction, One-Dimensional Transport,
 * and Inverse Geochemical Calculations
 */

class IPQ_DLL_EXPORT IPhreeqc : public PHRQ_io
{
public:
	/**
	 * Constructor.
	 *  @anchor IPhreeqc_cpp
	 *  @par Example:
	 *  @include IPhreeqc.cpp
	 */
	IPhreeqc(void);

	/**
	 * Destructor
	 */
	virtual ~IPhreeqc(void);

public:

	/**
	 *  Accumlulate line(s) for input to phreeqc.
	 *  @param line             The line(s) to add for input to phreeqc.
	 *  @retval VR_OK           Success
	 *  @retval VR_OUTOFMEMORY  Out of memory
	 *  @see                    ClearAccumulatedLines, OutputAccumulatedLines, RunAccumulated
	 */
	VRESULT                  AccumulateLine(const char *line);

	/**
	 *  Appends the given error message and increments the error count.
	 *  Internally used to create an error condition.
	 *  @param error_msg        The error message to display.
	 *  @return                 The current error count.
	 *  @see                    GetErrorString, GetErrorStringLine, GetErrorStringLineCount, OutputErrorString
	 */
	size_t                   AddError(const char* error_msg);

	/**
	 *  Appends the given warning message and increments the warning count.
	 *  Internally used to create a warning condition.
	 *  @param warning_msg      The warning message to display.
	 *  @return                 The current warning count.
	 *  @see                    GetWarningString, GetWarningStringLine, GetWarningStringLineCount, OutputWarningString
	 */
	size_t                   AddWarning(const char* warning_msg);

	/**
	 *  Clears the accumulated input buffer.  Input buffer is accumulated from calls to @ref AccumulateLine.
	 *  @see                    AccumulateLine, GetAccumulatedLines, OutputAccumulatedLines, RunAccumulated
	 */
	void                     ClearAccumulatedLines(void);

	/**
	 *  Retrieve the accumulated input string.  The accumulated input string can be run
	 *  with @ref RunAccumulated.
	 *  @return                 The accumulated input string.
	 *  @see                    AccumulateLine, ClearAccumulatedLines, OutputAccumulatedLines, RunAccumulated
	 */
	const std::string&       GetAccumulatedLines(void);

	/**
	 *  Retrieves the given component.
	 *  @param n                The zero-based index of the component to retrieve.
	 *  @return                 A null terminated string containing the given component.
	 *                          Returns an empty string if n is out of range.
	 *  @see                    GetComponentCount, ListComponents
	 */
	const char*              GetComponent(int n);

	/**
	 *  Retrieves the number of components in the current list of components.
	 *  @return                 The current count of components.
	 *  @see                    GetComponent, ListComponents
	 */
	size_t                   GetComponentCount(void);

	/**
	 *  Retrieves the current <B>SELECTED_OUTPUT</B> user number.  The initial setting is 1.
	 *  @return                 The current <b>SELECTED_OUTPUT</b> user number.
	 *  @see                    GetSelectedOutputColumnCount, GetSelectedOutputFileName, GetSelectedOutputRowCount, GetSelectedOutputString, GetSelectedOutputStringLine, GetSelectedOutputStringLineCount, GetSelectedOutputValue, SetCurrentSelectedOutputUserNumber
	 */
	int                      GetCurrentSelectedOutputUserNumber(void)const;

	/**
	 *  Retrieves the name of the dump file.  This file name is used if not specified within <B>DUMP</B> input.
	 *  The default value is <B><I>dump.id.out</I></B>, where id is obtained from @ref GetId.
	 *  @return filename        The name of the file to write <B>DUMP</B> output to.
	 *  @see                    GetDumpFileOn, GetDumpString, GetDumpStringOn, GetDumpStringLine, GetDumpStringLineCount, SetDumpFileName, SetDumpFileOn, SetDumpStringOn
	 */
	const char*              GetDumpFileName(void)const;

	/**
	 *  Retrieves the current value of the dump file switch.
	 *  @retval true            Output is written to the <B>DUMP</B> (<B><I>dump.id.out</I></B> if unspecified, where id is obtained from @ref GetId) file.
	 *  @retval false           No output is written.
	 *  @see                    GetDumpStringLine, GetDumpStringLineCount, GetDumpStringOn, GetDumpString, SetDumpFileOn, SetDumpStringOn
	 */
	bool                     GetDumpFileOn(void)const;

	/**
	 *  Retrieves the string buffer containing <b>DUMP</b> output.
	 *  @return                 A null terminated string containing <b>DUMP</b> output.
	 *  @pre                    
	 *      @ref SetDumpStringOn must have been set to true in order to receive <b>DUMP</b> output.
	 *  @see                    GetDumpStringLine, GetDumpFileOn, GetDumpStringLineCount, GetDumpStringOn, SetDumpFileOn, SetDumpStringOn
	 */
	const char*              GetDumpString(void)const;

	/**
	 *  Retrieves the given dump line.
	 *  @param n                The zero-based index of the line to retrieve.
	 *  @return                 A null terminated string containing the given line.
	 *                          Returns an empty string if n is out of range.
	 *  @pre                    @ref SetDumpStringOn must have been set to true.
	 *  @see                    GetDumpFileOn, GetDumpString, GetDumpStringLineCount, GetDumpStringOn, SetDumpFileOn, SetDumpStringOn
	 */
	const char*              GetDumpStringLine(int n);

	/**
	 *  Retrieves the number of lines in the current dump string buffer.
	 *  @return                 The number of lines.
	 *  @pre                    @ref SetDumpStringOn must have been set to true.
	 *  @see                    GetDumpFileOn, GetDumpString, GetDumpStringLine, GetDumpStringOn, SetDumpFileOn, SetDumpStringOn
	 */
	int                      GetDumpStringLineCount(void)const;

	/**
	 *  Retrieves the current value of the dump string switch.
	 *  @retval true            Output defined by the <B>DUMP</B> keyword is stored.
	 *  @retval false           No output is stored.
	 *  @see                    GetDumpFileOn, GetDumpString, GetDumpStringLine, GetDumpStringLineCount, SetDumpFileOn, SetDumpStringOn
	 */
	bool                     GetDumpStringOn(void)const;

	/**
	 *  Retrieves the name of the error file. The default value is <B><I>phreeqc.id.err</I></B>, where id is obtained from @ref GetId.
	 *  @return filename        The name of the file to write to.
	 *  @see                    GetErrorFileOn, GetErrorString, GetErrorStringOn, GetErrorStringLine, GetErrorStringLineCount, SetErrorFileName, SetErrorFileOn, SetErrorStringOn
	 */
	const char*              GetErrorFileName(void)const;

	/**
	 *  Retrieves the current value of the error file switch.
	 *  @retval true            Errors are written to the <B><I>phreeqc.id.err</I></B> (where id is obtained from @ref GetId) file.
	 *  @retval false           No errors are written.
	 *  @see                    SetErrorFileOn
	 */
	bool                     GetErrorFileOn(void)const;

	/**
	 *  Retrieves the current value of the error switch.
	 *  @retval true            Error messages are sent to the error file and to the string buffer
	 *  @retval false           No errors are sent.
	 *  @see                    SetErrorOn
	 */
	bool                     GetErrorOn(void)const;

	/**
	 *  Retrieves the error messages from the last call to @ref RunAccumulated, @ref RunFile, @ref RunString, @ref LoadDatabase, or @ref LoadDatabaseString.
	 *  @return                 A null terminated string containing error messages.
	 *  @see                    GetErrorStringLine, GetErrorStringLineCount, GetErrorFileOn, OutputErrorString, SetErrorFileOn
	 */
	const char*              GetErrorString(void);

	/**
	 *  Retrieves the given error line.
	 *  @return                 A null terminated string containing the given line of the error string buffer.
	 *  @param n                The zero-based index of the line to retrieve.
	 *  @see                    GetErrorStringLineCount, OutputErrorString
	 */
	const char*              GetErrorStringLine(int n);

	/**
	 *  Retrieves the number of lines in the current error string buffer.
	 *  @return                 The number of lines.
	 *  @see                    GetErrorStringLine, OutputErrorString
	 */
	int                      GetErrorStringLineCount(void)const;

	/**
	 *  Retrieves the current value of the error string switch.
	 *  @retval true            Error output is stored.
	 *  @retval false           No error output is stored.
	 *  @see                    GetErrorFileOn, GetErrorString, GetErrorStringLine, GetErrorStringLineCount, SetErrorFileOn, SetErrorStringOn
	 */
	bool                     GetErrorStringOn(void)const;

	/**
	 *  Retrieves the id of this object.  Each instance receives an id which is incremented for each instance
	 *  starting with the value zero.
	 *  @return                 The id.
	 */
	int                      GetId(void)const;

	/**
	 *  Retrieves the name of the log file. The default value is <B><I>phreeqc.id.log</I></B>, where id is obtained from @ref GetId.
	 *  @return filename        The name of the file to write to.
	 *  @see                    GetLogFileOn, GetLogString, GetLogStringOn, GetLogStringLine, GetLogStringLineCount, SetLogFileName, SetLogFileOn, SetLogStringOn
	 */
	const char*              GetLogFileName(void)const;

	/**
	 *  Retrieves the current value of the log file switch.
	 *  @retval true            Log messages are written to the <B><I>phreeqc.id.log</I></B> (where id is obtained from @ref GetId) file.
	 *  @retval false           No log messages are written.
	 *  @remarks
	 *      Logging must be enabled through the use of the KNOBS -logfile option in order to receive any log messages.
	 *  @see                    SetLogFileOn
	 */
	bool                     GetLogFileOn(void)const;

	/**
	 *  Retrieves the string buffer containing phreeqc log output.
	 *  @return                 A null terminated string containing log output.
	 *  @pre
	 *      @ref SetLogStringOn must have been set to true and enabled through the use of the KNOBS -logfile option in order to receive any log messages.
	 *  @see                    GetLogStringLine, GetLogFileOn, GetLogStringLineCount, GetLogStringOn, SetLogFileOn, SetLogStringOn
	 */
	const char*              GetLogString(void)const;

	/**
	 *  Retrieves the given log line.
	 *  @param n                The zero-based index of the line to retrieve.
	 *  @return                 A null terminated string containing the given line.
	 *                          Returns an empty string if n is out of range.
	 *  @pre                    @ref SetLogStringOn must have been set to true and enabled through the use of the KNOBS -logfile option in order to receive any log messages.
	 *  @see                    GetLogFileOn, GetLogString, GetLogStringLineCount, GetLogStringOn, SetLogFileOn, SetLogStringOn
	 */
	const char*              GetLogStringLine(int n)const;

	/**
	 *  Retrieves the number of lines in the current log string buffer.
	 *  @return                 The number of lines.
	 *  @pre                    @ref SetLogStringOn must have been set to true and enabled through the use of the KNOBS -logfile option in order to receive any log messages.
	 *  @see                    GetLogFileOn, GetLogString, GetLogStringLine, GetLogStringOn, SetLogFileOn, SetLogStringOn
	 */
	int                      GetLogStringLineCount(void)const;

	/**
	 *  Retrieves the current value of the log string switch.
	 *  @retval true            Log output is stored.
	 *  @retval false           No log output is stored.
	 *  @see                    GetLogFileOn, GetLogString, GetLogStringLine, GetLogStringLineCount, SetLogFileOn, SetLogStringOn
	 */
	bool                     GetLogStringOn(void)const;

	/**
	 *  Retrieves the nth user number of the currently defined <B>SELECTED_OUTPUT</B> blocks.
	 *  @param n                The zero-based index of the <B>SELECTED_OUTPUT</B> user number to retrieve.
	 *  @return                 The nth defined user number; a negative value indicates an error occurred.
	 *  @see                    GetCurrentSelectedOutputUserNumber, GetSelectedOutputCount, SetCurrentSelectedOutputUserNumber
	 */
	int                      GetNthSelectedOutputUserNumber(int n)const;

	/**
	 *  Retrieves the name of the output file. The default value is <B><I>phreeqc.id.out</I></B>, where id is obtained from @ref GetId.
	 *  @return filename        The name of the file to write phreeqc output to.
	 *  @see                    GetOutputFileOn, GetOutputString, GetOutputStringOn, GetOutputStringLine, GetOutputStringLineCount, SetOutputFileName, SetOutputFileOn, SetOutputStringOn
	 */
	const char*              GetOutputFileName(void)const;

	/**
	 *  Retrieves the current value of the output file switch.
	 *  @retval true            Output is written to the <B><I>phreeqc.id.out</I></B> (where id is obtained from @ref GetId) file.
	 *  @retval false           No output is written.
	 *  @see                    GetOutputFileOn, GetOutputString, GetOutputStringOn, GetOutputStringLine, GetOutputStringLineCount, SetOutputFileName, SetOutputFileOn, SetOutputStringOn
	 */
	bool                     GetOutputFileOn(void)const;

	/**
	 *  Retrieves the string buffer containing phreeqc output.
	 *  @return                 A null terminated string containing phreeqc output.
	 *  @pre
	 *      @ref SetOutputStringOn must have been set to true in order to receive output.
	 *  @see                    GetOutputStringLine, GetOutputFileOn, GetOutputStringLineCount, GetOutputStringOn, SetOutputFileOn, SetOutputStringOn
	 */
	const char*              GetOutputString(void)const;

	/**
	 *  Retrieves the given output line.
	 *  @param n                The zero-based index of the line to retrieve.
	 *  @return                 A null terminated string containing the given line.
	 *                          Returns an empty string if n is out of range.
	 *  @pre                    @ref SetOutputStringOn must have been set to true.
	 *  @see                    GetOutputFileOn, GetOutputString, GetOutputStringLineCount, GetOutputStringOn, SetOutputFileOn, SetOutputStringOn
	 */
	const char*              GetOutputStringLine(int n)const;

	/**
	 *  Retrieves the number of lines in the current output string buffer.
	 *  @return                 The number of lines.
	 *  @pre                    @ref SetOutputStringOn must have been set to true.
	 *  @see                    GetOutputFileOn, GetOutputString, GetOutputStringLine, GetOutputStringOn, SetOutputFileOn, SetOutputStringOn
	 */
	int                      GetOutputStringLineCount(void)const;

	/**
	 *  Retrieves the current value of the output string switch.
	 *  @retval true            Phreeqc output is stored.
	 *  @retval false           No phreeqc output is stored.
	 *  @see                    GetOutputFileOn, GetOutputString, GetOutputStringLine, GetOutputStringLineCount, SetOutputFileOn, SetOutputStringOn
	 */
	bool                     GetOutputStringOn(void)const;

	/**
	 *  Retrieves the number of columns in the current selected-output buffer (see @ref SetCurrentSelectedOutputUserNumber).
	 *  @return                 The number of columns.
	 *  @see                    GetCurrentSelectedOutputUserNumber, GetSelectedOutputRowCount, GetSelectedOutputValue, SetCurrentSelectedOutputUserNumber
	 */
	int                      GetSelectedOutputColumnCount(void)const;

	/**
	 *  Retrieves the count of <B>SELECTED_OUTPUT</B> blocks that are currently defined.
	 *  @return                 The number of <B>SELECTED_OUTPUT</B> blocks.
	 *  @see                    GetCurrentSelectedOutputUserNumber, GetNthSelectedOutputUserNumber, SetCurrentSelectedOutputUserNumber
	 */
	int                      GetSelectedOutputCount(void)const;

	/**
	 *  Retrieves the name of the current selected output file (see @ref SetCurrentSelectedOutputUserNumber).  This file name is used if not specified within <B>SELECTED_OUTPUT</B> input.
	 *  The default value is <B><I>selected_n.id.out</I></B>, where id is obtained from @ref GetId.
	 *  @return filename        The name of the file to write to.
	 *  @see                    GetCurrentSelectedOutputUserNumber, GetSelectedOutputFileOn, GetSelectedOutputString, GetSelectedOutputStringOn, GetSelectedOutputStringLine, GetSelectedOutputStringLineCount, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileName, SetSelectedOutputFileOn, SetSelectedOutputStringOn
	 */
	const char*              GetSelectedOutputFileName(void)const;

	/**
	 *  Retrieves the current selected-output file switch (see @ref SetCurrentSelectedOutputUserNumber).
	 *  @retval true            Output is written to the selected-output (<B><I>selected_n.id.out</I></B> if unspecified, where id is obtained from @ref GetId) file.
	 *  @retval false           No output is written.
	 *  @see                    GetSelectedOutputValue, GetSelectedOutputColumnCount, GetSelectedOutputRowCount, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn
	 */
	bool                     GetSelectedOutputFileOn(void)const;

	/**
	 *  Retrieves the number of rows in the current selected-output buffer (see @ref SetCurrentSelectedOutputUserNumber).
	 *  @return                 The number of rows.
	 *  @see                    GetCurrentSelectedOutputUserNumber, GetSelectedOutputColumnCount, GetSelectedOutputFileOn, GetSelectedOutputValue, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn
	 */
	int                      GetSelectedOutputRowCount(void)const;

	/**
	 *  Retrieves the string buffer containing <b>SELECTED_OUTPUT</b> for the currently selected user number (see @ref SetCurrentSelectedOutputUserNumber).
	 *  @return                 A null terminated string containing <b>SELECTED_OUTPUT</b>.
	 *  @pre
	 *      @ref SetSelectedOutputStringOn must have been set to true in order to receive <b>SELECTED_OUTPUT</b>.
	 *  @see                    GetCurrentSelectedOutputUserNumber, GetSelectedOutputStringLine, GetSelectedOutputFileOn, GetSelectedOutputStringLineCount, GetSelectedOutputStringOn, GetSelectedOutputString, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn, SetSelectedOutputStringOn
	 */
	const char*              GetSelectedOutputString(void)const;

	/**
	 *  Retrieves the given selected output line of the currently selected user number (see @ref SetCurrentSelectedOutputUserNumber).
	 *  @param n                The zero-based index of the line to retrieve.
	 *  @return                 A null terminated string containing the given line.
	 *                          Returns an empty string if n is out of range.
	 *  @pre                    @ref SetSelectedOutputStringOn must have been set to true.
	 *  @see                    GetCurrentSelectedOutputUserNumber, GetSelectedOutputFileOn, GetSelectedOutputString, GetSelectedOutputStringLineCount, GetSelectedOutputStringOn, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn, SetSelectedOutputStringOn
	 */
	const char*              GetSelectedOutputStringLine(int n);

	/**
	 *  Retrieves the number of lines in the current selected output string buffer (see @ref SetCurrentSelectedOutputUserNumber).
	 *  @return                 The number of lines.
	 *  @pre                    @ref SetSelectedOutputStringOn must have been set to true.
	 *  @see                    GetCurrentSelectedOutputUserNumber, GetSelectedOutputFileOn, GetSelectedOutputString, GetSelectedOutputStringLine, GetSelectedOutputStringOn, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn, SetSelectedOutputStringOn
	 */
	int                      GetSelectedOutputStringLineCount(void)const;

	/**
	 *  Retrieves the value of the current selected output string switch (see @ref SetCurrentSelectedOutputUserNumber).
	 *  @retval true            Output defined by the <B>SELECTED_OUTPUT</B> keyword is stored.
	 *  @retval false           No output is stored.
	 *  @see                    GetSelectedOutputFileOn, GetSelectedOutputString, GetSelectedOutputStringLine, GetSelectedOutputStringLineCount, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn, SetSelectedOutputStringOn
	 */
	bool                     GetSelectedOutputStringOn(void)const;

	/**
	 *  Returns the @c VAR associated with the specified row and column.  The current <b>SELECTED_OUTPUT</b> block is set using the @ref SetCurrentSelectedOutputUserNumber method.
	 *  @param row              The row index.
	 *  @param col              The column index.
	 *  @param pVAR             Pointer to the @c VAR to receive the requested data.
	 *  @retval VR_OK           Success.
	 *  @retval VR_INVALIDROW   The given row is out of range.
	 *  @retval VR_INVALIDCOL   The given column is out of range.
	 *  @retval VR_OUTOFMEMORY  Memory could not be allocated.
	 *  @retval VR_BADINSTANCE  The given id is invalid.
	 *  @see                    GetCurrentSelectedOutputUserNumber, GetSelectedOutputColumnCount, GetSelectedOutputFileOn, GetSelectedOutputRowCount, GetSelectedOutputValue2, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn
	 *  @remarks
	 *      Row 0 contains the column headings to the selected_ouput.
	 *  @par Examples:
	 *  The headings will include a suffix and/or prefix in order to differentiate the
	 *  columns.
	 *  @htmlonly
	<p>
	<table border=1>

	<TR VALIGN="top">
	<TH width=65%>
	Input
	</TH>
	<TH width=35%>
	Headings
	</TH>
	</TR>

	<TR VALIGN="top">
	<TD width=65%>
	<CODE><PRE>
	  SELECTED_OUTPUT
		-reset false
		-totals Ca Na
	</PRE></CODE>
	</TD>
	<TD width=35%>
	<CODE><PRE>
	  Ca(mol/kgw)  Na(mol/kgw)
	</PRE></CODE>
	</TD>
	</TR>

	<TR VALIGN="top">
	<TD width=65%>
	<CODE><PRE>
	  SELECTED_OUTPUT
		-reset false
		-molalities Fe+2 Hfo_sOZn+
	</PRE></CODE>
	</TD>
	<TD width=35%>
	<CODE><PRE>
	  m_Fe+2(mol/kgw)  m_Hfo_sOZn+(mol/kgw)
	</PRE></CODE>
	</TD>
	</TR>

	<TR VALIGN="top">
	<TD width=65%>
	<CODE><PRE>
	  SELECTED_OUTPUT
		-reset false
		-activities H+ Ca+2
	</PRE></CODE>
	</TD>
	<TD width=35%>
	<CODE><PRE>
	  la_H+  la_Ca+2
	</PRE></CODE>
	</TD>
	</TR>

	<TR VALIGN="top">
	<TD width=65%>
	<CODE><PRE>
	  SELECTED_OUTPUT
		-reset false
		-equilibrium_phases Calcite Dolomite
	</PRE></CODE>
	</TD>
	<TD width=35%>
	<CODE><PRE>
	  Calcite  d_Calcite  Dolomite  d_Dolomite
	</PRE></CODE>
	</TD>
	</TR>

	<TR VALIGN="top">
	<TD width=65%>
	<CODE><PRE>
	  SELECTED_OUTPUT
		-reset false
		-saturation_indices CO2(g) Siderite
	</PRE></CODE>
	</TD>
	<TD width=35%>
	<CODE><PRE>
	  si_CO2(g)  si_Siderite
	</PRE></CODE>
	</TD>
	</TR>

	<TR VALIGN="top">
	<TD width=65%>
	<CODE><PRE>
	  SELECTED_OUTPUT
		-reset false
		-gases CO2(g) N2(g)
	</PRE></CODE>
	</TD>
	<TD width=35%>
	<CODE><PRE>
	  pressure "total mol" volume g_CO2(g) g_N2(g)
	</PRE></CODE>
	</TD>
	</TR>

	<TR VALIGN="top">
	<TD width=65%>
	<CODE><PRE>
	  SELECTED_OUTPUT
		-reset false
		-kinetic_reactants CH2O Pyrite
	</PRE></CODE>
	</TD>
	<TD width=35%>
	<CODE><PRE>
	  k_CH2O dk_CH2O k_Pyrite dk_Pyrite
	</PRE></CODE>
	</TD>
	</TR>

	<TR VALIGN="top">
	<TD width=65%>
	<CODE><PRE>
	  SELECTED_OUTPUT
		-reset false
		-solid_solutions CaSO4 SrSO4
	</PRE></CODE>
	</TD>
	<TD width=35%>
	<CODE><PRE>
	  s_CaSO4 s_SrSO4
	</PRE></CODE>
	</TD>
	</TR>

	</table>
	 *  @endhtmlonly
	 */
	VRESULT                  GetSelectedOutputValue(int row, int col, VAR* pVAR);

	/**
	 *  Returns the associated data with the specified row and column.  The current <b>SELECTED_OUTPUT</b> block is set using the @ref SetCurrentSelectedOutputUserNumber method.
	 *  @param row               The row index.
	 *  @param col               The column index.
	 *  @param vtype             Receives the variable type.  See @ref VAR_TYPE.
	 *  @param dvalue            Receives the numeric value when (VTYPE=@ref TT_DOUBLE) or (VTYPE=@ref TT_LONG).
	 *  @param svalue            Receives the string variable when (VTYPE=@ref TT_STRING).  When (VTYPE=@ref TT_DOUBLE) or (VTYPE=@ref TT_LONG) this variable is filled with a string equivalent of DVALUE.
	 *  @param svalue_length     The length of the svalue buffer.
	 *  @retval IPQ_OK           Success.
	 *  @retval IPQ_INVALIDROW   The given row is out of range.
	 *  @retval IPQ_INVALIDCOL   The given column is out of range.
	 *  @retval IPQ_OUTOFMEMORY  Memory could not be allocated.
	 *  @retval IPQ_BADINSTANCE  The given id is invalid.
	 *  @see                     GetSelectedOutputFileOn, GetSelectedOutputColumnCount, GetSelectedOutputRowCount, GetSelectedOutputValue, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn
	 *  @remarks
	 *  Row 0 contains the column headings to the selected_ouput.
	 *  @par Examples:
	 *  See @ref GetSelectedOutputValue.
	 */
	VRESULT                  GetSelectedOutputValue2(int row, int col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length);

	/**
	 *  Retrieves the string buffer containing the version in the form of X.X.X-XXXX.
	 *  @return                 A null terminated string containing the IPhreeqc version number.
	 */
	static const char*       GetVersionString(void);

	/**
	 *  Retrieves the warning messages from the last call to @ref RunAccumulated, @ref RunFile, @ref RunString, @ref LoadDatabase, or @ref LoadDatabaseString.
	 *  @return                 A null terminated string containing warning messages.
	 *  @see                    GetWarningStringLine, GetWarningStringLineCount, OutputWarningString
	 */
	const char*              GetWarningString(void);

	/**
	 *  Retrieves the given warning line.
	 *  @param n                The zero-based index of the line to retrieve.
	 *  @return                 A null terminated string containing the given warning line message.
	 *  @see                    GetWarningStringLineCount, OutputWarningString
	 */
	const char*              GetWarningStringLine(int n);

	/**
	 *  Retrieves the number of lines in the current warning string buffer.
	 *  @return                 The number of lines.
	 *  @see                    GetWarningStringLine, GetWarningString, OutputWarningString
	 */
	int                      GetWarningStringLineCount(void)const;

	/**
	 *  Retrieves the current list of components.
	 *  @return                 The current list of components.
	 *  @see                    GetComponent, GetComponentCount
	 */
	std::list< std::string > ListComponents(void);

	/**
	 *  Load the specified database file into phreeqc.
	 *  @param filename         The name of the phreeqc database to load.
	 *                          The full path (or relative path with respect to the working directory) will be required if the file is not
	 *                          in the current working directory.
	 *  @return                 The number of errors encountered.
	 *  @see                    LoadDatabaseString
	 *  @remarks
	 *      All previous definitions are cleared.
	 */
	int                      LoadDatabase(const char* filename);

	/**
	 *  Load the specified string as a database into phreeqc.
	 *  @param input            String containing data to be used as the phreeqc database.
	 *  @return                 The number of errors encountered.
	 *  @see                    LoadDatabaseString
	 *  @remarks
	 *      All previous definitions are cleared.
	 */
	int                      LoadDatabaseString(const char* input);

	/**
	 *  Output the accumulated input buffer to stdout.  The input buffer can be run with a call to @ref RunAccumulated.
	 *  @see                    AccumulateLine, ClearAccumulatedLines, RunAccumulated
	 */
	void                     OutputAccumulatedLines(void);

	/**
	 *  Output the error messages normally stored in the <B><I>phreeqc.id.err</I></B> (where id is obtained from @ref GetId)
	 *  file to stdout.
	 *  @see                    GetErrorStringLine, GetErrorStringLineCount, GetErrorFileOn, SetErrorFileOn
	 */
	void                     OutputErrorString(void);

	/**
	 *  Output the warning messages to stdout.
	 *  @see                    GetWarningStringLine, GetWarningStringLineCount, GetWarningString
	 */
	void                     OutputWarningString(void);

	/**
	 *  Runs the input buffer as defined by calls to @ref AccumulateLine.
	 *  @return                 The number of errors encountered.
	 *  @see                    AccumulateLine, ClearAccumulatedLines, OutputAccumulatedLines, RunFile, RunString
	 *  @remarks
	 *      The accumulated input is cleared at the next call to @ref AccumulateLine.
	 *  @pre
	 *      @ref LoadDatabase/@ref LoadDatabaseString must have been called and returned 0 (zero) errors.
	 */
	int                      RunAccumulated(void);

	/**
	 *  Runs the specified phreeqc input file.
	 *  @param filename         The name of the phreeqc input file to run.
	 *  @return                 The number of errors encountered during the run.
	 *  @see                    RunAccumulated, RunString
	 *  @pre
	 *      @ref LoadDatabase/@ref LoadDatabaseString must have been called and returned 0 (zero) errors.
	 */
	int                      RunFile(const char* filename);

	/**
	 *  Runs the specified string as input to phreeqc.
	 *  @param input            String containing phreeqc input.
	 *  @return                 The number of errors encountered during the run.
	 *  @see                    RunAccumulated, RunFile
	 *  @pre
	 *      @ref LoadDatabase/@ref LoadDatabaseString must have been called and returned 0 (zero) errors.
	 */
	int                      RunString(const char* input);

	/**
	 *  Sets a C callback function for Basic programs. The syntax for the Basic command is
	 *  10 result = CALLBACK(x1, x2, string$)
	 *  The syntax for the C function is
	 *  double my_callback(double x1, double x2, const char * string)
	 *  @param fcn              The name of a user-defined function.
	 *  @param cookie1          A user defined value to be passed to the callback function.
	 *  @see                    SetBasicFortranCallback
	 */
	void                     SetBasicCallback(double (*fcn)(double x1, double x2, const char *str, void *cookie), void * cookie1);

	/**
	 *  Sets a Fortran callback function for Basic programs. The syntax for the Basic command is
	 *  10 result = CALLBACK(x1, x2, string$)
	 *  The syntax for the Fortran function is
	 *  real(kind=8) my_callback(x1, x2, string), where x1 and x2 are real(kind=8) and string is a character variable.
	 *  @param fcn              The name of a user-defined function.
	 *  @see                    SetBasicCallback
	 */
#ifdef IPHREEQC_NO_FORTRAN_MODULE
	void                     SetBasicFortranCallback(double (*fcn)(double *x1, double *x2, const char *str, size_t l));
#else
	void                     SetBasicFortranCallback(double (*fcn)(double *x1, double *x2, const char *str, int l));
#endif

	/**
	 *  Sets the current <B>SELECTED_OUTPUT</B> user number for use in subsequent calls to (@ref GetSelectedOutputColumnCount, 
     *  @ref GetSelectedOutputFileName, @ref GetSelectedOutputRowCount, @ref GetSelectedOutputString, @ref GetSelectedOutputStringLine, 
     *  @ref GetSelectedOutputStringLineCount, @ref GetSelectedOutputValue, @ref GetSelectedOutputValue2) routines.
	 *  The initial setting is 1.
	 *  @param n                The user number as specified in the <B>SELECTED_OUTPUT</B> block.
	 *  @retval VR_OK           Success
	 *  @retval VR_INVALIDARG   The given user number has not been defined.
	 *  @see                    GetCurrentSelectedOutputUserNumber, GetSelectedOutputColumnCount, GetSelectedOutputFileName, GetSelectedOutputRowCount, GetSelectedOutputString, GetSelectedOutputStringLine, GetSelectedOutputStringLineCount, GetSelectedOutputValue
	 */
	VRESULT                  SetCurrentSelectedOutputUserNumber(int n);

	/**
	 *  Sets the name of the dump file.  This file name is used if not specified within <B>DUMP</B> input.
	 *  The default value is <B><I>dump.id.out</I></B>, where id is obtained from @ref GetId.
	 *  @param filename         The name of the file to write <B>DUMP</B> output to.
	 *  @see                    GetDumpFileName, GetDumpFileOn, GetDumpString, GetDumpStringOn, GetDumpStringLine, GetDumpStringLineCount, SetDumpStringOn
	 */
	void                     SetDumpFileName(const char *filename);

	/**
	 *  Sets the dump file switch on or off.  This switch controls whether or not phreeqc writes to the <B>DUMP</B> (<B><I>dump.id.out</I></B>
	 *  if unspecified, where id is obtained from @ref GetId) file.
	 *  The initial setting is false.
	 *  @param bValue           If true, turns on output to the <B>DUMP</B> file;
	 *                          if false, turns off output to the <B>DUMP</B> file.
	 *  @see                    GetDumpFileOn, GetDumpString, GetDumpStringOn, GetDumpStringLine, GetDumpStringLineCount, SetDumpStringOn
	 */
	void                     SetDumpFileOn(bool bValue);

	/**
	 *  Sets the dump string switch on or off.  This switch controls whether or not the data normally sent
	 *  to the dump file are stored in a buffer for retrieval.  The initial setting is false.
	 *  @param bValue           If true, captures the output defined by the <B>DUMP</B> keyword into a string buffer;
	 *                          if false, output defined by the <B>DUMP</B> keyword is not captured to a string buffer.
	 *  @see                    GetDumpFileOn, GetDumpString, GetDumpStringOn, GetDumpStringLine, GetDumpStringLineCount, SetDumpFileOn
	 */
	void                     SetDumpStringOn(bool bValue);

	/**
	 *  Sets the name of the error file. The default value is <B><I>phreeqc.id.err</I></B>, where id is obtained from @ref GetId.
	 *  @param filename         The name of the file to write error output to.
	 *  @see                    GetErrorFileName, GetErrorFileOn, GetErrorString, GetErrorStringOn, GetErrorStringLine, GetErrorStringLineCount, SetErrorFileOn, SetErrorStringOn
	 */
	void                     SetErrorFileName(const char *filename);

	/**
	 *  Sets the error file switch on or off.  This switch controls whether or not
	 *  error messages are written to the <B><I>phreeqc.id.err</I></B> (where id is obtained from @ref GetId) file.
	 *  The initial setting is true.
	 *  @param bValue           If true, writes errors to the error file; if false, no errors are written to the error file.
	 *  @see                    GetErrorStringLine, GetErrorStringLineCount, GetErrorFileOn, OutputErrorString
	 */
	void                     SetErrorFileOn(bool bValue);

	/**
	 *  Sets the error switch on or off.  This switch controls whether
	 *  error messages are are generated and displayed.
	 *  The initial setting is true.
	 *  @param bValue           If true, error messages are sent to the error file and error string buffer; if false, no error messages are generated.
	 *  @see                    GetErrorOn, GetErrorStringLine, GetErrorStringLineCount, GetErrorFileOn, OutputErrorString
	 */
	void                     SetErrorOn(bool bValue);

	/**
	 *  Sets the error string switch on or off.  This switch controls whether or not the data normally sent
	 *  to the error file are stored in a buffer for retrieval.  The initial setting is true.
	 *  @param bValue           If true, captures error output into a string buffer; if false, error output is not captured to a string buffer.
	 *  @see                    GetErrorFileOn, GetErrorString, GetErrorStringOn, GetErrorStringLine, GetErrorStringLineCount, SetErrorFileOn
	 */
	void                     SetErrorStringOn(bool bValue);

	/**
	 *  Sets the name of the log file. The default value is <B><I>phreeqc.id.log</I></B>, where id is obtained from @ref GetId.
	 *  @param filename         The name of the file to write log output to.
	 *  @see                    GetLogFileName, GetLogFileOn, GetLogString, GetLogStringOn, GetLogStringLine, GetLogStringLineCount, SetLogFileOn, SetLogStringOn
	 */
	void                     SetLogFileName(const char *filename);

	/**
	 *  Sets the log file switch on or off.  This switch controls whether or not phreeqc
	 *  writes log messages to the <B><I>phreeqc.id.log</I></B> (where id is obtained from @ref GetId) file.  The initial setting is false.
	 *  @param bValue           If true, turns on output to the log file; if false, no log messages are written to the log file.
	 *  @remarks
	 *      Logging must be enabled through the use of the KNOBS -logfile option in order to receive any log messages.
	 *  @see                    GetLogFileOn
	 */
	void                     SetLogFileOn(bool bValue);

	/**
	 *  Sets the log string switch on or off.  This switch controls whether or not the data normally sent
	 *  to the log file are stored in a buffer for retrieval.  The initial setting is false.
	 *  @param bValue           If true, captures log output into a string buffer; if false, log output is not captured to a string buffer.
	 *  @see                    GetLogFileOn, GetLogString, GetLogStringOn, GetLogStringLine, GetLogStringLineCount, SetLogFileOn
	 */
	void                     SetLogStringOn(bool bValue);

	/**
	 *  Sets the name of the output file. The default value is <B><I>phreeqc.id.out</I></B>, where id is obtained from @ref GetId.
	 *  @param filename         The name of the file to write phreeqc output to.
	 *  @see                    GetOutputFileName, GetOutputFileOn, GetOutputString, GetOutputStringOn, GetOutputStringLine, GetOutputStringLineCount, SetOutputFileOn, SetOutputStringOn
	 */
	void                     SetOutputFileName(const char *filename);

	/**
	 *  Sets the output file switch on or off.  This switch controls whether or not phreeqc
	 *  writes to the <B><I>phreeqc.id.out</I></B> file (where id is obtained from @ref GetId).  This is the output that is normally generated
	 *  when phreeqc is run.  The initial setting is false.
	 *  @param bValue           If true, writes output to the output file; if false, no output is written to the output file.
	 *  @see                    GetOutputFileOn
	 */
	void                     SetOutputFileOn(bool bValue);

	/**
	 *  Sets the output string switch on or off.  This switch controls whether or not the data normally sent
	 *  to the output file are stored in a buffer for retrieval.  The initial setting is false.
	 *  @param bValue           If true, captures output into a string buffer; if false, output is not captured to a string buffer.
	 *  @see                    GetOutputFileOn, GetOutputString, GetOutputStringOn, GetOutputStringLine, GetOutputStringLineCount, SetOutputFileOn
	 */
	void                     SetOutputStringOn(bool bValue);

	/**
	 *  Sets the name of the current selected output file (see @ref SetCurrentSelectedOutputUserNumber).  This file name is used if not specified within <B>SELECTED_OUTPUT</B> input.
	 *  The default value is <B><I>selected_n.id.out</I></B>, where id is obtained from @ref GetId.
	 *  @param filename         The name of the file to write <B>SELECTED_OUTPUT</B> output to.
	 *  @see                    GetSelectedOutputFileName, GetSelectedOutputFileOn, GetSelectedOutputString, GetSelectedOutputStringOn, GetSelectedOutputStringLine, GetSelectedOutputStringLineCount, SetCurrentSelectedOutputUserNumber, SetSelectedOutputStringOn
	 */
	void                     SetSelectedOutputFileName(const char *filename);

	/**
	 *  Sets the selected-output file switch on or off.  This switch controls whether or not phreeqc writes output to
	 *  the current <B>SELECTED_OUTPUT</B> (<B><I>selected_n.id.out</I></B> if unspecified, where id is obtained from @ref GetId) file.
	 *  The initial setting is false.
	 *  @param bValue           If true, writes output to the selected-output file; if false, no output is written to the selected-output file.
	 *  @see                    GetSelectedOutputColumnCount, GetSelectedOutputFileOn, GetSelectedOutputRowCount, GetSelectedOutputValue, SetCurrentSelectedOutputUserNumber
	 */
	void                     SetSelectedOutputFileOn(bool bValue);

	/**
	 *  Sets the selected output string switch on or off.  This switch controls whether or not the data normally sent
	 *  to the current <B>SELECTED_OUTPUT</B> file (see @ref SetCurrentSelectedOutputUserNumber) are stored in a buffer for retrieval.
	 *  The initial setting is false.
	 *  @param bValue           If true, captures the output defined by the <B>SELECTED_OUTPUT</B> keyword into a string buffer;
	 *                          if false, output defined by the <B>SELECTED_OUTPUT</B> keyword is not captured to a string buffer.
	 *  @see                    GetSelectedOutputFileOn, GetSelectedOutputString, GetSelectedOutputStringOn, GetSelectedOutputStringLine, GetSelectedOutputStringLineCount, SetCurrentSelectedOutputUserNumber, SetSelectedOutputFileOn
	 */
	void                     SetSelectedOutputStringOn(bool bValue);

public:
	// overrides
	virtual void error_msg(const char *str, bool stop=false);
	virtual void log_msg(const char * str);
	virtual void output_msg(const char *str);
	virtual void punch_msg(const char *str);
	virtual void screen_msg(const char *str);
	virtual void warning_msg(const char *str);

	virtual void fpunchf(const char *name, const char *format, double d);
	virtual void fpunchf(const char *name, const char *format, char * d);
	virtual void fpunchf(const char *name, const char *format, int d);
	virtual void fpunchf_end_row(const char *format);

	virtual bool output_open(const char *file_name, std::ios_base::openmode mode = std::ios_base::out);
	virtual bool punch_open(const char *file_name, std::ios_base::openmode mode = std::ios_base::out, int n_user = 1);

protected:
	int EndRow(void);
	void AddSelectedOutput(const char* name, const char* format, va_list argptr);
	void UnLoadDatabase(void);

	void check_database(const char* sz_routine);
	int close_input_files(void);
	int close_output_files(void);
	void open_output_files(const char* sz_routine);

	void do_run(const char* sz_routine, std::istream* pis, PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie);

	void update_errors(void);

	int load_db(const char* filename);
	int load_db_str(const char* filename);
	int test_db(void);

	bool get_sel_out_file_on(int n)const;
	std::string sel_file_name(int n_user);

	std::string create_file_name(const char *prefix, const char *suffix);

	bool get_sel_out_string_on(int n)const;

protected:
#if defined(_MSC_VER)
/* disable warning C4251: 'identifier' : class 'type' needs to have dll-interface to be used by clients of class 'type2' */
#pragma warning(disable:4251)
#endif

	bool                       DatabaseLoaded;
	bool                       ClearAccumulated;
	bool                       UpdateComponents;
	std::map< int, bool >      SelectedOutputFileOnMap;

	bool                       OutputFileOn;

	bool                       LogFileOn;
	bool                       ErrorFileOn;
	bool                       DumpOn;

	bool                       DumpStringOn;

	bool                       OutputStringOn;
	std::string                OutputString;
	std::vector< std::string > OutputLines;

	bool                       LogStringOn;
	std::string                LogString;
	std::vector< std::string > LogLines;

	bool                       ErrorStringOn;
	IErrorReporter            *ErrorReporter;
	std::string                ErrorString;
	std::vector< std::string > ErrorLines;

	bool                       WarningStringOn;
	IErrorReporter            *WarningReporter;
	std::string                WarningString;
	std::vector< std::string > WarningLines;

	int                                           CurrentSelectedOutputUserNumber;
	std::map< int, CSelectedOutput* >             SelectedOutputMap;
	std::string                                   StringInput;

	std::string                DumpString;
	std::vector< std::string > DumpLines;

	std::list< std::string >   Components;
	std::list< std::string >   EquilibriumPhasesList;
	const std::list<std::string> &GetEquilibriumPhasesList() { return this->EquilibriumPhasesList; };
	std::list< std::string >   GasComponentsList;
	const std::list<std::string> &GetGasComponentsList() { return this->GasComponentsList; };
	std::list< std::string >   KineticReactionsList;
	const std::list<std::string> &GetKineticReactionsList() { return this->KineticReactionsList; };
	std::list< std::string >   SolidSolutionComponentsList;
	const std::list<std::string> &GetSolidSolutionComponentsList() { return this->SolidSolutionComponentsList; };
	std::list< std::string >   SolidSolutionNamesList;
	const std::list<std::string> &GetSolidSolutionNamesList() { return this->SolidSolutionNamesList; };
	//std::list< std::string >   SurfaceSpeciesList;
	//const std::list<std::string> &GetSurfaceSpeciesList() { return this->SurfaceSpeciesList; };
	std::list< std::string >   SurfaceTypeList;
	const std::list<std::string> &GetSurfaceTypeList() { return this->SurfaceTypeList; };
	std::list< std::string >   SurfaceNamesList;
	const std::list<std::string> &GetSurfaceNamesList() { return this->SurfaceNamesList; };
	//std::list< std::string >   ExchangeSpeciesList;
	//const std::list<std::string> &GetExchangeSpeciesList() { return this->ExchangeSpeciesList; };
	std::list< std::string >   ExchangeNamesList;
	const std::list<std::string> &GetExchangeNamesList() { return this->ExchangeNamesList; };

	std::map< int, std::string > SelectedOutputFileNameMap;

	std::string                OutputFileName;
	std::string                ErrorFileName;
	std::string                LogFileName;
	std::string                DumpFileName;

	std::map< int, bool >                         SelectedOutputStringOn;
	std::map< int, std::string >                  SelectedOutputStringMap;
	std::map< int, std::vector< std::string > >   SelectedOutputLinesMap;

protected:
	Phreeqc* PhreeqcPtr;
	FILE *input_file;
	FILE *database_file;

	friend class IPhreeqcLib;
	static std::map<size_t, IPhreeqc*> Instances;
	static size_t InstancesIndex;
	size_t Index;

	static std::string Version;

#if defined(_MSC_VER)
/* reset warning C4251 */
#pragma warning(default:4251)
#endif

#if defined(CPPUNIT)
	friend class TestIPhreeqc;
	friend class TestSelectedOutput;
#endif

private:
	/**
	 *  Copy constructor not supported
	 */
	IPhreeqc(const IPhreeqc&);

	/**
	 *  operator= not supported
	 */
	IPhreeqc& operator=(const IPhreeqc&);
};

#endif // INC_IPHREEQC_HPP
