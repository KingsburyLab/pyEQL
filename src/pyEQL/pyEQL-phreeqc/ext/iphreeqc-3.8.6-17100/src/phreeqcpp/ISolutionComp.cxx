#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <cassert>

#include "Utils.h"
#include "ISolutionComp.h"
#include "Parser.h"
#include "Solution.h"
#include "phqalloc.h"

#if defined(PHREEQCI_GUI)
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#endif

cxxISolutionComp::cxxISolutionComp(PHRQ_io *io):
PHRQ_base(io),
moles(0.0), 
input_conc(0.0), 
phase_si(0.0),
gfw(0.0)
{
}
cxxISolutionComp::~cxxISolutionComp(void)
{
}

/* ---------------------------------------------------------------------- */
CParser::STATUS_TYPE cxxISolutionComp::
read(const char *line_in, cxxSolution *solution_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Remove space between "kg" and "solution" or "water" in units
 */
	std::string line = line_in;
	Utilities::replace("Kg", "kg", line);
	Utilities::replace("KG", "kg", line);
	while (Utilities::replace("kg ", "kg", line) == TRUE);
/*
 *   Read master species list for mass balance equation
 */
	std::string master_list;
	std::string token;
	std::string::iterator b = line.begin(); 
	std::string::iterator e = line.end(); 
	{
		int j;
		while (((j = CParser::copy_token(token, b, e)) == CParser::TT_UPPER) ||
			(token[0] == '[') ||
			(Utilities::strcmp_nocase(token.c_str(), "ph") == 0) ||
			(Utilities::strcmp_nocase(token.c_str(), "pe") == 0))
		{
			Utilities::replace("(+", "(", token);
			if (master_list.size() > 0)
				master_list.append(" ");
			master_list.append(token);
		}
	}
	if (master_list.size() == 0)
	{
		error_msg
			("No element or master species given for concentration input.",
			  PHRQ_io::OT_CONTINUE);
		return (CParser::PARSER_ERROR);
	}
	this->Set_description(master_list.c_str());
/*
 *   Determine if reading alkalinity, allow equivalents for units
 */
	bool alk;
	Utilities::str_tolower(master_list);
	if (strstr(master_list.c_str(), "alk") == master_list.c_str())
	{
		alk = true;
	}
	else
	{
		alk = false;
	}
/*
 *   Read concentration
 */
	{
		LDBLE dummy;
		int j = sscanf(token.c_str(), SCANFORMAT, &dummy);
		if (j == 0)
		{
			std::ostringstream errstr;
			errstr << "Concentration data error for " << master_list << " in solution input.";
			error_msg(errstr.str().c_str(),  PHRQ_io::OT_CONTINUE);
			return (CParser::PARSER_ERROR);
		}
		else
		{
			this->Set_input_conc(dummy);
		}
		if ((j = CParser::copy_token(token, b, e)) == CParser::TT_EMPTY)
			return (CParser::PARSER_OK);
	}
/*
 *   Read optional data
 */
	std::string token1 = token;

/*
 *   Check for units info
 */
	CParser parser(this->io);
	if (solution_ptr->Get_initial_data() == NULL)
	{
		error_msg("Initial_data instance not defined in cxxISolutionComp::read", 1);
	}
	if (parser.check_units(token1, alk, false, solution_ptr->Get_initial_data()->Get_units().c_str(), false) == CParser::PARSER_OK)
	{
		if (parser.check_units(token1, alk, false, solution_ptr->Get_initial_data()->Get_units().c_str(), true) == CParser::PARSER_OK)
		{
			this->units = token1;
			if ((CParser::copy_token(token, b, e)) == CParser::TT_EMPTY)
				return (CParser::PARSER_OK);
		}
		else
		{
			return (CParser::PARSER_ERROR);
		}
	}
/*
 *   Check for "as" followed by formula to be used for gfw
 */
	token1 = token;
	Utilities::str_tolower(token1);
	if (strcmp(token1.c_str(), "as") == 0)
	{
		CParser::copy_token(token, b, e);
		this->as = token;
		if ((CParser::copy_token(token, b, e)) == CParser::TT_EMPTY)
			return (CParser::PARSER_OK);
/*
 *   Check for "gfw" followed by gram formula weight
 */
	}
	else if (strcmp(token1.c_str(), "gfw") == 0 || strcmp(token1.c_str(), "gfm") == 0)
	{
		if (CParser::copy_token(token, b, e) != DIGIT)
		{
			error_msg("Expecting gram formula weight.",  PHRQ_io::OT_CONTINUE);
			return (CParser::PARSER_ERROR);
		}
		else
		{
			(void)sscanf(token.c_str(), SCANFORMAT, &this->gfw);
			if ((CParser::copy_token(token, b, e)) == CParser::TT_EMPTY)
				return (CParser::PARSER_OK);
		}
	}
/*
 *   Check for redox couple for pe
 */
	if (Utilities::strcmp_nocase(token.c_str(), "pe") == 0)
	{
		this->pe_reaction = token;
		if ((CParser::copy_token(token, b, e)) == CParser::TT_EMPTY)
			return (CParser::PARSER_OK);
	}
	else if (strstr(token.c_str(), "/") != NULL)
	{
		if (parser.parse_couple(token) == CParser::PARSER_OK)
		{
			this->pe_reaction = token;
			if ((CParser::copy_token(token, b, e)) == CParser::TT_EMPTY)
				return (CParser::PARSER_OK);
		}
		else
		{
			return (CParser::PARSER_ERROR);
		}
	}
/*
 *   Must have phase
 */
	this->equation_name = token;
	if (CParser::copy_token(token, b, e) == CParser::TT_EMPTY)
		return (CParser::PARSER_OK);
/*
 *   Check for saturation index
 */
	{
		int j = sscanf(token.c_str(), SCANFORMAT,
			&(this->phase_si));
		if (j != 1)
		{
			error_msg("Expected saturation index.",  PHRQ_io::OT_CONTINUE);
			return (CParser::PARSER_ERROR);
		}
	}
	return (CParser::PARSER_OK);

}
