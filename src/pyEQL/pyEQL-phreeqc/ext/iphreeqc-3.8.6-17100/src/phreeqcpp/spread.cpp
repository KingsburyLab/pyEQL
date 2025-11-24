#ifndef boolean
typedef unsigned char boolean;
#endif
#include "Phreeqc.h"
#include "phqalloc.h"
#include "Solution.h"
#include "Utils.h"

#define STRING 11
#define NUMBER 12
#define MIXED 13
#define EOL 14

#define OPTION_EOF -1
#define OPTION_KEYWORD -2
#define OPTION_ERROR -3
#define OPTION_DEFAULT -4
#define OPT_1 -5

#if defined(PHREEQCI_GUI)
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#endif

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_solution_spread(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads solution data
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
	class spread_row *heading, *row_ptr, *units;
	int count, strings, numbers;
	int spread_lines;
	const char* cptr;
	class defaults soln_defaults;
	int return_value, opt;
	const char* next_char;
	const char *opt_list[] = {
		"temp",					/* 0 */
		"temperature",			/* 1 */
		"dens",					/* 2 */
		"density",				/* 3 */
		"units",				/* 4 */
		"redox",				/* 5 */
		"ph",					/* 6 */
		"pe",					/* 7 */
		"unit",					/* 8 */
		"isotope",				/* 9 */
		"water",				/* 10 */
		"isotope_uncertainty",	/* 11 */
		"uncertainty",			/* 12 */
		"uncertainties",		/* 13 */
		"pressure",				/* 14 */
		"press"		  		    /* 15 */
	};
	int count_opt_list = 16;
/*
 * Initialize defaults
 */
	soln_defaults.temp = 25;
	soln_defaults.density = 1.0;
	soln_defaults.calc_density = false;
	soln_defaults.units = string_hsave("mmol/kgw");
	soln_defaults.redox = string_hsave("pe");
	soln_defaults.ph = 7.0;
	soln_defaults.pe = 4.0;
	soln_defaults.water = 1.0;
	soln_defaults.pressure = 1.0;

#ifdef PHREEQCI_GUI
	free_spread();
#endif


	/* fill in soln_defaults.iso */
	soln_defaults.iso.resize(count_iso_defaults);

	/* all iso[i].name is hsave'd, so no conflicts */
	// memcpy(&soln_defaults.iso[0], iso_defaults,
	// 	   soln_defaults.iso.size() * sizeof(class iso));
	for (size_t i = 0; i < count_iso_defaults; ++i) {
		soln_defaults.iso[i].name        = iso_defaults[i].name;
		soln_defaults.iso[i].value       = iso_defaults[i].value;
		soln_defaults.iso[i].uncertainty = iso_defaults[i].uncertainty;
	}

	heading = NULL;
	units = NULL;
	return_value = UNKNOWN;
	spread_lines = 0;
	CParser parser(this->phrq_io);
/*
 *   Loop on solutions
 */
	for (;;)
	{
		std::string token, token1;
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (spread_lines == 0 && opt != OPTION_DEFAULT)
		{
			row_ptr = string_to_spread_row(line);
			cptr = line;
			count = numbers = strings = 0;
			int j;
			while (((j = copy_token(token, &cptr)) != EMPTY))
			{
				count++;
				if (j == UPPER || j == LOWER)
					strings++;
				if (j == DIGIT)
					numbers++;
			}
			/*
			 * Is 2nd token all number
			 */
			cptr = line;
			copy_token(token, &cptr);
			j = copy_token(token, &cptr);
			bool num = false;
			if (j == DIGIT)
			{
				char* ptr;
				(void)strtod(token.c_str(), &ptr);
				cptr = ptr;
				int j1 = copy_token(token1, &cptr);
				if (j1 != EMPTY)
				{
					num = FALSE;
				}
				else
				{
					num = TRUE;
				}
			}
			/*
			 *   Starts with hyphen
			 */
			cptr = line;
			copy_token(token, &cptr);
			if (token[0] == '-')
			{
				/* opt = opt; */
			}
			else
			{
				switch (opt)
				{
				case 0:		/* temp */
				case 1:		/* temperature */
				case 10:		/* water */
					if ((count == 2 || count == 3) && num == TRUE)
					{
						/* opt = opt; */
					}
					else
					{
						opt = OPTION_DEFAULT;
					}
					break;
				case 2:		/* dens */
				case 3:		/* density */
					copy_token(token, &cptr);
					if (count == 2 || count == 3 && (num == TRUE || token[0] == 'c' || token[0] == 'C'))
					{
						/* opt = opt; */
					}
					else
					{
						opt = OPTION_DEFAULT;
					}
					break;
				case 6:		/* ph */
				case 7:		/* pe */
					if ((count == 2 || count == 3 || count == 4)
						&& num == TRUE)
					{
						/* opt = opt; */
					}
					else
					{
						opt = OPTION_DEFAULT;
					}
					break;
				case 5:		/* redox */
				case 4:		/* units */
				case 8:		/* unit */
					if (count == 2)
					{
						/* opt = opt; */
					}
					else
					{
						opt = OPTION_DEFAULT;
					}
					break;
				case 9:		/* isotope */
					if (row_ptr->count > 4)
					{
						opt = OPTION_DEFAULT;
					}
					else
					{
						/* opt = opt; */
					}
					break;
				case 11:		/* isotope_uncertainty */
				case 12:		/* uncertainty */
				case 13:		/* uncertainties */
					if (row_ptr->count > 3)
					{
						opt = OPTION_DEFAULT;
					}
					else
					{
						/* opt = opt; */
					}
					break;
				case 14: /* pressure */
				case 15: /* press */
					(void)sscanf(next_char, SCANFORMAT, &(soln_defaults.pressure));
					break;
				}
			}
			spread_row_free(row_ptr);
		}
		if (opt == OPTION_DEFAULT)
		{
			if (spread_lines == 0)
			{
				opt = 100;
			}
			spread_lines++;
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in SOLUTION keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case OPTION_DEFAULT:	/* solution definition */
			row_ptr = string_to_spread_row(line);
			if (spread_lines == 2)
			{
				numbers = 0;
				strings = 0;
				//for (int i = 0; i < heading->count; i++)
				for (int i = 0; i < row_ptr->count; i++)
				{
					if (row_ptr->type_vector[i] == STRING)
					{
						strings++;
					}
					else if (row_ptr->type_vector[i] == NUMBER)
					{
						numbers++;
					}
				}
				if (numbers == 0)
				{
					units = row_ptr;
					break;
				}
			}
			spread_row_to_solution(heading, units, row_ptr, soln_defaults);
#ifdef PHREEQCI_GUI
			add_row(row_ptr);
#endif
			spread_row_free(row_ptr);
			break;
		case 0:				/* temperature */
		case 1:
			(void)sscanf(next_char, SCANFORMAT, &(soln_defaults.temp));
			break;
		case 2:				/* density */
		case 3:
			{
				int j = copy_token(token, &next_char);
				if (j == DIGIT)
				{
					if (sscanf(token.c_str(), SCANFORMAT, &dummy) != 1)
					{
						error_msg("Expecting numeric value for density.", PHRQ_io::OT_CONTINUE);
						error_msg(line_save, PHRQ_io::OT_CONTINUE);
						input_error++;
					}
					else
					{
						soln_defaults.density = dummy;
						copy_token(token, &next_char);
						if (token[0] == 'c' || token[0] == 'C')
							soln_defaults.calc_density = true;
					}
				}
				else if (j != EMPTY)
				{
					if (token[0] != 'c' && token[0] != 'C')
					{
						error_msg("Options following density are numeric value or c[alculate].", PHRQ_io::OT_CONTINUE);
						error_msg(line_save, PHRQ_io::OT_CONTINUE);
						input_error++;
					}
					else
						soln_defaults.calc_density = true;
				}
			}
			break;
		case 4:				/* units */
		case 8:				/* unit */
			if (copy_token(token, &next_char) == EMPTY)
				break;
			{
				if (check_units(token, FALSE, FALSE, NULL, TRUE) == OK)
				{
					soln_defaults.units = string_hsave(token.c_str());
				}
				else
				{
					input_error++;
				}
			}
			break;
		case 5:				/* redox */
			if (copy_token(token, &next_char) == EMPTY)
				break;
			if (parser.parse_couple(token) == OK)
			{
				soln_defaults.redox = string_hsave(token.c_str());
			}
			else
			{
				input_error++;
			}
			break;
		case 6:				/* ph */
			copy_token(token, &next_char);
			(void)sscanf(token.c_str(), SCANFORMAT, &(soln_defaults.ph));
			if (copy_token(token, &next_char) != EMPTY)
			{
				warning_msg
					("Not possible to use phase name or saturation index in definition of default pH in SOLUTION_SPREAD.");
			}
			break;
		case 7:				/* pe */
			copy_token(token, &next_char);
			(void)sscanf(token.c_str(), SCANFORMAT, &(soln_defaults.pe));
			if (copy_token(token, &next_char) != EMPTY)
			{
				warning_msg
					("Not possible to use phase name or saturation index in definition of default pe in SOLUTION_SPREAD.");
			}
			break;
		case 11:				/* isotope_uncertainty */
		case 12:				/* uncertainty */
		case 13:				/* uncertainties */
			{
				if (copy_token(token, &next_char) != DIGIT)
				{
					input_error++;
					error_string = sformatf( "Expected isotope name to"
						" begin with an isotopic number.");
					error_msg(error_string, CONTINUE);
					error_string = sformatf( "In read_solution_spread isotope_uncertainty\n");
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					error_string = sformatf( "\t%s\t%s\n", "token:     ", token.c_str());
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					error_string = sformatf( "\t%s\t%s\n", "next_char: ", next_char);
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					error_string = sformatf( "\t%s\t%s\n", "line_save: ", line_save);
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					continue;
				}
				size_t i;
				for (i = 0; i < soln_defaults.iso.size(); i++)
				{
					if (strcmp(token.c_str(), soln_defaults.iso[i].name) == 0)
					{
						break;
					}
				}
				if (i == soln_defaults.iso.size())
				{
					soln_defaults.iso.resize((size_t)i + 1);
					soln_defaults.iso[i].name = string_hsave(token.c_str());
					soln_defaults.iso[i].value = NAN;
					soln_defaults.iso[i].uncertainty = NAN;
				}

				/* read and store isotope ratio uncertainty */
				int j;
				if ((j = copy_token(token, &next_char)) != EMPTY)
				{
					if (j != DIGIT)
					{
						input_error++;
						error_string = sformatf(
							"Expected numeric value for uncertainty in isotope ratio.");
						error_msg(error_string, CONTINUE);
						continue;
					}
					else
					{
						(void)sscanf(token.c_str(), SCANFORMAT,
							&(soln_defaults.iso[i].uncertainty));
					}
				}
				else
				{
					soln_defaults.iso[i].uncertainty = NAN;
				}
			}
			break;
		case 10:				/* water */
			{
				int j = copy_token(token, &next_char);
				if (j != DIGIT)
				{
					input_error++;
					error_string = sformatf(
						"Expected numeric value for mass of water in solution.");
					error_msg(error_string, CONTINUE);
				}
				else
				{
					(void)sscanf(token.c_str(), SCANFORMAT, &(soln_defaults.water));
				}
			}
			break;
		case 9:				/* isotope */
			{
				if (copy_token(token, &next_char) != DIGIT)
				{
					input_error++;
					error_string = sformatf( "Expected isotope name to"
						" begin with an isotopic number.");
					error_msg(error_string, CONTINUE);
					error_string = sformatf( "In read_solution_spread isotope\n");
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					error_string = sformatf( "\t%s\t%s\n", "token:     ", token.c_str());
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					error_string = sformatf( "\t%s\t%s\n", "next_char: ", next_char);
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					error_string = sformatf( "\t%s\t%s\n", "line_save: ", line_save);
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					continue;
				}
				size_t i;
				for (i = 0; i < soln_defaults.iso.size(); i++)
				{
					if (strcmp(token.c_str(), soln_defaults.iso[i].name) == 0)
					{
						break;
					}
				}
				if (i == soln_defaults.iso.size())
				{
					soln_defaults.iso.resize(i + 1);
					soln_defaults.iso[i].name = string_hsave(token.c_str());
					soln_defaults.iso[i].value = NAN;
					soln_defaults.iso[i].uncertainty = NAN;
				}
				/* read and store isotope ratio */
				if (copy_token(token, &next_char) != DIGIT)
				{
					input_error++;
					error_string = sformatf(
						"Expected numeric value for default isotope ratio.");
					error_msg(error_string, CONTINUE);
					break;
				}
				(void)sscanf(token.c_str(), SCANFORMAT, &(soln_defaults.iso[i].value));
				/* read and store isotope ratio uncertainty */
				int j;
				if ((j = copy_token(token, &next_char)) != EMPTY)
				{
					if (j != DIGIT)
					{
						input_error++;
						error_string = sformatf(
							"Expected numeric value for uncertainty in isotope ratio.");
						error_msg(error_string, CONTINUE);
						continue;
					}
					else
					{
						(void)sscanf(token.c_str(), SCANFORMAT,
							&(soln_defaults.iso[i].uncertainty));
					}
				}
			}
			break;
		case 14: /* pressure */
			(void)sscanf(next_char, SCANFORMAT, &(soln_defaults.pressure));
			break;
		case 100:				/* read headings */
			heading = string_to_spread_row(line);
			{
				int i;
				for (i = 0; i < heading->count; i++)
				{
					while (replace(" ", "", heading->str_vector[i]) == TRUE);
					while (replace(",", "_", heading->str_vector[i]) == TRUE);
				}
			}
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
#ifdef PHREEQCI_GUI
	if (heading)
	{
		assert(g_spread_sheet.heading == NULL);
		g_spread_sheet.heading = copy_row(heading);
	}
	if (units)
	{
		assert(g_spread_sheet.units == NULL);
		g_spread_sheet.units = copy_row(units);
	}
	g_spread_sheet.defaults = soln_defaults;
#endif
	spread_row_free(heading);
	spread_row_free(units);

	soln_defaults.iso.clear();
	return (return_value);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
spread_row_to_solution(class spread_row *heading, class spread_row *units,
					   class spread_row *data, class defaults defaults)
/* ---------------------------------------------------------------------- */
{
	Keywords::KEYWORDS next_keyword_save;
	int n_user, n_user_end;
	std::string description;
	std::string token, token1, string;
	CParser parser(this->phrq_io);

	int return_value, opt;
	const char* next_char;
	const char *opt_list[] = {
		"temp",					/* 0 */
		"temperature",			/* 1 */
		"dens",					/* 2 */
		"density",				/* 3 */
		"units",				/* 4 */
		"redox",				/* 5 */
		"ph",					/* 6 */
		"pe",					/* 7 */
		"unit",					/* 8 */
		"isotope",				/* 9 */
		"water",				/* 10 */
		"description",			/* 11 */
		"desc",					/* 12 */
		"descriptor",			/* 13 */
		"pressure",				/* 14 */
		"press",			    /* 15 */
		"potential"			    /* 16 */
	};
	int count_opt_list = 17;

/*
 *      look for solution number
 */
	n_user = n_user_end = -1;
	{
		int i; 
		for (i = 0; i < heading->count; i++)
		{
			if (strcmp_nocase(heading->str_vector[i].c_str(), "number") == 0)
			{
				break;
			}
		}
		if (i == heading->count || data->type_vector[i] == EMPTY
			|| data->count <= i)
		{
			n_user = -1;
		}
		else if (data->type_vector[i] == STRING)
		{
			input_error++;
			error_string = sformatf(
				"Expected solution number or number range in 'number' column, found:  %s.",
				data->str_vector[i].c_str());
			error_msg(error_string, CONTINUE);
		}
		else
		{
			string = "solution_s ";
			string.append(data->str_vector[i]);
			next_keyword_save = next_keyword;
			next_keyword = Keywords::KEY_SOLUTION_SPREAD;
			cxxNumKeyword nk;
			nk.read_number_description(string);
			n_user = nk.Get_n_user();
			n_user_end = nk.Get_n_user_end();
			description = nk.Get_description();
			next_keyword = next_keyword_save;
			Rxn_new_solution.insert(n_user);
		}
	}
/*
 *   set up solution
 */
	
	cxxSolution temp_solution;
	temp_solution.Set_n_user(n_user);
	temp_solution.Set_n_user_end(n_user_end);
	temp_solution.Set_new_def(true);
	temp_solution.Create_initial_data();
	cxxISolution * initial_data_ptr = temp_solution.Get_initial_data();
	if (use.Get_solution_in() == FALSE)
	{
		use.Set_solution_in(true);
		use.Set_n_solution_user(n_user);
	}
/*
 *   Set default ph, temp, density, pe, units
 */
	temp_solution.Set_description(description);
	temp_solution.Set_tc(defaults.temp);
	temp_solution.Set_patm(defaults.pressure);
	temp_solution.Set_potV(0);
	temp_solution.Set_ph(defaults.ph);
	temp_solution.Set_density(defaults.density);
	initial_data_ptr->Set_calc_density(defaults.calc_density);
	temp_solution.Set_pe(defaults.pe);
	temp_solution.Set_mass_water(defaults.water);
	temp_solution.Set_ah2o(1.0);
	temp_solution.Set_mu(1e-7);
	initial_data_ptr->Set_units(defaults.units);
	initial_data_ptr->Set_default_pe(defaults.redox);
	{
		CReaction temp_chem_reaction;
		initial_data_ptr->Get_pe_reactions()[defaults.redox] = temp_chem_reaction;
	}
/*
 *   Read concentration data
 */
	return_value = UNKNOWN;
	for (int i = 0; i < heading->count; i++)
	{
		if (strcmp_nocase(heading->str_vector[i].c_str(), "number") == 0)
			continue;
		if (strcmp_nocase(heading->str_vector[i].c_str(), "uncertainty") == 0)
			continue;
		if (strcmp_nocase(heading->str_vector[i].c_str(), "uncertainties") == 0)
			continue;
		if (strcmp_nocase(heading->str_vector[i].c_str(), "isotope_uncertainty") ==
			0)
			continue;
		/*
		 *  Copy in element name
		 */
		if (heading->type_vector[i] == EMPTY)
			continue;
		string = heading->str_vector[i];
		string.append(" ");
		/*
		 *  Copy in concentration data
		 */
		if (i >= data->count || data->type_vector[i] == EMPTY)
			continue;
		string.append(data->str_vector[i]);
		string.append(" ");
		/*
		 *  Copy in concentration data
		 */
		if (units != NULL && i < units->count
			&& units->type_vector[i] != EMPTY)
		{
			string.append(units->str_vector[i]);
		}
/*
 *   Parse string just like read_solution input 
 */
		char * char_string = string_duplicate(string.c_str());
		next_char = char_string;
		opt = get_option_string(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT && heading->type_vector[i] == NUMBER)
		{
			opt = 9;
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in SOLUTION keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* temperature */
		case 1:
			(void)sscanf(next_char, SCANFORMAT, &dummy);
			temp_solution.Set_tc(dummy);
			break;
		case 2:				/* density */
		case 3:
			{
				int j = copy_token(token, &next_char);
				(void)sscanf(token.c_str(), SCANFORMAT, &dummy);
				temp_solution.Set_density(dummy);
				j = copy_token(token, &next_char);
				if (j != EMPTY)
				{
					if (token[0] != 'c' && token[0] != 'C')
					{
						error_msg("Only option following density is c[alculate].", PHRQ_io::OT_CONTINUE);
						error_msg(line_save, PHRQ_io::OT_CONTINUE);
						input_error++;
					}
					else
					{
						initial_data_ptr->Set_calc_density(true);
					}
				}
			}
			break;
		case 4:				/* units */
		case 8:				/* unit */
			if (copy_token(token, &next_char) == EMPTY)
				break;
			{
				if (check_units(token, false, false, initial_data_ptr->Get_units().c_str(), true) == OK)
				{
					initial_data_ptr->Set_units(token);
				}
				else
				{
					input_error++;
				}
			}
			break;
		case 5:				/* redox */
			if (copy_token(token, &next_char) == EMPTY)
				break;
			if (parser.parse_couple(token) == OK)
			{
				const char * pe_str = string_hsave(token.c_str());
				initial_data_ptr->Set_default_pe(pe_str);
				CReaction temp_chem_reaction;
				initial_data_ptr->Get_pe_reactions()[token] = temp_chem_reaction;
			}
			else
			{
				input_error++;
			}
			break;
		case 6:				/* ph */
			{
				cxxISolutionComp temp_comp(this->phrq_io);
				if (temp_comp.read(char_string, &temp_solution) == CParser::PARSER_ERROR)
				{
					input_error++;
					break;
				}
				
				temp_solution.Set_ph(temp_comp.Get_input_conc());
				if (temp_comp.Get_equation_name().size() == 0)
				{
					break;
					
				}
				temp_comp.Set_description("H(1)");
				initial_data_ptr->Get_comps()[temp_comp.Get_description()] = temp_comp;
			}
			break;
		case 7:				/* pe */
			{
				cxxISolutionComp temp_comp(this->phrq_io);
				if (temp_comp.read(char_string, &temp_solution) == CParser::PARSER_ERROR)
				{
					input_error++;
					break;
				}
				temp_solution.Set_pe(temp_comp.Get_input_conc());
				if (temp_comp.Get_equation_name().size() == 0)
				{
					break;
				}
				temp_comp.Set_description("E");
				initial_data_ptr->Get_comps()[temp_comp.Get_description()] = temp_comp;
			}
			break;
		case 9:				/* isotope */
			{
				next_char = char_string;
				cxxSolutionIsotope temp_isotope;
				if (copy_token(token, &next_char) !=  CParser::TT_DIGIT)
				{
					input_error++;
					error_string = sformatf( "Expected isotope name to"
						" begin with an isotopic number.");
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					error_string = sformatf( "In spread_row_to_solution isotope\n");
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					error_string = sformatf( "\t%s\t%s\n", "token:       ", token.c_str());
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					error_string = sformatf( "\t%s\t%s\n", "next_char:   ", next_char);
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					error_string = sformatf( "\t%s\t%s\n", "char_string: ", char_string);
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					error_string = sformatf( "\t%s\t%s\n", "line_save:   ", line_save);
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
//class spread_row
//{
//	int count;
//	int empty, string, number;
//	char **char_vector;
//	LDBLE *d_vector;
//	int *type_vector;
//};
					error_string = sformatf("Heading spread_row\n");
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					if (heading)
					{
						for (int ii = 0; ii < heading->count; ii++)
						{
							error_string = sformatf("%d\t%s\n",ii,heading->str_vector[ii].c_str());
							error_msg(error_string, PHRQ_io::OT_CONTINUE);
						}
					}
					else
					{
						error_msg("heading is null", PHRQ_io::OT_CONTINUE);
					}

					error_string = sformatf("Data spread_row\n");
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					if (data)
					{
						for (int ii = 0; ii < data->count; ii++)
						{
							error_string = sformatf("%d\t%s\t%d\n",ii,data->str_vector[ii].c_str(),data->type_vector[ii]);
							error_msg(error_string, PHRQ_io::OT_CONTINUE);
						}
					}
					else
					{
						error_msg("Data is null", PHRQ_io::OT_CONTINUE);
					}
					error_string = sformatf("Units spread_row\n");
					error_msg(error_string, PHRQ_io::OT_CONTINUE);
					if (units)
					{
						for (int ii = 0; ii < units->count; ii++)
						{
							error_string = sformatf("%d\t%s\n",ii,units->str_vector[ii].c_str());
							error_msg(error_string, PHRQ_io::OT_CONTINUE);
						}
					}
					else
					{
						error_msg("Units is null", PHRQ_io::OT_CONTINUE);
					}
					free_check_null(char_string);
					continue;
				}
				temp_isotope.Set_isotope_name(token.c_str());

				/* read and save element name */
				{
					char *temp_iso_name = string_duplicate(token.c_str());
					const char* cptr1 = temp_iso_name;
					get_num(&cptr1, &dummy);
					temp_isotope.Set_isotope_number(dummy);
					if (cptr1[0] == '\0' || isupper((int)cptr1[0]) == FALSE)
					{
						error_msg("Expecting element name.", PHRQ_io::OT_CONTINUE);
						error_msg(line_save, PHRQ_io::OT_CONTINUE);
						input_error++;
						temp_iso_name = (char*)free_check_null(temp_iso_name);
						char_string = (char*)free_check_null(char_string);
						return (CParser::PARSER_ERROR);
					}
					temp_isotope.Set_elt_name(cptr1);
					temp_iso_name = (char*)free_check_null(temp_iso_name);
				}
				/* read and store isotope ratio */
				if (copy_token(token, &next_char) != CParser::TT_DIGIT)
				{
					input_error++;
					error_string = sformatf(
						"Expected numeric value for isotope ratio.");
					error_msg(error_string, CONTINUE);
					free_check_null(char_string);
					continue;
				}
				(void)sscanf(token.c_str(), SCANFORMAT, &dummy);
				temp_isotope.Set_ratio(dummy);
				temp_isotope.Set_ratio_uncertainty(NAN);

				/* read and store isotope ratio uncertainty */
				int j;
				if ((j = copy_token(token, &next_char)) != CParser::TT_EMPTY)
				{
					if (j != DIGIT)
					{
						input_error++;
						error_string = sformatf(
							"Expected numeric value for uncertainty in isotope ratio.");
						error_msg(error_string, PHRQ_io::OT_CONTINUE);
						free_check_null(char_string);
						continue;
					}
					(void)sscanf(token.c_str(), SCANFORMAT, &dummy);
					temp_isotope.Set_ratio_uncertainty(dummy);
				}
				temp_solution.Get_isotopes()[temp_isotope.Get_isotope_name()] = temp_isotope;
			}
			break;
		case 10:				/* water */
			{
				int j = copy_token(token, &next_char);
				if (j == EMPTY)
				{
					temp_solution.Set_mass_water(1.0);
				}
				else if (j != DIGIT)
				{
					input_error++;
					error_string = sformatf(
						"Expected numeric value for mass of water in solution.");
					error_msg(error_string, CONTINUE);
				}
				else
				{
					(void)sscanf(token.c_str(), SCANFORMAT, &dummy);
					temp_solution.Set_mass_water(dummy);
				}
			}
			break;
		case 11:				/* description */
		case 12:				/* desc */
		case 13:				/* descriptor */
			{
				temp_solution.Set_description(next_char);
			}
			break;
		case 14:				/* pressure */
		case 15:				/* press */
			{
				if (sscanf(next_char, SCANFORMAT, &dummy) == 1)
				{
					temp_solution.Set_patm(dummy);
				}
			}
			break;
		case 16:				/* pote, V */
			{
				if (sscanf(next_char, SCANFORMAT, &dummy) == 1)
				{
					temp_solution.Set_potV(dummy);
				}
			}
			break;
		case OPTION_DEFAULT:
/*
 *   Read concentration
 */
			{
				next_char = char_string;
				if (copy_token(token, &next_char) == LOWER)
				{
					free_check_null(char_string);
					continue;
				}
				cxxISolutionComp temp_comp(this->phrq_io);
				if (temp_comp.read(char_string, &temp_solution) == CParser::PARSER_ERROR)
				{
#ifdef SKIP
					input_error++;
					break;
#endif
				}
				initial_data_ptr->Get_comps()[temp_comp.Get_description()] = temp_comp;
				if (temp_comp.Get_pe_reaction().size() > 0)
				{
					CReaction temp_chem_reaction;
					initial_data_ptr->Get_pe_reactions()[temp_comp.Get_pe_reaction()] = temp_chem_reaction;
				}
			}
			break;
		}
		free_check_null(char_string);
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
/*
 *   fix up default units and default pe
 */
	std::map < std::string, cxxISolutionComp >::iterator it;
	for (it = initial_data_ptr->Get_comps().begin(); it != initial_data_ptr->Get_comps().end(); it++)
	{
		token = it->first;
		Utilities::str_tolower(token);
		if (it->second.Get_units().size() == 0)
		{
			it->second.Set_units(initial_data_ptr->Get_units().c_str());
		}
		else
		{
			bool alk = false;
			if (strstr(token.c_str(), "alk") == token.c_str())
				alk = true;
			std::string token1 = it->second.Get_units();
			if (check_units(token1, alk, true, initial_data_ptr->Get_units().c_str(), true) ==	CParser::PARSER_ERROR)
			{
				input_error++;
			}
			else
			{
				it->second.Set_units(token1.c_str());
			}
		}
		if (it->second.Get_pe_reaction().size() == 0)
		{
			it->second.Set_pe_reaction(initial_data_ptr->Get_default_pe());
		}
	}
	if (n_user >= 0)
	{
		Rxn_solution_map[n_user] = temp_solution;
	}
	else
	{
		unnumbered_solutions.push_back(temp_solution);
	}
	return (return_value);
}
/* ---------------------------------------------------------------------- */
class spread_row * Phreeqc::
string_to_spread_row(char *string)
/* ---------------------------------------------------------------------- */
{
	int j;
	std::string token;
	const char* cptr;
/*
 *   Allocate space
 */
	class spread_row* spread_row_ptr = new class spread_row;
	if (spread_row_ptr == NULL)
	{
		malloc_error();
		return spread_row_ptr;
	}
	spread_row_ptr->count = 0;
	spread_row_ptr->empty = 0;
	spread_row_ptr->string = 0;
	spread_row_ptr->number = 0;
	cptr = string;
/*
 *   Split by tabs, reallocate space
 */
	for (;;)
	{
		j = copy_token_tab(token, &cptr);
		if (j == EOL)
			break;
		spread_row_ptr->str_vector.push_back(token);
		if (j == EMPTY || token.size() == 0)
		{
			spread_row_ptr->empty++;
			spread_row_ptr->type_vector.push_back(EMPTY);
		}
		else if (j == UPPER || j == LOWER)
		{
			spread_row_ptr->string++;
			spread_row_ptr->type_vector.push_back(STRING);
		}
		else if (j == DIGIT)
		{
			spread_row_ptr->number++;
			spread_row_ptr->type_vector.push_back(NUMBER);
		}
		else
		{
			input_error++;
			error_msg("Unknown input in string_to_spread_row keyword.", CONTINUE);
			error_string = sformatf("\tcopy_token j: %d, token: %s\n", j, token.c_str());
			error_msg(error_string, CONTINUE);
			error_msg(line_save, CONTINUE);
		}
		spread_row_ptr->count++;
	}
	assert(spread_row_ptr->count == spread_row_ptr->str_vector.size());
	assert(spread_row_ptr->count == spread_row_ptr->type_vector.size());
	assert(spread_row_ptr->count == spread_row_ptr->empty + spread_row_ptr->string + spread_row_ptr->number);
	return (spread_row_ptr);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
spread_row_free(class spread_row *spread_row_ptr)
/* ---------------------------------------------------------------------- */
{
	if (spread_row_ptr == NULL)
		return (OK);
	spread_row_ptr->str_vector.clear();
	spread_row_ptr->type_vector.clear();
	delete spread_row_ptr;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
copy_token_tab(std::string& token, const char **cptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies from **cptr to *token until first tab is encountered.
 *
 *   Arguments:
 *      *token       output, place to store token
 *
 *     **cptr        input, character string to read token from
 *                   output, next position after token
 *   Returns:
 *      UPPER,
 *      LOWER,
 *      DIGIT,
 *      EMPTY,
 *      EOL,
 *      UNKNOWN.
 */
	//int i, return_value;
	int return_value;
	char c;
/*
 *   Strip leading spaces
 */
	token.clear();
	while ((c = **cptr) == ' ')
		(*cptr)++;
/*
 *   Check what we have
 */
	if (isupper((int) c) || c == '[')
	{
		return_value = UPPER;
	}
	else if (islower((int) c))
	{
		return_value = LOWER;
	}
	else if (isdigit((int) c) || c == '.' || c == '-')
	{
		return_value = DIGIT;
	}
	else if (c == '\0')
	{
		return_value = EOL;
		return (return_value);
	}
	else if (c == '\t')
	{
		return_value = EMPTY;
	}
	else
	{
		return_value = UNKNOWN;
	}
/*
 *   Begin copying to token
 */
	//i = 0;
	for (;;)
	{
		c = **cptr;
		if (c == '\t')
		{
			(*cptr)++;
			break;
		}
		else if (c == '\0')
		{
			break;
		}
		else
		{
			token.push_back(c);
			(*cptr)++;
			//i++;
		}
	}
	return (return_value);
}

/* ---------------------------------------------------------------------- */
 int Phreeqc::
get_option_string(const char **opt_list, int count_opt_list, const char **next_char)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read a line and check for options
 */
	int j;
	int opt_l, opt;
	const char *opt_ptr;
	char option[MAX_LENGTH];

	opt_ptr = *next_char;
	if (opt_ptr[0] == '-')
	{
		opt_ptr++;
		copy_token(option, &opt_ptr, &opt_l);
		if (find_option(&(option[1]), &opt, opt_list, count_opt_list, FALSE)
			== OK)
		{
			j = opt;
			*next_char = opt_ptr;
		}
		else
		{
			error_msg("Unknown option.", CONTINUE);
			error_msg(*next_char, CONTINUE);
			input_error++;
			j = OPTION_ERROR;
		}
	}
	else
	{
		copy_token(option, &opt_ptr, &opt_l);
		if (find_option(&(option[0]), &opt, opt_list, count_opt_list, TRUE)
			== OK)
		{
			j = opt;
			*next_char = opt_ptr;
		}
		else
		{
			j = OPTION_DEFAULT;
		}
	}
	return (j);
}

#if defined(PHREEQCI_GUI)
/* ---------------------------------------------------------------------- */
void Phreeqc::
free_spread(void)
/* ---------------------------------------------------------------------- */
{
	int i;
	spread_row_free(g_spread_sheet.heading);
	spread_row_free(g_spread_sheet.units);
	for (i = 0; i < g_spread_sheet.rows.size(); ++i)
	{
		spread_row_free(g_spread_sheet.rows[i]);
	}
	g_spread_sheet.rows.clear();
	g_spread_sheet.defaults.iso.clear();
	g_spread_sheet.defaults.redox = NULL;
	g_spread_sheet.defaults.units = NULL;

	g_spread_sheet.heading = NULL;
	g_spread_sheet.units   = NULL;
	g_spread_sheet.defaults.iso.clear();
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
add_row(class spread_row *spread_row_ptr)
/* ---------------------------------------------------------------------- */
{
	g_spread_sheet.rows.push_back(copy_row(spread_row_ptr));
}

/* ---------------------------------------------------------------------- */
class spread_row * Phreeqc::
copy_row(class spread_row *spread_row_ptr)
/* ---------------------------------------------------------------------- */
{
	spread_row *copy = new spread_row(*spread_row_ptr);
	if (copy == NULL)
		malloc_error();
	return copy;
}
#endif
