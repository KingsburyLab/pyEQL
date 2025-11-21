#include "Utils.h"
#include "Phreeqc.h"
#include "phqalloc.h"
#include "NameDouble.h"
#include "Exchange.h"
#include "Solution.h"
#include <time.h>

#if defined(PHREEQCI_GUI)
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#endif

/* ---------------------------------------------------------------------- */
double Phreeqc::
calc_alk(CReaction& rxn_ref)
/* ---------------------------------------------------------------------- */
{
	LDBLE return_value;
	class master* master_ptr;

	return_value = 0.0;
	class rxn_token* r_token = &rxn_ref.token[1];
	while (r_token->s != NULL)
	{
		master_ptr = r_token->s->secondary;
		if (master_ptr == NULL)
		{
			master_ptr = r_token->s->primary;
		}
		if (master_ptr == NULL)
		{
			error_string = sformatf(
				"Non-master species in secondary reaction, %s.",
				rxn_ref.token[0].s->name);
			error_msg(error_string, CONTINUE);
			input_error++;
			break;
		}
		return_value += r_token->coef * master_ptr->alk;
		//if (strcmp(r_token->name, "e-") == 0 && strcmp(rxn_ref.token[0].name,"e-") != 0)
		//{
		//	std::cerr << rxn_ref.token[0].name << "   Non-master species has e- in reaction.\n";
		//}
		r_token++;
	}
	return (return_value);
}/* ---------------------------------------------------------------------- */
double Phreeqc::
calc_delta_v(CReaction& r_ref, bool phase)
/* ---------------------------------------------------------------------- */
{
	/* calculate delta_v from molar volumes */
	double d_v = 0.0;
	if (phase)
	{
		/* for phases: reactants have coef's < 0, products have coef's > 0, v.v. for species */
		for (size_t i = 1; r_ref.Get_tokens()[i].s; i++)
		{
			if (!r_ref.Get_tokens()[i].s)
				continue;
			d_v += r_ref.Get_tokens()[i].coef * r_ref.Get_tokens()[i].s->logk[vm_tc];
		}
	}
	else
	{
		for (size_t i = 0; r_ref.token[i].name /*|| r_ptr->token[i].s*/; i++)
		{
			if (!r_ref.Get_tokens()[i].s)
				continue;
			d_v -= r_ref.Get_tokens()[i].coef * r_ref.Get_tokens()[i].s->logk[vm_tc];
		}
	}
	return d_v;
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
calc_dielectrics(LDBLE tc, LDBLE pa)
/* ---------------------------------------------------------------------- */
{
	/* Relative dielectric constant of pure water, eps as a function of (P, T)
	   Bradley and Pitzer, 1979, JPC 83, 1599.
	   (newer data in Fernandez et al., 1995, JPCRD 24, 33,
				  and Fernandez et al., 1997, JPCRD 26, 1125, show its correctness)
	   + d(eps)/d(P), Debye-Hueckel A and B, and Av (for Av, see Pitzer et al., 1984, JPCRD 13, p. 4)
	*/
	if (llnl_temp.size() > 0) return OK;
	if (tc > 350.)
	{
		tc = 350.;
	}
	LDBLE T = tc + 273.15;
	LDBLE u1 = 3.4279e2, u2 = -5.0866e-3, u3 = 9.469e-7, u4 = -2.0525,
		u5 = 3.1159e3, u6 = -1.8289e2, u7 = -8.0325e3, u8 = 4.2142e6,
		u9 = 2.1417;
	LDBLE d1000 = u1 * exp(T * (u2 + T * u3)); // relative dielectric constant at 1000 bar
	LDBLE c = u4 + u5 / (u6 + T);
	LDBLE b = u7 + u8 / T + u9 * T;
	LDBLE pb = pa * 1.01325; // pa in bar
	eps_r = d1000 + c * log((b + pb) / (b + 1e3)); // relative dielectric constant
	if (eps_r <= 0)
	{
		eps_r = 10.;
		warning_msg("Relative dielectric constant is negative.\nTemperature is out of range of parameterization.");
	}

	/* qe^2 / (eps_r * kB * T) = 4.803204e-10**2 / 1.38065e-16 / (eps_r * T)
							   = 1.671008e-3 (esu^2 / (erg/K)) / (eps_r * T) */
	LDBLE e2_DkT = 1.671008e-3 / (eps_r * T);

	DH_B = sqrt(8 * pi * AVOGADRO * e2_DkT * rho_0 / 1e3);  // Debye length parameter, 1/cm(mol/kg)^-0.5

	DH_A = DH_B * e2_DkT / (2. * LOG_10); //(mol/kg)^-0.5

	/* A0 in pitzer */
	if (pitzer_model || sit_model)
	{
		A0 = DH_B * e2_DkT / 6.0;
		if (pitzer_model && aphi != NULL)
		{
			calc_pitz_param(aphi, T, 298.15);
			A0 = aphi->p;
		}
	}

	/* Debye-Hueckel limiting slope = DH_B *  e2_DkT * RT * (d(ln(eps_r)) / d(P) - compressibility) */
	DH_Av = DH_B * e2_DkT * R_LITER_ATM * 1e3 * T * (c / (b + pb) * 1.01325 / eps_r - kappa_0 / 3.); // (cm3/mol)(mol/kg)^-0.5

	DH_B /= 1e8; // kappa, 1/Angstrom(mol/kg)^-0.5

	/* the Born functions, * 41.84 to give molal volumes in cm3/mol... */
	ZBrn = (-1 / eps_r + 1.0) * 41.84004;
	QBrn = c / (b + pb) / eps_r / eps_r * 41.84004;
	/* dgdP from subroutine gShok2 in supcrt92, g is neglected here (at tc < 300)...
	   and, dgdP is small. Better, adapt Wref to experimental Vm's */
	dgdP = 0;
	//if (tc > 150 && rho_0 < 1.0)
	//{
	//	LDBLE sc[7] = {1, -0.2037662e+01,  0.5747000e-02, -0.6557892e-05,
	//				0.6107361e+01, -0.1074377e-01,  0.1268348e-04};
	//	LDBLE csc[4] = {1, 0.3666666e+02, -0.1504956e-9,   0.5017997e-13};
	//	LDBLE sa = sc[1] + tc * (sc[2] + tc * sc[3]);
	//	LDBLE sb = sc[4] + tc * (sc[5] + tc * sc[6]);

	//	dgdP = - sa * sb * pow(1.0 - rho_0, sb - 1.0) * rho_0 * kappa_0 / 1.01325;

	//	LDBLE ft = pow((tc - 155.0)/300.0, 4.8) + csc[1] * pow((tc - 155.0)/300.0, 16.0);
	//	LDBLE dfdP   = ft * (-3.0 * csc[2] * pow(1000.0 - pb, 2) - 4.0 * csc[3] * pow(1000.0 - pb, 3)); 
	//	dgdP -= dfdP;
	//}

	return (OK);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
calc_rho_0(LDBLE tc, LDBLE pa)
/* ---------------------------------------------------------------------- */
{
	/* Density of pure water
        Wagner and Pruss, 2002, JPCRD 31, 387, eqn. 2.6, along the saturation pressure line +
		interpolation 0 - 300 oC, 0.006 - 1000 atm...
    */
	if (llnl_temp.size() > 0) return OK;
	if (tc > 350.)
	{
		if (need_temp_msg < 1)
		{
			std::ostringstream w_msg;
			w_msg << "Fitting range for dielectric constant of pure water is 0-350 C.\n";
			w_msg << "Fitting range for density along the saturation pressure line is 0-374 C,\n";
			w_msg << "                         for higher pressures up to 1000 atm    0-300 C.\n";
			w_msg << "Using temperature of 350 C for dielectric and density calculation.";
			warning_msg(w_msg.str().c_str());
			need_temp_msg++;
		}
		tc = 350.;
	}
	LDBLE T = tc + 273.15;
	//eqn. 2.6...
	LDBLE Tc = 647.096, th = 1 - T / Tc;
	LDBLE b1 = 1.99274064, b2 = 1.09965342, b3 = -0.510839303,
		b4 = -1.75493479, b5 = -45.5170352, b6 = -6.7469445e5;
    rho_0_sat = 322.0 * (1.0 + b1 * pow(th, (LDBLE) 1./3.) + b2 * pow(th, (LDBLE) 2./3.) + b3 * pow(th, (LDBLE) 5./3.) +\
               b4 * pow(th, (LDBLE) 16./3.) + b5 * pow(th, (LDBLE) 43./3.) + b6 * pow(th, (LDBLE) 110./3));
	//pressure...
	LDBLE p0 =  5.1880000E-02 + tc * (-4.1885519E-04 + tc * ( 6.6780748E-06 + tc * (-3.6648699E-08 + tc *  8.3501912E-11)));
	LDBLE p1 = -6.0251348E-06 + tc * ( 3.6696407E-07 + tc * (-9.2056269E-09 + tc * ( 6.7024182E-11 + tc * -1.5947241E-13)));
	LDBLE p2 = -2.2983596E-09 + tc * (-4.0133819E-10 + tc * ( 1.2619821E-11 + tc * (-9.8952363E-14 + tc *  2.3363281E-16)));
	LDBLE p3 =  7.0517647E-11 + tc * ( 6.8566831E-12 + tc * (-2.2829750E-13 + tc * ( 1.8113313E-15 + tc * -4.2475324E-18)));
	/* The minimal pressure equals the saturation pressure... */
	if (ah2o_x <= 1.0)
		p_sat = exp(11.6702 - 3816.44 / (T - 46.13)) * ah2o_x;
	else
		p_sat = exp(11.6702 - 3816.44 / (T - 46.13));
	//ah2o_x0 = ah2o_x; // for updating rho in model(): compare with new ah2o_x
	if (pa < p_sat || (use.Get_solution_ptr() && use.Get_solution_ptr()->Get_patm() < p_sat))
	{
		pa = p_sat;
	}
	if (!use.Get_gas_phase_in())
		patm_x = pa;
	pa -= (p_sat - 1e-6);
	rho_0 = rho_0_sat + pa * (p0 + pa * (p1 + pa * (p2 + sqrt(pa) * p3)));
	if (rho_0 < 0.01)
		rho_0 = 0.01;

	/* compressibility, d(ln(rho)) / d(P), 1/atm... */
	kappa_0 = (p0 + pa * (2 * p1 + pa * (3 * p2 + sqrt(pa) * 3.5 * p3))) / rho_0;

	return (rho_0 / 1e3);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
compute_gfw(const char *string, LDBLE * gfw)
/* ---------------------------------------------------------------------- */
{
/*
 *    Input:  string contains a chemical formula
 *    Output:  gfw contains the calculated gfw
 */
	std::string str(string);
	std::map<std::string, double>::iterator it;
	it = gfw_map.find(str);
	if (it != gfw_map.end())
	{
		*gfw = it->second;
		return OK;
	}

	int i;
	char token[MAX_LENGTH];
	const char* cptr;

	count_elts = 0;
	paren_count = 0;
	Utilities::strcpy_safe(token, MAX_LENGTH, string);
	cptr = token;
	if (get_elts_in_species(&cptr, 1.0) == ERROR)
	{
		return (ERROR);
	}
	*gfw = 0.0;
	for (i = 0; i < count_elts; i++)
	{
		if (elt_list[i].elt->gfw <= 0.0)
		{
			return (ERROR);
		}
		*gfw += elt_list[i].coef * (elt_list[i].elt)->gfw;
	}
	gfw_map[str] = *gfw;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
copy_token(char *token_ptr, const char **cptr, int *length)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies from **cptr to *token_ptr until first space is encountered.
 *
 *   Arguments:
 *      *token_ptr  output, place to store token
 *
 *     **cptr        input, character string to read token from
 *                  output, next position after token
 *
 *       length     output, length of token
 *
 *   Returns:
 *      UPPER,
 *      LOWER,
 *      DIGIT,
 *      EMPTY,
 *      UNKNOWN.
 */
	int i, return_value;
	char c;

/*
 *   Read to end of whitespace
 */
	while (isspace((int) (c = **cptr)))
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
		return_value = EMPTY;
	}
	else
	{
		return_value = UNKNOWN;
	}
/*
 *   Begin copying to token
 */
	i = 0;
	while ((!isspace((int) (c = **cptr))) &&
		   /*              c != ',' && */
		   c != ';' && c != '\0')
	{
		token_ptr[i] = c;
		(*cptr)++;
		i++;
	}
	token_ptr[i] = '\0';
	*length = i;
	return (return_value);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
copy_token(std::string &token, const char **cptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies from **cptr to *token until first space is encountered.
 *
 *   Arguments:
 *      &token_ptr  output, place to store token
 *
 *     **cptr        input, character string to read token from
 *                  output, next position after token
 *
 *   Returns:
 *      UPPER,
 *      LOWER,
 *      DIGIT,
 *      EMPTY,
 *      UNKNOWN.
 */
	int return_value;
	char c;

/*
 *   Read to end of whitespace
 */
	token.clear();
	while (isspace((int) (c = **cptr)))
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
		return_value = EMPTY;
	}
	else
	{
		return_value = UNKNOWN;
	}
/*
 *   Begin copying to token
 */
	char c_char[2];
	c_char[1] = '\0';
	while ((!isspace((int) (c = **cptr))) &&
		   /*              c != ',' && */
		   c != ';' && c != '\0')
	{
		c_char[0] = c;
		token.append(c_char);
		(*cptr)++;
	}
	return (return_value);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
dup_print(const char* cptr, int emphasis)
/* ---------------------------------------------------------------------- */
{
/*
 *   print character string to output and logfile
 *   if emphasis == TRUE the print is set off by
 *   a row of dashes before and after the character string.
 *
 */
	int l;

	if (pr.headings == FALSE)
		return (OK);
	std::string save_in(cptr);
	l = (int) strlen(cptr);
	if (emphasis == TRUE)
	{
		std::string dash;
		dash.resize(l, '-');
		output_msg(sformatf("%s\n%s\n%s\n\n", dash.c_str(), save_in.c_str(), dash.c_str()));
		log_msg(sformatf("%s\n%s\n%s\n\n", dash.c_str(), save_in.c_str(), dash.c_str()));
	}
	else
	{
		output_msg(sformatf("%s\n\n", save_in.c_str()));
		log_msg(sformatf("%s\n\n", save_in.c_str()));
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
equal(LDBLE a, LDBLE b, LDBLE eps)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks equality between two LDBLE precision numbers
 */
	if (fabs(a - b) <= eps)
		return (TRUE);
	return (FALSE);
}

/* ---------------------------------------------------------------------- */
void * Phreeqc::
free_check_null(void *ptr)
/* ---------------------------------------------------------------------- */
{
	if (ptr != NULL)
	{
		PHRQ_free(ptr);
	}
	return (NULL);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
get_token(const char** eqnaddr, std::string& string, LDBLE* l_z, int* l)
/* ---------------------------------------------------------------------- */
/*
 *   Function finds next species in equation, coefficient has already
 *   been removed.  Also determines charge on the species.
 *
 *   Arguments:
 *      *eqnaddr   input, pointer to position in eqn to start parsing
 *                 output, pointer to a pointer to next position in eqn to start
 *                         parsing.
 *      *string    input pointer to place to store token
 *      *l_z        charge on token
 *      *l         length of token
 *
 *   Returns:
 *      ERROR,
 *      OK.
 */
{
	int i, j;
	int ltoken, lcharge;
	char c;
	const char* cptr, * ptr1, * rest;
	char charge[MAX_LENGTH];

	string.clear();
	rest = *eqnaddr;
	cptr = *eqnaddr;
	i = 0;
	/*
	 *   Find end of token or beginning of charge
	 */
	while (((c = *cptr) != '+') && (c != '-') && (c != '=') && (c != '\0'))
	{
		string.push_back(c);
		i++;
		if (c == '[')
		{
			cptr++;
			while ((c = *cptr) != ']')
			{
				if (c == '\0')
				{
					error_string = sformatf(
						"No final bracket \"]\" for element name, %s.",
						string.c_str());
					error_msg(error_string, CONTINUE);
					return (ERROR);
				}
				string.push_back(c);
				i++;
				cptr++;
			}
			string.push_back(c);
			i++;
		}

		cptr++;
	}
	ltoken = i;
	/*
	 *   Check for an empty string
	 */
	if (i == 0)
	{
		error_string = sformatf("NULL string detected in get_token, %s.", rest);
		error_msg(error_string, CONTINUE);
		return (ERROR);
	}
	/*
	 *   End of token is = or \0, charge is zero
	 */
	if (c == '=' || c == '\0')
	{
		*eqnaddr = cptr;
		lcharge = 0;
		*l_z = 0.0;
	}
	else
	{
		/*
		 *   Copy characters into charge until next species or end is detected
		 */
		j = 0;
		ptr1 = cptr;
		while ((isalpha((int)(c = *ptr1)) == FALSE) &&
			(c != '(') &&
			(c != ')') &&
			(c != ']') && (c != '[') && (c != '=') && (c != '\0'))
		{
			charge[j++] = c;
			/* error if no more space */
			if (j >= MAX_LENGTH)
			{
				error_msg("The charge on a species has exceeded MAX_LENGTH characters.",
					CONTINUE);
				return (ERROR);
			}
			ptr1++;
		}
		/*
		 *   Go back to last + or - if not end of side,
		 *   everything before the last + or - in charge is part of the charge
		 */
		if ((c != '=') && (c != '\0'))
		{
			while (((c = *ptr1) != '+') && (c != '-'))
			{
				j--;
				ptr1--;
			}
		}
		charge[j] = '\0';
		lcharge = j;
		*eqnaddr = ptr1;
		/*
		 *   Charge has been written, now need to check if charge has legal format
		 */
		if (get_charge(charge, MAX_LENGTH, l_z) == OK)
		{
			string.append(charge);
		}
		else
		{
			return (ERROR);
		}
	}
	*l = ltoken + lcharge;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
isamong(char c, const char *s_l)
/* ---------------------------------------------------------------------- */
/*
 *   Function checks if c is among the characters in the string s
 *
 *   Arguments:
 *      c     input, character to check
 *     *s     string of characters
 *
 *   Returns:
 *      TRUE  if c is in set,
 *      FALSE if c in not in set.
 */
{
	int i;

	for (i = 0; s_l[i] != '\0'; i++)
	{
		if (c == s_l[i])
		{
			return (TRUE);
		}
	}
	return (FALSE);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
islegit(const char c)
/* ---------------------------------------------------------------------- */
/*
 *   Function checks for legal characters for chemical equations
 *
 *   Argument:
 *      c     input, character to check
 *
 *   Returns:
 *      TRUE  if c is in set,
 *      FALSE if c in not in set.
 */
{
	if (isalpha((int) c) || isdigit((int) c) || isamong(c, "+-=().:_[]"))
	{
		return (TRUE);
	}
	return (FALSE);
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
malloc_error(void)
/* ---------------------------------------------------------------------- */
{
	error_msg("NULL pointer returned from malloc or realloc.", CONTINUE);
	error_msg("Program terminating.", STOP);
	return;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
print_centered(const char *string)
/* ---------------------------------------------------------------------- */
{
	int i, l, l1, l2;
	char token[MAX_LENGTH];

	l = (int) strlen(string);
	l1 = (79 - l) / 2;
	l2 = 79 - l - l1;
	for (i = 0; i < l1; i++)
		token[i] = '-';
	token[i] = '\0';
	Utilities::strcat_safe(token, MAX_LENGTH, string);
	for (i = 0; i < l2; i++)
		token[i + l1 + l] = '-';
	token[79] = '\0';
	output_msg(sformatf("%s\n\n", token));
	return (OK);
}
bool Phreeqc::
replace(const char *str1, const char *str2, std::string & str)
{
	size_t pos = str.find(str1);
	if (pos != std::string::npos)
	{
		size_t l = strlen(str1);
		str.replace(pos, l, str2);
		return true;
	}
	return false;
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
replace(const char *str1, const char *str2, char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function replaces str1 with str2 in str
 *
 *   Arguments:
 *      str1     search str for str1
 *      str2     replace str1 if str1 found in str
 *      str      string to be searched
 *
 *   Returns
 *      TRUE     if string was replaced
 *      FALSE    if string was not replaced
 */
	int l, l1, l2;
	char* ptr_start;

	ptr_start = strstr(str, str1);
/*
 *   Str1 not found, return
 */
	if (ptr_start == NULL)
		return (FALSE);
/*
 *   Str1 found, replace Str1 with Str2
 */
	l = (int) strlen(str);
	l1 = (int) strlen(str1);
	l2 = (int) strlen(str2);
/*
 *   Make gap in str long enough for str2
 */
	/* The plus one includes the terminating NULL */
	memmove(ptr_start + l2, ptr_start + l1, l - (ptr_start - str + l1) + 1);
/*
 *   Copy str2 into str
 */
	memcpy(ptr_start, str2, l2);
	return (TRUE);
}
void Phreeqc::
replace(std::string &stds, const char* str1, const char* str2)
{
	size_t pos, l;
	l = strlen(str1);
	while ((pos = stds.find(str1)) != std::string::npos) {
		stds.replace(pos, l, str2);
	}
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
space(void **ptr, int i, int *max, int struct_size)
/* ---------------------------------------------------------------------- */
{
/*
 *   Routine has 4 functions, allocate space, reallocate space, test to
 *   determine whether space is available, and free space.
 *
 *   Arguments:
 *      ptr          pointer to malloced space
 *      i            value for test
 *         i = INIT,          allocates space
 *         i >= 0 && i < max, space is available, return
 *         i >= max,          reallocate space
 *         i = FREE,          free space.
 *      max          maximum value for i with available space
 *      struct size  size of structure to be allocated
 */
/*
 *   Return if space exists
 */
	if ((i >= 0) && (i + 1 < *max))
	{
		return;
	}
/*
 *   Realloc space
 */
	if (i + 1 >= *max)
	{
		if (*max > 1000)
		{
			*max += 1000;
		}
		else
		{
			*max *= 2;
		}
		if (i + 1 > *max)
			*max = i + 1;
		*ptr = PHRQ_realloc(*ptr, (size_t) (*max) * struct_size);
		if (*ptr == NULL)
			malloc_error();
		return;
	}
/*
 *   Allocate space
 */
	if (i == INIT)
	{
/*		free(*ptr); */
		*ptr = PHRQ_malloc((size_t) (*max) * struct_size);
		if (*ptr == NULL)
			malloc_error();
		return;
	}
/*
 *   Free space
 */
/*
	if ( i == FREE ) {
		free(*ptr);
		return;
	}
 */
/*
 *   Error
 */
	error_msg("Illegal argument to function space.", CONTINUE);
	error_msg("Program terminating.", STOP);
	return;
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
squeeze_white(char *s_l)
/* ---------------------------------------------------------------------- */
/*
 *   Delete all white space from string s
 *
 *   Argument:
 *      *s_l input, character string, possibly containing white space
 *           output, character string with all white space removed
 *
 *   Return: void
 */
{
	int i, j;

	for (i = j = 0; s_l[i] != '\0'; i++)
	{
		if (!isspace((int) s_l[i]))
			s_l[j++] = s_l[i];
	}
	s_l[j] = '\0';
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
str_tolower(char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Replaces string, str, with same string, lower case
 */
	char* ptr;
	ptr = str;
	while (*ptr != '\0')
	{
		*ptr = (char) tolower(*ptr);
		ptr++;
	}
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
str_toupper(char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Replaces string, str, with same string, lower case
 */
	char* ptr;
	ptr = str;
	while (*ptr != '\0')
	{
		*ptr = (char) toupper(*ptr);
		ptr++;
	}
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
strcmp_nocase(const char *str1, const char *str2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare two strings disregarding case
 */
	int c1, c2;
	while ((c1 = tolower(*str1++)) == (c2 = tolower(*str2++)))
	{
		if (c1 == '\0')
			return (0);
	}
	if (c1 < c2)
		return (-1);
	return (1);
}
void Phreeqc::str_tolower(std::string &name)
{
	std::transform(name.begin(), name.end(), name.begin(), tolower);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
strcmp_nocase_arg1(const char *str1, const char *str2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare two strings disregarding case
 */
	int c1, c2;
	while ((c1 = tolower(*str1++)) == (c2 = *str2++))
	{
		if (c1 == '\0')
			return (0);
	}
	if (c1 < c2)
		return (-1);
	return (1);
}

/* ---------------------------------------------------------------------- */
char * Phreeqc::
#if !defined(NDEBUG) && defined(WIN32_MEMORY_DEBUG)
_string_duplicate(const char *token, const char *szFileName, int nLine)
#else
string_duplicate(const char *token)
#endif
/* ---------------------------------------------------------------------- */
{
	int l;
	char *str;

	if (token == NULL)
		return NULL;
	l = (int) strlen(token);
#if !defined(NDEBUG) && defined(WIN32_MEMORY_DEBUG)
	str = (char *) _malloc_dbg((size_t) (l + 1) * sizeof(char), _NORMAL_BLOCK, szFileName, nLine);
#else
	str = (char *) PHRQ_malloc(((size_t)l + 1) * sizeof(char));
#endif

	if (str == NULL)
		malloc_error();
	strcpy(str, token);
	return (str);
}

/* ---------------------------------------------------------------------- */
const char * Phreeqc::
string_hsave(const char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *      Save character string str
 *
 *      Arguments:
 *         str   input string to save.
 *
 *      Returns:
 *         starting address of saved string (str)
 */
	if (str == NULL) return (NULL);
	std::map<std::string, std::string *>::const_iterator it;
	it = strings_map.find(str);
	if (it != strings_map.end())
	{
		return (it->second->c_str());
	}

	std::string *stdstr = new std::string(str);
	strings_map[*stdstr] = stdstr;
	return(stdstr->c_str());
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
strings_map_clear()
/* ---------------------------------------------------------------------- */
{
/*
 *      Save character string str
 *
 *      Arguments:
 *         str   input string to save.
 *
 *      Returns:
 *         starting address of saved string (str)
 */
	std::map<std::string, std::string *>::iterator it;
	for (it = strings_map.begin(); it != strings_map.end(); it++)
	{
		delete it->second;
	}
	strings_map.clear();
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
under(LDBLE xval)
/* ---------------------------------------------------------------------- */
{
/*
 *   Exponentiate a number, x, but censor large and small numbers
 *   log values less than MIN_LM are set to 0
 *   log values greater than MAX_LM are set to 10**MAX_LM
 */
/*	if (xval < MIN_LM) { */
	if (xval < -40.)
	{
		return (0.0);
	}
	if (xval > MAX_LM)  
	{
		return ( MAX_M );
	}
	return (pow ((LDBLE) 10.0, xval));
}
#ifndef PHREEQCI_GUI
/* ---------------------------------------------------------------------- */
int Phreeqc::
status(int count, const char *str, bool rk_string)
/* ---------------------------------------------------------------------- */
{
	char sim_str[20];
	char state_str[45];
	char spin_str[2];
	clock_t t2;

	if (pr.status == FALSE || phast == TRUE)
		return (OK);

	if (state == INITIALIZE)
	{
		screen_string = sformatf("\n%-80s", "Initializing...");
		screen_msg(screen_string.c_str());
		status_on = true;
		return (OK);
	}
#ifdef NPP
	t2 = clock();
	if (((state < ADVECTION && reaction_step < count_total_steps) || 
		(state == ADVECTION && (advection_step < count_ad_shifts || cell_no < count_cells)) || 
		(state == TRANSPORT && (transport_step < count_shifts || (mixrun < nmix /*&&
		                        cell_no < count_cells*/))))
		&& (int) (1e3 / CLOCKS_PER_SEC * (t2 - status_timer)) < status_interval)
		return (OK);
#endif

	switch (state)
	{
	case INITIALIZE:
		break;
	case TRANSPORT:
		if (str != NULL)
		{
			if (rk_string)
			{

				screen_string = screen_string.substr(0, 43);
				screen_string.append(str);
				status_string = screen_string;
			}
			else
			{
				screen_string = "\r";
				screen_string.append(str);
				status_string = screen_string;
			}
			status_on = true;
		}
	case PHAST:
		break;
	default:
		// if str not NULL, print it
		if (str != NULL && !rk_string)
		{
			screen_string = "\r";
			screen_string.append(str);
			status_string = screen_string;
		}
		else
		// print state
		{
			std::string stdstr;
			if (str != NULL && rk_string)
			{
				stdstr = str;
			}
			snprintf(sim_str, sizeof(sim_str), "\rSimulation %d.", simulation);
			snprintf(state_str, sizeof(state_str), " ");
			snprintf(spin_str, sizeof(spin_str), " ");
			switch (state)
			{
			default:
				break;
			case INITIAL_SOLUTION:
				snprintf(state_str, sizeof(state_str), "Initial solution %d.", use.Get_solution_ptr()->Get_n_user());
				break;
			case INITIAL_EXCHANGE:
				snprintf(state_str, sizeof(state_str), "Initial exchange %d.", use.Get_exchange_ptr()->Get_n_user());
				break;
			case INITIAL_SURFACE:
				snprintf(state_str, sizeof(state_str), "Initial surface %d.", use.Get_surface_ptr()->Get_n_user());
				break;
			case INVERSE:
				snprintf(state_str, sizeof(state_str), "Inverse %d. Models = %d.", use.Get_inverse_ptr()->n_user, count);
				break;
			case REACTION:
				if (use.Get_kinetics_in() == TRUE)
				{
					snprintf(state_str, sizeof(state_str), "Kinetic step %d.", reaction_step);
				}
				else
				{
					snprintf(state_str, sizeof(state_str), "Reaction step %d.", reaction_step);
				}
				break;
			case ADVECTION:
				snprintf(state_str, sizeof(state_str), "Advection, shift %d.", advection_step);
				break;
			}
			spinner++;
			if (spinner == 1)
			{
				spin_str[0] = '/';
			}
			else if (spinner == 2)
			{
				spin_str[0] = '-';
			}
			else
			{
				spin_str[0] = '\\';
				spinner = 0;
			}
			if (use.Get_kinetics_in() == TRUE)
			{
				screen_string = sformatf("%-15s%-27s%38s", sim_str, state_str, stdstr.c_str());
				status_string = screen_string;
			}
			else
			{
				screen_string = sformatf("%-15s%-27s%1s%45s", sim_str, state_str, spin_str, stdstr.c_str());
				status_string = screen_string;
			}
		}
		status_on = true;
		break;
	}

#ifndef NPP
	t2 = clock();
	if ((int) (1e3 / CLOCKS_PER_SEC * (t2 - status_timer)) > status_interval)
#endif
	{
		status_timer = t2;
		screen_msg(status_string.c_str());
		status_string.clear();
	}
	return (OK);
}
#endif /*PHREEQCI_GUI */
# include	<assert.h>

/* ---------------------------------------------------------------------- */
int Phreeqc::
string_trim(char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function trims white space from left and right of string
 *
 *   Arguments:
 *      str      string to trime
 *
 *   Returns
 *      TRUE     if string was changed
 *      FALSE    if string was not changed
 *      EMPTY    if string is all whitespace
 */
	int i, l, start, end, length;
	char* ptr_start;

	l = (int) strlen(str);
	/*
	 *   leading whitespace
	 */
	for (i = 0; i < l; i++)
	{
		if (isspace((int) str[i]))
			continue;
		break;
	}
	if (i == l)
		return (EMPTY);
	start = i;
	ptr_start = &(str[i]);
	/*
	 *   trailing whitespace
	 */
	for (i = l - 1; i >= 0; i--)
	{
		if (isspace((int) str[i]))
			continue;
		break;
	}
	end = i;
	if (start == 0 && end == l)
		return (FALSE);
	length = end - start + 1;
	memmove((void *) str, (void *) ptr_start, (size_t) length);
	str[length] = '\0';

	return (TRUE);
}
void Phreeqc::string_trim_left(std::string& str)
{
	const std::string& chars = "\t\n ";
	str.erase(0, str.find_first_not_of(chars));
}
void Phreeqc::string_trim_right(std::string& str)
{
	const std::string& chars = "\t\n ";
	str.erase(str.find_last_not_of(chars) + 1);
}
void Phreeqc::string_trim(std::string& str)
{
	const std::string& chars = "\t\n ";
	str.erase(0, str.find_first_not_of(chars));
	str.erase(str.find_last_not_of(chars) + 1);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
string_trim_right(char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function trims white space from right of string
 *
 *   Arguments:
 *      str      string to trime
 *
 *   Returns
 *      TRUE     if string was changed
 *      FALSE    if string was not changed
 *      EMPTY    if string is all whitespace
 */
	int i, l, end, length;

	l = (int) strlen(str);
	for (i = l - 1; i >= 0; i--)
	{
		if (isspace((int) str[i]))
			continue;
		break;
	}
	end = i;
	length = end + 1;
	str[length] = '\0';
	if (end == 0)
		return (EMPTY);
	if (end == l)
		return (FALSE);
	return (TRUE);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
string_trim_left(char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function trims white space from left of string
 *
 *   Arguments:
 *      str      string to trime
 *
 *   Returns
 *      TRUE     if string was changed
 *      FALSE    if string was not changed
 *      EMPTY    if string is all whitespace
 */
	int i, l, start, end, length;
	char* ptr_start;

	l = (int) strlen(str);
	/*
	 *   leading whitespace
	 */
	for (i = 0; i < l; i++)
	{
		if (isspace((int) str[i]))
			continue;
		break;
	}
	if (i == l)
		return (EMPTY);
	start = i;
	ptr_start = &(str[i]);
	end = l;
	if (start == 0 && end == l)
		return (FALSE);
	length = end - start + 1;
	memmove((void *) str, (void *) ptr_start, (size_t) length);
	str[length] = '\0';

	return (TRUE);
}

/* ---------------------------------------------------------------------- */
char * Phreeqc::
string_pad(const char *str, int i)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function returns new string padded to width i
 *       or returns str if longer
 *   Arguments:
 *      str      string to pad
 *
 *   Returns
 *      new string of with i
 */
	int j, l, max;
	char *str_ptr;

	l = (int) strlen(str);
	max = l;
	if (l < i)
		max = i;
	str_ptr = (char *) PHRQ_malloc((((size_t)max + 1) * sizeof(char)));
	if (str_ptr == NULL)
		malloc_error();
	else
		strcpy(str_ptr, str);
	if (i > l)
	{
		for (j = l; j < i; j++)
		{
			str_ptr[j] = ' ';
		}
		str_ptr[i] = '\0';
	}
	return (str_ptr);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
get_input_errors()
/* ---------------------------------------------------------------------- */
{
	if (input_error == 0)
	{
		return phrq_io->Get_io_error_count();
	}
	return input_error;
}
