#include "Phreeqc.h"
#include "phqalloc.h"
#include "Exchange.h"
#include "Solution.h"

#if defined(PHREEQCI_GUI)
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#endif

/* ---------------------------------------------------------------------- */
int Phreeqc::
sit_init(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Initialization for SIT
 */
	sit_model = FALSE;
	sit_params.clear();
	OTEMP = -100.;
	OPRESS = -100.;
	return OK;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
sit_tidy(void)
/* ---------------------------------------------------------------------- */
{
	/*
	 *      Make lists of species for cations, anions, neutral
	 */
	int i, j;
	/*
	* Ensure new parameters are calculated
	*/
	OTEMP = -100.;
	OPRESS = -100.;
	/*
	 *  allocate pointers to species structures
	 */
	spec.clear();
	spec.resize(3 * s.size(), NULL);

	cations = &spec[0];
	neutrals = &(spec[s.size()]);
	anions = &(spec[2 * s.size()]);
	sit_MAXCATIONS = (int)s.size();
	sit_FIRSTANION = 2 * (int)s.size();
	sit_MAXNEUTRAL = (int)s.size();
	sit_count_cations = 0;
	sit_count_anions = 0;
	sit_count_neutrals = 0;
	if (itmax < 200) itmax = 200;
	/*
	 *  allocate other arrays for SIT
	 */
	sit_IPRSNT.resize(3 * s.size());
	sit_M.resize(3 * s.size());
	sit_LGAMMA.resize(3 * s.size());

	for (i = 0; i < (int)s.size(); i++)
	{
		if (s[i] == s_eminus)
			continue;
		if (s[i] == s_h2o)
			continue;
		if (s[i]->type == EX || s[i]->type == SURF) 
			continue;
		if (s[i]->z < -.001)
		{
			anions[sit_count_anions++] = s[i];
		}
		else if (s[i]->z > .001)
		{
			cations[sit_count_cations++] = s[i];
		}
		else
		{
			neutrals[sit_count_neutrals++] = s[i];
		}
	}
	/*
	 * no ethetas
	 */
	/*
	 *  put species numbers in sit_params
	 */
	for (i = 0; i < (int)sit_params.size(); i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (sit_params[i]->species[j] == NULL)
				continue;
			sit_params[i]->ispec[j] = sit_ISPEC(sit_params[i]->species[j]);
			if ((j < 2 && sit_params[i]->ispec[j] == -1) ||
				(j == 3
				 && (sit_params[i]->type == TYPE_PSI
					 || sit_params[i]->type == TYPE_ZETA)
				 && sit_params[i]->ispec[j] == -1))
			{
				input_error++;
				error_string = sformatf(
						"Species for Pitzer parameter not defined in SOLUTION_SPECIES, %s",
						sit_params[i]->species[j]);
				error_msg(error_string, CONTINUE);
			}
		}
	}	/* remake map */
	{
		sit_param_map.clear();
		for (int j = 0; j < (int)sit_params.size(); j++)
		{	
			std::set< std::string > header;
			for (int i = 0; i < 3; i++)
			{
				if (sit_params[j]->species[i] != NULL) header.insert(sit_params[j]->species[i]);
			}
			std::ostringstream key_str;
			key_str << sit_params[j]->type << " ";
			std::set< std::string >::iterator it = header.begin();
			for(; it != header.end(); ++it)
			{
				key_str << *it << " ";
			}
			std::string key = key_str.str().c_str();
			sit_param_map[key] = j;
		}
		assert ((int) sit_param_map.size() == (int)sit_params.size());
	}
	if (get_input_errors() > 0) return (ERROR);
	return OK;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
sit_ISPEC(const char *name)
/* ---------------------------------------------------------------------- */
/*
 *      Find species number in spec for character string species name
 */
{
	int i;
	for (i = 0; i < 3 * (int)s.size(); i++)
	{
		if (spec[i] == NULL)
			continue;
		if (name == spec[i]->name)
		{
			return (i);
		}
	}
	return (-1);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_sit(void)
/* ---------------------------------------------------------------------- */
{
  /*
   *      Reads advection information
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
  /*
   *   Read advection parameters: 
   *        number of cells;
   *        number of shifts;
   */
  int n;
  class pitz_param *pzp_ptr;
  pitz_param_type pzp_type;

  int return_value, opt, opt_save;
  const char* next_char;
  const char *opt_list[] = {
    "epsilon",					/* 0 */
    "epsilon1"					/* 1 */
  };
  int count_opt_list = 2;
  /*
   *   Read lines
   */
  opt_save = OPTION_ERROR;
  return_value = UNKNOWN;
  n = -1;
  pzp_type = TYPE_Other;
  pitzer_pe = TRUE;
  for (;;)
  {
    opt = get_option(opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      opt = opt_save;
    }
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_DEFAULT:
      pzp_ptr = pitz_param_read(line, n);
	  if (pzp_ptr != NULL)
	  {
		  pzp_ptr->type = pzp_type;
		  sit_param_store(pzp_ptr);
      }
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg("Unknown input in SIT keyword.", CONTINUE);
      error_msg(line_save, CONTINUE);
      break;
    case 0:				/* epsilon */
      pzp_type = TYPE_SIT_EPSILON;
      n = 2;
      opt_save = OPTION_DEFAULT;
      break;
    case 1:				/* epsilon1 */
      pzp_type = TYPE_SIT_EPSILON_MU;
      n = 2;
      opt_save = OPTION_DEFAULT;
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  sit_model = TRUE;
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
calc_sit_param(class pitz_param *pz_ptr, LDBLE TK, LDBLE TR)
/* ---------------------------------------------------------------------- */
{
	LDBLE param;
	/*
	 */

	if (fabs(TK - TR) < 0.01)
	{
		param = pz_ptr->a[0];
	}
	else
	{
		param = (pz_ptr->a[0] +
			 pz_ptr->a[1] * (1.e0 / TK - 1.e0 / TR) +
			 pz_ptr->a[2] * log(TK / TR) +
			 pz_ptr->a[3] * (TK - TR) + 
			 pz_ptr->a[4] * (TK * TK - TR * TR));
	}
	pz_ptr->p = param;
	switch (pz_ptr->type)
	{
	case TYPE_SIT_EPSILON:
		pz_ptr->U.eps = param;
		break;
	case TYPE_SIT_EPSILON_MU:
		pz_ptr->U.eps1 = param;
		break;
	case TYPE_Other:
	default:
		error_msg("Should not be TYPE_Other in function calc_sit_param",
				  STOP);
		break;
	}
	return OK;
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
sit(void)
/* ---------------------------------------------------------------------- */
{
  int i, i0, i1;
  LDBLE param, z0, z1;
  LDBLE A, AGAMMA, T;
	/*
	   LDBLE CONV, XI, XX, OSUM, BIGZ, DI, F, XXX, GAMCLM, 
	   CSUM, PHIMAC, OSMOT, BMXP, ETHEAP, CMX, BMX, PHI,
	   BMXPHI, PHIPHI, AW, A, B;
	 */
/*
	LDBLE CONV, XI, XX, OSUM, BIGZ, DI, F, XXX, GAMCLM, CSUM, PHIMAC, OSMOT,
		B;
*/
	LDBLE XI, XX, OSUM, DI, F, OSMOT, B;
	LDBLE I, TK;
	/*
	   C
	   C     INITIALIZE
	   C
	 */
	//CONV = 1.0 / LOG_10;
	XI = 0.0e0;
	XX = 0.0e0;
	OSUM = 0.0e0;
	/*n
	   I = *I_X;
	   TK = *TK_X;
	 */
	I = mu_x;
	TK = tk_x;
	/*      DH_AB(TK, &A, &B); */
	/*
	   C
	   C     TRANSFER DATA FROM TO sit_M
	   C
	 */
	double log_min = log10(MIN_TOTAL);
	for (size_t j = 0; j < s_list.size(); j++)
	{
		i = s_list[j];
		if (spec[i]->lm > log_min)
		{
			sit_M[i] = under(spec[i]->lm);
		}
		else
		{
			sit_M[i] = 0.0;
		}
	}
	//for (i = 0; i < 3 * (int)s.size(); i++)
	//{
	//	sit_IPRSNT[i] = FALSE;
	//	sit_M[i] = 0.0;
	//	if (spec[i] != NULL && spec[i]->in == TRUE)
	//	{
	//		if (spec[i]->type == EX ||
	//			spec[i]->type == SURF || spec[i]->type == SURF_PSI)
	//			continue;
	//		sit_M[i] = under(spec[i]->lm);
	//		if (sit_M[i] > MIN_TOTAL)
	//			sit_IPRSNT[i] = TRUE;
	//	}
	//}
	/*
	   C
	   C     COMPUTE SIT COEFFICIENTS' TEMPERATURE DEPENDENCE
	   C
	 */
	PTEMP_SIT(TK);
	for (size_t j = 0; j < s_list.size(); j++)
	{
		int i = s_list[j];
		sit_LGAMMA[i] = 0.0;
		XX = XX + sit_M[i] * fabs(spec[i]->z);
		XI = XI + sit_M[i] * spec[i]->z * spec[i]->z;
		OSUM = OSUM + sit_M[i];
	}
	//for (i = 0; i < 2 * (int)s.size() + sit_count_anions; i++)
	//{
	//	sit_LGAMMA[i] = 0.0;
	//	if (sit_IPRSNT[i] == TRUE)
	//	{
	//		XX = XX + sit_M[i] * fabs(spec[i]->z);
	//		XI = XI + sit_M[i] * spec[i]->z * spec[i]->z;
	//		OSUM = OSUM + sit_M[i];
	//	}
	//}
	I = XI / 2.0e0;
	I = mu_x;   // Added equation for MU
	DI = sqrt(I);
	/*
	   C
	   C     CALCULATE F & GAMCLM
	   C
	 */
	AGAMMA = 3*sit_A0; /* Grenthe p 379 */
	A = AGAMMA / LOG_10;
	/*
	*  F is now for log10 gamma
	*/

	B = 1.5;
	F = -A * (DI / (1.0e0 + B * DI));


	/*OSMOT = -(sit_A0) * pow(I, 1.5e0) / (1.0e0 + B * DI);*/
	T = 1.0 + B*DI;
	OSMOT = -2.0*A/(B*B*B)*(T - 2.0*log(T) - 1.0/T);
	/*
	 *  Sums for sit_LGAMMA, and OSMOT
	 *  epsilons are tabulated for log10 gamma (not ln gamma)
	 */
	for (size_t j = 0; j < param_list.size(); j++)
	{
		int i = param_list[j];
		i0 = sit_params[i]->ispec[0];
		i1 = sit_params[i]->ispec[1];
		//if (sit_IPRSNT[i0] == FALSE || sit_IPRSNT[i1] == FALSE) continue;
		z0 = spec[i0]->z;
		z1 = spec[i1]->z;
		param = sit_params[i]->p;
		switch (sit_params[i]->type)
		{
		case TYPE_SIT_EPSILON:
			sit_LGAMMA[i0] += sit_M[i1] * param;
			sit_LGAMMA[i1] += sit_M[i0] * param;
			if (z0 == 0.0 && z1 == 0.0)
			{
				OSMOT += sit_M[i0] * sit_M[i1] * param / 2.0;
			}
			else
			{
				OSMOT += sit_M[i0] * sit_M[i1] * param;
			}
			break;
		case TYPE_SIT_EPSILON_MU:
			sit_LGAMMA[i0] += sit_M[i1] * I * param;
			sit_LGAMMA[i1] += sit_M[i0] * I * param;
			OSMOT += sit_M[i0] * sit_M[i1] * param;
			if (z0 == 0.0 && z1 == 0.0)
			{
				OSMOT += sit_M[i0] * sit_M[i1] * param * I / 2.0;
			}
			else
			{
				OSMOT += sit_M[i0] * sit_M[i1] * param * I;
			}
			break;
		default:
		case TYPE_Other:
			error_msg("TYPE_Other in pitz_param list.", STOP);
			break;
		}
	}

	/*
	 *  Add F and CSUM terms to sit_LGAMMA
	 */
	for (size_t j = 0; j < ion_list.size(); j++)
	{
		int i = ion_list[j];
		z0 = spec[i]->z;
		sit_LGAMMA[i] += z0 * z0 * F;
	}
	//for (i = 0; i < sit_count_cations; i++)
	//{
	//	z0 = spec[i]->z;
	//	sit_LGAMMA[i] += z0 * z0 * F;
	//}
	//for (i = 2 * (int)s.size(); i < 2 * (int)s.size() + sit_count_anions; i++)
	//{
	//	z0 = spec[i]->z;
	//	sit_LGAMMA[i] += z0 * z0 * F;
	//}
	/*
	   C
	   C     CONVERT TO MACINNES CONVENTION
	   C
	 */
	/*COSMOT = 1.0e0 + 2.0e0 * OSMOT / OSUM;*/
	COSMOT = 1.0e0 + OSMOT*LOG_10 / OSUM;
	/*
	   C
	   C     CALCULATE THE ACTIVITY OF WATER
	   C
	 */
	AW = exp(-OSUM * COSMOT / 55.50837e0);
	/*if (AW > 1.0) AW = 1.0;*/
	/*s_h2o->la=log10(AW); */
	mu_x = I;
	for (size_t j = 0; j < s_list.size(); j++)
	{
		int i = s_list[j];
		spec[i]->lg_pitzer = sit_LGAMMA[i];
	}
//	for (i = 0; i < 2 * (int)s.size() + sit_count_anions; i++)
//	{
//		if (sit_IPRSNT[i] == FALSE)	continue;
//		spec[i]->lg_pitzer = sit_LGAMMA[i];
///*
//		   output_msg(sformatf( "%d %s:\t%e\t%e\t%e\t%e \n", i, spec[i]->name, sit_M[i], spec[i]->la, spec[i]->lg_pitzer, spec[i]->lg));
//*/
//	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
sit_clean_up(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free all allocated memory, except strings
 */
	int i;

	for (i = 0; i < (int)sit_params.size(); i++)
	{
		delete sit_params[i]; 
	}
	sit_params.clear();
	sit_param_map.clear();
	sit_LGAMMA.clear();
	sit_IPRSNT.clear();
	spec.clear();
	//delete aphi; 
	sit_M.clear(); 

	return OK;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
set_sit(int initial)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sets initial guesses for unknowns if initial == TRUE
 *   Revises guesses whether initial is true or not
 */
	int i;
	cxxSolution *solution_ptr;
/*
 *   Set initial log concentrations to zero
 */
	iterations = -1;
	solution_ptr = use.Get_solution_ptr();
	for (i = 0; i < (int)this->s_x.size(); i++)
	{
		s_x[i]->lm = LOG_ZERO_MOLALITY;
		s_x[i]->lg_pitzer = 0.0;
	}
	if (initial == TRUE || set_and_run_attempt > 0)
	{
		for (i = 0; i < (int)this->s_x.size(); i++)
		{
			s_x[i]->lg = 0.0;
		}
	}
/*
 *   Set master species activities
 */
	tc_x = solution_ptr->Get_tc();
	tk_x = tc_x + 273.15;

	patm_x = solution_ptr->Get_patm(); // done in calc_rho_0(tc, pa)
	potV_x = solution_ptr->Get_potV();

/*
 *   H+, e-, H2O
 */
	mass_water_aq_x = solution_ptr->Get_mass_water();
	mu_x = solution_ptr->Get_mu();
	s_h2o->moles = mass_water_aq_x / gfw_water;
	s_h2o->la = log10(solution_ptr->Get_ah2o());
	AW = pow((LDBLE) 10.0E0, s_h2o->la);
	s_hplus->la = -solution_ptr->Get_ph();
	s_hplus->lm = s_hplus->la;
	s_hplus->moles = exp(s_hplus->lm * LOG_10) * mass_water_aq_x;
	s_eminus->la = -solution_ptr->Get_pe();
	if (initial == TRUE) sit_initial_guesses();
	if (dl_type_x != cxxSurface::NO_DL)	initial_surface_water();
	sit_revise_guesses();
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
sit_initial_guesses(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Make initial guesses for activities of master species and
 *   ionic strength
 */
	int i;
	cxxSolution *solution_ptr;

	solution_ptr = use.Get_solution_ptr();
	mu_x =
		s_hplus->moles +
		exp((solution_ptr->Get_ph() - 14.) * LOG_10) * mass_water_aq_x;
	mu_x /= mass_water_aq_x;
	s_h2o->la = 0.0;
	for (i = 0; i < count_unknowns; i++)
	{
		if (x[i] == ph_unknown || x[i] == pe_unknown)
			continue;
		if (x[i]->type < CB)
		{
			mu_x +=
				x[i]->moles / mass_water_aq_x * 0.5 * x[i]->master[0]->s->z *
				x[i]->master[0]->s->z;
			x[i]->master[0]->s->la = log10(x[i]->moles / mass_water_aq_x);
		}
		else if (x[i]->type == CB)
		{
			x[i]->master[0]->s->la =
				log10(0.001 * x[i]->moles / mass_water_aq_x);
		}
		else if (x[i]->type == SOLUTION_PHASE_BOUNDARY)
		{
			x[i]->master[0]->s->la =
				log10(0.001 * x[i]->moles / mass_water_aq_x);
		}
		else if (x[i]->type == EXCH)
		{
			if (x[i]->moles <= 0)
			{
				x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
			}
			else
			{
				x[i]->master[0]->s->la = log10(x[i]->moles);
			}
		}
		else if (x[i]->type == SURFACE)
		{
			if (x[i]->moles <= 0)
			{
				x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
			}
			else
			{
				x[i]->master[0]->s->la = log10(0.1 * x[i]->moles);
			}
		}
		else if (x[i]->type == SURFACE_CB)
		{
			x[i]->master[0]->s->la = 0.0;
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
sit_revise_guesses(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Revise molalities species
 */
	int i;
	int l_iter, max_iter, repeat, fail;
	LDBLE weight, f;

	max_iter = 100;
	/* gammas(mu_x); */
	l_iter = 0;
	repeat = TRUE;
	fail = FALSE;
	double d = 2;
	double logd = log10(d);
	while (repeat == TRUE && fail == FALSE)
	{
		l_iter++;
		if (debug_set == TRUE)
		{
			output_msg(sformatf( "\nBeginning set iteration %d.\n",
					   l_iter));
		}
		if (l_iter == max_iter + 1)
		{
			log_msg(sformatf(
					   "Did not converge in set, iteration %d.\n",
					   iterations));
			fail = TRUE;
		}
		if (l_iter > 2 * max_iter)
		{
			log_msg(sformatf(
					   "Did not converge with relaxed criteria in set.\n"));
			return (OK);
		}
		molalities(TRUE);
		/*pitzer(); */
		/*s_h2o->la = 0.0; */
		/*molalities(TRUE); */
		mb_sums();
		if (state < REACTION)
		{
			sum_species();
		}
		else
		{
			for (i = 0; i < count_unknowns; i++)
			{
				x[i]->sum = x[i]->f;
			}
		}
		/*n
		   if (debug_set == TRUE) {
		   pr.species = TRUE;
		   pr.all = TRUE;
		   print_species();
		   }
		 */
		repeat = FALSE;
		for (i = 0; i < count_unknowns; i++)
		{
			if (x[i] == ph_unknown || x[i] == pe_unknown)
				continue;
			if (x[i]->type == MB ||
				/*			    x[i]->type == ALK || */
					x[i]->type == CB ||
					x[i]->type == SOLUTION_PHASE_BOUNDARY ||
					x[i]->type == EXCH || x[i]->type == SURFACE)
			{

				if (debug_set == TRUE)
				{
					output_msg(sformatf(
						"\n\t%5s  at beginning of set %d: %e\t%e\t%e\n",
						x[i]->description, l_iter, (double) x[i]->sum,
						(double) x[i]->moles,
						(double) x[i]->master[0]->s->la));
				}
				if (fabs(x[i]->moles) < 1e-30)
					x[i]->moles = 0;
				f = fabs(x[i]->sum);
				if (f == 0 && x[i]->moles == 0)
				{
					x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
					continue;
				}
				else if (f == 0)
				{
					repeat = TRUE;
					x[i]->master[0]->s->la += logd;
/*!!!!*/ if (x[i]->master[0]->s->la < -999.)
						x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
				}
				else if (f > d * fabs(x[i]->moles)
					|| f < 1.0/d * fabs(x[i]->moles))
				{
					weight = (f < 1.0/d * fabs(x[i]->moles)) ? 0.3 : 1.0;
					if (x[i]->moles <= 0)
					{
						x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
					}
					else
					{
						repeat = TRUE;
						x[i]->master[0]->s->la +=
							weight * log10(fabs(x[i]->moles / x[i]->sum));
					}
					if (debug_set == TRUE)
					{
						output_msg(sformatf(
							"\t%5s not converged in set %d: %e\t%e\t%e\n",
							x[i]->description, l_iter,
							(double) x[i]->sum, (double) x[i]->moles,
							(double) x[i]->master[0]->s->la));
					}
				}
			}
			else if (x[i]->type == ALK)
			{
				f = total_co2;
				if (fail == TRUE && f < 1.5 * fabs(x[i]->moles))
				{
					continue;
				}
				if (f > 1.5 * fabs(x[i]->moles)
					|| f < 1.0/d * fabs(x[i]->moles))
				{
					repeat = TRUE;
					weight = (f < 1.0/d * fabs(x[i]->moles)) ? 0.3 : 1.0;
					x[i]->master[0]->s->la += weight *
						log10(fabs(x[i]->moles / x[i]->sum));
					if (debug_set == TRUE)
					{
						output_msg(sformatf(
								   "%s not converged in set. %e\t%e\t%e\n",
								   x[i]->description, (double) x[i]->sum,
								   (double) x[i]->moles,
								   (double) x[i]->master[0]->s->la));
					}
				}
			}
		}
	}
	log_msg(sformatf( "Iterations in sit_revise_guesses: %d\n", l_iter));
	/*mu_x = mu_unknown->f * 0.5 / mass_water_aq_x; */
	if (mu_x <= 1e-8)
	{
		mu_x = 1e-8;
	}
	/*gammas(mu_x); */
	return (OK);
}
//#define ORIGINAL
#ifdef ORIGINAL
/* ---------------------------------------------------------------------- */
int Phreeqc::
jacobian_sit(void)
/* ---------------------------------------------------------------------- */
{
	std::vector<double> base;
	LDBLE d, d1, d2;
	int i, j;
Restart:
	size_t pz_max_unknowns = max_unknowns;
	//k_temp(tc_x, patm_x);
	if (full_pitzer == TRUE)
	{
		molalities(TRUE);
		sit();
		residuals();
	}
	base = residual; // std::vectors
	d = 0.0001;
	d1 = d * LOG_10;
	d2 = 0;
	for (i = 0; i < count_unknowns; i++)
	{
		switch (x[i]->type)
		{
		case MB:
		case ALK:
		case CB:
		case SOLUTION_PHASE_BOUNDARY:
		case EXCH:
		case SURFACE:
		case SURFACE_CB:
		case SURFACE_CB1:
		case SURFACE_CB2:
			x[i]->master[0]->s->la += d;
			d2 = d1;
			break;
		case AH2O:
			x[i]->master[0]->s->la += d;
			d2 = d1;
			break;
		case PITZER_GAMMA:
			if (!full_pitzer) 
				continue;
			x[i]->s->lg += d;
			d2 = d;
			break;
		case MH2O:
			mass_water_aq_x *= (1.0 + d);
			x[i]->master[0]->s->moles = mass_water_aq_x / gfw_water;
			d2 = log(1.0 + d);
			break;
		case MH:
			s_eminus->la += d;
			d2 = d1;
			break;
			/*
			if (pitzer_pe == TRUE)
			{
				s_eminus->la += d;
				d2 = d1;
				break;
			}
			else
			{
				continue;
			}
			*/
		case GAS_MOLES:
			if (gas_in == FALSE)
				continue;
			d2 = (x[i]->moles > 1 ? 1 : 20);
			d2 *= d * x[i]->moles;
			if (d2 < 1e-14)
				d2 = 1e-14;
			x[i]->moles += d2;
			break;
		case MU:
			//continue;
			d2 = d * mu_x;
			mu_x += d2;
			//k_temp(tc_x, patm_x);
			gammas(mu_x);
			break;
		case PP:
		case SS_MOLES:
			continue;
			break;
		}
		molalities(TRUE);
		if (max_unknowns > pz_max_unknowns) 
		{
			gammas_sit();
			jacobian_sums();
			goto Restart;
		}
		if (full_pitzer == TRUE)
			sit();
		mb_sums();
		residuals();
		for (j = 0; j < count_unknowns; j++)
		{
			my_array[(size_t)j * (count_unknowns + 1) + (size_t)i] =
				-(residual[j] - base[j]) / d2;
		}
		switch (x[i]->type)
		{
		case MB:
		case ALK:
		case CB:
		case SOLUTION_PHASE_BOUNDARY:
		case EXCH:
		case SURFACE:
		case SURFACE_CB:
		case SURFACE_CB1:
		case SURFACE_CB2:
		case AH2O:
			x[i]->master[0]->s->la -= d;
			break;
		case MH:
			s_eminus->la -= d;
			if (my_array[(size_t)i * (count_unknowns + 1) + (size_t)i] == 0)
			{
				my_array[(size_t)i * (count_unknowns + 1) + (size_t)i] =
					exp(s_h2->lm * LOG_10) * 2;
			}
			break;
		case PITZER_GAMMA:
			x[i]->s->lg -= d;
			break;
		case MH2O:
			mass_water_aq_x /= (1 + d);
			x[i]->master[0]->s->moles = mass_water_aq_x / gfw_water;
			break;
		case MU:
			mu_x -= d2;
			//k_temp(tc_x, patm_x);
			gammas(mu_x);
			break;
		case GAS_MOLES:
			if (gas_in == FALSE)
				continue;
			x[i]->moles -= d2;
			break;
		}
	}
	molalities(TRUE);
	if (full_pitzer == TRUE)
		sit();
	mb_sums();
	residuals();
	return OK;
}
#else
/* ---------------------------------------------------------------------- */
int Phreeqc::
jacobian_sit(void)
/* ---------------------------------------------------------------------- */
{
	std::vector<double> base;
	LDBLE d, d1, d2;
	int i, j;
	std::vector<class phase*> phase_ptrs;
	std::vector<class phase> base_phases;
	cxxGasPhase base_gas_phase;
	cxxSurface base_surface;
Restart:
	if (use.Get_surface_ptr() != NULL)
	{
		base_surface = *use.Get_surface_ptr();
	}
	if (use.Get_gas_phase_ptr() != NULL)
	{
		cxxGasPhase* gas_phase_ptr = use.Get_gas_phase_ptr();
		base_gas_phase = *gas_phase_ptr;
		base_phases.resize(gas_phase_ptr->Get_gas_comps().size());
		for (size_t i = 0; i < gas_phase_ptr->Get_gas_comps().size(); i++)
		{
			const cxxGasComp* gas_comp_ptr = &(gas_phase_ptr->Get_gas_comps()[i]);
			class phase* phase_ptr = phase_bsearch(gas_comp_ptr->Get_phase_name().c_str(), &j, FALSE);
			phase_ptrs.push_back(phase_ptr);
			base_phases[i] = *phase_ptr;
		}
	}
	calculating_deriv = 1;
	size_t pz_max_unknowns = max_unknowns;
	//k_temp(tc_x, patm_x);
	molalities(TRUE);
	if (full_pitzer == TRUE)
	{
		
		sit();
	}
	mb_sums();
	residuals();
	base = residual; // std::vectors
	d = 0.0001;
	d1 = d * LOG_10;
	d2 = 0;
	for (i = 0; i < count_unknowns; i++)
	{
		switch (x[i]->type)
		{
		case MB:
		case ALK:
		case CB:
		case SOLUTION_PHASE_BOUNDARY:
		case EXCH:
		case SURFACE:
		case SURFACE_CB:
		case SURFACE_CB1:
		case SURFACE_CB2:
			x[i]->master[0]->s->la += d;
			d2 = d1;
			break;
		case AH2O:
			x[i]->master[0]->s->la += d;
			d2 = d1;
			break;
		case PITZER_GAMMA:
			if (!full_pitzer)
				continue;
			x[i]->s->lg += d;
			d2 = d;
			break;
		case MH2O:
			mass_water_aq_x *= (1.0 + d);
			x[i]->master[0]->s->moles = mass_water_aq_x / gfw_water;
			d2 = log(1.0 + d);
			break;
		case MH:
			s_eminus->la += d;
			d2 = d1;
			break;
			/*
			if (pitzer_pe == TRUE)
			{
				s_eminus->la += d;
				d2 = d1;
				break;
			}
			else
			{
				continue;
			}
			*/
		case GAS_MOLES:
			if (gas_in == FALSE)
				continue;
			d2 = (x[i]->moles > 1 ? 1 : 20);
			d2 *= d * x[i]->moles;
			if (d2 < 1e-14)
				d2 = 1e-14;
			x[i]->moles += d2;
			break;
		case MU:
			//continue;
			d2 = d * mu_x;
			mu_x += d2;
			//k_temp(tc_x, patm_x);
			gammas_sit();
			break;
		case PP:
		case SS_MOLES:
			continue;
			break;
		}
		molalities(TRUE);
		if (max_unknowns > pz_max_unknowns)
		{
			gammas_sit();
			jacobian_sums();
			goto Restart;
		}
		if (full_pitzer == TRUE)
			sit();
		mb_sums();
		residuals();
		for (j = 0; j < count_unknowns; j++)
		{
			my_array[(size_t)j * (count_unknowns + 1) + (size_t)i] =
				-(residual[j] - base[j]) / d2;
		}
		switch (x[i]->type)
		{
		case MB:
		case ALK:
		case CB:
		case SOLUTION_PHASE_BOUNDARY:
		case EXCH:
		case SURFACE:
		case SURFACE_CB:
		case SURFACE_CB1:
		case SURFACE_CB2:
		case AH2O:
			x[i]->master[0]->s->la -= d;
			break;
		case MH:
			s_eminus->la -= d;
			if (my_array[(size_t)i * (count_unknowns + 1) + (size_t)i] == 0)
			{
				my_array[(size_t)i * (count_unknowns + 1) + (size_t)i] =
					exp(s_h2->lm * LOG_10) * 2;
			}
			break;
		case PITZER_GAMMA:
			x[i]->s->lg -= d;
			break;
		case MH2O:
			mass_water_aq_x /= (1 + d);
			x[i]->master[0]->s->moles = mass_water_aq_x / gfw_water;
			break;
		case MU:
			mu_x -= d2;
			//k_temp(tc_x, patm_x);
			gammas_sit();
			break;
		case GAS_MOLES:
			if (gas_in == FALSE)
				continue;
			x[i]->moles -= d2;
			break;
		}
		if (use.Get_surface_ptr() != NULL)
		{
			*use.Get_surface_ptr() = base_surface;
		}
		if (use.Get_gas_phase_ptr() != NULL)
		{
			*use.Get_gas_phase_ptr() = base_gas_phase;
			for (size_t g = 0; g < base_phases.size(); g++)
			{
				*phase_ptrs[g] = base_phases[g];
			}
		}
	}
	molalities(TRUE);
	if (full_pitzer == TRUE)
		sit();
	mb_sums();
	residuals();
	//for (i = 0; i < count_unknowns; i++)
	//{
	//	//Debugging
	//	if (fabs(2.0 * (residual[i] - base[i]) / (residual[i] + base[i])) > 1e-2 &&
	//		fabs(residual[i]) + fabs(base[i]) > 1e-8)
	//	{
	//		std::cerr << i << ": " << x[i]->description << "  " << residual[i] << "  " << base[i] << std::endl;
	//	}
	//}
	calculating_deriv = 0;
	return OK;
}
#endif
/* ---------------------------------------------------------------------- */
int Phreeqc::
model_sit(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   model is called after the equations have been set up by prep
 *   and initial guesses have been made in set.
 * 
 *   Here is the outline of the calculation sequence:
 *      residuals--residuals are calculated, if small we are done
 *      sum_jacobian--jacobian is calculated 
 *      ineq--inequality solver is called
 *      reset--estimates of unknowns revised, if changes are small solution
 *         has been found, usually convergence is found in residuals.
 *      gammas--new activity coefficients
 *      molalities--calculate molalities
 *      mb_sums--calculate mass-balance sums
 *      mb_gases--decide if gas_phase exists
 *      mb_ss--decide if solid_solutions exists
 *      switch_bases--check to see if new basis species is needed
 *         reprep--rewrite equations with new basis species if needed
 *         sit_revise_guesses--revise unknowns to get initial mole balance
 *      check_residuals--check convergence one last time
 *         sum_species--calculate sums of elements from species concentrations
 *
 *      An additional pass through may be needed if unstable phases still exist
 *         in the phase assemblage. 
 */
	int l_kode, return_kode;
	int r;
	int count_infeasible, count_basis_change;
	int debug_model_save;
	int mass_water_switch_save;

/*	debug_model = TRUE; */
/*	debug_prep = TRUE; */
/*	debug_set = TRUE; */
	/* mass_water_switch == TRUE, mass of water is constant */
	mass_water_switch_save = mass_water_switch;
	if (mass_water_switch_save == FALSE && delay_mass_water == TRUE)
	{
		mass_water_switch = TRUE;
	}
	debug_model_save = debug_model;
	pe_step_size_now = pe_step_size;
	step_size_now = step_size;
#ifdef NPP
	if (!use.Get_kinetics_in()) status(0, NULL);
#else
	status(0, NULL);
#endif
	iterations = 0;
	gamma_iterations = 0;
	count_basis_change = count_infeasible = 0;
	stop_program = FALSE;
	remove_unstable_phases = FALSE;
	if (always_full_pitzer == TRUE)
	{
		full_pitzer = TRUE;
	}
	else
	{
		full_pitzer = FALSE;
	}
	//sit_make_lists();
	for (;;)
	{
		mb_gases();
		mb_ss();
		l_kode = 1;
		while ((r = residuals()) != CONVERGED
			   || remove_unstable_phases == TRUE)
		{
#if defined(PHREEQCI_GUI)
			PhreeqcIWait(this);
#endif
			iterations++;
			overall_iterations++;
			if (iterations > itmax - 1 && debug_model == FALSE
				&& pr.logfile == TRUE)
			{
				set_forward_output_to_log(TRUE);
				debug_model = TRUE;
			}
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "\nIteration %d\tStep_size = %f\n", iterations,
						   (double) step_size_now));
				output_msg(sformatf( "\t\tPe_step_size = %f\n\n",
						   (double) pe_step_size_now));
			}
			/*
			 *   Iterations exceeded
			 */
			if (iterations > itmax)
			{
				error_string = sformatf( "Maximum iterations exceeded, %d\n",
						itmax);
				warning_msg(error_string);
				stop_program = TRUE;
				break;
			}
			/*
			 *   Calculate jacobian
			 */
			gammas_sit();
			jacobian_sums();
			jacobian_sit();
			/*
			 *   Full matrix with pure phases
			 */
			if (r == OK || remove_unstable_phases == TRUE)
			{
				return_kode = ineq(l_kode);
				if (return_kode != OK)
				{
					if (debug_model == TRUE)
					{
						output_msg(sformatf(
								   "Ineq had infeasible solution, "
								   "kode %d, iteration %d\n", return_kode,
								   iterations));
					}
					log_msg(sformatf( "Ineq had infeasible solution, "
							   "kode %d, iteration %d\n", return_kode,
							   iterations));
					count_infeasible++;
				}
				if (return_kode == 2)
				{
					ineq(0);
				}
				reset();
			}
			gammas_sit();
			if (full_pitzer == TRUE)
				sit();
			if (always_full_pitzer == TRUE)
			{
				full_pitzer = TRUE;
			}
			else
			{
				full_pitzer = FALSE;
			}
			if (molalities(FALSE) == ERROR)
			{
				sit_revise_guesses();
			}
			if (use.Get_surface_ptr() != NULL &&
				use.Get_surface_ptr()->Get_dl_type() != cxxSurface::NO_DL &&
				use.Get_surface_ptr()->Get_related_phases() == TRUE)
				initial_surface_water();
			mb_sums();
			mb_gases();
			mb_ss();
/*
 *   Switch bases if necessary
 */
			if (switch_bases() == TRUE)
			{
				
				count_basis_change++;
				//count_unknowns -= (int)this->s_x.size();
				count_unknowns -= sit_aqueous_unknowns;
				reprep();
				full_pitzer = false;
			}
			/* debug
			   species_list_sort();
			   sum_species();
			   print_species();
			   print_exchange();
			   print_surface();
			 */
			if (stop_program == TRUE)
			{
				break;
			}
		}
/*
 *   Check for stop_program
 */

		if (stop_program == TRUE)
		{
			break;
		}
		if (check_residuals() == ERROR)
		{
			stop_program = TRUE;
			break;
		}
		/* remove_unstable_phases is set in check_residuals */
		if (remove_unstable_phases == FALSE && mass_water_switch_save == FALSE
			&& mass_water_switch == TRUE)
		{
			log_msg(sformatf(
					   "\nChanging water switch to FALSE. Iteration %d.\n",
					   iterations));
			mass_water_switch = FALSE;
			continue;
		}
		gamma_iterations++;
		if (gamma_iterations > itmax)
		{
			error_string = sformatf( "Maximum gamma iterations exceeded, %d\n",
					itmax);
			warning_msg(error_string);
			stop_program = TRUE;
			break;
		}
		if (check_gammas_sit() != TRUE)
		{
			full_pitzer = TRUE;
			continue;
		}
		if (remove_unstable_phases == FALSE)
			break;
		if (debug_model == TRUE)
		{
			output_msg(sformatf(
					   "\nRemoving unstable phases. Iteration %d.\n",
					   iterations));
		}
		log_msg(sformatf( "\nRemoving unstable phases. Iteration %d.\n",
				   iterations));
	}
	log_msg(sformatf( "\nNumber of infeasible solutions: %d\n",
			   count_infeasible));
	log_msg(sformatf( "Number of basis changes: %d\n\n",
			   count_basis_change));
	log_msg(sformatf( "Number of iterations: %d\n", iterations));
	log_msg(sformatf( "Number of gamma iterations: %d\n\n", gamma_iterations));
	debug_model = debug_model_save;
	set_forward_output_to_log(FALSE);
	if (stop_program == TRUE)
	{
		return (ERROR);
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
check_gammas_sit(void)
/* ---------------------------------------------------------------------- */
{
	LDBLE old_mu, tol;
	int converge, i;

	old_mu = mu_x;
	sit();
	molalities(TRUE);
	mb_sums();
	converge = TRUE;
	tol = convergence_tolerance * 10.;
	for (i = 0; i < count_unknowns; i++)
	{
		if (x[i]->type != PITZER_GAMMA)
			continue;
		if (fabs(x[i]->s->lg - x[i]->s->lg_pitzer) > tol)
		{
			converge = FALSE;
		}
	}
	if (fabs(old_mu - mu_x) > tol)
	{
		converge = FALSE;
	}
	if ((pow((LDBLE) 10.0, s_h2o->la) - AW) > tol)
	{
		converge = FALSE;
	}
	return converge;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
gammas_sit()
/* ---------------------------------------------------------------------- */
{
/*
 *   Need exchange gammas for pitzer
 */
	int i, j;
	LDBLE coef;
	/* Initialize */
	k_temp(tc_x, patm_x);
/*
 *   Calculate activity coefficients
 */
	for (i = 0; i < (int)this->s_x.size(); i++)
	{
		switch (s_x[i]->gflag)
		{
		case 0:				/* uncharged */
		case 1:				/* Davies */
		case 2:				/* Extended D-H, WATEQ D-H */
		case 3:				/* Always 1.0 */
			break;
		case 4:				/* Exchange */
			/* Now calculated in next loop */
			break;
		case 5:				/* Always 1.0 */
			break;
		case 6:				/* Surface */
/*
 *   Find moles of sites. 
 *   s_x[i]->equiv is stoichiometric coefficient of sites in species
 */
			for (j = 1; s_x[i]->rxn_x.token[j].s != NULL; j++)
			{
				if (s_x[i]->rxn_x.token[j].s->type == SURF)
				{
					s_x[i]->alk =
						s_x[i]->rxn_x.token[j].s->primary->unknown->moles;
					break;
				}
			}
			if (s_x[i]->alk > 0)
			{
				s_x[i]->lg = log10(s_x[i]->equiv / s_x[i]->alk);
				s_x[i]->dg = 0.0;
			}
			else
			{
				s_x[i]->lg = 0.0;
				s_x[i]->dg = 0.0;
			}
			break;
		case 7:				/* LLNL */
			break;
		case 8:				/* LLNL CO2 */
			break;
		case 9:				/* activity water */
			s_x[i]->lg = log10(exp(s_h2o->la * LOG_10) * gfw_water);
			s_x[i]->dg = 0.0;
			break;
		}
/*
		if (mu_unknown != NULL) {
			if (fabs(residual[mu_unknown->number]) > 0.1 &&
			    fabs(residual[mu_unknown->number])/mu_x > 0.5) {
				s_x[i]->dg = 0.0;
			}
		}
 */
	}
	/*
	 *  calculate exchange gammas 
	 */

	if (use.Get_exchange_ptr() != NULL)
	{
		for (i = 0; i < (int)this->s_x.size(); i++)
		{
			switch (s_x[i]->gflag)
			{
			case 0:			/* uncharged */
			case 1:			/* Davies */
			case 2:			/* Extended D-H, WATEQ D-H */
			case 3:			/* Always 1.0 */
			case 5:			/* Always 1.0 */
			case 6:			/* Surface */
			case 7:			/* LLNL */
			case 8:			/* LLNL CO2 */
			case 9:			/* activity water */
				break;
			case 4:			/* Exchange */

				/*
				 *   Find CEC
				 *   z contains valence of cation for exchange species, alk contains cec
				 */
				/* !!!!! */
				for (j = 1; s_x[i]->rxn_x.token[j].s != NULL; j++)
				{
					if (s_x[i]->rxn_x.token[j].s->type == EX)
					{
						s_x[i]->alk =
							s_x[i]->rxn_x.token[j].s->primary->unknown->
							moles;
						break;
					}
				}
				/*
				 *   Master species is a dummy variable with meaningless activity and mass
				 */
				s_x[i]->lg = 0.0;
				s_x[i]->dg = 0.0;
				if (s_x[i]->primary != NULL)
				{
					break;
				}
				/*
				 *   All other species
				 */

				/* modific 29 july 2005... */
				if (s_x[i]->equiv != 0 && s_x[i]->alk > 0)
				{
					s_x[i]->lg = log10(fabs(s_x[i]->equiv) / s_x[i]->alk);
				}
				if (use.Get_exchange_ptr()->Get_pitzer_exchange_gammas())
				{
					/* Assume equal gamma's of solute and exchangeable species...  */
					for (j = 1; s_x[i]->rxn_x.token[j].s != NULL; j++)
					{
						if (s_x[i]->rxn_x.token[j].s->type == EX)
							continue;
						coef = s_x[i]->rxn_x.token[j].coef;
						s_x[i]->lg += coef * s_x[i]->rxn_x.token[j].s->lg;
						s_x[i]->dg += coef * s_x[i]->rxn_x.token[j].s->dg;
					}
				}
			}
		}
	}
/* ...end modific 29 july 2005 */

	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
PTEMP_SIT(LDBLE TK)
/* ---------------------------------------------------------------------- */
{
/*
C
C     SUBROUTINE TO CALCULATE TEMPERATURE DEPENDENCE OF PITZER PARAMETER
C
*/
	LDBLE TR = 298.15;

	if (fabs(TK - OTEMP) < 0.001 && fabs(patm_x - OPRESS) < 0.1)	return OK; 
/*
C     Set DW0
*/
	DW0 = rho_0 = calc_rho_0(TK - 273.15, patm_x);
	VP = patm_x;
	for (size_t j = 0; j < param_list.size(); j++)
	{
		int i = param_list[j];
		calc_sit_param(sit_params[i], TK, TR);
	}
	calc_dielectrics(TK - 273.15, patm_x);
	sit_A0 = A0;
	OTEMP = TK;
	OPRESS = patm_x;
	return OK;
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
sit_make_lists(void)
/* ---------------------------------------------------------------------- */
{
	double log_min = log10(MIN_TOTAL);
	s_list.clear();
	cation_list.clear();
	neutral_list.clear();
	anion_list.clear();
	ion_list.clear();
	param_list.clear();
	OTEMP = -100.0;	
	for (int j = 0; j < 3; j++)
	{
		int min, max;
		switch (j)
		{
		case 0:
			min = 0;
			max = sit_count_cations;
			break;
		case 1:
			min = (int)s.size();
			max = (int)s.size() + sit_count_neutrals;
			break;
		case 2:
			min = 2*(int)s.size();
			max = 2*(int)s.size() + sit_count_anions;
			break;
		}
		for (int i = min; i < max; i++)
		{
			sit_IPRSNT[i] = FALSE;
			sit_M[i] = 0.0;
			if (spec[i] != NULL && spec[i]->in == TRUE)
			{
				if (spec[i]->type == EX ||
					spec[i]->type == SURF || spec[i]->type == SURF_PSI)
					continue;	
				sit_IPRSNT[i] = TRUE;	
				s_list.push_back(i);	
				if (i < (int)s.size())
				{
					cation_list.push_back(i);
				}
				if (i >= (int)s.size() && i < 2*(int)s.size())
				{
					neutral_list.push_back(i);
				}
				if (i >= 2*(int)s.size())
				{
					anion_list.push_back(i);
				}
				if (i < (int)s.size() || i >= 2*(int)s.size())
				{
					ion_list.push_back(i);
				}
				if (spec[i]->lm > log_min)
				{
					sit_M[i] = under(spec[i]->lm);
				}
			}
		}
	}
	for (int i = 0; i < (int)sit_params.size(); i++)
	{
		int i0 = sit_params[i]->ispec[0];
		int i1 = sit_params[i]->ispec[1];
		if (sit_IPRSNT[i0] == FALSE || sit_IPRSNT[i1] == FALSE) continue;
		param_list.push_back(i);
	}
}
