#include "Phreeqc.h"
#include <algorithm>			// std::replace

#include "NameDouble.h"
#include "Solution.h"
#include "Reaction.h"
#include "PPassemblage.h"
#include "Exchange.h"
#include "Surface.h"
#include "GasPhase.h"
#include "SSassemblage.h"
#include "cxxKinetics.h"
#include "phqalloc.h"
#include "PBasic.h"
#include "Temperature.h"
#include "SSassemblage.h"
#include "Utils.h"

#if defined(PHREEQCI_GUI)
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#endif

//const const_iso Phreeqc::iso_defaults[] = {
//	{"13C", -10, 1},
//	{"13C(4)", -10, 1},
//	{"13C(-4)", -50, 5},
//	{"34S", 10, 1},
//	{"34S(6)", 10, 1},
//	{"34S(-2)", -30, 5},
//	{"2H", -28, 1},
//	{"2H(1)", -28, 1},
//	{"2H(0)", -28, 1},
//	{"18O", -5, .1},
//	{"18O(-2)", -5, .1},
//	{"18O(0)", -5, .1},
//	{"87Sr", .71, .01},
//	{"11B", 20, 5}
//};
const const_iso Phreeqc::iso_defaults[] = {
	const_iso("13C", -10, 1),
	const_iso("13C(4)", -10, 1),
	const_iso("13C(-4)", -50, 5),
	const_iso("34S", 10, 1),
	const_iso("34S(6)", 10, 1),
	const_iso("34S(-2)", -30, 5),
	const_iso("2H", -28, 1),
	const_iso("2H(1)", -28, 1),
	const_iso("2H(0)", -28, 1),
	const_iso("18O", -5, .1),
	const_iso("18O(-2)", -5, .1),
	const_iso("18O(0)", -5, .1),
	const_iso("87Sr", .71, .01),
	const_iso("11B", 20, 5)
};

const int Phreeqc::count_iso_defaults = (sizeof(iso_defaults) / sizeof(class const_iso));

Phreeqc::~Phreeqc(void)
{

	clean_up();
	
	PHRQ_free_all();
	if (phrq_io == &ioInstance)
	{
		this->phrq_io->clear_istream();
		this->phrq_io->close_ostreams();
	}
}

void Phreeqc::set_phast(int tf)
{
	this->phast = tf;
}
size_t Phreeqc::list_components(std::list<std::string> &list_c)
/*
 *	 Find all elements in any class definition
 */
{
	cxxNameDouble accumulator;
	//accumulator.add("H", 1);
	//accumulator.add("O", 1);

	// solutions
	{
		std::map<int, cxxSolution>::const_iterator cit = Rxn_solution_map.begin();
		for (; cit !=  Rxn_solution_map.end(); cit++)
		{
			cxxSolution entity(cit->second);
			accumulator.add_extensive(entity.Get_totals(), 1.0);
		}
	}

	// irreversible reactions
	{
		std::map<int, cxxReaction>::const_iterator cit = Rxn_reaction_map.begin();
		for (; cit !=  Rxn_reaction_map.end(); cit++)
		{
			cxxReaction r_ptr(cit->second);
			reaction_calc(&r_ptr);
			accumulator.add_extensive(r_ptr.Get_elementList(), 1.0);
		}
	}

	// pure phases
	{
		std::map<int, cxxPPassemblage>::const_iterator cit = Rxn_pp_assemblage_map.begin();
		for (; cit !=  Rxn_pp_assemblage_map.end(); cit++)
		{
			cxxPPassemblage entity = cit->second;
			entity.totalize(this);
			accumulator.add_extensive(entity.Get_eltList(), 1.0);
		}
	}
	// exchangers
	{
		std::map<int, cxxExchange>::const_iterator cit = Rxn_exchange_map.begin();
		for (; cit !=  Rxn_exchange_map.end(); cit++)
		{
			cxxExchange entity = cit->second;
			entity.totalize();
			accumulator.add_extensive(entity.Get_totals(), 1.0);
		}
	}

	// surfaces
	{
		std::map<int, cxxSurface>::const_iterator cit = Rxn_surface_map.begin();
		for (; cit !=  Rxn_surface_map.end(); cit++)
		{
			cxxSurface entity = cit->second;
			entity.totalize();
			accumulator.add_extensive(entity.Get_totals(), 1.0);
		}
	}
	// gas phases
	{
		std::map<int, cxxGasPhase>::const_iterator cit = Rxn_gas_phase_map.begin();
		for (; cit !=  Rxn_gas_phase_map.end(); cit++)
		{
			cxxGasPhase entity = cit->second;
			entity.totalize(this);
			accumulator.add_extensive(entity.Get_totals(), 1.0);
		}
	}

	// solid-solutions
	{
		std::map<int, cxxSSassemblage>::const_iterator cit = Rxn_ss_assemblage_map.begin();
		for (; cit !=  Rxn_ss_assemblage_map.end(); cit++)
		{
			cxxSSassemblage entity = cit->second;
			entity.totalize(this);
			accumulator.add_extensive(entity.Get_totals(), 1.0);
		}
	}
	// kinetics
	{
		std::map<int, cxxKinetics>::iterator it = Rxn_kinetics_map.begin();
		for (; it !=  Rxn_kinetics_map.end(); it++)
		{
			calc_dummy_kinetic_reaction_tally(&(it->second));
			cxxKinetics entity = it->second;
			accumulator.add_extensive(entity.Get_totals(), 1.0);
		}
	}
	// Put in all primaries
	cxxNameDouble::iterator it;
	for (it = accumulator.begin(); it != accumulator.end(); it++)
	{
		if (it->first == "Charge") continue;
		char string[MAX_LENGTH];
		Utilities::strcpy_safe(string, MAX_LENGTH, it->first.c_str());
		class master *master_ptr = master_bsearch_primary(string);
		if (master_ptr == NULL) continue;
		if (master_ptr->type != AQ) continue;
		accumulator.add(master_ptr->elt->name, 1);
	}
	// print list
	for (it = accumulator.begin(); it != accumulator.end(); it++)
	{
		class master *master_ptr = master_bsearch(it->first.c_str());
		if (master_ptr == NULL) continue;
		if (master_ptr->type != AQ) continue;
		if (master_ptr->primary == 0) continue;
		if (it->first == "Charge") continue;
		if (it->first == "O") continue;
		if (it->first == "H") continue;
		list_c.push_back(it->first);
	}
	return(list_c.size());
}
size_t Phreeqc::list_EquilibriumPhases(std::list<std::string> &list_pp)
/*
*	 Find all elements in any class definition
*/
{
	std::set<std::string> accumulator;
	// pure phases
	{
		std::map<int, cxxPPassemblage>::const_iterator cit = Rxn_pp_assemblage_map.begin();
		for (; cit != Rxn_pp_assemblage_map.end(); cit++)
		{
			cxxPPassemblage entity = cit->second;
			std::set<std::string> pp = entity.GetPhases(this);
			std::set<std::string>::iterator ppit = pp.begin();
			for (; ppit != pp.end(); ppit++)
			{ 
				accumulator.insert(*ppit);
			}
		}
	}
	list_pp.clear();
	std::set<std::string>::iterator it = accumulator.begin();
	for (; it != accumulator.end(); it++)
	{
		list_pp.insert(list_pp.end(),*it);
	}
	return(list_pp.size());
}
size_t Phreeqc::list_GasComponents(std::list<std::string> &list_gc)
/*
*	 Find all elements in any class definition
*/
{
	std::set<std::string> accumulator;
	// pure phases
	{
		std::map<int, cxxGasPhase>::const_iterator cit = Rxn_gas_phase_map.begin();
		for (; cit != Rxn_gas_phase_map.end(); cit++)
		{
			cxxGasPhase entity = cit->second;
			std::vector<cxxGasComp> &gc = entity.Get_gas_comps();
			for (size_t i = 0; i < gc.size(); i++)
			{
				int j;
				phase * p = phase_bsearch(gc[i].Get_phase_name().c_str(), &j, 0);
				accumulator.insert(p->name);
			}
		}
	}
	list_gc.clear();
	std::set<std::string>::iterator it = accumulator.begin();
	for (; it != accumulator.end(); it++)
	{
		list_gc.insert(list_gc.end(), *it);
	}
	return(list_gc.size());
}
size_t Phreeqc::list_KineticReactions(std::list<std::string> &list_kr)
/*
*	 Find all kinetic reactions
*/
{
	std::set<std::string> accumulator;
	// Kinetics
	{
		std::map<int, cxxKinetics>::const_iterator cit = Rxn_kinetics_map.begin();
		for (; cit != Rxn_kinetics_map.end(); cit++)
		{
			cxxKinetics entity = cit->second;
			for (size_t i = 0; i < entity.Get_kinetics_comps().size(); i++)
			{
				std::string ratename = entity.Get_kinetics_comps()[i].Get_rate_name();
				int j;
				rate *r = rate_search(ratename.c_str(), &j);
				if (r != NULL)
				{
					accumulator.insert(r->name);
				}
			}
		}
	}
	list_kr.clear();
	std::set<std::string>::iterator it = accumulator.begin();
	for (; it != accumulator.end(); it++)
	{
		list_kr.insert(list_kr.end(), *it);
	}
	return(list_kr.size());
}
size_t Phreeqc::list_SolidSolutions(std::list<std::string> &list_comps, std::list<std::string> &list_names)
/*
*	 Find all elements in any class definition
*/
{
	std::vector< std::set<std::string> > ss_sets;
	std::vector<std::string> ss_names;
	// solid solutions
	std::map<int, cxxSSassemblage>::const_iterator cit = Rxn_ss_assemblage_map.begin();
	// Fill vectors, ss names and related set of component names
	for (; cit != Rxn_ss_assemblage_map.end(); cit++)
	{
		cxxSSassemblage entity = cit->second;
		std::map<std::string, cxxSS> &SSs = entity.Get_SSs();
		std::map<std::string, cxxSS>::iterator ssit = SSs.begin();
		for (; ssit != SSs.end(); ssit++)
		{
			std::string ssname = ssit->second.Get_name();
			std::set<std::string> accumulator_phases;
			for (size_t i = 0; i < ssit->second.Get_ss_comps().size(); i++)
			{
				std::string pname = ssit->second.Get_ss_comps()[i].Get_name();
				int j;
				phase * p = phase_bsearch(pname.c_str(), &j, 0);
				accumulator_phases.insert(p->name);
			}
			ss_names.push_back(ssname);
			ss_sets.push_back(accumulator_phases);
		}
	}
	// need to merge into exclusive sets of solid solution components
	bool repeat = true;
	while (repeat)
	{
		repeat = false;
		for (int i = 0; i < (int) ss_sets.size() - 1; i++)
		{
			for (int j = i + 1; j < (int) ss_sets.size(); j++)
			{
				// locate any common component
				std::set<std::string>::iterator it = ss_sets[j].begin();
				for (; it != ss_sets[j].end(); it++)
				{
					if (ss_sets[i].find(*it) != ss_sets[i].end())
					{
						repeat = true;
						break;
					}
				}
				// merge sets and clear second set
				if (repeat)
				{
					for (it = ss_sets[j].begin(); it != ss_sets[j].end(); it++)
					{
						ss_sets[i].insert(*it);
					}
					ss_sets[j].clear();
					break;
				}
			}
			if (repeat) break;
		}
	}
	list_comps.clear();
	list_names.clear();
	// Write lists
	for (size_t i = 0; i < ss_sets.size(); i++)
	{
		std::set<std::string>::iterator it = ss_sets[i].begin();
		for (; it != ss_sets[i].end(); it++)
		{
			list_names.push_back(ss_names[i]);
			list_comps.push_back(*it);
		}
	}
	return(list_comps.size());
}
size_t Phreeqc::list_Surfaces(std::list<std::string> &list_surftype, std::list<std::string> &list_surfname)
/*
*	 Find all surface types and surfaces
*/
{
	std::set<std::pair<std::string,std::string> > accumulator;
	// Surfaces
	{
		std::map<int, cxxSurface>::const_iterator cit = Rxn_surface_map.begin();
		for (; cit != Rxn_surface_map.end(); cit++)
		{
			cxxSurface entity = cit->second;
			std::vector<cxxSurfaceComp> &scomps = entity.Get_surface_comps();
			//std::vector<cxxSurfaceCharge> &scharges = entity.Get_surface_charges();
			for (size_t i = 0; i < scomps.size(); i++)
			{
				std::pair<std::string, std::string> p(scomps[i].Get_master_element(), scomps[i].Get_charge_name());
				accumulator.insert(p);
			}
		}
	}
	list_surftype.clear();
	list_surfname.clear();
	std::set<std::pair<std::string, std::string> >::iterator it = accumulator.begin();
	for (; it != accumulator.end(); it++)
	{
		list_surftype.push_back(it->first);
		list_surfname.push_back(it->second);
	}
	return(list_surfname.size());
}
size_t Phreeqc::list_Exchangers(std::list<std::string> &list_exname)
/*
*	 Find all exchangers
*/
{
	std::set<std::string> accumulator;
	// Exchangers
	std::map<int, cxxExchange>::const_iterator cit = Rxn_exchange_map.begin();
	for (; cit != Rxn_exchange_map.end(); cit++)
	{
		cxxExchange entity = cit->second;
		std::vector<cxxExchComp> &ecomps = entity.Get_exchange_comps();
		for (size_t i = 0; i < ecomps.size(); i++)
		{
			std::string exname = "";
			cxxNameDouble nd = ecomps[i].Get_totals();
			cxxNameDouble::iterator it = nd.begin();
			for (; it != nd.end(); it++)
			{
				class master *m = master_bsearch(it->first.c_str());
				if (m != NULL)
				{
					if (m->type == EX)
					{
						exname = it->first;
						break;
					}
				}
			}
			if (exname != "")
			{
				accumulator.insert(exname);
			}
		}
	}
	list_exname.clear();
	std::set< std::string>::iterator it = accumulator.begin();
	for (; it != accumulator.end(); it++)
	{
		list_exname.push_back(*it);
	}
	return(list_exname.size());
}
Phreeqc::Phreeqc(PHRQ_io *io)
{
	user_print = NULL;
	sformatf_buffer = NULL;
	basic_interpreter = NULL;
	count_elts = 0;
	aphi = NULL;
	// phrq_io
	if (io)
	{
		this->phrq_io = io;
	}
	else
	{
		this->phrq_io = &this->ioInstance;
	}
	// auto PHRQ_io ioInstance;

	// initialize data members
	init();

#if defined(SWIG) || defined(SWIG_IPHREEQC)
	basicCallback = NULL;
#endif
}
void Phreeqc::init(void)
{
	same_model                      = FALSE;
	current_tc                      = NAN;
	current_pa                      = NAN;
	current_mu                      = NAN;
	mu_terms_in_logk                = true;
	current_A                       = 0.0;
	current_x                       = 0.0;
	fix_current                     = 0.0;
	/* ----------------------------------------------------------------------
	*   STRUCTURES
	* ---------------------------------------------------------------------- */
/*
 *	 last model
 */
	last_model.force_prep           = true;
	last_model.gas_phase_type       = cxxGasPhase::GP_UNKNOWN;
	last_model.gas_phase.clear();
	last_model.ss_assemblage.clear();
	last_model.pp_assemblage.clear();
	last_model.add_formula.clear();
	last_model.si.clear();
	last_model.dl_type              = cxxSurface::NO_DL;
	last_model.surface_type         = cxxSurface::UNKNOWN_DL;

	current_selected_output         = NULL;
	current_user_punch              = NULL;
	high_precision                  = false;
	MIN_LM = -30.0;			    /* minimum log molality allowed before molality set to zero */
	LOG_ZERO_MOLALITY = -30;	/* molalities <= LOG_ZERO_MOLALITY are considered equal to zero */
	MIN_RELATED_LOG_ACTIVITY = -30;
	MIN_TOTAL = 1e-25;
	MIN_TOTAL_SS = MIN_TOTAL/100;
	MIN_RELATED_SURFACE = MIN_TOTAL*100;
	// auto Rxn_temperature_map;
	// auto Rxn_pressure_map;

	/* ----------------------------------------------------------------------
	*   Surface
	* --------------------------------------------------------------------- */
	g_iterations               = -1;
	G_TOL                      = 1e-8;
	// auto Rxn_surface_map;
	// auto charge_group_map;
	change_surf_count          = 0;
	change_surf                = NULL;
	/* ----------------------------------------------------------------------
	*   Exchange
	* ---------------------------------------------------------------------- */
	// auto Rxn_exchange_map;

	/* ----------------------------------------------------------------------
	*   Kinetics
	* ---------------------------------------------------------------------- */
	// auto Rxn_kinetics_map;

	/*----------------------------------------------------------------------
	*   Save
	*---------------------------------------------------------------------- */
	save_init(-1);             // set initial save values

	/*----------------------------------------------------------------------
	*   Inverse
	*---------------------------------------------------------------------- */
	count_inverse			= 0;
	/*----------------------------------------------------------------------
	*   Mix
	*---------------------------------------------------------------------- */
	// auto Rxn_mix_map;
	// auto Dispersion_mix_map;
	// auto Rxn_solution_mix_map;
	// auto Rxn_exchange_mix_map;
	// auto Rxn_gas_phase_mix_map;
	// auto Rxn_kinetics_mix_map;
	// auto Rxn_pp_assemblage_mix_map;
	// auto Rxn_ss_assemblage_mix_map;
	// auto Rxn_surface_mix_map;
	/*----------------------------------------------------------------------
	*   Irreversible reaction
	*---------------------------------------------------------------------- */
	run_cells_one_step = false;
	// auto Rxn_reaction_map;
	/*----------------------------------------------------------------------
	*   Gas phase
	*---------------------------------------------------------------------- */
	// auto Rxn_gas_phase_map;
	/*----------------------------------------------------------------------
	*   Solid solution
	*---------------------------------------------------------------------- */
	// auto Rxn_ss_assemblage_map;
	/*----------------------------------------------------------------------
	*   Pure-phase assemblage
	*---------------------------------------------------------------------- */
	// auto Rxn_pp_assemblage_map;
	/*----------------------------------------------------------------------
	*   Species_list
	*---------------------------------------------------------------------- */
	/*----------------------------------------------------------------------
	*   Jacobian and Mass balance lists
	*---------------------------------------------------------------------- */

	/*----------------------------------------------------------------------
	*   Solution
	*---------------------------------------------------------------------- */
	// auto Rxn_solution_map;
	// auto unnumbered_solutions;
	save_species = false;
	/*----------------------------------------------------------------------
	*   Global solution
	*---------------------------------------------------------------------- */
	new_x                   = FALSE;
	tc_x                    = 0;
	tk_x                    = 0;
	patm_x                  = 1;
	last_patm_x             = 1;
	potV_x                  = 0;
	numerical_fixed_volume  = false;
	force_numerical_fixed_volume = false;
	//switch_numerical        = false;
	ph_x                    = 0;
	solution_pe_x           = 0;
	mu_x                    = 0;
	ah2o_x                  = 1.0;
	total_h_x               = 0;
	total_o_x               = 0;
	cb_x                    = 0;
	total_ions_x            = 0;
	mass_water_aq_x         = 0;
	mass_water_surfaces_x   = 0;
	mass_water_bulk_x       = 0;
	// auto pe_x
	// auto isotopes_x
	// auto default_pe_x
	dl_type_x                = cxxSurface::NO_DL;
	total_carbon             = 0;
	total_co2                = 0;
	total_alkalinity         = 0;
	gfw_water                = 0;
	step_x                   = 0;
	kin_time_x               = 0;
	/*----------------------------------------------------------------------
	*   Transport data
	*---------------------------------------------------------------------- */
	count_cells              = 1;
	count_shifts             = 1;
	ishift                   = 1;
	bcon_first = bcon_last   = 3;
	correct_disp             = FALSE;
	tempr                    = 2.0;
	timest                   = 0.0;
	simul_tr                 = 0;
	diffc                    = 0.3e-9;
	heat_diffc               = -0.1;
	cell                     = 0;
	mcd_substeps             = 1.0;
	print_modulus            = 1;
	punch_modulus            = 1;
	dump_in                  = FALSE;
	dump_modulus             = 0;
	transport_warnings       = TRUE;
	old_cells                = 0;
	max_cells                = 0;
	all_cells                = 0;
	multi_Dflag              = FALSE;
	interlayer_Dflag         = FALSE;
	implicit                 = FALSE;
	max_mixf                 = 1.0;
	min_dif_LM               = -30.0;
	default_Dw               = 0;
	correct_Dw               = 0;
	multi_Dpor               = 0;
	interlayer_Dpor          = 0.1;
	multi_Dpor_lim           = 0;
	interlayer_Dpor_lim      = 0;
	multi_Dn                 = 0;
	interlayer_tortf         = 100.0;
	cell_no                  = 0;
	fix_current              = 0.0;
	/*----------------------------------------------------------------------
	*   Advection data
	*---------------------------------------------------------------------- */
	count_ad_cells           = 1;
	count_ad_shifts          = 1;
	print_ad_modulus         = 1;
	punch_ad_modulus         = 1;
	advection_kin_time       = 0.0;
	advection_kin_time_defined = FALSE;
	advection_warnings       = TRUE;
	/*----------------------------------------------------------------------
	*   Tidy data
	*---------------------------------------------------------------------- */
	new_model                = TRUE;
	new_exchange             = FALSE;
	new_pp_assemblage        = FALSE;
	new_surface              = FALSE;
	new_reaction             = FALSE;
	new_temperature          = FALSE;
	new_mix                  = FALSE;
	new_solution             = FALSE;
	new_gas_phase            = FALSE;
	new_inverse              = FALSE;
	new_punch                = FALSE;
	new_ss_assemblage        = FALSE;
	new_kinetics             = FALSE;
	new_copy                 = FALSE;
	new_pitzer               = FALSE;
	/*----------------------------------------------------------------------
	*   Elements
	*---------------------------------------------------------------------- */
	element_h_one            = NULL;
	/*----------------------------------------------------------------------
	*   Element List
	*---------------------------------------------------------------------- */
	count_elts               = 0;
	/*----------------------------------------------------------------------
	*   Species
	*---------------------------------------------------------------------- */
	s_h2o					= NULL;
	s_hplus					= NULL;
	s_h3oplus				= NULL;
	s_eminus				= NULL;
	s_co3					= NULL;
	s_h2					= NULL;
	s_o2					= NULL;
	/*----------------------------------------------------------------------
	*   Phases
	*---------------------------------------------------------------------- */

	/*----------------------------------------------------------------------
	*   Master species
	*---------------------------------------------------------------------- */

	/*----------------------------------------------------------------------
	*   Unknowns
	*---------------------------------------------------------------------- */
	count_unknowns          = 0;
	max_unknowns            = 0;
	ah2o_unknown            = NULL;
	alkalinity_unknown      = NULL;
	carbon_unknown          = NULL;
	charge_balance_unknown  = NULL;
	exchange_unknown        = NULL;
	mass_hydrogen_unknown   = NULL;
	mass_oxygen_unknown     = NULL;
	mb_unknown              = NULL;
	mu_unknown              = NULL;
	pe_unknown              = NULL;
	ph_unknown              = NULL;
	pure_phase_unknown      = NULL;
	solution_phase_boundary_unknown = NULL;
	surface_unknown         = NULL;
	gas_unknown             = NULL;
	ss_unknown              = NULL;
	// auto gas_unknowns;
	/*----------------------------------------------------------------------
	*   Reaction work space
	*---------------------------------------------------------------------- */
	// struct trxn;	
	for (int i = 0; i < MAX_LOG_K_INDICES; i++)
	{
		trxn.logk[i] = 0;
	}
	for (int i = 0; i < 3; i++)
	{
		trxn.dz[i] = 0;
	}
	count_trxn              = 0;

	/* ----------------------------------------------------------------------
	*   Print
	* ---------------------------------------------------------------------- */
	pr.all                  = TRUE;
	pr.initial_solutions    = TRUE;
	pr.initial_exchangers   = TRUE;
	pr.reactions            = TRUE;
	pr.gas_phase            = TRUE;
	pr.ss_assemblage        = TRUE;
	pr.pp_assemblage        = TRUE;
	pr.surface              = TRUE;
	pr.exchange             = TRUE;
	pr.kinetics             = TRUE;
	pr.totals               = TRUE;
	pr.eh                   = TRUE;
	pr.species              = TRUE;
	pr.saturation_indices   = TRUE;
	pr.irrev                = TRUE;
	pr.mix                  = TRUE;
	pr.reaction             = TRUE;
	pr.use                  = TRUE;
	pr.logfile              = FALSE;
	pr.punch                = TRUE;
	pr.status               = TRUE;
	pr.inverse              = TRUE;
	pr.dump                 = TRUE;
	pr.user_print           = TRUE;
	pr.headings             = TRUE;
	pr.user_graph           = TRUE;
	pr.echo_input           = TRUE;
	pr.warnings             = 100;
	pr.initial_isotopes     = TRUE;
	pr.isotope_ratios       = TRUE;
	pr.isotope_alphas       = TRUE;
	pr.hdf                  = FALSE;
	pr.alkalinity           = FALSE;
	status_on               = true;
#ifdef NPP
	status_interval         = 40;
#else
	status_interval         = 250;
#endif
	status_timer            = clock();
	count_warnings          = 0;
	/* ----------------------------------------------------------------------
	*   RATES
	* ---------------------------------------------------------------------- */
	rate_m					= 0;
	rate_m0					= 0;
	rate_time				= 0;
	rate_kin_time           = 1.0;
	rate_sim_time_start		= 0;
	rate_sim_time_end		= 0;
	rate_sim_time			= 0;
	rate_moles				= 0;
	initial_total_time		= 0;
	// auto rate_p
	count_rate_p            = 0;
	/* ----------------------------------------------------------------------
	*   USER PRINT COMMANDS
	* ---------------------------------------------------------------------- */
	//user_print				= NULL;
	n_user_punch_index      = 0;
	fpunchf_user_s_warning  = 0;
	fpunchf_user_buffer[0]  = 0;

#if defined PHREEQ98 
	class rate *user_graph;
	char **user_graph_headings;
	int user_graph_count_headings;
#endif
#if defined MULTICHART
	// auto chart_handler;
	chart_handler.Set_io(phrq_io);
#endif
	/* ----------------------------------------------------------------------
	*   GLOBAL DECLARATIONS
	* ---------------------------------------------------------------------- */
	error_string            = NULL;
	simulation				= 0;
	state                   = INITIALIZE;
	reaction_step           = 0;
	transport_step          = 0;
	transport_start         = 0;
	advection_step          = 0;
	stop_program            = FALSE;
	incremental_reactions   = FALSE;
	count_strings           = 0;
	input_error             = 0;
	next_keyword            = Keywords::KEY_NONE;
	parse_error             = 0;
	paren_count             = 0;
	iterations              = 0;
	gamma_iterations        = 0;
	density_iterations = 0;
	run_reactions_iterations= 0;
	overall_iterations      = 0;
	max_line				= MAX_LINE;
	line                    = NULL;
	line_save				= NULL;
	LOG_10                  = log(10.0);
	debug_model             = FALSE;
	debug_prep              = FALSE;
	debug_set               = FALSE;
	debug_mass_action       = FALSE;
	debug_mass_balance      = FALSE;
	debug_diffuse_layer     = FALSE;
	debug_inverse           = FALSE;
#ifdef USE_LONG_DOUBLE
	/* from float.h, sets tolerance for cl1 routine */
	inv_tol_default         = pow((long double) 10, (long double) -LDBL_DIG + 5);
#else
	inv_tol_default         = pow((double) 10, (double) -DBL_DIG + 5);
#endif
	itmax                   = 100;
	max_tries               = 1000;
#ifdef USE_LONG_DOUBLE
	/* from float.h, sets tolerance for cl1 routine */
	ineq_tol                = pow((long double) 10, (long double) -LDBL_DIG);
#elif NPP
// appt:
	ineq_tol                = pow((double) 10, (double) -DBL_DIG + 2);
#else
	ineq_tol                = pow((double) 10, (double) -DBL_DIG);
#endif
	convergence_tolerance   = 1e-8;	
	step_size				= 100.;
	pe_step_size			= 10.;
	step_size_now           = step_size;
	pe_step_size_now        = pe_step_size;
	pp_scale				= 1.0;
	pp_column_scale			= 1.0;
	diagonal_scale			= FALSE;
	mass_water_switch		= FALSE;
	delay_mass_water		= FALSE;
	equi_delay      		= 0;
	dampen_ah2o             = false;
	censor					= 0.0;
	aqueous_only			= 0;
	negative_concentrations = FALSE;
	calculating_deriv		= FALSE;
	numerical_deriv			= FALSE;
	count_total_steps       = 0;
	phast                   = FALSE;
	output_newline          = true;
	//selected_output_file_name = NULL;
	dump_file_name			= NULL;
	remove_unstable_phases  = FALSE;
	// auto screen_string;
	spread_length           = 10;
	/* ----------------------------------------------------------------------
	*   ISOTOPES
	* ---------------------------------------------------------------------- */
	initial_solution_isotopes = FALSE;

	phreeqc_mpi_myself		= 0;
	first_read_input		= TRUE;
	print_density		    = 0;
	print_viscosity		    = 0;
	cell_pore_volume	    = 0;
	cell_volume			    = 0;
	cell_porosity		    = 0;
	cell_saturation		    = 0;
	sys_tot                 = 0;

	V_solutes               = 0.0;
	viscos                  = 0.0;
	viscos_0                = 0.0;
	viscos_0_25             = 0.0;
	density_x               = 0.0;
	rho_0                   = 0.0;
	kappa_0                 = 0.0;
	p_sat                   = 0.0;
	eps_r                   = EPSILON;
	DH_A                    = 0.0;
	DH_B                    = 0.0;
	DH_Av                   = 0.0;
	QBrn                    = 0.0;
	ZBrn                    = 0.0;
	dgdP                    = 0.0;

	need_temp_msg           = 0;
	solution_mass           = 0;
	solution_volume         = 0;
	/* phqalloc.cpp ------------------------------- */
	s_pTail                 = NULL;
	/* Basic */
	//basic_interpreter       = NULL;
	basic_callback_ptr      = NULL;
	basic_callback_cookie   = NULL;
	basic_fortran_callback_ptr  = NULL;

#ifdef SKIP
	/* dw.cpp ------------------------------- */
	/* COMMON /QQQQ/ */	
	Q0                      = 0;
	Q5                      = 0;
	GASCON                  = 0.461522e0;
	TZ                      = 647.073e0;
	AA                      = 1.e0;
	Z                       = 0;
	DZ                      = 0;
	Y                       = 0;
	G1                      = 11.e0;
	G2                      = 44.333333333333e0;
	GF                      = 3.5e0;
	B1                      = 0;
	B2                      = 0;
	B1T                     = 0;
	B2T                     = 0;
	B1TT                    = 0;
	B2TT                    = 0;
#endif
	/* gases.cpp ------------------------------- */
	a_aa_sum                = 0;
	b2                      = 0;
	b_sum                   = 0;
	R_TK                    = 0;
	/* input.cpp ------------------------------- */
	check_line_return       = 0;  
	reading_db              = FALSE;
	/* integrate.cpp ------------------------------- */
	midpoint_sv             = 0;
	z_global                = 0;
	xd_global               = 0;
	alpha_global            = 0;
	/* integrate.cpp ------------------------------- */
	max_row_count           = 50;
	max_column_count        = 50;
	carbon                  = FALSE;
	count_rows              = 0;
	count_optimize          = 0;
	col_phases              = 0;
	col_redox               = 0;
	col_epsilon             = 0;
	col_ph                  = 0;
	col_water               = 0;
	col_isotopes            = 0;
	col_phase_isotopes      = 0;
	row_mb                  = 0;
	row_fract               = 0;
	row_charge              = 0;
	row_carbon              = 0;
	row_isotopes            = 0;
	row_epsilon             = 0;
	row_isotope_epsilon     = 0;
	row_water               = 0;
	klmd                    = 0;
	nklmd                   = 0;
	n2d                     = 0;
	kode                    = 0;
	iter                    = 0;
	toler                   = 0;
	error                   = 0;
	max_pct                 = 0;
	scaled_error            = 0;
	master_alk              = NULL;
	max_good                = 0;
	max_bad                 = 0;
	max_minimal             = 0;
	count_good              = 0;
	count_bad               = 0;
	count_minimal           = 0;
	count_calls             = 0;
	soln_bits               = 0;
	phase_bits              = 0;
	current_bits            = 0;
	temp_bits               = 0;
	netpath_file            = NULL;
	count_inverse_models    = 0;
	count_pat_solutions     = 0;
	for (int i = 0; i < 32; i++)
	{
		min_position[i]     = 0;
		max_position[i]     = 0;
		now[i]              = 0;
	}
	/* kinetics.cpp ------------------------------- */
	count_pp = count_pg = count_ss = 0; 
	cvode_kinetics_ptr      = NULL;
	cvode_test              = FALSE;
	cvode_error             = FALSE;
	cvode_n_user            = -99;
	cvode_n_reactions       = -99;
	cvode_step_fraction     = 0.0;
	cvode_rate_sim_time     = 0.0;
	cvode_rate_sim_time_start = 0.0;
	cvode_last_good_time    = 0.0;
	cvode_prev_good_time    = 0.0;
	cvode_last_good_y       = NULL;
	cvode_prev_good_y       = NULL;
	kinetics_machEnv        = NULL;
	kinetics_y              = NULL;
	kinetics_abstol         = NULL;
	kinetics_cvode_mem      = NULL;
	cvode_pp_assemblage_save= NULL;
	cvode_ss_assemblage_save= NULL;
	set_and_run_attempt     = 0;
	/* model.cpp ------------------------------- */
	gas_in                  = FALSE;
	min_value               = 1e-10;

	/* phrq_io_output.cpp ------------------------------- */
	forward_output_to_log   = 0;
	/* phreeqc_files.cpp ------------------------------- */
#ifdef NPP
	default_data_base = "c:\\phreeqc\\database\\phreeqc.dat";
#else
	default_data_base = "phreeqc.dat";
#endif
	/* Pitzer  */	
	pitzer_model			= FALSE;
	sit_model				= FALSE;
	pitzer_pe				= FALSE;
	full_pitzer = FALSE;
	always_full_pitzer = FALSE;
	ICON					= TRUE;
	IC                      = -1;
	COSMOT                  = 0;
	AW                      = 0;
	VP                      = 0;
	DW0                     = 0;
	// auto pitz_param_map
	use_etheta				= TRUE;
	OTEMP					= -100.;
	OPRESS					= -100.;
	A0                      = 0;
	cations                 = NULL;
	anions                  = NULL;
	neutrals                = NULL;
	count_cations           = 0;
	count_anions            = 0;
	count_neutrals          = 0;
	MAXCATIONS              = 0;
	FIRSTANION              = 0;
	MAXNEUTRAL              = 0;
	mcb0                    = NULL;
	mcb1                    = NULL;
	mcc0                    = NULL;
	for (int i = 0; i < 23; i++)
	{
		BK[i]				= 0.0;
		DK[i]				= 0.0;
	}
	dummy                   = 0;
	/* print.cpp ------------------------------- */
	if (sformatf_buffer != NULL)
	{
		sformatf_buffer = (char*)free_check_null(sformatf_buffer);
	}
	sformatf_buffer = (char *) PHRQ_calloc(256 , sizeof(char));
	if (sformatf_buffer == NULL) 
			malloc_error();
	sformatf_buffer_size = 256;
	/* read.cpp */
	prev_next_char          = NULL;
#if defined PHREEQ98 
	int shifts_as_points;
#endif
	/* read_class.cxx */
	// auto dump_info
	// auto delete_info
	// auto run_info
	run_info.Set_io(phrq_io);
	/* readtr.cpp */
	// auto dump_file_name_cpp;
	/* sit.cpp ------------------------------- */

	sit_A0                  = 0;
	sit_count_cations       = 0;
	sit_count_anions        = 0;
	sit_count_neutrals      = 0;
	sit_MAXCATIONS          = 0;
	sit_FIRSTANION          = 0;
	sit_MAXNEUTRAL          = 0;
	/* tidy.cpp ------------------------------- */
	a0                      = 0;
	a1                      = 0;
	kc                      = 0;
	kb                      = 0;
	/* tally.cpp ------------------------------- */
	t_buffer                = NULL;
	tally_count_component   = 0;
	//tally_table             = NULL;
	count_tally_table_columns = 0;
	count_tally_table_rows  = 0;
	/* transport.cpp ------------------------------- */
	sol_D                   = NULL;
	sol_D_dbg               = NULL;
	J_ij                    = NULL;
	J_ij_il                 = NULL;
	J_ij_count_spec         = 0;
	m_s                     = NULL;
	count_m_s               = 0;
	tot1_h                  = 0;
	tot1_o                  = 0;
	tot2_h                  = 0;
	tot2_o                  = 0;
	diffc_max               = 0;
	diffc_tr                = 0;
	J_ij_sum                = 0;
	transp_surf             = FALSE;
	heat_mix_array          = NULL;
	temp1                   = NULL;
	temp2                   = NULL;
	nmix                    = 0;
	heat_nmix               = 0;
	heat_mix_f_imm          = 0;
	heat_mix_f_m            = 0;
	warn_MCD_X              = 0;
	warn_fixed_Surf         = 0;
	/* utilities.cpp ------------------------------- */
	spinner                 = 0;
	// keycount;
	keycount.resize(Keywords::KEY_COUNT_KEYWORDS);
	for (int i = 0; i < Keywords::KEY_COUNT_KEYWORDS; i++)
	{
		keycount[i] = 0;
	}

	return;
}
/*-----------------------------------------------------*/
Phreeqc::Phreeqc(const Phreeqc &src)
{
	user_print = NULL;
	sformatf_buffer = NULL;
	basic_interpreter = NULL;
	count_elts = 0;
	aphi = NULL;
	//this->phrq_io = src.phrq_io;
	this->phrq_io = &this->ioInstance;
	this->init();
	this->initialize();
	InternalCopy(&src);
}
void
Phreeqc::InternalCopy(const Phreeqc* pSrc)
{
	// phrq_io
	//this->phrq_io = new PHRQ_io;
	same_model = FALSE;
	current_tc = pSrc->current_tc;
	current_pa = pSrc->current_pa;
	current_mu = pSrc->current_mu;
	mu_terms_in_logk = pSrc->mu_terms_in_logk;

	/* ----------------------------------------------------------------------
	*   STRUCTURES
	* ---------------------------------------------------------------------- */
	//last_model, accept init
	high_precision = pSrc->high_precision;
	// Maps
	Rxn_temperature_map = pSrc->Rxn_temperature_map;
	Rxn_pressure_map = pSrc->Rxn_pressure_map;
	g_iterations = -1;
	G_TOL = pSrc->G_TOL;
	Rxn_surface_map = pSrc->Rxn_surface_map;
	change_surf_count = pSrc->change_surf_count;
	change_surf = change_surf_alloc(change_surf_count + 1);
	for (int ii = 0; ii < change_surf_count; ii++)
	{
		change_surf[ii].comp_name = string_hsave(pSrc->change_surf[ii].comp_name);
		change_surf[ii].fraction = pSrc->change_surf[ii].fraction;
		change_surf[ii].new_comp_name = string_hsave(pSrc->change_surf[ii].new_comp_name);
		change_surf[ii].new_Dw = pSrc->change_surf[ii].new_Dw;
		change_surf[ii].cell_no = pSrc->change_surf[ii].cell_no;
		change_surf[ii].next = pSrc->change_surf[ii].next;
	}
	Rxn_exchange_map = pSrc->Rxn_exchange_map;
	Rxn_kinetics_map = pSrc->Rxn_kinetics_map;
	use_kinetics_limiter = pSrc->use_kinetics_limiter;
	save_values = pSrc->save_values;
	save_strings = pSrc->save_strings;
	save = pSrc->save;
	//class copier copy_solution;
	//class copier copy_pp_assemblage;
	//class copier copy_exchange;
	//class copier copy_surface;
	//class copier copy_ss_assemblage;
	//class copier copy_gas_phase;
	//class copier copy_kinetics;
	//class copier copy_mix;
	//class copier copy_reaction;
	//class copier copy_temperature;
	//class copier copy_pressure;
	//	Inverse not implemented
	//std::vector<class inverse> inverse;
	count_inverse = 0;
	/* rate parameters */
	rate_parameters_pk = pSrc->rate_parameters_pk;
	rate_parameters_svd = pSrc->rate_parameters_svd;
	rate_parameters_hermanska = pSrc->rate_parameters_hermanska;
	// Mean gammas
	mean_gammas = pSrc->mean_gammas;
	//   Mix
	Rxn_mix_map = pSrc->Rxn_mix_map;
	Dispersion_mix_map = pSrc->Dispersion_mix_map;
	Rxn_solution_mix_map = pSrc->Rxn_solution_mix_map;
	Rxn_exchange_mix_map = pSrc->Rxn_exchange_mix_map;
	Rxn_gas_phase_mix_map = pSrc->Rxn_gas_phase_mix_map;
	Rxn_kinetics_mix_map = pSrc->Rxn_kinetics_mix_map;
	Rxn_pp_assemblage_mix_map = pSrc->Rxn_pp_assemblage_mix_map;
	Rxn_ss_assemblage_mix_map = pSrc->Rxn_ss_assemblage_mix_map;
	Rxn_surface_mix_map = pSrc->Rxn_surface_mix_map;
	//List new definitions
	//std::set<int> Rxn_new_exchange;
	//std::set<int> Rxn_new_gas_phase;
	//std::set<int> Rxn_new_kinetics;     // not used
	//std::set<int> Rxn_new_mix;          // not used
	//std::set<int> Rxn_new_pp_assemblage;
	//std::set<int> Rxn_new_pressure;     // not used
	//std::set<int> Rxn_new_reaction;     // not used
	//std::set<int> Rxn_new_solution;
	//std::set<int> Rxn_new_ss_assemblage;
	//std::set<int> Rxn_new_surface;
	//std::set<int> Rxn_new_temperature;  // not used
	Rxn_reaction_map = pSrc->Rxn_reaction_map;
	Rxn_gas_phase_map = pSrc->Rxn_gas_phase_map;
	Rxn_ss_assemblage_map = pSrc->Rxn_ss_assemblage_map;
	Rxn_pp_assemblage_map = pSrc->Rxn_pp_assemblage_map;

	std::vector<class species_list> species_list;
	// will be rebuilt
	//std::vector<class list0> sum_jacob0;	
	//std::vector<class list1> sum_mb1; 
	//std::vector<class list1> sum_jacob1;	
	//std::vector<class list2> sum_mb2; 
	//std::vector<class list2> sum_jacob2; 
	//std::vector<class list2> sum_delta; 
	// Solution
	Rxn_solution_map = pSrc->Rxn_solution_map;
	unnumbered_solutions = pSrc->unnumbered_solutions;
	save_species = pSrc->save_species;
	// Global solution
	title_x = pSrc->title_x;
	last_title_x = pSrc->last_title_x;
	//new_x                   = FALSE;
	description_x = pSrc->description_x;
	//new_x                   = FALSE;
	description_x = pSrc->description_x;
	//tc_x                    = 0;
	//tk_x                    = 0;
	//patm_x                  = 1;
	//last_patm_x             = 1;
	//numerical_fixed_volume  = false;
	//force_numerical_fixed_volume = false;
	//ph_x                    = 0;
	//solution_pe_x           = 0;
	//mu_x                    = 0;
	//ah2o_x                  = 1.0;
	//density_x               = 0;
	//total_h_x               = 0;
	//total_o_x               = 0;
	//cb_x                    = 0;
	//total_ions_x            = 0;
	//mass_water_aq_x         = 0;
	//mass_water_surfaces_x   = 0;
	//mass_water_bulk_x       = 0;
	//units_x
	//pe_x
	//isotopes_x
	//default_pe_x
	//dl_type_x                = cxxSurface::NO_DL;
	//total_carbon             = 0;
	//total_co2                = 0;
	//total_alkalinity         = 0;
	gfw_water = pSrc->gfw_water;
	//step_x                   = 0;
	//kin_time_x               = 0;
	//   Transport data
	count_cells = pSrc->count_cells;
	count_shifts = pSrc->count_shifts;
	ishift = pSrc->ishift;
	bcon_first = pSrc->bcon_first;
	bcon_last = pSrc->bcon_last;
	correct_disp = pSrc->correct_disp;
	tempr = pSrc->tempr;
	timest = pSrc->timest;
	simul_tr = pSrc->simul_tr;
	diffc = pSrc->diffc;
	heat_diffc = pSrc->heat_diffc;
	cell = pSrc->cell;
	mcd_substeps = pSrc->mcd_substeps;
	stag_data = pSrc->stag_data;
	print_modulus = pSrc->print_modulus;
	punch_modulus = pSrc->punch_modulus;
	dump_in = pSrc->dump_in;
	dump_modulus = pSrc->dump_modulus;
	transport_warnings = pSrc->transport_warnings;
	// cell_data 
	cell_data = pSrc->cell_data;
	old_cells = pSrc->old_cells;
	max_cells = pSrc->max_cells;
	if (stag_data.count_stag > 0)
	{
		max_cells = (max_cells - 2) / (1 + stag_data.count_stag);
	}
	all_cells = pSrc->all_cells;
	max_cells = pSrc->max_cells;
	multi_Dflag = pSrc->multi_Dflag;
	interlayer_Dflag = pSrc->interlayer_Dflag;
	implicit = pSrc->implicit;
	max_mixf = pSrc->max_mixf;
	min_dif_LM = pSrc->min_dif_LM;
	default_Dw = pSrc->default_Dw;
	correct_Dw = pSrc->correct_Dw;
	multi_Dpor = pSrc->multi_Dpor;
	interlayer_Dpor = pSrc->interlayer_Dpor;
	multi_Dpor_lim = pSrc->multi_Dpor_lim;
	interlayer_Dpor_lim = pSrc->interlayer_Dpor_lim;
	multi_Dn = pSrc->multi_Dn;
	interlayer_tortf = pSrc->interlayer_tortf;
	cell_no = pSrc->cell_no;
	mixrun = pSrc->mixrun;
	//  Advection data
	count_ad_cells = pSrc->count_ad_cells;
	count_ad_shifts = pSrc->count_ad_shifts;
	print_ad_modulus = pSrc->print_ad_modulus;
	punch_ad_modulus = pSrc->punch_ad_modulus;
	advection_punch = pSrc->advection_punch;
	advection_print = pSrc->advection_print;
	advection_kin_time = pSrc->advection_kin_time;
	advection_kin_time_defined = pSrc->advection_kin_time_defined;
	advection_warnings = pSrc->advection_warnings;
	// Tidy data
	new_model = TRUE;
	new_exchange = FALSE;
	new_pp_assemblage = FALSE;
	new_surface = FALSE;
	new_reaction = FALSE;
	new_temperature = FALSE;
	new_mix = FALSE;
	new_solution = FALSE;
	new_gas_phase = FALSE;
	new_inverse = FALSE;
	new_punch = FALSE;
	new_ss_assemblage = FALSE;
	new_kinetics = FALSE;
	new_copy = FALSE;
	new_pitzer = FALSE;
	// Elements
	for (int i = 0; i < (int)pSrc->elements.size(); i++)
	{
		const char* ptr = string_hsave(pSrc->elements[i]->name);
		class element* elt_ptr = element_store(ptr);
		elt_ptr->gfw = pSrc->elements[i]->gfw;
	}
	element_h_one = element_store("H(1)");
	// Element List
	count_elts = 0;
	// Reaction
	run_cells_one_step = pSrc->run_cells_one_step;
	//// logk
	//logk.clear();
	//for (size_t i = 0; i < pSrc->logk.size(); i++)
	//{
	//	class logk* tlk = new class logk;
	//	*tlk = *pSrc->logk[i];
	//	tlk->name = string_hsave(pSrc->logk[i]->name);
	//	logk.push_back(tlk);
	//}
	for (int i = 0; i < (int)pSrc->logk.size(); i++)
	{
		class logk* logk_ptr = logk_store(pSrc->logk[i]->name, FALSE);
		//memcpy(logk_ptr, pSrc->logk[i], sizeof(class logk));
		*logk_ptr = *pSrc->logk[i];
		logk_ptr->name = string_hsave(pSrc->logk[i]->name);
		logk_ptr->add_logk.resize(pSrc->logk[i]->add_logk.size());
		for (size_t j = 0; j < logk_ptr->add_logk.size(); j++)
		{
			logk_ptr->add_logk[j].coef = pSrc->logk[i]->add_logk[j].coef;
			logk_ptr->add_logk[j].name = string_hsave(pSrc->logk[i]->add_logk[j].name);
		}
	}
	// s, species
	for (int i = 0; i < (int)pSrc->s.size(); i++)
	{
		class species* s_ptr = s_store(pSrc->s[i]->name, pSrc->s[i]->z, FALSE);
		//memcpy(s_ptr, pSrc->s[i], sizeof(class species));
		*s_ptr = *pSrc->s[i];
		// fix up all pointers
		s_ptr->name = string_hsave(pSrc->s[i]->name);
		s_ptr->mole_balance = NULL;
		if (pSrc->s[i]->mole_balance != NULL)
		{
			s_ptr->mole_balance = string_hsave(pSrc->s[i]->mole_balance);
		}
		s_ptr->primary = NULL;
		s_ptr->secondary = NULL;
		s_ptr->add_logk.resize(pSrc->s[i]->add_logk.size());
		for (size_t j = 0; j < s_ptr->add_logk.size(); j++)
		{
			s_ptr->add_logk[j].coef = pSrc->s[i]->add_logk[j].coef;
			s_ptr->add_logk[j].name = string_hsave(pSrc->s[i]->add_logk[j].name);
		}
		//next_elt
		s_ptr->next_elt = elt_list_internal_copy(pSrc->s[i]->next_elt);
		s_ptr->next_secondary = elt_list_internal_copy(pSrc->s[i]->next_secondary);
		s_ptr->next_sys_total = elt_list_internal_copy(pSrc->s[i]->next_sys_total);
		//rxn
		s_ptr->rxn = CReaction_internal_copy(pSrc->s[i]->rxn);
		s_ptr->rxn_s = CReaction_internal_copy(pSrc->s[i]->rxn_s);
		s_ptr->rxn_x = CReaction_internal_copy(pSrc->s[i]->rxn_x);
	}
	s_diff_layer = pSrc->s_diff_layer;
	//s_x  will be built
	s_h2o = s_search("H2O");
	s_hplus = s_search("H+");
	s_h3oplus = s_search("H3O+");
	s_eminus = s_search("e-");
	s_co3 = s_search("CO3-2");
	s_h2 = s_search("H2");
	s_o2 = s_search("O2");
	/*----------------------------------------------------------------------
	*   Phases
	*---------------------------------------------------------------------- */
	for (int i = 0; i < (int)pSrc->phases.size(); i++)
	{
		class phase* phase_ptr = phase_store(pSrc->phases[i]->name);
		//memcpy(phase_ptr, pSrc->phases[i], sizeof(class phase));
		*phase_ptr = *pSrc->phases[i];
		// clean up pointers
		phase_ptr->name = string_hsave(pSrc->phases[i]->name);
		phase_ptr->formula = string_hsave(pSrc->phases[i]->formula);
		//add_logk
		phase_ptr->add_logk.resize(pSrc->phases[i]->add_logk.size());
		for (size_t j = 0; j < phase_ptr->add_logk.size(); j++)
		{
			phase_ptr->add_logk[j].coef = pSrc->phases[i]->add_logk[j].coef;
			phase_ptr->add_logk[j].name = string_hsave(pSrc->phases[i]->add_logk[j].name);
		}
		//next_elt
		phase_ptr->next_elt = elt_list_internal_copy(pSrc->phases[i]->next_elt);
		phase_ptr->next_sys_total = elt_list_internal_copy(pSrc->phases[i]->next_sys_total);
		//rxn
		phase_ptr->rxn = CReaction_internal_copy(pSrc->phases[i]->rxn);
		phase_ptr->rxn_s = CReaction_internal_copy(pSrc->phases[i]->rxn_s);
		phase_ptr->rxn_x = CReaction_internal_copy(pSrc->phases[i]->rxn_x);
	}
	// Master species
	for (size_t i = 0; i < pSrc->master.size(); i++)
	{
		master.resize(i + 1);
		master[i] = new class master;
		//memcpy(master[i], pSrc->master[i], sizeof(class master));
		*master[i] = *pSrc->master[i];
		// clean up pointers
		master[i]->gfw_formula = string_hsave(pSrc->master[i]->gfw_formula);
		master[i]->elt = element_store(pSrc->master[i]->elt->name);
		master[i]->unknown = NULL;
		master[i]->s = s_store(pSrc->master[i]->s->name, pSrc->master[i]->s->z, FALSE);
		//rxn_primary
		master[i]->rxn_primary = CReaction_internal_copy(pSrc->master[i]->rxn_primary);
		master[i]->rxn_secondary = CReaction_internal_copy(pSrc->master[i]->rxn_secondary);
	}
	// Unknowns will be built
	//x                       = NULL;
	//count_unknowns          = 0;
	//max_unknowns            = 0;
	//ah2o_unknown            = NULL;
	//alkalinity_unknown      = NULL;
	//carbon_unknown          = NULL;
	//charge_balance_unknown  = NULL;
	//exchange_unknown        = NULL;
	//mass_hydrogen_unknown   = NULL;
	//mass_oxygen_unknown     = NULL;
	//mb_unknown              = NULL;
	//mu_unknown              = NULL;
	//pe_unknown              = NULL;
	//ph_unknown              = NULL;
	//pure_phase_unknown      = NULL;
	//solution_phase_boundary_unknown = NULL;
	//surface_unknown         = NULL;
	//gas_unknown             = NULL;
	//ss_unknown              = NULL;
	//gas_unknowns;
	//mb_unknowns
	//   Reaction work space
	// class reaction_temp trxn;
	count_trxn = 0;
	// Print
	pr = pSrc->pr;
	status_on = pSrc->status_on;
	status_interval = pSrc->status_interval;
	status_timer = clock();
	status_string.clear();
	count_warnings = 0;
	// RATES
	//rates = pSrc->rates;
	for (size_t i = 0; i < pSrc->rates.size(); i++)
	{
		rates.push_back(*rate_copy(&pSrc->rates[i]));
	}
	//rate_m					= 0;
	//rate_m0					= 0;
	//rate_time				= 0;
	//rate_kin_time           = 1.0;
	//rate_sim_time_start		= 0;
	//rate_sim_time_end		= 0;
	//rate_sim_time			= 0;
	//rate_moles				= 0;
	initial_total_time = pSrc->initial_total_time;
	//rate_p
	count_rate_p = 0;
	// User print
	user_print = rate_copy(pSrc->user_print);
	// For now, User Punch is NOT copied
	n_user_punch_index = pSrc->n_user_punch_index;
	fpunchf_user_s_warning = pSrc->fpunchf_user_s_warning;
	//fpunchf_user_buffer[0]  = 0;
#if defined MULTICHART
	// auto chart_handler;
	chart_handler.Set_io(phrq_io);
#endif
	/* ----------------------------------------------------------------------
	*   GLOBAL DECLARATIONS
	* ---------------------------------------------------------------------- */
	error_string = NULL;
	simulation = pSrc->simulation;
	//state                   = INITIALIZE;
	//reaction_step           = 0;
	//transport_step          = 0;
	//transport_start         = 0;
	//advection_step          = 0;
	//stop_program            = FALSE;
	incremental_reactions = pSrc->incremental_reactions;
	// Constants
	MIN_LM = pSrc->MIN_LM;			    /* minimum log molality allowed before molality set to zero */
	LOG_ZERO_MOLALITY = pSrc->LOG_ZERO_MOLALITY;	/* molalities <= LOG_ZERO_MOLALITY are considered equal to zero */
	MIN_RELATED_LOG_ACTIVITY = pSrc->MIN_RELATED_LOG_ACTIVITY;
	MIN_TOTAL = pSrc->MIN_TOTAL;
	MIN_TOTAL_SS = pSrc->MIN_TOTAL_SS;
	MIN_RELATED_SURFACE = pSrc->MIN_RELATED_SURFACE;
	simulation = pSrc->simulation;
	//my_array, 
	//delta, 
	//residual
	input_error = 0;
	next_keyword = Keywords::KEY_NONE;
	parse_error = 0;
	paren_count = 0;
	iterations = 0;
	gamma_iterations = 0;
	density_iterations = 0;
	run_reactions_iterations = 0;
	overall_iterations = 0;
	free_check_null(line);
	free_check_null(line_save);
	max_line = pSrc->max_line;
	line = (char*)PHRQ_malloc(max_line * sizeof(char));
	line_save = (char*)PHRQ_malloc(max_line * sizeof(char));
	LOG_10 = pSrc->LOG_10;
	// Debug
	debug_model = pSrc->debug_model;
	debug_prep = pSrc->debug_prep;
	debug_set = pSrc->debug_set;
	debug_diffuse_layer = pSrc->debug_diffuse_layer;
	debug_inverse = pSrc->debug_inverse;
	//
	inv_tol_default = pSrc->inv_tol_default;
	itmax = pSrc->itmax;
	max_tries = pSrc->max_tries;
	ineq_tol = pSrc->ineq_tol;
	convergence_tolerance = pSrc->convergence_tolerance;
	step_size = pSrc->step_size;
	pe_step_size = pSrc->pe_step_size;
	step_size_now = step_size;
	pe_step_size_now = pe_step_size;
	pp_scale = pSrc->pp_scale;
	pp_column_scale = pSrc->pp_column_scale;
	diagonal_scale = pSrc->diagonal_scale;
	mass_water_switch = pSrc->mass_water_switch;
	delay_mass_water = pSrc->delay_mass_water;
	equi_delay = pSrc->equi_delay;
	dampen_ah2o = pSrc->dampen_ah2o;
	censor = pSrc->censor;
	aqueous_only = pSrc->aqueous_only;
	negative_concentrations = pSrc->negative_concentrations;
	calculating_deriv = pSrc->calculating_deriv;
	numerical_deriv = pSrc->numerical_deriv;
	count_total_steps = 0;
	phast = FALSE;
	output_newline = true;
	// llnl
	a_llnl = pSrc->a_llnl;
	b_llnl = pSrc->b_llnl;
	bdot_llnl = pSrc->bdot_llnl;
	llnl_temp = pSrc->llnl_temp;
	llnl_adh = pSrc->llnl_adh;
	llnl_bdh = pSrc->llnl_bdh;
	llnl_bdot = pSrc->llnl_bdot;
	llnl_co2_coefs = pSrc->llnl_co2_coefs;

	// Not implemented for now
	SelectedOutput_map = pSrc->SelectedOutput_map;
	{
		std::map<int, SelectedOutput>::iterator it = SelectedOutput_map.begin();
		for (; it != SelectedOutput_map.end(); it++)
		{
			//phrq_io->punch_open(it->second.Get_file_name().c_str());
			//it->second.Set_punch_ostream(phrq_io->Get_punch_ostream());
			//phrq_io->Set_punch_ostream(NULL);
			it->second.Set_punch_ostream(NULL);
		}
	}
	SelectedOutput_map.clear();

	UserPunch_map = pSrc->UserPunch_map;
	std::map<int, UserPunch>::iterator it = UserPunch_map.begin();
	for (; it != UserPunch_map.end(); it++)
	{
		class rate* rate_new = new class rate;
		rate_new = rate_copy(it->second.Get_rate());
		it->second.Set_rate(rate_new);
		it->second.Set_PhreeqcPtr(this);
	}

	remove_unstable_phases = FALSE;
	//screen_string;
	spread_length = pSrc->spread_length;
	//maps set by store below
	//std::map<std::string, std::string*> strings_map;
	//std::map<std::string, class element*> elements_map;
	//std::map<std::string, class species*> species_map;
	//std::map<std::string, class phase*> phases_map;
	//std::map<std::string, class logk*> logk_map;
	//std::map<std::string, class master_isotope*> master_isotope_map;
	/* ----------------------------------------------------------------------
	*   ISOTOPES
	* ---------------------------------------------------------------------- */
	for (int i = 0; i < (int)pSrc->master_isotope.size(); i++)
	{
		class master_isotope* master_isotope_ptr = master_isotope_store(pSrc->master_isotope[i]->name, FALSE);
		// memcpy(master_isotope_ptr, pSrc->master_isotope[i], sizeof(class master_isotope));
		*master_isotope_ptr = *pSrc->master_isotope[i];
		master_isotope_ptr->name = string_hsave(pSrc->master_isotope[i]->name);
		int n;
		master_isotope_ptr->master = NULL;
		if (pSrc->master_isotope[i]->master)
		{
			master_isotope_ptr->master = master_search(pSrc->master_isotope[i]->master->elt->name, &n);
		}
		if (master_isotope_ptr->master == NULL)
		{
			//error_msg("Error in copy constructor for master_isotope.", STOP);
		}
		master_isotope_ptr->elt = NULL;
		if (pSrc->master_isotope[i]->elt)
		{
			master_isotope_ptr->elt = element_store(pSrc->master_isotope[i]->elt->name);
		}
		master_isotope_ptr->units = NULL;
		if (pSrc->master_isotope[i]->units)
		{
			master_isotope_ptr->units = string_hsave(pSrc->master_isotope[i]->units);
		}
	}
	initial_solution_isotopes = pSrc->initial_solution_isotopes;
	// Calculate values
	for (int i = 0; i < pSrc->calculate_value.size(); i++)
	{
		class calculate_value* calculate_value_ptr = calculate_value_store(pSrc->calculate_value[i]->name, FALSE);
		calculate_value_ptr->value = pSrc->calculate_value[i]->value;
		calculate_value[i]->commands = pSrc->calculate_value[i]->commands;
	}
	// More isotopes
	for (int i = 0; i < (int)pSrc->isotope_ratio.size(); i++)
	{
		class isotope_ratio* isotope_ratio_ptr = isotope_ratio_store(pSrc->isotope_ratio[i]->name, FALSE);
		isotope_ratio_ptr->name = string_hsave(pSrc->isotope_ratio[i]->name);
		isotope_ratio_ptr->isotope_name = string_hsave(pSrc->isotope_ratio[i]->isotope_name);
		isotope_ratio_ptr->ratio = pSrc->isotope_ratio[i]->ratio;
		isotope_ratio_ptr->converted_ratio = pSrc->isotope_ratio[i]->converted_ratio;
	}
	//std::map<std::string, class isotope_ratio*> isotope_ratio_map;
	for (int i = 0; i < (int)pSrc->isotope_alpha.size(); i++)
	{
		class isotope_alpha* isotope_alpha_ptr = isotope_alpha_store(pSrc->isotope_alpha[i]->name, FALSE);
		isotope_alpha_ptr->named_logk = string_hsave(pSrc->isotope_alpha[i]->named_logk);
		isotope_alpha_ptr->value = pSrc->isotope_alpha[i]->value;
	}
	//std::map<std::string, class isotope_alpha*> isotope_alpha_map;
	// Misc
	phreeqc_mpi_myself = 0;
	first_read_input = pSrc->first_read_input;
	user_database = pSrc->user_database;
	//have_punch_name			= pSrc->have_punch_name;
	print_density = pSrc->print_density;
	print_viscosity = pSrc->print_viscosity;
	viscos = pSrc->viscos;
	viscos_0 = pSrc->viscos_0;
	viscos_0_25 = pSrc->viscos_0_25; // viscosity of the solution, of pure water, of pure water at 25 C
	density_x = pSrc->density_x;
	solution_volume_x = pSrc->solution_volume_x;
	solution_mass_x = pSrc->solution_mass_x;
	kgw_kgs = pSrc->kgw_kgs;
	cell_pore_volume = pSrc->cell_pore_volume;
	cell_porosity = pSrc->cell_porosity;
	cell_volume = pSrc->cell_volume;
	cell_saturation = pSrc->cell_saturation;
	sys.clear();
	sys_tot = pSrc->sys_tot;
	// solution properties
	V_solutes = pSrc->V_solutes;
	rho_0 = pSrc->rho_0;
	kappa_0 = pSrc->kappa_0;
	p_sat = pSrc->p_sat;
	eps_r = pSrc->eps_r;
	DH_A = pSrc->DH_A;
	DH_B = pSrc->DH_B;
	DH_Av = pSrc->DH_Av;
	QBrn = pSrc->QBrn;
	ZBrn = pSrc->ZBrn;
	dgdP = pSrc->dgdP;
	//
	need_temp_msg = pSrc->need_temp_msg;
	solution_mass = pSrc->solution_mass;
	solution_volume = pSrc->solution_volume;
	s_pTail = NULL;
	//basic_interpreter = NULL;
	/* cl1.cpp ------------------------------- */
	//std::vector<double> x_arg, res_arg, scratch;
	// gases.cpp 
	a_aa_sum = pSrc->a_aa_sum;
	b2 = pSrc->b2;
	b_sum = pSrc->b_sum;
	R_TK = pSrc->R_TK;
	gas_binary_parameters = pSrc->gas_binary_parameters;
	/* input.cpp ------------------------------- */
	check_line_return = 0;
	reading_db = FALSE;
	/* integrate.cpp ------------------------------- */
	midpoint_sv = pSrc->midpoint_sv;
	z_global = pSrc->z_global;
	xd_global = pSrc->xd_global;
	alpha_global = pSrc->alpha_global;
	/* inverse.cpp ------------------------------- */	/* integrate.cpp ------------------------------- */
	max_row_count = pSrc->max_row_count;
	max_column_count = pSrc->max_column_count;
	carbon = pSrc->carbon;
	//std::vector<const char*> col_name, row_name;
	count_rows = pSrc->count_rows;
	count_optimize = pSrc->count_optimize;
	col_phases = pSrc->col_phases;
	col_redox = pSrc->col_redox;
	col_epsilon = pSrc->col_epsilon;
	col_ph = pSrc->col_ph;
	col_water = pSrc->col_water;
	col_isotopes = pSrc->col_isotopes;
	col_phase_isotopes = pSrc->col_phase_isotopes;
	row_mb = pSrc->row_mb;
	row_fract = pSrc->row_fract;
	row_charge = pSrc->row_charge;
	row_carbon = pSrc->row_carbon;
	row_isotopes = pSrc->row_isotopes;
	row_epsilon = pSrc->row_epsilon;
	row_isotope_epsilon = pSrc->row_isotope_epsilon;
	row_water = pSrc->row_water;
	//std::vector<double> inv_zero, array1, inv_res, inv_delta1, delta2,
	//	delta3, inv_cu, delta_save;
	//std::vector<double> min_delta, max_delta;
	//std::vector<int> inv_iu, inv_is;
	klmd = pSrc->klmd;
	nklmd = pSrc->nklmd;
	n2d = pSrc->n2d;
	kode = pSrc->kode;
	iter = pSrc->iter;
	toler = pSrc->toler;
	error = pSrc->error;
	max_pct = pSrc->max_pct;
	scaled_error = pSrc->scaled_error;
	master_alk = NULL;
	//std::vector<int> row_back, col_back;
	//std::vector<unsigned long> good, bad, minimal;
	max_good = pSrc->max_good;
	max_bad = pSrc->max_bad;
	max_minimal = pSrc->max_minimal;
	count_good = pSrc->count_good;
	count_bad = pSrc->count_bad;
	count_minimal = pSrc->count_minimal;
	count_calls = pSrc->count_calls;
	soln_bits = pSrc->soln_bits;
	phase_bits = pSrc->phase_bits;
	current_bits = pSrc->current_bits;
	temp_bits = pSrc->temp_bits;
	netpath_file = NULL;
	count_inverse_models = pSrc->count_inverse_models;
	count_pat_solutions = pSrc->count_pat_solutions;
	for (int i = 0; i < 32; i++)
	{
		min_position[i] = pSrc->min_position[i];
		max_position[i] = pSrc->max_position[i];
		now[i] = pSrc->now[i];
	}
	//std::vector <std::string> inverse_heading_names;
	/* kinetics.cpp ------------------------------- */
	count_pp = count_pg = count_ss = 0;
	cvode_kinetics_ptr = NULL;
	cvode_test = FALSE;
	cvode_error = FALSE;
	cvode_n_user = -99;
	cvode_n_reactions = -99;
	cvode_step_fraction = 0.0;
	cvode_rate_sim_time = 0.0;
	cvode_rate_sim_time_start = 0.0;
	cvode_last_good_time = 0.0;
	cvode_prev_good_time = 0.0;
	cvode_last_good_y = NULL;
	cvode_prev_good_y = NULL;
	kinetics_machEnv = NULL;
	kinetics_y = NULL;
	kinetics_abstol = NULL;
	kinetics_cvode_mem = NULL;
	cvode_pp_assemblage_save = NULL;
	cvode_ss_assemblage_save = NULL;
	//std::vector<double> m_temp, m_original, rk_moles, x0_moles;
	set_and_run_attempt = 0;
	/* model.cpp ------------------------------- */
	gas_in = FALSE;
	min_value = 1e-10;
	//std::vector<double> normal, ineq_array, res, cu, zero, delta1;
	//std::vector<int> iu, is, back_eq;									 
	/* phrq_io_output.cpp ------------------------------- */
	forward_output_to_log = pSrc->forward_output_to_log;
	/* phreeqc_files.cpp ------------------------------- */
	default_data_base = pSrc->default_data_base;
	// Pitzer 
	pitzer_model = pSrc->pitzer_model;
	sit_model = pSrc->sit_model;
	pitzer_pe = pSrc->pitzer_pe;


	full_pitzer = pSrc->full_pitzer;
	always_full_pitzer = pSrc->always_full_pitzer;
	ICON = pSrc->ICON;
	IC = pSrc->IC;
	COSMOT = pSrc->COSMOT;
	AW = pSrc->AW;
	VP = pSrc->VP;
	DW0 = pSrc->DW0;
	for (int i = 0; i < (int)pSrc->pitz_params.size(); i++)
	{
		pitz_param_store(pSrc->pitz_params[i]);
	}

	//pitz_param_map = pSrc->pitz_param_map; created by store
	for (int i = 0; i < (int)pSrc->theta_params.size(); i++)
	{
		size_t count_theta_params = theta_params.size();
		theta_params.resize(count_theta_params + 1);
		theta_params[count_theta_params] = new class theta_param;
		*theta_params[count_theta_params] = *pSrc->theta_params[i];
	}
	use_etheta = pSrc->use_etheta;
	OTEMP = pSrc->OTEMP;
	OPRESS = pSrc->OPRESS;
	A0 = pSrc->A0;
	aphi = pitz_param_copy(pSrc->aphi);
	// will be rebuilt
	spec = pSrc->spec;
	cations = pSrc->cations;
	anions = pSrc->anions;
	neutrals = pSrc->neutrals;
	count_cations = pSrc->count_cations;
	count_anions = pSrc->count_anions;
	count_neutrals = pSrc->count_neutrals;
	MAXCATIONS = pSrc->MAXCATIONS;
	FIRSTANION = pSrc->FIRSTANION;
	MAXNEUTRAL = pSrc->MAXNEUTRAL;
	mcb0 = pSrc->mcb0;
	mcb1 = pSrc->mcb1;
	mcc0 = pSrc->mcc0;
	IPRSNT = pSrc->IPRSNT;
	M = pSrc->M;
	LGAMMA = pSrc->LGAMMA;
	for (int i = 0; i < 23; i++)
	{
		BK[i] = pSrc->BK[i];
		DK[i] = pSrc->DK[i];
	}
	dummy = 0;
	/* print.cpp ------------------------------- */
	/*
	sformatf_buffer = (char *) PHRQ_malloc(256 * sizeof(char));
	if (sformatf_buffer == NULL)
		malloc_error();
	sformatf_buffer_size = 256;
	*/
	/* read.cpp */
	prev_next_char = NULL;
#if defined PHREEQ98 
	int shifts_as_points;
#endif
	/* read_class.cxx */
	// auto dump_info
	// auto delete_info
	// auto run_info
	/*
	run_info.Set_io(phrq_io);
	*/
	/* readtr.cpp */
	// auto dump_file_name_cpp;
	/* sit.cpp ------------------------------- */
	for (int i = 0; i < (int)pSrc->sit_params.size(); i++)
	{
		sit_param_store(pSrc->sit_params[i]);
	}
	//sit_param_map = pSrc->sit_param_map; // filled by store
	sit_A0 = pSrc->sit_A0;
	sit_count_cations = pSrc->sit_count_cations;
	sit_count_anions = pSrc->sit_count_anions;
	sit_count_neutrals = pSrc->sit_count_neutrals;
	sit_MAXCATIONS = pSrc->sit_MAXCATIONS;
	sit_FIRSTANION = pSrc->sit_FIRSTANION;
	sit_MAXNEUTRAL = pSrc->sit_MAXNEUTRAL;
	sit_IPRSNT = pSrc->sit_IPRSNT;
	sit_M = pSrc->sit_M;
	sit_LGAMMA = pSrc->sit_LGAMMA;
	s_list = pSrc->s_list;
	cation_list = pSrc->cation_list;
	neutral_list = pSrc->neutral_list;
	anion_list = pSrc->anion_list;
	ion_list = pSrc->ion_list;
	param_list = pSrc->param_list;

	/* tidy.cpp ------------------------------- */
	//a0                      = 0;
	//a1                      = 0;
	//kc                      = 0;
	//kb                      = 0;
	/* tally.cpp ------------------------------- */
	//t_buffer                = NULL;
	//tally_count_component   = 0;
	//tally_table             = NULL;
	//count_tally_table_columns = 0;
	//count_tally_table_rows  = 0;

	/* transport.cpp ------------------------------- */
	/* storage is created and freed in transport.cpp */
	sol_D = NULL;
	sol_D_dbg = NULL;
	J_ij = NULL;
	J_ij_il = NULL;
	J_ij_count_spec = pSrc->J_ij_count_spec;
	m_s = NULL;
	count_m_s = pSrc->count_m_s;
	tot1_h = pSrc->tot1_h;
	tot1_o = pSrc->tot1_o;
	tot2_h = pSrc->tot2_h;
	tot2_o = pSrc->tot2_o;
	diffc_max = pSrc->diffc_max;
	diffc_tr = pSrc->diffc_tr;
	J_ij_sum = pSrc->J_ij_sum;
	transp_surf = pSrc->transp_surf;
	heat_mix_array = NULL;
	temp1 = NULL;
	temp2 = NULL;
	nmix = pSrc->nmix;
	heat_nmix = pSrc->heat_nmix;
	heat_mix_f_imm = pSrc->heat_mix_f_imm;
	heat_mix_f_m = pSrc->heat_mix_f_m;
	warn_MCD_X = pSrc->warn_MCD_X;
	warn_fixed_Surf = pSrc->warn_fixed_Surf;
	current_x = pSrc->current_x;
	current_A = pSrc->current_A;
	fix_current = pSrc->fix_current;

	/* utilities.cpp ------------------------------- */
	//spinner                 = 0;
	//// keycount;
	//for (int i = 0; i < Keywords::KEY_COUNT_KEYWORDS; i++)
	//{
	//	keycount.push_back(0);
	//}
	spinner = pSrc->spinner;
	gfw_map = pSrc->gfw_map;
	//rates_map = pSrc->rates_map;
	sum_species_map = pSrc->sum_species_map;
	sum_species_map_db = pSrc->sum_species_map_db;

	// make sure new_model gets set
	this->keycount[Keywords::KEY_SOLUTION_SPECIES] = 1;
	this->tidy_model();
	return;
}
// Operator overloaded using a member function
Phreeqc &Phreeqc::operator=(const Phreeqc &rhs) 
{
	if (this == &rhs)      // Same object?
		return *this; 

	// clean up this here
	this->clean_up();

	this->PHRQ_free_all();
	if (this->phrq_io == &this->ioInstance)
	{
		this->phrq_io->clear_istream();
		this->phrq_io->close_ostreams();
	}

	// copy Phreeqc object to this
	//this->phrq_io = rhs.phrq_io;
	//this->phrq_io = new PHRQ_io;
#if !defined(R_SO)
	this->phrq_io->Set_output_ostream(&std::cout);
	this->phrq_io->Set_error_ostream(&std::cerr);
#endif	
	this->init();
	this->initialize();
	this->InternalCopy(&rhs);
	return *this;
}

int Phreeqc::next_user_number(Keywords::KEYWORDS key)
{
	switch (key)
	{
	case Keywords::KEY_REACTION_TEMPERATURE:
		return Utilities::Rxn_next_user_number(Rxn_temperature_map);
		break;
	case Keywords::KEY_REACTION_PRESSURE:
		return Utilities::Rxn_next_user_number(Rxn_pressure_map);
		break;
	case Keywords::KEY_SURFACE:
		return Utilities::Rxn_next_user_number(Rxn_surface_map);
		break;
	case Keywords::KEY_EXCHANGE:
		return Utilities::Rxn_next_user_number(Rxn_exchange_map);
		break;
	case Keywords::KEY_KINETICS:
		return Utilities::Rxn_next_user_number(Rxn_kinetics_map);
		break;
	case Keywords::KEY_MIX:
		return Utilities::Rxn_next_user_number(Rxn_mix_map);
		break;
	case Keywords::KEY_REACTION:
		return Utilities::Rxn_next_user_number(Rxn_reaction_map);
		break;
	case Keywords::KEY_GAS_PHASE:
		return Utilities::Rxn_next_user_number(Rxn_gas_phase_map);
		break;
	case Keywords::KEY_SOLID_SOLUTIONS:
		return Utilities::Rxn_next_user_number(Rxn_ss_assemblage_map);
		break;
	case Keywords::KEY_EQUILIBRIUM_PHASES:
		return Utilities::Rxn_next_user_number(Rxn_pp_assemblage_map);
		break;
	case Keywords::KEY_SOLUTION:
		return Utilities::Rxn_next_user_number(Rxn_solution_map);
		break;
	default:
		assert(false);
		return -999;
	}
}
