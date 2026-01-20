#ifndef _INC_GLOBAL_STRUCTURES_H
#define _INC_GLOBAL_STRUCTURES_H
#include "Surface.h"
#include "GasPhase.h"
/* ----------------------------------------------------------------------
 *   #define DEFINITIONS
 * ---------------------------------------------------------------------- */
#if !defined(NAN)
#  if defined(_MSC_VER) && (_MSC_VER <= 1700) // VS2012
//   https://learn.microsoft.com/en-us/cpp/preprocessor/predefined-macros?view=msvc-170     
#    include <limits>
#    define NAN std::numeric_limits<double>::signaling_NaN()
#  else
#    define NAN nan("1")
#  endif
#endif
#define MISSING -9999.999            
#include "NA.h"   /* NA = not available */

#define F_C_MOL 96493.5			/* C/mol or joule/volt-eq */
#define F_KJ_V_EQ  96.4935		/* kJ/volt-eq */
#define F_KCAL_V_EQ 23.0623		/* kcal/volt-eq */
#define R_LITER_ATM 0.0820597	/* L-atm/deg-mol */
#define R_KCAL_DEG_MOL 0.00198726	/* kcal/deg-mol */
#define R_KJ_DEG_MOL 0.00831470	/* kJ/deg-mol */
#define EPSILON 78.5			/* dialectric constant, dimensionless. Is calculated as eps_r(P, T) in calc_dielectrics. Update the code?? */
#define EPSILON_ZERO 8.854e-12	/* permittivity of free space, C/V-m = C**2/m-J */
#define JOULES_PER_CALORIE 4.1840
#define PASCAL_PER_ATM 1.01325E5 /* conversion from atm to Pa */
#define AVOGADRO 6.02252e23		/* atoms / mole */
#define pi 3.14159265358979
#define AH2O_FACTOR 0.017

#define TRUE 1
#define FALSE 0
#define OK 1
#define ERROR 0
#define STOP 1
#define CONTINUE 0

#define OPTION_EOF -1
#define OPTION_KEYWORD -2
#define OPTION_ERROR -3
#define OPTION_DEFAULT -4
#define OPT_1 -5

#define DISP 2
#define STAG 3
#define NOMIX 4
#define MIX_BS 5 // mix boundary solutions in electromigration

#define CONVERGED 2
#define MASS_BALANCE 3

#define REWRITE 2
#define INIT -1

/* check_line values, plus EMPTY, EOF, OK */
#define KEYWORD 3

/* copy_token values */
#define EMPTY 2
#define UPPER 4
#define LOWER 5
#define DIGIT 6
#define UNKNOWN 7
#define OPTION 8

/* species types */
#define AQ 0
#define HPLUS 1
#define H2O 2
#define EMINUS 3
#define SOLID 4
#define EX 5
#define SURF 6
#define SURF_PSI 7
#define SURF_PSI1 8
#define SURF_PSI2 9

/* unknown types */
#define MB 10
#define ALK 11
#define CB 12
#define SOLUTION_PHASE_BOUNDARY 13
#define MU 14
#define AH2O 15
#define MH 16
#define MH2O 17
#define PP 18
#define EXCH 19
#define SURFACE 20
#define SURFACE_CB 21
#define SURFACE_CB1 22
#define SURFACE_CB2 23
#define GAS_MOLES 24
#define SS_MOLES 25
#define PITZER_GAMMA 26
#define SLACK 28
/* state */
#define INITIALIZE	       0
#define INITIAL_SOLUTION   1
#define INITIAL_EXCHANGE   2
#define INITIAL_SURFACE 3
#define INITIAL_GAS_PHASE  4
#define REACTION		   5
#define INVERSE		 6
#define ADVECTION		 7
#define TRANSPORT		 8
#define PHAST		     9

/* constraints in mass balance */
#define EITHER 0
#define DISSOLVE 1
#define PRECIPITATE -1

/* gas phase type */
#define PRESSURE 1
#define VOLUME 2

#define MAX_PP_ASSEMBLAGE 10	/* default estimate of the number of phase assemblages */
#define MAX_ADD_EQUATIONS 20	/* maximum number of equations added together to reduce eqn to
								   master species */
#define MAX_ELEMENTS 50			/* default estimate of the number of elements */
#define MAX_LENGTH 256			/* maximum number of characters component name */
#define MAX_LINE 4096			/* estimate of maximum line length */
#define MAX_MASS_BALANCE 10		/* initial guess of number mass balance equations for a solution */
#define MAX_MASTER 50			/* default estimate of the number of master species */
#define MAX_ELTS 15				/* default estimate for maximum number of times elements occur in
								   an equation */
#define MAX_PHASES 500			/* initial guess of number of phases defined */
#define MAX_S 500				/* default estimate for maximum number of species in aqueous model */
#define MAX_SUM_JACOB0 50		/* list used to calculate jacobian */
#define MAX_SUM_JACOB1 500		/* list used to calculate jacobian */
#define MAX_SUM_JACOB2 500		/* list used to calculate jacobian */
#define MAX_SUM_MB 500			/* list used to calculate mass balance sums */
#define MAX_TRXN 16				/* default estimate for maximum number of components in an eqn */
#define MAX_UNKNOWNS 15			/* default estimate for maximum number of unknowns in model */
#define TOL 1e-9				/* tolerance for comparisons of double numbers */
#define MAX_LM 3.0				/* maximum log molality allowed in intermediate iterations */
#define MAX_M 1000.0
#ifdef USE_DECIMAL128
// #define MIN_LM -80.0			/* minimum log molality allowed before molality set to zero */
// #define LOG_ZERO_MOLALITY -80	/* molalities <= LOG_ZERO_MOLALITY are considered equal to zero */
// #define MIN_TOTAL 1e-60
// #define MIN_TOTAL_SS MIN_TOTAL/100
// #define MIN_RELATED_SURFACE MIN_TOTAL*100
// #define MIN_RELATED_LOG_ACTIVITY -60
#else
// #define MIN_LM -30.0			/* minimum log molality allowed before molality set to zero */
// #define LOG_ZERO_MOLALITY -30	/* molalities <= LOG_ZERO_MOLALITY are considered equal to zero */
// #define MIN_TOTAL 1e-25
// #define MIN_TOTAL_SS MIN_TOTAL/100
// #define MIN_RELATED_SURFACE MIN_TOTAL*100
// #define MIN_RELATED_LOG_ACTIVITY -30
#endif
#define REF_PRES_PASCAL 1.01325E5   /* Reference pressure: 1 atm */
#define MAX_P_NONLLNL 1500.0

//
// Typedefs and structure definitions
//
typedef enum { kcal, cal, kjoules, joules } DELTA_H_UNIT;
typedef enum { cm3_per_mol, dm3_per_mol, m3_per_mol } DELTA_V_UNIT;
enum entity_type
{ Solution, Reaction, Exchange, Surface, Gas_phase, Pure_phase, Ss_phase,
	Kinetics, Mix, Temperature, Pressure, UnKnown
};

typedef enum {
	logK_T0,
	delta_h,
	T_A1,
	T_A2,
	T_A3,
	T_A4,
	T_A5,
	T_A6,
	delta_v,	/* set in calc_delta_v: calculated molar volume-change of the reaction */
	vm_tc,		/* set in calc_vm: calculated molal volume of the species at tc */
	vm0,		/* read: molar volume of a phase */
	vma1, vma2, vma3, vma4, /* read: a1..a4 from supcrt, see calc_vm */
	wref,       /* from supcrt */
	b_Av,		/* b in z^2 * A_v * log(1 + b * I^0.5) / (2 * b) */
	vmi1, vmi2, vmi3, vmi4, /* ionic strength terms: (i1 + i2/(TK - 228) + i3 * (TK - 228) ) * I^i4 */
	MAX_LOG_K_INDICES	/* Keep this definition at the end of the enum */
} LOG_K_INDICES;

typedef struct PHRQMemHeader
{
	struct PHRQMemHeader *pNext;	/* memory allocated just after this one */
	struct PHRQMemHeader *pPrev;	/* memory allocated just prior to this one */
	size_t size;				/* memory request + sizeof(PHRQMemHeader) */
#if !defined(NDEBUG)
	char *szFileName;			/* file name */
	int nLine;					/* line number */
	int dummy;					/* alignment */
#endif
} PHRQMemHeader;

struct Change_Surf
{
	const char *comp_name;
	LDBLE fraction;
	const char *new_comp_name;
	LDBLE new_Dw;
	int cell_no;
	int next;
};
/*----------------------------------------------------------------------
 *   CReaction
 *---------------------------------------------------------------------- */
class rxn_token
{
public:
	~rxn_token() {};
	rxn_token()
	{
		s = NULL;
		coef = 0.0;
		name = NULL;
	}
	class species* s;
	LDBLE coef;
	const char* name;
};
class CReaction
{
public:
	CReaction(void);
	CReaction(size_t ntoken);
	~CReaction(void) {}
	double* Get_logk(void) { return this->logk; }
	void   Set_logk(double* d);
	double* Get_dz(void) { return this->dz; }
	void   Set_dz(double* d);
	size_t size() { return token.size(); }
	std::vector<class rxn_token>& Get_tokens(void) { return this->token; }
	void Set_tokens(const std::vector<class rxn_token>& t) { this->token = t; }

public:
	double logk[MAX_LOG_K_INDICES];
	double dz[3];
	std::vector<class rxn_token> token;
};
class save
{
public:
	~save() {};
	save()
	{
		solution = 0;
		n_solution_user = 0;
		n_solution_user_end = 0;
		mix = 0;
		n_mix_user = 0;
		n_mix_user_end = 0;
		reaction = 0;
		n_reaction_user = 0;
		n_reaction_user_end = 0;
		pp_assemblage = 0;
		n_pp_assemblage_user = 0;
		n_pp_assemblage_user_end = 0;
		exchange = 0;
		n_exchange_user = 0;
		n_exchange_user_end = 0;
		kinetics = 0;
		n_kinetics_user = 0;
		n_kinetics_user_end = 0;
		surface = 0;
		n_surface_user = 0;
		n_surface_user_end = 0;
		gas_phase = 0;
		n_gas_phase_user = 0;
		n_gas_phase_user_end = 0;
		ss_assemblage = 0;
		n_ss_assemblage_user = 0;
		n_ss_assemblage_user_end = 0;
	}
	int solution;
	int n_solution_user;
	int n_solution_user_end;
	int mix;
	int n_mix_user;
	int n_mix_user_end;
	int reaction;
	int n_reaction_user;
	int n_reaction_user_end;
	int pp_assemblage;
	int n_pp_assemblage_user;
	int n_pp_assemblage_user_end;
	int exchange;
	int n_exchange_user;
	int n_exchange_user_end;
	int kinetics;
	int n_kinetics_user;
	int n_kinetics_user_end;
	int surface;
	int n_surface_user;
	int n_surface_user_end;
	int gas_phase;
	int n_gas_phase_user;
	int n_gas_phase_user_end;
	int ss_assemblage;
	int n_ss_assemblage_user;
	int n_ss_assemblage_user_end;
};

/*----------------------------------------------------------------------
 *   Copy
 *---------------------------------------------------------------------- */
class copier
{
public:
	~copier() {};
	copier()
	{}
	std::vector<int> n_user;
	std::vector<int> start;
	std::vector<int> end;
};
/*----------------------------------------------------------------------
 *   Inverse
 *---------------------------------------------------------------------- */
class inv_elts
{
public:
	~inv_elts() {};
	inv_elts()
	{
		name = NULL;
		master = NULL;
		row = 0;
		//uncertainties.clear();
	}
	const char* name;
	class master* master;
	size_t row;
	std::vector<double> uncertainties;
};
class isotope
{
public:
	~isotope() {};
	isotope()
	{
		isotope_number = 0;
		elt_name = NULL;
		isotope_name = NULL;
		total = 0;
		ratio = 0;
		ratio_uncertainty = 0;
		x_ratio_uncertainty = 0;
		master = NULL;
		primary = NULL;
		coef = 0;					/* coefficient of element in phase */
	}
	LDBLE isotope_number;
	const char* elt_name;
	const char* isotope_name;
	LDBLE total;
	LDBLE ratio;
	LDBLE ratio_uncertainty;
	LDBLE x_ratio_uncertainty;
	class master* master;
	class master* primary;
	LDBLE coef;
};
class inv_isotope
{
public:
	~inv_isotope() {};
	inv_isotope()
	{
		isotope_name = NULL;
		isotope_number = 0;
		elt_name = NULL;
		//uncertainties.clear();
	}
	const char* isotope_name;
	LDBLE isotope_number;
	const char* elt_name;
	std::vector<double> uncertainties;
};
class inv_phases
{
public:
	~inv_phases() {};
	inv_phases()
	{
		name = NULL;
		phase = NULL;
		column = 0;
		constraint = EITHER;
		force = FALSE;
		//isotopes.clear();
	}
	const char* name;
	class phase* phase;
	int column;
	int constraint;
	int force;
	std::vector<class isotope> isotopes;
};
class inverse
{
public:
	~inverse() {};
	inverse()
	{
		n_user = -1;
		description = NULL;
		new_def = FALSE;
		minimal = FALSE;
		range = FALSE;
		mp = FALSE;
		mp_censor = 1e-20;
		range_max = 1000.0;
		tolerance = 1e-10;
		mp_tolerance = 1e-12;
		//uncertainties.clear();
		//ph_uncertainties.clear();
		water_uncertainty = 0.0;
		mineral_water = TRUE;
		carbon = TRUE;
		//dalk_dph.clear();
		//dalk_dc.clear();
		count_solns = 0;
		//solns.clear();
		//force_solns.clear();
		//elts.clear();
		//phases.clear();
		count_redox_rxns = 0;
		//isotopes.clear();
		//i_u.clear();
		//isotope_unknowns.clear();
		netpath = NULL;
		pat = NULL;
	}
	int n_user;
	char* description;
	int new_def;
	int minimal;
	int range;
	int mp;
	LDBLE mp_censor;
	LDBLE range_max;
	LDBLE tolerance;
	LDBLE mp_tolerance;
	std::vector<double> uncertainties;
	std::vector<double> ph_uncertainties;
	LDBLE water_uncertainty;
	int mineral_water;
	int carbon;
	std::vector<double> dalk_dph;
	std::vector<double> dalk_dc;
	size_t count_solns;
	std::vector<int> solns;
	std::vector<bool> force_solns;
	std::vector<class inv_elts> elts;
	std::vector<class inv_phases> phases;
	size_t count_redox_rxns;
	std::vector<class inv_isotope> isotopes;
	std::vector<class inv_isotope> i_u;
	std::vector<class isotope> isotope_unknowns;
	const char* netpath;
	const char* pat;
};
/*----------------------------------------------------------------------
 *   Jacobian and Mass balance lists
 *---------------------------------------------------------------------- */
class Model
{
public:
	Model()
	{
		force_prep = true;
		gas_phase_type = cxxGasPhase::GP_UNKNOWN;
		numerical_fixed_volume = false;
		dl_type = cxxSurface::NO_DL;
		surface_type = cxxSurface::UNKNOWN_DL;
	};
	~Model()
	{
	};
	bool force_prep;
	bool numerical_fixed_volume;
	cxxGasPhase::GP_TYPE gas_phase_type;
	std::vector<class phase*> gas_phase;
	std::vector<const char*> ss_assemblage;
	std::vector<class phase*> pp_assemblage;
	std::vector<double> si;
	std::vector<const char*> add_formula;
	cxxSurface::DIFFUSE_LAYER_TYPE dl_type;
	cxxSurface::SURFACE_TYPE surface_type;
	std::vector<const char*> surface_comp;
	std::vector<const char*> surface_charge;
};
class name_coef
{
public:
	~name_coef() {};
	name_coef()
	{
		name = NULL;
		coef = 0;
	}
	const char* name;
	LDBLE coef;
};
/*----------------------------------------------------------------------
 *   Species_list
 *---------------------------------------------------------------------- */
class species_list
{
public:
	~species_list() {};
	species_list()
	{
		master_s = NULL;
		s = NULL;
		coef = 0;
	}
	class species* master_s;
	class species* s;
	LDBLE coef;
};
/*----------------------------------------------------------------------
 *   Jacobian and Mass balance lists
 *---------------------------------------------------------------------- */
class list0
{
public:
	~list0() {};
	list0()
	{
		target = NULL;
		coef = 0;
	}
	LDBLE* target;
	LDBLE coef;
};
class list1
{
public:
	~list1() {};
	list1()
	{
		source = NULL;
		target = NULL;
	}
	LDBLE* source;
	LDBLE* target;
};
class list2
{
public:
	~list2() {};
	list2()
	{
		source = NULL;
		target = NULL;
		coef = 0;
	}
	LDBLE* source;
	LDBLE* target;
	LDBLE coef;
};
class iso
{
public:
	~iso() {};
	iso()
	{
		name = NULL;
		value = 0;
		uncertainty = 0.05;
	}
	const char* name;
	LDBLE value;
	LDBLE uncertainty;
};
/*----------------------------------------------------------------------
 *   Transport data
 *---------------------------------------------------------------------- */
class stag_data
{
public:
	~stag_data() {};
	stag_data()
	{
		count_stag = 0;
		exch_f = 0;
		th_m = 0;
		th_im = 0;
	}
	int count_stag;
	LDBLE exch_f;
	LDBLE th_m;
	LDBLE th_im;
};
class cell_data
{
public:
	~cell_data() {};
	cell_data()
	{
		length = 1;
		mid_cell_x = 1.;
		disp = 1.0;
		temp = 25.;
		// free (uncharged) porewater porosities
		por = 0.1;
		// interlayer water porosities
		por_il = 0.01;
		// potential (V)
		potV = 0;
		punch = FALSE;
		print = FALSE;
		same_model = FALSE;
	}
	LDBLE length;
	LDBLE mid_cell_x;
	LDBLE disp;
	LDBLE temp;
	LDBLE por;
	LDBLE por_il;
	LDBLE potV;
	int punch;
	int print;
	int same_model;
};
/*----------------------------------------------------------------------
 *   Keywords
 *---------------------------------------------------------------------- */
//class key
//{
//public:
//	~key() {};
//	key()
//	{
//		name = NULL;
//		keycount = 0;
//	}
//	char* name;
//	int keycount = 0;
// };
//class const_key
//{
//public:
//	~const_key() {};
//	const_key()
//	{
//		name = NULL;
//		keycount = 0;
//	}
//	const char* name;
//	int keycount;
//};
/*----------------------------------------------------------------------
 *   Elements
 *---------------------------------------------------------------------- */
class element
{
public:
	~element() {};
	element()
	{
		// element name
		name = NULL;
		/*    int in; */
		master = NULL;
		primary = NULL;
		gfw = 0;
	}
	const char* name;
	class master* master;
	class master* primary;
	LDBLE gfw;
};
/*----------------------------------------------------------------------
 *   Element List
 *---------------------------------------------------------------------- */
class elt_list
{
public:
	~elt_list() {};
	elt_list()
	{	/* list of name and number of elements in an equation */
		elt = NULL;	/* pointer to element structure */
		coef = 0.0;			/* number of element e's in eqn */
	}
	class element* elt;
	LDBLE coef;
};
/*----------------------------------------------------------------------
 *   Species
 *---------------------------------------------------------------------- */
class species
{
public:
	~species() {};
	species()
	{	/* all data pertinent to an aqueous species */
		name = NULL;          // name of species 
		mole_balance = NULL;  // formula for mole balance 
		in = FALSE;           // set internally if species in model
		number = 0;
		// points to master species list, NULL if not primary master
		primary = NULL;
		// points to master species list, NULL if not secondary master
		secondary = NULL;
		gfw = 0;              // gram formula wt of species
		z = 0;                // charge of species
		dw = 0;		// tracer diffusion coefficient in water at 25oC, m2/s
		dw_t = 0;	// correct Dw for temperature: Dw(TK) = Dw(298.15) * exp(dw_t / TK - dw_t / 298.15)
		// parms for calc'ng SC = SC0 * exp(-dw_a * z * mu^0.5 / (1 + DH_B * dw_a2 * mu^0.5) / (1 + mu^dw_a3))
		// with DHO: ka = DH_B * dw_a * (1 + DD(V_apparent)^dw_a2 * sqrt_mu, dw_a3 is a switch, see calc_SC in PBasic
		dw_a = 0;
		dw_a2 = 0;
		dw_a3 = 0;
		dw_a_visc = 0;   // exponent in viscosity correction of SC
		dw_a_v_dif = 0;  // exponent in viscosity correction of D, the diffusion coefficient of the species
		dw_t_SC = 0;     // contribution to SC, for calc'ng transport number with BASIC
		dw_t_visc = 0;   // contribution to viscosity
		dw_corr = 0;	 // dw corrected for mu and TK
		erm_ddl = 0;     // enrichment factor in DDL
		equiv = 0;       // equivalents in exchange species
		alk = 0;	     // alkalinity of species, used for cec in exchange
		carbon = 0;      // stoichiometric coefficient of carbon in species
		co2 = 0;         // stoichiometric coefficient of C(4) in species
		h = 0;           // stoichiometric coefficient of H in species
		// stoichiometric coefficient of O in species
		o = 0;
		// WATEQ Debye Huckel a and b-dot; active_fraction coef for exchange species
		dha = 0, dhb = 0, a_f = 0;
		lk = 0;           // log10 k at working temperature
		// log kt0, delh, 6 coefficients analytical expression + volume terms
		for (size_t i = 0; i < MAX_LOG_K_INDICES; i++) logk[i] = 0;
		// 7 coefficients analytical expression for B, D, anion terms and pressure in Jones_Dole viscosity eqn
		for (size_t i = 0; i < 10; i++) Jones_Dole[i] = 0;
		// regression coefficients to calculate temperature dependent phi_0and b_v of Millero density model
		for (size_t i = 0; i < 7; i++) millero[i] = 0;
		original_units = kjoules;  // enum with original delta H units
		//add_logk.clear();
		lg = 0;            // log10 activity coefficient, gamma
		lg_pitzer = 0;     // log10 activity coefficient, from pitzer calculation
		lm = 0;            // log10 molality
		la = 0;		       // log10 activity
		dg = 0;		       // gamma term for jacobian
		dg_total_g = 0;
		moles = 0;		   // moles in solution; moles/mass_water = molality
		type = 0;          // flag indicating presence in model and types of equations
		gflag = 0;		   // flag for preferred activity coef eqn
		exch_gflag = 0;    // flag for preferred activity coef eqn
		// vector of elements
		//next_elt.clear();
		//next_secondary.clear();
		//next_sys_total.clear();
		// switch to check equation for charge and element balance
		check_equation = TRUE;
		//CReaction rxn;   // data base reaction
		//CReaction rxn_s; // reaction converted to secondary and primary master species
		//CReaction rxn_x; // reaction to be used in model
		// (1 + sum(g)) * moles
		tot_g_moles = 0;
		// sum(moles*g*Ws/Waq)
		tot_dh2o_moles = 0;
		for (size_t i = 0; i < 5; i++) cd_music[i] = 0;
		for (size_t i = 0; i < 3; i++) dz[i] = 0;
		original_deltav_units = cm3_per_mol;
	}
	const char* name;
	const char* mole_balance; 
	int in;
	int number;
	class master* primary;
	class master* secondary;
	LDBLE gfw;
	LDBLE z;
	LDBLE dw;
	LDBLE dw_t;
	LDBLE dw_a;
	LDBLE dw_a2;
	LDBLE dw_a3;
	LDBLE dw_a_visc;
	LDBLE dw_a_v_dif;
	LDBLE dw_t_SC;
	LDBLE dw_t_visc;
	LDBLE dw_corr;
	LDBLE erm_ddl;
	LDBLE equiv;
	LDBLE alk;
	LDBLE carbon;
	LDBLE co2;
	LDBLE h;
	LDBLE o;	
	LDBLE dha, dhb, a_f;
	LDBLE lk;
	LDBLE logk[MAX_LOG_K_INDICES];
	LDBLE Jones_Dole[10];
	LDBLE millero[7];
	DELTA_H_UNIT original_units;
	std::vector<class name_coef> add_logk;
	LDBLE lg;
	LDBLE lg_pitzer;
	LDBLE lm;
	LDBLE la;
	LDBLE dg;
	LDBLE dg_total_g;
	LDBLE moles;
	int type;
	int gflag;
	int exch_gflag;
	std::vector<class elt_list> next_elt;
	std::vector<class elt_list> next_secondary;
	std::vector<class elt_list> next_sys_total;
	int check_equation;
	CReaction rxn;
	CReaction rxn_s;
	CReaction rxn_x;
	LDBLE tot_g_moles;
	LDBLE tot_dh2o_moles;
	LDBLE cd_music[5];
	LDBLE dz[3];
	DELTA_V_UNIT original_deltav_units;
};
class logk
{
public:
	~logk() {};
	logk()
	{	/* Named log K's */
		name = NULL;		 // name of species 
		lk = 0.0;	         // log10 k at working temperature                   
		// log kt0, delh, 6 coefficients analalytical expression 
		for (size_t i = 0; i < MAX_LOG_K_INDICES; i++) log_k[i] = 0;
		// enum with original delta H units 
		original_units = kjoules;	   
		done = FALSE;
		//add_logk.clear();
		// log kt0, delh, 5 coefficients analalytical expression 
		for (size_t i = 0; i < MAX_LOG_K_INDICES; i++) log_k_original[i] = 0;
		original_deltav_units = cm3_per_mol;
	}
	const char* name;
	LDBLE lk;
	LDBLE log_k[MAX_LOG_K_INDICES];
	DELTA_H_UNIT original_units;
	int done;
	std::vector<class name_coef> add_logk;
	LDBLE log_k_original[MAX_LOG_K_INDICES];
	DELTA_V_UNIT original_deltav_units;
};
/*----------------------------------------------------------------------
 *   Phases
 *---------------------------------------------------------------------- */
class phase
{
public:
	~phase() {};
	phase()
	{	/* all data pertinent to a pure solid phase */
		name = NULL;                //name of species 
		formula = NULL;				// chemical formula 
		in = FALSE;					// species used in model if TRUE 
		lk = 0;					    // log10 k at working temperature 
		// log kt0, delh, 6 coefficients analalytical expression
		for (size_t i = 0; i < MAX_LOG_K_INDICES; i++) logk[i] = 0;
		// enum with original delta H units 
		original_units = kjoules;	
		original_deltav_units = cm3_per_mol;
		//add_logk.clear();
		moles_x = 0;
		delta_max = 0;
		p_soln_x = 0;
		fraction_x = 0;
		log10_lambda = 0;
		log10_fraction_x = 0;
		dn = 0, dnb = 0, dnc = 0;
		gn = 0, gntot = 0;
		gn_n = 0, gntot_n = 0;
		// gas: critical TK, critical P(atm), Pitzer acentric coeff 
		t_c = 0, p_c = 0, omega = 0;  
		// Peng-Robinson parm's
		pr_a = 0, pr_b = 0, pr_alpha = 0;	
		// Temperature (K), Pressure (atm)
		pr_tk = 0, pr_p = 0;
		// fugacity coefficient (-) 
		pr_phi = 0;		
		// for calculating multicomponent phi 
		pr_aa_sum2 = 0;
		// delta_v[0] = [1] + [2]*T + [3]/T + [4]*log10(T) + [5]/T^2 + [6]*T^2 + [7]*P
		for (size_t i = 0; i < 9; i++) delta_v[i] = 0;
		// si adapter: log10(phi) - delta_v[0] * (P - 1) /RT
		pr_si_f = 0;	
		// Peng-Robinson in the calc's, or not
		pr_in = false;
		// flag indicating presence in model and types of equations
		type = SOLID;	
		// list of elements in phase 
		//next_elt.clear();	    
		//next_sys_total.clear();
		// switch to check equation for charge and element balance
		check_equation = TRUE;
		// data base reaction
		//CReaction rxn;
		// reaction converted to secondary and primary master species
		//CReaction rxn_s;
		// reaction to be used in model
		//CReaction rxn_x;
		// equation contains solids or gases                    
		replaced = FALSE;                       
		in_system = FALSE;
	}
	const char* name;
	const char* formula;
	int in;
	LDBLE lk;
	LDBLE logk[MAX_LOG_K_INDICES];
	DELTA_H_UNIT original_units;
	DELTA_V_UNIT original_deltav_units;
	std::vector<class name_coef> add_logk;
	LDBLE moles_x;
	LDBLE delta_max;
	LDBLE p_soln_x;
	LDBLE fraction_x;
	LDBLE log10_lambda;
	LDBLE log10_fraction_x;
	LDBLE dn, dnb, dnc;
	LDBLE gn, gntot;
	LDBLE gn_n, gntot_n;
	LDBLE t_c, p_c, omega; 
	LDBLE pr_a, pr_b, pr_alpha;
	LDBLE pr_tk, pr_p;
	LDBLE pr_phi;
	LDBLE pr_aa_sum2;
	LDBLE delta_v[9];
	LDBLE pr_si_f;
	bool pr_in;
	int type;	
	std::vector<class elt_list> next_elt;
	std::vector<class elt_list> next_sys_total;
	int check_equation;
	CReaction rxn;
	CReaction rxn_s;
	CReaction rxn_x;
	int replaced;
	int in_system;
};
/*----------------------------------------------------------------------
 *   Master species
 *---------------------------------------------------------------------- */
class master
{
public:
	~master() {};
	master()
	{	
		// TRUE if in model, FALSE if out, REWRITE if other mb eq
		in = 0;
		// sequence number in list of masters
		number = 0;
		// saved to determine if model has changed
		last_model = FALSE;
		// AQ or EX
		type = 0;
		// TRUE if master species is primary
		primary = FALSE;
		// coefficient of element in master species
		coef = 0;
		// total concentration for element or valence state
		total = 0;
		isotope_ratio = 0;
		isotope_ratio_uncertainty = 0;
		isotope = 0;
		total_primary = 0;
		// element structure
		elt = NULL;
		// alkalinity of species
		alk = 0;
		// default gfw for species
		gfw = 1;
		// formula from which to calculate gfw
		gfw_formula = NULL;
		// pointer to unknown structure
		unknown = NULL;
		// pointer to species structure
		s = NULL;
		// reaction writes master species in terms of primary  master species
		//CReaction rxn_primary;
		// reaction writes master species in terms of secondary master species
		//CReaction rxn_secondary;
		pe_rxn = NULL;
		minor_isotope = FALSE;
	}
	int in;
	size_t number;
	int last_model;
	int type;
	int primary;
	LDBLE coef;
	LDBLE total;
	LDBLE isotope_ratio;
	LDBLE isotope_ratio_uncertainty;
	int isotope;
	LDBLE total_primary;
	class element* elt;
	LDBLE alk;
	LDBLE gfw;
	const char* gfw_formula;
	class unknown* unknown;
	class species* s;
	CReaction rxn_primary;
	CReaction rxn_secondary;
	const char* pe_rxn;
	int minor_isotope;
};
/*----------------------------------------------------------------------
 *   Unknowns
 *---------------------------------------------------------------------- */
class unknown
{
public:
	~unknown() {};
	unknown()
	{
		type = 0;
		moles = 0;
		ln_moles = 0;
		f = 0;
		sum = 0;
		delta = 0;
		la = 0;
		number = 0;
		description = NULL;
		//master.clear();
		phase = NULL;
		si = 0;
		n_gas_phase_user = 0;
		s = NULL;
		exch_comp = NULL;
		pp_assemblage_comp_name = NULL;
		pp_assemblage_comp_ptr = NULL;
		ss_name = NULL;
		ss_ptr = NULL;
		ss_comp_name = NULL;
		ss_comp_ptr = NULL;
		ss_comp_number = 0;
		ss_in = FALSE;
		surface_comp = NULL;
		surface_charge = NULL;
		related_moles = 0;
		potential_unknown = NULL;
		potential_unknown1 = NULL;
		potential_unknown2 = NULL;
		// list for CD_MUSIC of comps that contribute to 0 plane mass-balance term
		//comp_unknowns.clear();
		phase_unknown = NULL;
		mass_water = 1;
		dissolve_only = FALSE;
		inert_moles = 0;
		V_m = 0;
		pressure = 1;
		mb_number = 0;
		iteration = 0;
	}
	int type;
	LDBLE moles;
	LDBLE ln_moles;
	LDBLE f;
	LDBLE sum;
	LDBLE delta;
	LDBLE la;
	size_t number;
	const char* description;
	std::vector<class master*> master;
	class phase* phase;
	LDBLE si;
	int n_gas_phase_user;
	class species* s;
	const char* exch_comp;
	const char* pp_assemblage_comp_name;
	void* pp_assemblage_comp_ptr;
	const char* ss_name;
	void* ss_ptr;
	const char* ss_comp_name;
	void* ss_comp_ptr;
	int ss_comp_number;
	int ss_in;
	const char* surface_comp;
	const char* surface_charge;
	LDBLE related_moles;
	class unknown* potential_unknown;
	class unknown* potential_unknown1;
	class unknown* potential_unknown2;
	std::vector<class unknown*> comp_unknowns;
	class unknown* phase_unknown;
	LDBLE mass_water;
	int dissolve_only;
	LDBLE inert_moles;
	LDBLE V_m;
	LDBLE pressure;
	int mb_number;
	int iteration;
};
/*----------------------------------------------------------------------
 *   Reaction work space
 *---------------------------------------------------------------------- */
class rxn_token_temp
{
public:
	~rxn_token_temp() {};
	rxn_token_temp()
	{	// data for equations, aq. species or minerals
		name = NULL;		// pointer to a species name (formula)
		z = 0;		// charge on species 
		s = NULL;
		unknown = NULL;
		coef = 0;			// coefficient of species name 
	}
	const char* name;
	LDBLE z;
	class species* s;
	class unknown* unknown;
	LDBLE coef;
};
class reaction_temp
{
public:
	~reaction_temp() {};
	reaction_temp()
	{
		for (size_t i = 0; i < MAX_LOG_K_INDICES; i++) logk[i] = 0;
		for (size_t i = 0; i < 3; i++) dz[i] = 0;
		//token.clear();
	}
	LDBLE logk[MAX_LOG_K_INDICES];
	LDBLE dz[3];
	std::vector<class rxn_token_temp> token;
};
class unknown_list
{
public:
	~unknown_list() {};
	unknown_list()
	{
		unknown = NULL;
		source = NULL;
		gamma_source = NULL;
		coef = 0;
	}
	class unknown* unknown;
	LDBLE* source;
	LDBLE* gamma_source;
	LDBLE coef;
};
/* ----------------------------------------------------------------------
 *   Print
 * ---------------------------------------------------------------------- */
class prints
{
public:
	~prints() {};
	prints()
	{
		all = 0;
		initial_solutions = 0;
		initial_exchangers = 0;
		reactions = 0;
		gas_phase = 0;
		ss_assemblage = 0;
		pp_assemblage = 0;
		surface = 0;
		exchange = 0;
		kinetics = 0;
		totals = 0;
		eh = 0;
		species = 0;
		saturation_indices = 0;
		irrev = 0;
		mix = 0;
		reaction = 0;
		use = 0;
		logfile = 0;
		punch = 0;
		status = 0;
		inverse = 0;
		dump = 0;
		user_print = 0;
		headings = 0;
		user_graph = 0;
		echo_input = 0;
		warnings = 0;
		initial_isotopes = 0;
		isotope_ratios = 0;
		isotope_alphas = 0;
		hdf = 0;
		alkalinity = 0;
	}
	int all;
	int initial_solutions;
	int initial_exchangers;
	int reactions;
	int gas_phase;
	int ss_assemblage;
	int pp_assemblage;
	int surface;
	int exchange;
	int kinetics;
	int totals;
	int eh;
	int species;
	int saturation_indices;
	int irrev;
	int mix;
	int reaction;
	int use;
	int logfile;
	int punch;
	int status;
	int inverse;
	int dump;
	int user_print;
	int headings;
	int user_graph;
	int echo_input;
	int warnings;
	int initial_isotopes;
	int isotope_ratios;
	int isotope_alphas;
	int hdf;
	int alkalinity;
};
/* ----------------------------------------------------------------------
 *   RATES
 * ---------------------------------------------------------------------- */
class rate
{
public:
	~rate() {};
	rate()
	{
		name = NULL;
		//std::string commands;
		new_def = 0;
		linebase = NULL;
		varbase = NULL;
		loopbase = NULL;
	}
	const char* name;
	std::string commands;
	int new_def;
	void* linebase;
	void* varbase;
	void* loopbase;
};
/* ----------------------------------------------------------------------
 *   GLOBAL DECLARATIONS
 * ---------------------------------------------------------------------- */
class spread_row
{
public:
	~spread_row() {};
	spread_row()
	{
		count = 0;
		empty = 0, string = 0, number = 0;
		//char_vector.clear();
		//d_vector.clear();
		//type_vector.clear();
	}
	size_t count;
	size_t empty, string, number;
	std::vector<std::string> str_vector;
	std::vector<int> type_vector;
};
class defaults
{
public:
	~defaults() {};
	defaults()
	{
		temp = 25;
		density = 1;
		calc_density = false;
		units = NULL;
		redox = NULL;
		ph = 7;
		pe = 4;
		water = 1;
		//iso.clear();
		pressure = 1;	/* pressure in atm */
	}
	LDBLE temp;
	LDBLE density;
	bool calc_density;
	const char* units;
	const char* redox;
	LDBLE ph;
	LDBLE pe;
	LDBLE water;
	std::vector<class iso> iso;
	LDBLE pressure;
};
class spread_sheet
{
public:
	~spread_sheet() {};
	spread_sheet()
	{
		heading = NULL;
		units = NULL;
		//class defaults defaults;
	}
	class spread_row* heading;
	class spread_row* units;
	std::vector<class spread_row*> rows;
	class defaults defaults;
};
/* ----------------------------------------------------------------------
 *   ISOTOPES
 * ---------------------------------------------------------------------- */
class master_isotope
{
public:
	~master_isotope() {};
	master_isotope()
	{
		name = NULL;
		master = NULL;
		elt = NULL;
		units = NULL;
		standard = 0;
		ratio = 0;
		moles = 0;
		total_is_major = 0;
		minor_isotope = 0;
	}
	const char* name;
	class master* master;
	class element* elt;
	const char* units;
	LDBLE standard;
	LDBLE ratio;
	LDBLE moles;
	int total_is_major;
	int minor_isotope;
};
class calculate_value
{
public:
	~calculate_value() {};
	calculate_value()
	{
		name = NULL;
		value = 0;
		//commands.clear();
		new_def = 0;
		calculated = 0;
		linebase = NULL;
		varbase = NULL;
		loopbase = NULL;
	}
	const char* name;
	LDBLE value;
	std::string commands;
	int new_def;
	int calculated;
	void* linebase;
	void* varbase;
	void* loopbase;
};
class isotope_ratio
{
public:
	isotope_ratio()
	{
		name = NULL;
		isotope_name = NULL;
		ratio = 0;
		converted_ratio = 0;
	}
	~isotope_ratio() {};

	const char* name;
	const char* isotope_name;
	LDBLE ratio;
	LDBLE converted_ratio;
};
class isotope_alpha
{
public:
	isotope_alpha()
	{
		name = NULL;
		named_logk = NULL;
		value = 0;
	}
	~isotope_alpha() {};
	const char* name;
	const char* named_logk;
	LDBLE value;
};
class system_species
{
public:
	~system_species() {};
	system_species()
	{
		name = NULL;
		type = NULL;
		moles = 0;
	}
	char* name;
	char* type;
	LDBLE moles;
};
/* tally.c ------------------------------- */
class tally_buffer
{
public:
	~tally_buffer() {};
	tally_buffer()
	{
		name = NULL;
		master = NULL;
		moles = 0;
		gfw = 0;
	}
	const char* name;
	class master* master;
	LDBLE moles;
	LDBLE gfw;
};
class tally
{
public:
	~tally() {};
	tally()
	{
		name = NULL;
		type = UnKnown;
		add_formula = NULL;
		moles = 0;
		//formula.clear();
		/*
		 * first total is initial
		 * second total is final
		 * third total is difference (final - initial)
		 */
		for(size_t i = 0; i < 3; i++) total[i]= NULL;
	}
	const char* name;
	enum entity_type type;
	const char* add_formula;
	LDBLE moles;
	std::vector<class elt_list> formula;
	/*
	 * first total is initial
	 * second total is final
	 * third total is difference (final - initial)
	 */
	class tally_buffer* total[3];
};
/* transport.c ------------------------------- */
class spec
{
public:
	~spec() {};
	spec()
	{
		// name of species
		name = NULL;
		// name of aqueous species in EX species
		aq_name = NULL;
		// type: AQ or EX
		type = 0;
		// activity
		a = 0;
		// log(concentration)
		lm = 0;
		// log(gamma)
		lg = 0;
		// concentration for AQ, equivalent fraction for EX
		c = 0;
		// charge number
		z = 0;
		// free water diffusion coefficient, m2/s
		Dw = 0;
		// temperature and viscosity corrected free water diffusion coefficient, m2/s
		Dwt = 0;
		// temperature factor for Dw
		dw_t = 0;
		// viscosity factor for Dw
		dw_a_v_dif = 0;
		// enrichment factor in ddl
		erm_ddl = 0;
	}
	const char* name;
	const char* aq_name;
	int type;
	LDBLE a;
	LDBLE lm;
	LDBLE lg;
	LDBLE c;
	LDBLE z;
	LDBLE Dw;
	LDBLE Dwt;
	LDBLE dw_t;
	LDBLE dw_a_v_dif;
	LDBLE erm_ddl;
};

class sol_D
{
public:
	~sol_D() {};
	sol_D()
	{
		// number of aqueous + exchange species
		count_spec = 0;
		// number of exchange species
		count_exch_spec = 0;
		// total moles of X-, max X- in transport step in sol_D[1], tk
		exch_total = 0, x_max = 0, tk_x = 0;
		// viscos_0 at I = 0
		viscos_0 = 0;
		// viscosity of solution
		viscos = 0;
		spec = NULL;
		spec_size = 0;
	}
	int count_spec;
	int count_exch_spec;
	LDBLE exch_total, x_max, tk_x;
	LDBLE viscos_0, viscos;
	class spec* spec;
	int spec_size;
};
class J_ij
{
public:
	~J_ij() {};
	J_ij()
	{
		name = NULL;
		// species change in cells i and j
		tot1 = 0;
		tot2 = 0;
		tot_stag = 0;
		charge = 0;
	}
	const char* name;
	LDBLE tot1, tot2, tot_stag, charge;
};
class J_ij_save
{
public:
	~J_ij_save() {};
	J_ij_save()
	{
		// species change in cells i and j
		flux_t = 0;
		flux_c = 0;
	}
	double flux_t, flux_c;
};
class M_S
{
public:
	~M_S() {};
	M_S()
	{
		name = NULL;
		// master species transport in cells i and j 
		tot1 = 0;
		tot2 = 0;
		tot_stag = 0;
		charge = 0;
	}
	const char* name;
	LDBLE tot1, tot2, tot_stag, charge;
};
// Pitzer definitions
typedef enum
{ TYPE_B0, TYPE_B1, TYPE_B2, TYPE_C0, TYPE_THETA, TYPE_LAMBDA, TYPE_ZETA,
  TYPE_PSI, TYPE_ETHETA, TYPE_ALPHAS, TYPE_MU, TYPE_ETA, TYPE_Other,
  TYPE_SIT_EPSILON, TYPE_SIT_EPSILON_MU, TYPE_APHI
} pitz_param_type;
class pitz_param
{
public:
	~pitz_param() {};
	pitz_param()
	{
		for(size_t i = 0; i < 3; i++) species[i] = NULL;
		for (size_t i = 0; i < 3; i++) ispec[i] = -1;
		type = TYPE_Other;
		p = 0;
		U.b0 = 0;
		for (size_t i = 0; i < 6; i++) a[i] = 0;
		alpha = 0;
		os_coef = 0;
		for (size_t i = 0; i < 3; i++) ln_coef[i] = 0;
		thetas = NULL;
	}
	const char* species[3];
	int ispec[3];
	pitz_param_type type;
	LDBLE p;
	union
	{
		LDBLE b0;
		LDBLE b1;
		LDBLE b2;
		LDBLE c0;
		LDBLE theta;
		LDBLE lambda;
		LDBLE zeta;
		LDBLE psi;
		LDBLE alphas;
		LDBLE mu;
		LDBLE eta;
		LDBLE eps;
		LDBLE eps1;
		LDBLE aphi;
	} U;
	LDBLE a[6];
	LDBLE alpha;
	LDBLE os_coef;
	LDBLE ln_coef[3];
	class theta_param* thetas;
};
class theta_param
{
public:
	~theta_param() {};
	theta_param()
	{
		zj = 0;
		zk = 0;
		etheta = 0;
		ethetap = 0;
	}
	LDBLE zj;
	LDBLE zk;
	LDBLE etheta;
	LDBLE ethetap;
};
class const_iso
{
public:
	~const_iso() {};
	const_iso()
	{
		name = NULL;
		value = 0;
		uncertainty = 0;
	}
	const_iso(const char *n, LDBLE v, LDBLE u)
	{
		name = n;
		value = v;
		uncertainty = u;
	}
	const char* name;
	LDBLE value;
	LDBLE uncertainty;
};

#endif /* _INC_GLOBAL_STRUCTURES_H  */
