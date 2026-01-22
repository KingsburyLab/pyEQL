#ifndef _INC_PHREEQC_H
#define _INC_PHREEQC_H
#if defined(WIN32)
#  if defined(PHREEQCI_GUI)
#    ifndef WINVER
#      define WINVER 0x0400
#    endif
#    include <afx.h>
#  endif
#  include <windows.h>
#  if defined(PHREEQCI_GUI)
#    include "../../resource.h"
#  endif
#endif
#if defined(WIN32_MEMORY_DEBUG)
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif
#ifndef boolean
typedef unsigned char boolean;
#endif
/* ----------------------------------------------------------------------
*   INCLUDE FILES
* ---------------------------------------------------------------------- */
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <cmath>
#include <errno.h>
#include <cfloat>
#include "phrqtype.h"
#include "cvdense.h"	
#include "runner.h"
#include "dumper.h"
#include "PHRQ_io.h"
#include "SelectedOutput.h"
#include "UserPunch.h"
#ifdef MULTICHART
#include "ChartHandler.h"
#endif
#include "Keywords.h"
#include "Pressure.h"
#include "cxxMix.h"
#include "Use.h"
#include "Surface.h"
#ifdef SWIG_SHARED_OBJ
#include "thread.h"
#endif

class cxxNameDouble;
class cxxKinetics;
class cxxKineticsComp;
class cxxExchange;
class cxxExchComp;
class cxxGasPhase;
class cxxTemperature;
class cxxPPassemblage;
class cxxPPassemblageComp;
class cxxReaction;
class cxxSolution;
class cxxISolutionComp;
class cxxSolutionIsotope;
class cxxSSassemblage;
class cxxSS;
class cxxStorageBin;
class PBasic;

#include "global_structures.h"

class Phreeqc
{
public:
	Phreeqc(PHRQ_io* io = NULL);
	Phreeqc(const Phreeqc& src);
	void InternalCopy(const Phreeqc* pSrc);
	Phreeqc& operator=(const Phreeqc& rhs);
	~Phreeqc(void);

public:
	//
	// Phreeqc class methods
	//

	// advection.cpp -------------------------------
	int advection(void);

	// basicsubs.cpp -------------------------------
	int basic_compile(const char* commands, void** lnbase, void** vbase, void** lpbase);
	int basic_run(char* commands, void* lnbase, void* vbase, void* lpbase);
	void basic_free(void);
#ifdef IPHREEQC_NO_FORTRAN_MODULE
	double basic_callback(double x1, double x2, const char* str);
#else
	double basic_callback(double x1, double x2, const char* str);
#endif
	void register_basic_callback(double (*fcn)(double x1, double x2, const char* str, void* cookie), void* cookie1);
#ifdef IPHREEQC_NO_FORTRAN_MODULE
	void register_fortran_basic_callback(double (*fcn)(double* x1, double* x2, const char* str, size_t l));
#else
	void register_fortran_basic_callback(double (*fcn)(double* x1, double* x2, const char* str, int l));
#endif

	LDBLE activity(const char* species_name);
	LDBLE activity_coefficient(const char* species_name);
	LDBLE log_activity_coefficient(const char* species_name);
	LDBLE aqueous_vm(const char* species_name);
	LDBLE phase_vm(const char* phase_name);
	LDBLE diff_c(const char* species_name);
	LDBLE setdiff_c(const char * species_name, double d, double d_v_d);
	LDBLE flux_mcd(const char* species_name, int option);
	LDBLE sa_declercq(double type, double sa, double d, double m, double m0, double gfw);
	LDBLE calc_SC(void);
	/* VP: Density Start */
	LDBLE calc_dens(void);
	/* VP: Density End */
	LDBLE calc_logk_n(const char* name);
	LDBLE calc_logk_p(const char* name);
	LDBLE calc_logk_s(const char* name);
	LDBLE calc_deltah_s(const char* name);
	LDBLE calc_deltah_p(const char* name);
	LDBLE dh_a0(const char* name);
	LDBLE dh_bdot(const char* name);
	LDBLE calc_surface_charge(const char* surface_name);
	LDBLE calc_t_sc(const char* name);
	LDBLE diff_layer_total(const char* total_name, const char* surface_name);
	LDBLE edl_species(const char* surf_name, LDBLE* count, char*** names, LDBLE** moles, LDBLE* area, LDBLE* thickness);
	int get_edl_species(cxxSurfaceCharge& charge_ref);
	LDBLE equi_phase(const char* phase_name);
	LDBLE equi_phase_delta(const char* phase_name);
	LDBLE equivalent_fraction(const char* name, LDBLE* eq, std::string& elt_name);
	LDBLE find_gas_comp(const char* gas_comp_name);
	LDBLE find_gas_p(void);
	LDBLE find_gas_vm(void);
	LDBLE find_misc1(const char* ss_name);
	LDBLE find_misc2(const char* ss_name);
	LDBLE find_ss_comp(const char* ss_comp_name);
	LDBLE get_calculate_value(const char* name);
	char* iso_unit(const char* total_name);
	LDBLE iso_value(const char* total_name);
	LDBLE kinetics_moles(const char* kinetics_name);
	LDBLE kinetics_moles_delta(const char* kinetics_name);
	LDBLE log_activity(const char* species_name);
	LDBLE log_molality(const char* species_name);
	LDBLE molality(const char* species_name);
	LDBLE pressure(void);
	LDBLE pr_pressure(const char* phase_name);
	LDBLE pr_phi(const char* phase_name);
	LDBLE saturation_ratio(const char* phase_name);
	int saturation_index(const char* phase_name, LDBLE* iap, LDBLE* si);
	int solution_number(void);
	LDBLE solution_sum_secondary(const char* total_name);
	LDBLE sum_match_gases(const char* stemplate, const char* name);
	LDBLE sum_match_species(const char* stemplate, const char* name);
	LDBLE sum_match_ss(const char* stemplate, const char* name);
	int match_elts_in_species(const char* name, const char* stemplate);
	int extract_bracket(const char** string, char* bracket_string);
	LDBLE surf_total(const char* total_name, const char* surface_name);
	LDBLE surf_total_no_redox(const char* total_name, const char* surface_name);
	static int system_species_compare(const void* ptr1, const void* ptr2);
	static int system_species_compare_name(const void* ptr1, const void* ptr2);
	LDBLE system_total(const char* total_name, LDBLE* count, char*** names,
		char*** types, LDBLE** moles, int i);
	std::string kinetics_formula(std::string kinetics_name, cxxNameDouble& stoichiometry);
	std::string phase_formula(std::string phase_name, cxxNameDouble& stoichiometry);
	std::string species_formula(std::string phase_name, cxxNameDouble& stoichiometry);
	std::string phase_equation(std::string phase_name, std::vector<std::pair<std::string, double> >& stoichiometry);
	std::string species_equation(std::string species_name, std::vector<std::pair<std::string, double> >& stoichiometry);
	LDBLE list_ss(std::string ss_name, cxxNameDouble& composition);
	int system_total_elements(void);
	int system_total_si(void);
	int system_total_aq(void);
	int system_total_ex(void);
	int system_total_surf(void);
	int system_total_gas(void);
	int system_total_equi(void);
	int system_total_kin(void);
	int system_total_ss(void);
	int system_total_elt(const char* total_name);
	int system_total_elt_secondary(const char* total_name);
	LDBLE total(const char* total_name);
	LDBLE total_mole(const char* total_name);
	int system_total_solids(cxxExchange* exchange_ptr,
		cxxPPassemblage* pp_assemblage_ptr,
		cxxGasPhase* gas_phase_ptr,
		cxxSSassemblage* ss_assemblage_ptr,
		cxxSurface* surface_ptr);

	static LDBLE f_rho(LDBLE rho_old, void* cookie);
	static LDBLE f_Vm(LDBLE v1, void* cookie);
	LDBLE calc_solution_volume(void);

	// cl1.cpp -------------------------------
	int cl1(int k, int l, int m, int n,
		int nklmd, int n2d,
		LDBLE* q,
		int* kode, LDBLE toler,
		int* iter, LDBLE* x, LDBLE* res, LDBLE* error,
		LDBLE* cu, int* iu, int* s, int check);
	void cl1_space(int check, int n2d, int klm, int nklmd);

	// cl1mp.cpp -------------------------------
	int cl1mp(int k, int l, int m, int n,
		int nklmd, int n2d,
		LDBLE* q_arg,
		int* kode, LDBLE toler,
		int* iter, LDBLE* x_arg, LDBLE* res_arg, LDBLE* error,
		LDBLE* cu_arg, int* iu, int* s, int check, LDBLE censor_arg);

	// class_main.cpp -------------------------------
	int write_banner(void);

	/* default.cpp */
public:
	//int close_input_files(void);
	//int close_output_files(void);
	//static int istream_getc(void* cookie);
	int process_file_names(int argc, char* argv[], std::istream** db_cookie,
		std::istream** input_cookie, int log);

	/* PHRQ_io_output.cpp */
	void screen_msg(const char* str);

	void echo_msg(const char* err_str);
	int warning_msg(const char* err_str);
	void set_forward_output_to_log(int value);
	int get_forward_output_to_log(void);

	// dump_ostream
	bool dump_open(const char* file_name);
	void dump_flush(void);
	void dump_close(void);
	void dump_msg(const char* str);

	// log_ostream
	bool log_open(const char* file_name);
	void log_flush(void);
	void log_close(void);
	void log_msg(const char* str);

	// error_ostream
	bool error_open(const char* file_name);
	void error_flush(void);
	void error_close(void);
	void error_msg(const char* str, bool stop = false);

	// output_ostream
	bool output_open(const char* file_name);
	void output_flush(void);
	void output_close(void);
	void output_msg(const char* str);

	// punch_ostream
	bool punch_open(const char* file_name, int n_user);
	void punch_flush(void);
	void punch_close(void);
	void punch_msg(const char* str);

	void fpunchf_heading(const char* name);
	void fpunchf(const char* name, const char* format, double d);
	void fpunchf(const char* name, const char* format, char* d);
	void fpunchf(const char* name, const char* format, int d);
	void fpunchf_user(int user_index, const char* format, double d);
	void fpunchf_user(int user_index, const char* format, char* d);
	int fpunchf_end_row(const char* format);
	// input.cpp -------------------------------
	int reading_database(void);
	void set_reading_database(int reading_database);
	int check_line(const char* string, int allow_empty, int allow_eof,
		int allow_keyword, int print);
	int check_line_impl(const char* string, int allow_empty,
		int allow_eof, int allow_keyword, int print);
	int get_line(void);
	//int get_logical_line(void* cookie, int* l);
	int read_database(void);
	int run_simulations(void);

	// integrate.cpp -------------------------------
	int calc_all_g(void);
	int calc_init_g(void);
	int initial_surface_water(void);
	int sum_diffuse_layer(cxxSurfaceCharge* surface_charge_ptr1);
	int calc_all_donnan(void);
	int calc_init_donnan(void);
	LDBLE calc_psi_avg(cxxSurfaceCharge * charge_ptr, LDBLE surf_chrg_eq, LDBLE nDbl, LDBLE f_free, std::vector<LDBLE> &zcorr);
	LDBLE g_function(LDBLE x_value);
	LDBLE midpnt(LDBLE x1, LDBLE x2, int n);
	void polint(LDBLE* xa, LDBLE* ya, int n, LDBLE xv, LDBLE* yv,
		LDBLE* dy);
	LDBLE qromb_midpnt(cxxSurfaceCharge* charge_ptr, LDBLE x1, LDBLE x2);

	// inverse.cpp -------------------------------
	int inverse_models(void);
	int add_to_file(const char* filename, const char* string);
	int bit_print(unsigned long bits, int l);
	int carbon_derivs(class inverse* inv_ptr);
	int check_isotopes(class inverse* inv_ptr);
	int check_solns(class inverse* inv_ptr);
	bool set_isotope_unknowns(class inverse* inv_ptrs);
	cxxSolutionIsotope* get_isotope(cxxSolution* solution_ptr, const char* elt);
	LDBLE get_inv_total(cxxSolution* solution_ptr, const char* elt);
	int isotope_balance_equation(class inverse* inv_ptr, int row, int n);
	int post_mortem(void);
	bool test_cl1_solution(void);
	unsigned long get_bits(unsigned long bits, int position, int number);
	unsigned long minimal_solve(class inverse* inv_ptr,
		unsigned long minimal_bits);
	void dump_netpath(class inverse* inv_ptr);
	int dump_netpath_pat(class inverse* inv_ptr);
	int next_set_phases(class inverse* inv_ptr, int first_of_model_size,
		int model_size);
	int phase_isotope_inequalities(class inverse* inv_ptr);
	int print_model(class inverse* inv_ptr);
	int punch_model_heading(class inverse* inv_ptr);
	int punch_model(class inverse* inv_ptr);
	void print_isotope(FILE* netpath_file, cxxSolution* solution_ptr,
		const char* elt, const char* string);
	void print_total(FILE* netpath_file, cxxSolution* solution_ptr,
		const char* elt, const char* string);
	void print_total_multi(FILE* netpath_file, cxxSolution* solution_ptr,
		const char* string, const char* elt0,
		const char* elt1, const char* elt2, const char* elt3,
		const char* elt4);

	void print_total_pat(FILE* netpath_file, const char* elt,
		const char* string);
	int range(class inverse* inv_ptr, unsigned long cur_bits);
	int save_bad(unsigned long bits);
	int save_good(unsigned long bits);
	int save_minimal(unsigned long bits);
	unsigned long set_bit(unsigned long bits, int position, int value);
	int setup_inverse(class inverse* inv_ptr);
	int set_initial_solution(int n_user_old, int n_user_new);
	int set_ph_c(class inverse* inv_ptr,
		int i, cxxSolution* soln_ptr_orig, int n_user_new,
		LDBLE d_alk, LDBLE ph_factor, LDBLE alk_factor);
	int shrink(class inverse* inv_ptr, LDBLE* array_in,
		LDBLE* array_out, int* k, int* l, int* m, int* n,
		unsigned long cur_bits, LDBLE* delta_l, int* col_back_l,
		int* row_back_l);
	int solve_inverse(class inverse* inv_ptr);
	int solve_with_mask(class inverse* inv_ptr, unsigned long cur_bits);
	int subset_bad(unsigned long bits);
	int subset_minimal(unsigned long bits);
	int superset_minimal(unsigned long bits);
	int write_optimize_names(class inverse* inv_ptr);

	// isotopes.cpp -------------------------------
	int add_isotopes(cxxSolution& solution_ptr);
	int calculate_values(void);
	int calculate_isotope_moles(class element* elt_ptr,
		cxxSolution* solution_ptr, LDBLE total_moles);
	LDBLE convert_isotope(class master_isotope* master_isotope_ptr, LDBLE ratio);
	int from_pcil(class master_isotope* master_isotope_ptr);
	int from_permil(class master_isotope* master_isotope_ptr, LDBLE major_total);
	int from_pct(class master_isotope* master_isotope_ptr, LDBLE major_total);
	int from_tu(class master_isotope* master_isotope_ptr);
	class calculate_value* calculate_value_alloc(void);
	int calculate_value_free(class calculate_value* calculate_value_ptr);
	class calculate_value* calculate_value_search(const char* name);
	class calculate_value* calculate_value_store(const char* name,
		int replace_if_found);
	class isotope_alpha* isotope_alpha_alloc(void);
	class isotope_alpha* isotope_alpha_search(const char* name);
	class isotope_alpha* isotope_alpha_store(const char* name,
		int replace_if_found);
	class isotope_ratio* isotope_ratio_alloc(void);
	class isotope_ratio* isotope_ratio_search(const char* name);
	class isotope_ratio* isotope_ratio_store(const char* name,
		int replace_if_found);
	class master_isotope* master_isotope_store(const char* name,
		int replace_if_found);
	class master_isotope* master_isotope_alloc(void);
	class master_isotope* master_isotope_search(const char* name);
	int print_initial_solution_isotopes(void);
	int print_isotope_ratios(void);
	int print_isotope_alphas(void);
	int punch_isotopes(void);
	int punch_calculate_values(void);
	int read_calculate_values(void);
	int read_isotopes(void);
	int read_isotope_ratios(void);
	int read_isotope_alphas(void);
	int calculate_value_init(class calculate_value* calculate_value_ptr);
	int isotope_alpha_init(class isotope_alpha* isotope_alpha_ptr);
	int isotope_ratio_init(class isotope_ratio* isotope_ratio_ptr);
	int master_isotope_init(class master_isotope* master_isotope_ptr);

	// kinetics.cpp -------------------------------
	void cvode_init(void);
	bool cvode_update_reactants(int i, int nsaver, bool save_it);
	int run_reactions(int i, LDBLE kin_time, int use_mix, LDBLE step_fraction);
	int set_and_run(int i, int use_mix, int use_kinetics, int nsaver,
		LDBLE step_fraction);
	int set_and_run_wrapper(int i, int use_mix, int use_kinetics, int nsaver,
		LDBLE step_fraction);
	int set_advection(int i, int use_mix, int use_kinetics, int nsaver);
	int free_cvode(void);
public:
	static void f(integertype N, realtype t, N_Vector y, N_Vector ydot,
		void* f_data);
	static void Jac(integertype N, DenseMat J, RhsFn f, void* f_data, realtype t,
		N_Vector y, N_Vector fy, N_Vector ewt, realtype h,
		realtype uround, void* jac_data, long int* nfePtr,
		N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

	int calc_final_kinetic_reaction(cxxKinetics* kinetics_ptr);
	int calc_kinetic_reaction(cxxKinetics* kinetics_ptr,
		LDBLE time_step);
	bool limit_rates(cxxKinetics* kinetics_ptr);
	int rk_kinetics(int i, LDBLE kin_time, int use_mix, int nsaver,
		LDBLE step_fraction);
	int set_reaction(int i, int use_mix, int use_kinetics);
	int set_transport(int i, int use_mix, int use_kinetics, int nsaver);
	int store_get_equi_reactants(int k, int kin_end);

	// mainsubs.cpp  -------------------------------
	std::ifstream* open_input_stream(std::string query, std::string& default_name, std::ios_base::openmode mode, bool batch);
	std::ofstream* open_output_stream(std::string query, std::string& default_name, std::ios_base::openmode mode, bool batch);
	int copy_entities(void);
	void do_mixes(void);
	void initialize(void);
	int initial_exchangers(int print);
	int initial_gas_phases(int print);
	int initial_solutions(int print);
	int step_save_exch(int n_user);
	int step_save_surf(int n_user);
	int initial_surfaces(int print);
	int reactions(void);
	int saver(void);
	int xsolution_save(int k_user);
	int xexchange_save(int n_user);
	int xgas_save(int n_user);
	int xpp_assemblage_save(int n_user);
	int xss_assemblage_save(int n_user);
	int xsurface_save(int n_user);
	int do_initialize(void);
	int do_status(void);
	void save_init(int i);
	int copy_use(int i);
	int set_use(void);

	// model.cpp -------------------------------
	int check_residuals(void);
	int free_model_allocs(void);
	int ineq(int kode);
	int model(void);
	int jacobian_sums(void);
	int mb_gases(void);
	int mb_ss(void);
	int mb_sums(void);
	int molalities(int allow_overflow);
	int reset(void);
	int residuals(void);
	int set(int initial);
	int sum_species(void);
	int surface_model(void);
	LDBLE ss_root(LDBLE a0, LDBLE a1, LDBLE kc, LDBLE kb, LDBLE xcaq,
		LDBLE xbaq);
	LDBLE ss_halve(LDBLE a0, LDBLE a1, LDBLE x0, LDBLE x1, LDBLE kc,
		LDBLE kb, LDBLE xcaq, LDBLE xbaq);
	LDBLE ss_f(LDBLE xb, LDBLE a0, LDBLE a1, LDBLE kc, LDBLE kb,
		LDBLE xcaq, LDBLE xbaq);
	int numerical_jacobian(void);
	void set_inert_moles(void);
	void unset_inert_moles(void);
#ifdef SLNQ
	int add_trivial_eqns(int rows, int cols, LDBLE* matrix);
	//int slnq(int n, LDBLE * a, LDBLE * delta, int ncols, int print);
#endif
	int calc_gas_pressures(void);
	int calc_fixed_volume_gas_pressures(void);
	double calc_gas_binary_parameter(std::string name1, std::string name2) const;
	int calc_ss_fractions(void);
	int gammas(LDBLE mu);
	int gammas_a_f(int i);
	int initial_guesses(void);
	int revise_guesses(void);
	int ss_binary(cxxSS* ss_ptr);
	int ss_ideal(cxxSS* ss_ptr);

	// parse.cpp -------------------------------
	int check_eqn(int association);
	int get_charge(char* charge, size_t charge_size, LDBLE* z);
	int get_elt(const char** t_ptr, std::string& element, int* i);
	int get_elts_in_species(const char** t_ptr, LDBLE coef);
	int get_num(const char** t_ptr, LDBLE* num);
	int get_secondary_in_species(const char** t_ptr, LDBLE coef);
	int parse_eq(char* eqn, std::vector<class elt_list>& new_elt_list, int association);
	int get_coef(LDBLE* coef, const char** eqnaddr);
	int get_secondary(const char** t_ptr, char* element, int* i);
	int get_species(const char** ptr);

	// phqalloc.cpp -------------------------------
public:
#if !defined(NDEBUG)
	void* PHRQ_malloc(size_t, const char*, int);
	void* PHRQ_calloc(size_t, size_t, const char*, int);
	void* PHRQ_realloc(void*, size_t, const char*, int);
#else
	void* PHRQ_malloc(size_t);
	void* PHRQ_calloc(size_t, size_t);
	void* PHRQ_realloc(void*, size_t);
#endif
	void PHRQ_free(void* ptr);
	void PHRQ_free_all(void);

public:

	// pitzer.cpp -------------------------------
	class pitz_param* pitz_param_read(char* string, int n);
	void pitz_param_store(class pitz_param* pzp_ptr);
	void sit_param_store(class pitz_param* pzp_ptr);
	class pitz_param* pitz_param_copy(const class pitz_param* src);
	class theta_param* theta_param_search(LDBLE zj, LDBLE zk);
	void pitzer_make_lists(void);
	int gammas_pz(bool exch_a_f);
	int model_pz(void);
	int pitzer(void);
	int pitzer_clean_up(void);
	int pitzer_init(void);
	int pitzer_tidy(void);
	int read_pitzer(void);
	int set_pz(int initial);
	int calc_pitz_param(class pitz_param* pz_ptr, LDBLE TK, LDBLE TR);
	int check_gammas_pz(void);
	int ISPEC(const char* name);
	LDBLE G(LDBLE Y);
	LDBLE GP(LDBLE Y);
	int ETHETAS(LDBLE ZJ, LDBLE ZK, LDBLE I, LDBLE* etheta,
		LDBLE* ethetap);
	void ETHETA_PARAMS(LDBLE X, LDBLE& JAY, LDBLE& JPRIME);
	int pitzer_initial_guesses(void);
	int pitzer_revise_guesses(void);
	int PTEMP(LDBLE TK);
	int jacobian_pz(void);

	// prep.cpp -------------------------------
	int add_potential_factor(void);
	int add_cd_music_factors(int n);
	int add_surface_charge_balance(void);
	int add_cd_music_charge_balances(int i);
	int build_gas_phase(void);
	int build_fixed_volume_gas(void);
	int build_jacobian_sums(int k);
	int build_mb_sums(void);
	int build_min_exch(void);
	int build_model(void);
	int build_pure_phases(void);
	int build_ss_assemblage(void);
	int build_solution_phase_boundaries(void);
	int build_species_list(int n);
	int build_min_surface(void);
	LDBLE calc_lk_phase(phase* p_ptr, LDBLE TK, LDBLE pa);
	LDBLE calc_PR(std::vector<class phase*> phase_ptrs, LDBLE P, LDBLE TK, LDBLE V_m);
	LDBLE calc_PR();
	int calc_vm(LDBLE tc, LDBLE pa);
	LDBLE calc_vm0(const char *species_name, LDBLE tc, LDBLE pa, LDBLE mu);
	int clear(void);
	int convert_units(cxxSolution* solution_ptr);
	class unknown* find_surface_charge_unknown(std::string& str_ptr, int plane);
	std::vector<class master*> get_list_master_ptrs(const char* cptr, class master* master_ptr);
	int inout(void);
	int is_special(class species* spec);
	int mb_for_species_aq(int n);
	int mb_for_species_ex(int n);
	int mb_for_species_surf(int n);
	int quick_setup(void);
	int resetup_master(void);
	int save_model(void);
	int setup_exchange(void);
	int setup_gas_phase(void);
	int setup_fixed_volume_gas(void);
	int setup_master_rxn(const std::vector<class master*>& master_ptr_list,
		const std::string& pe_rxn);
	int setup_pure_phases(void);
	int adjust_setup_pure_phases(void);
	int setup_related_surface(void);
	int setup_ss_assemblage(void);
	int setup_solution(void);
	int adjust_setup_solution(void);
	int setup_surface(void);
	int setup_unknowns(void);
	int store_dn(int k, LDBLE* source, int row, LDBLE coef_in,
		LDBLE* gamma_source);
	int store_jacob(LDBLE* source, LDBLE* target, LDBLE coef);
	int store_jacob0(int row, int column, LDBLE coef);
	int store_mb(LDBLE* source, LDBLE* target, LDBLE coef);
	int store_mb_unknowns(class unknown* unknown_ptr, LDBLE* LDBLE_ptr,
		LDBLE coef, LDBLE* gamma_ptr);
	int store_sum_deltas(LDBLE* source, LDBLE* target, LDBLE coef);
	int tidy_redox(void);
	int write_mb_eqn_x(void);
	int write_mb_for_species_list(int n);
	int write_mass_action_eqn_x(int stop);

	int check_same_model(void);
	int k_temp(LDBLE tc, LDBLE pa);
	LDBLE k_calc(LDBLE* logk, LDBLE tempk, LDBLE presPa);
	int prep(void);
	int reprep(void);
	int rewrite_master_to_secondary(class master* master_ptr1,
		class master* master_ptr2);
	int switch_bases(void);
	int write_phase_sys_total(int n);

	// print.cpp -------------------------------
	char* sformatf(const char* format, ...);
	int array_print(LDBLE* array_l, int row_count, int column_count,
		int max_column_count);
	int set_pr_in_false(void);
	int print_all(void);
	int print_exchange(void);
	int print_gas_phase(void);
	int print_master_reactions(void);
	int print_species(void);
	int print_surface(void);
	int print_user_print(void);
	int punch_all(void);
	int print_alkalinity(void);
	int print_diffuse_layer(cxxSurfaceCharge* surface_charge_ptr);
	int print_eh(void);
	int print_reaction(void);
	int print_kinetics(void);
	int print_mix(void);
	int print_pp_assemblage(void);
	int print_ss_assemblage(void);
	int print_saturation_indices(void);
	int print_surface_cd_music(void);
	int print_totals(void);
	int print_using(void);
	int punch_gas_phase(void);
	int punch_identifiers(void);
	int punch_kinetics(void);
	int punch_molalities(void);
	int punch_activities(void);
	int punch_pp_assemblage(void);
	int punch_ss_assemblage(void);
	int punch_saturation_indices(void);
	int punch_totals(void);
	int punch_user_punch(void);
#if defined MULTICHART
	int punch_user_graph(void);
#endif

	// read.cpp -------------------------------
	int read_input(void);
	int* read_list_ints_range(const char** ptr, int* count_ints, int positive,
		int* int_list);
	int read_list_ints_range(const char** cptr, bool positive, std::vector<int>& int_list);

	int read_log_k_only(const char* cptr, LDBLE* log_k);
	int read_t_c_only(const char* cptr, LDBLE* t_c);
	int read_p_c_only(const char* cptr, LDBLE* p_c);
	int read_omega_only(const char* cptr, LDBLE* omega);
	int read_number_description(const char* cptr, int* n_user, int* n_user_end,
		char** description, int allow_negative = FALSE);
	int check_key(const char* str);
	int check_units(std::string& tot_units, bool alkalinity, bool check_compatibility,
		const char* default_units, bool print);
	int find_option(const char* item, int* n, const char** list, int count_list,
		int exact);
	int get_option(const char** opt_list, int count_opt_list, const char** next_char);
	int get_true_false(const char* string, int default_value);

	int add_psi_master_species(char* token);
	int read_advection(void);
	int read_analytical_expression_only(const char* cptr, LDBLE* log_k);
	/* VP: Density Start */
	int read_millero_abcdef(const char* cptr, LDBLE* abcdef);
	/* VP: Density End */
	int read_viscosity_parms(const char* cptr, LDBLE* Jones_Dole);
	int read_copy(void);
	int read_debug(void);
	int read_delta_h_only(const char* cptr, LDBLE* delta_h,
		DELTA_H_UNIT* units);
	int read_aq_species_vm_parms(const char* cptr, LDBLE* delta_v);
	int read_vm_only(const char* cptr, LDBLE* delta_v,
		DELTA_V_UNIT* units);
	int read_phase_vm(const char* cptr, LDBLE* delta_v,
		DELTA_V_UNIT* units);
	int read_llnl_aqueous_model_parameters(void);
	int read_exchange(void);
	int read_exchange_master_species(void);
	int read_exchange_species(void);
	int read_gas_phase(void);
	int read_incremental_reactions(void);
	int read_inverse(void);
	int read_inv_balances(class inverse* inverse_ptr, const char* next_char);
	int read_inv_isotopes(class inverse* inverse_ptr, const char* cptr);
	int read_inv_phases(class inverse* inverse_ptr, const char* next_char);
	int read_kinetics(void);
	bool read_vector_doubles(const char** ptr, std::vector<double>& v);
	bool read_vector_ints(const char** cptr, std::vector<int>& v, int positive);
	bool read_vector_t_f(const char** ptr, std::vector<bool>& v);
	int read_master_species(void);
	int read_rate_parameters_pk(void);
	int read_rate_parameters_svd(void);
	int read_rate_parameters_hermanska(void);
	int read_mean_gammas(void);
	int read_gas_binary_parameters(void);
	int read_mix(void);
	int read_entity_mix(std::map<int, cxxMix>& mix_map);
	//int read_solution_mix(void);
	int read_named_logk(void);
	int read_phases(void);
	int read_print(void);
	int read_pp_assemblage(void);
	int read_rates(void);
	int read_reaction(void);
	int read_reaction_reactants(cxxReaction* reaction_ptr);
	int read_reaction_steps(cxxReaction* reaction_ptr);
	int read_solid_solutions(void);
	int read_temperature(void);
	//int read_reaction_temps(struct temperature* temperature_ptr);
	int read_reaction_pressure(void);
	int read_reaction_pressure_raw(void);
	int read_save(void);
	int read_selected_output(void);
	int read_solution(void);
	int read_species(void);
	int read_surface(void);
	int read_surface_master_species(void);
	int read_surface_species(void);
	int read_use(void);
	int read_title(void);
	int read_user_print(void);
	int read_user_punch(void);
#if defined MULTICHART
	int read_user_graph_handler();
#endif
	int next_keyword_or_option(const char** opt_list, int count_opt_list);
	int cleanup_after_parser(CParser& parser);

	// ReadClass.cxx
	int read_dump(void);
	int read_delete(void);
	int read_run_cells(void);
	int streamify_to_next_keyword(std::istringstream& lines);
	int dump_entities(void);
	int delete_entities(void);
	int run_as_cells(void);
	void dump_ostream(std::ostream& os);

	// readtr.cpp -------------------------------
	int read_transport(void);
	int dump(void);
	//int dump_exchange(int k);
	//int dump_gas_phase(int k);
	//int dump_kinetics(int k);
	//int dump_mix(int k);
	//int dump_pp_assemblage(int k);
	//int dump_reaction(int k);
	//int dump_ss_assemblage(int k);
	//int dump_solution(int k);
	//int dump_surface(int k);
	int dump_cpp(void);
	int read_line_LDBLEs(const char* next_char, LDBLE** d, int* count_d,
		int* count_alloc);

	// sit.cpp -------------------------------
	int gammas_sit(void);
	int model_sit(void);
	int sit(void);
	int sit_clean_up(void);
	int sit_init(void);
	int sit_tidy(void);
	int read_sit(void);
	int set_sit(int initial);
	int calc_sit_param(class pitz_param* pz_ptr, LDBLE TK, LDBLE TR);
	int check_gammas_sit(void);
	int sit_ISPEC(const char* name);
	/*int DH_AB (LDBLE TK, LDBLE *A, LDBLE *B);*/
	int sit_initial_guesses(void);
	int sit_revise_guesses(void);
	int PTEMP_SIT(LDBLE tk);
	void sit_make_lists(void);
	int jacobian_sit(void);

	// spread.cpp -------------------------------
	int read_solution_spread(void);
	int copy_token_tab(std::string& token, const char** cptr);
	int get_option_string(const char** opt_list, int count_opt_list,
		const char** next_char);
	int spread_row_free(class spread_row* spread_row_ptr);
	int spread_row_to_solution(class spread_row* heading,
		class spread_row* units,
		class spread_row* data,
		class defaults defaults);
	class spread_row* string_to_spread_row(char* string);
#ifdef PHREEQCI_GUI
	void add_row(class spread_row* spread_row_ptr);
	void free_spread(void);
	class spread_row* copy_row(class spread_row* spread_row_ptr);
#endif

	// step.cpp -------------------------------
	int step(LDBLE step_fraction);
	int xsolution_zero(void);
	int add_exchange(cxxExchange* exchange_ptr);
	int add_gas_phase(cxxGasPhase* gas_phase_ptr);
	int add_kinetics(cxxKinetics* kinetics_ptr);
	int add_mix(cxxMix* mix_ptr);
	int add_pp_assemblage(cxxPPassemblage* pp_assemblage_ptr);
	int add_reaction(cxxReaction* reaction_ptr, int step_number, LDBLE step_fraction);
	int add_ss_assemblage(cxxSSassemblage* ss_assemblage_ptr);
	int add_solution(cxxSolution* solution_ptr, LDBLE extensive,
		LDBLE intensive);
	int add_surface(cxxSurface* surface_ptr);
	int check_pp_assemblage(cxxPPassemblage* pp_assemblage_ptr);
	int gas_phase_check(cxxGasPhase* gas_phase_ptr);
	int pp_assemblage_check(cxxPPassemblage* pp_assemblage_ptr);
	int reaction_calc(cxxReaction* reaction_ptr);
	int solution_check(void);
	int ss_assemblage_check(cxxSSassemblage* ss_assemblage_ptr);

	// structures.cpp -------------------------------
	int clean_up(void);
	int reinitialize(void);

	int copier_add(class copier* copier_ptr, int n_user, int start, int end);
	int copier_clear(class copier* copier_ptr);
	//
	CReaction CReaction_internal_copy(CReaction& rxn_ref);
	double rxn_find_coef(CReaction& r_ptr, const char* str);
	//
	static int element_compare(const void* ptr1, const void* ptr2);
	class element* element_store(const char* element);
	//
	int add_elt_list(const cxxNameDouble& nd, LDBLE coef);
	int add_elt_list(const std::vector<class elt_list>& el, double coef);
	int change_hydrogen_in_elt_list(LDBLE charge);
	int elt_list_combine(void);
	static int elt_list_compare(const void* ptr1, const void* ptr2);
	std::vector<class elt_list> elt_list_internal_copy(const std::vector<class elt_list>& el);
	std::vector<class elt_list> elt_list_vsave(void);
	cxxNameDouble elt_list_NameDouble(void);
	//
	enum entity_type get_entity_enum(char* name);
	//
	class inverse* inverse_alloc(void);
	int inverse_delete(int i);
	static int inverse_isotope_compare(const void* ptr1, const void* ptr2);
	class inverse* inverse_search(int n_user, int* n);
	int inverse_sort(void);
	//
	class logk* logk_alloc(void);
	int logk_copy2orig(class logk* logk_ptr);
	class logk* logk_store(const char* name, int replace_if_found);
	class logk* logk_search(const char* name);
	//
	class master* master_alloc(void);
	static int master_compare(const void* ptr1, const void* ptr2);
	int master_delete(const char* cptr);
	class master* master_bsearch(const char* cptr);
	class master* master_bsearch_primary(const char* cptr);
	class master* master_bsearch_secondary(const char* cptr);
	class master* master_search(const char* cptr, int* n);
	class master* surface_get_psi_master(const char* name, int plane);
	//
	class phase* phase_bsearch(const char* cptr, int* j, int print);
#ifdef OBSOLETE
	static int phase_compare(const void* ptr1, const void* ptr2);
#endif
	int phase_delete(int i);
	class phase* phase_store(const char* name);
	//
	class rate* rate_bsearch(const char* cptr, int* j);
	int rate_free(class rate* rate_ptr);
	class rate* rate_copy(const class rate* rate_ptr);
	class rate* rate_search(const char* name, int* n);
	int rate_sort(void);
	//
	static int s_compare(const void* ptr1, const void* ptr2);
	int s_delete(int i);
	class species* s_search(const char* name);
	class species* s_store(const char* name, LDBLE z, int replace_if_found);
	//
	static int isotope_compare(const void* ptr1, const void* ptr2);
	//
	static int species_list_compare_alk(const void* ptr1, const void* ptr2);
	static int species_list_compare_master(const void* ptr1, const void* ptr2);
	int species_list_sort(void);
	//
	struct Change_Surf* change_surf_alloc(int count);
	//
	int system_duplicate(int i, int save_old);
	//
	//
	bool phase_rxn_to_trxn(class phase* phase_ptr, CReaction& rxn_ptr);
	bool trxn_add(CReaction& r_ptr, double coef, bool combine);
	bool trxn_add_phase(CReaction& r_ref, double coef, bool combine);
	int trxn_combine(void);
	static int trxn_compare(const void* ptr1, const void* ptr2);
	bool trxn_copy(CReaction& rxn_ref);
	LDBLE trxn_find_coef(const char* str, int start);
	int trxn_multiply(LDBLE coef);
	int trxn_print(void);
	int trxn_reverse_k(void);
	int trxn_sort(void);
	int trxn_swap(const char* token);

	class unknown* unknown_alloc(void);
	int unknown_delete(int i);
	int unknown_free(class unknown* unknown_ptr);
	int entity_exists(const char* name, int n_user);
	static int inverse_compare(const void* ptr1, const void* ptr2);
	int inverse_free(class inverse* inverse_ptr);
	int logk_init(class logk* logk_ptr);
	static int master_compare_string(const void* ptr1, const void* ptr2);
	int master_free(class master* master_ptr);
	class phase* phase_alloc(void);
	static int phase_compare_string(const void* ptr1, const void* ptr2);
	int phase_free(class phase* phase_ptr);
	int phase_init(class phase* phase_ptr);
	static int rate_compare(const void* ptr1, const void* ptr2);
	static int rate_compare_string(const void* ptr1, const void* ptr2);
	class species* s_alloc(void);
	int s_free(class species* s_ptr);
	int s_init(class species* s_ptr);
	static int species_list_compare(const void* ptr1, const void* ptr2);

	void Use2cxxStorageBin(cxxStorageBin& sb);
	void phreeqc2cxxStorageBin(cxxStorageBin& sb);
	void phreeqc2cxxStorageBin(cxxStorageBin& sb, int n);
	void cxxStorageBin2phreeqc(cxxStorageBin& sb, int n);
	void cxxStorageBin2phreeqc(cxxStorageBin& sb);

	/* tally.cpp */
	void add_all_components_tally(void);
	int build_tally_table(void);
	int calc_dummy_kinetic_reaction_tally(cxxKinetics* kinetics_ptr);
	int diff_tally_table(void);
	int extend_tally_table(void);
	int free_tally_table(void);
	int fill_tally_table(int* n_user, int index_conservative, int n_buffer);
	int get_tally_table_rows_columns(int* rows, int* columns);
	int get_tally_table_column_heading(int column, int* type, char* string);
	int get_tally_table_row_heading(int column, char* string);
	int store_tally_table(LDBLE* array, int row_dim, int col_dim,
		LDBLE fill_factor);
	int zero_tally_table(void);
	int elt_list_to_tally_table(class tally_buffer* buffer_ptr);
	int master_to_tally_table(class tally_buffer* buffer_ptr);
	int get_all_components(void);
	int print_tally_table(void);
	int set_reaction_moles(int n_user, LDBLE moles);
	int set_reaction_temperature(int n_user, LDBLE tc);
	int set_kinetics_time(int n_user, LDBLE step);

	// tidy.cpp -------------------------------
	int add_other_logk(LDBLE* source_k, std::vector<class name_coef>& add_logk);
	int add_logks(class logk* logk_ptr, int repeats);
	LDBLE halve(LDBLE f(LDBLE x, void*), LDBLE x0, LDBLE x1, LDBLE tol);
	int replace_solids_gases(void);
	int ss_prep(LDBLE t, cxxSS* ss_ptr, int print);
	int select_log_k_expression(LDBLE* source_k, LDBLE* target_k);
	int slnq(int n, LDBLE* a, LDBLE* delta, int ncols, int print);
public:
	int tidy_punch(void);
	int tidy_model(void);
	int check_species_input(void);
	LDBLE coef_in_master(class master* master_ptr);
	int reset_last_model(void);
	int rewrite_eqn_to_primary(void);
	int rewrite_eqn_to_secondary(void);
	int species_rxn_to_trxn(class species* s_ptr);
	int tidy_logk(void);
	int tidy_exchange(void);
	int tidy_min_exchange(void);
	int update_min_exchange(void);
	int tidy_kin_exchange(void);
	int update_kin_exchange(void);
	int tidy_gas_phase(void);
	int tidy_inverse(void);
	int tidy_isotopes(void);
	int tidy_isotope_ratios(void);
	int tidy_isotope_alphas(void);
	int tidy_kin_surface(void);
	int update_kin_surface(void);
	int tidy_master_isotope(void);
	int tidy_min_surface(void);
	int update_min_surface(void);
	int tidy_phases(void);
	int tidy_pp_assemblage(void);
	int tidy_solutions(void);
	int tidy_ss_assemblage(void);
	int tidy_species(void);
	int tidy_surface(void);
	int scan(LDBLE f(LDBLE x, void*), LDBLE* xx0, LDBLE* xx1);
	static LDBLE f_spinodal(LDBLE x, void*);
	int solve_misc(LDBLE* xxc1, LDBLE* xxc2, LDBLE tol);
	int ss_calc_a0_a1(cxxSS* ss_ptr);

	// transport.cpp -------------------------------
	int transport(void);
	void print_punch(int i, boolean active);
	int set_initial_moles(int i);
	cxxSurface sum_surface_comp(cxxSurface* source1, LDBLE f1,
		cxxSurface* source2, std::string charge_name, LDBLE f2,
		LDBLE new_Dw);
	int reformat_surf(const char* comp_name, LDBLE fraction, const char* new_comp_name,
		LDBLE new_Dw, int cell);
	LDBLE viscosity(cxxSurface *surf_ptr);
	LDBLE calc_f_visc(const char *name);
	LDBLE calc_vm_Cl(void);
	int multi_D(LDBLE DDt, int mobile_cell, int stagnant);
	LDBLE find_J(int icell, int jcell, LDBLE mixf, LDBLE DDt, int stagnant);
	void calc_b_ij(int icell, int jcell, int k, LDBLE b_i, LDBLE b_j, LDBLE g_i, LDBLE g_j, LDBLE free_i, LDBLE free_j, int stagnant);
	void diffuse_implicit(LDBLE DDt, int stagnant);
	int fill_spec(int cell_no, int ref_cell);
	LDBLE moles_from_redox_states(cxxSolution* sptr, const char* name);
	LDBLE moles_from_donnan_layer(cxxSurface* sptr, const char* name, LDBLE moles_needed);
	LDBLE add_MCD_moles(LDBLE moles, LDBLE min_mol, int i, cxxSolution* sptr, const char* name);
	int fill_m_s(class J_ij* J_ij, int J_ij_count_spec, int i, int stagnant);
	static int sort_species_name(const void* ptr1, const void* ptr2);
	int disp_surf(LDBLE stagkin_time);
	int diff_stag_surf(int mobile_cell);
	int check_surfaces(cxxSurface* surface_ptr1, cxxSurface* surface_ptr2);
	cxxSurface mobile_surface_copy(cxxSurface* surface_old_ptr,
		int n_user_new,
		bool move_old);
	void transport_cleanup(void);
	int init_mix(void);
	int init_heat_mix(int nmix);
	int heat_mix(int heat_nmix);
	int mix_stag(int i, LDBLE stagkin_time, int punch,
		LDBLE step_fraction_kin);

	// utilities.cpp -------------------------------
public:
	double calc_alk(CReaction& rxn_ptr);
	double calc_delta_v(CReaction& r_ref, bool phase);
	LDBLE calc_dielectrics(LDBLE tc, LDBLE pa);
	LDBLE calc_rho_0(LDBLE tc, LDBLE pa);
	int compute_gfw(const char* string, LDBLE* gfw);
	static int copy_token(char* token_ptr, const char** ptr, int* length);
	static int copy_token(std::string& token, const char** ptr);
	int dup_print(const char* cptr, int emphasis);
	int equal(LDBLE a, LDBLE b, LDBLE eps);
	void* free_check_null(void* ptr);
	int get_token(const char** eqnaddr, std::string& string, LDBLE* z, int* l);
	int islegit(const char c);
	void malloc_error(void);
	int print_centered(const char* string);
	static int replace(const char* str1, const char* str2, char* str);
	static void replace(std::string &stds, const char* str1, const char* str2);
	static bool replace(const char* str1, const char* str2, std::string& str);
	static int strcmp_nocase(const char* str1, const char* str2);
	static int strcmp_nocase_arg1(const char* str1, const char* str2);
	static void str_tolower(std::string& name);
	void space(void** ptr, int i, int* max, int struct_size);
	void squeeze_white(char* s_l);
	int status(int count, const char* str, bool kinetics = false);
	void str_tolower(char* str);
	void str_toupper(char* str);
#if !defined(NDEBUG) && defined(WIN32_MEMORY_DEBUG)
	char* _string_duplicate(const char* token, const char* szFileName, int nLine);
#else
	char* string_duplicate(const char* token);
#endif
	const char* string_hsave(const char* str);
	void strings_map_clear();
protected:
	char* string_pad(const char* str, int i);
	static int string_trim(char* str);
	static int string_trim_right(char* str);
	static int string_trim_left(char* str);
	static void string_trim(std::string& str);
	static void string_trim_left(std::string& str);
	static void string_trim_right(std::string& str);
	static LDBLE under(LDBLE xval);
	int get_input_errors(void);
	int isamong(char c, const char* s_l);
public:
	int main_method(int argc, char* argv[]);
	void set_phast(int);
	int next_user_number(Keywords::KEYWORDS key);
	size_t list_components(std::list<std::string>& list_c);
	size_t list_EquilibriumPhases(std::list<std::string>& list_pp);
	size_t list_GasComponents(std::list<std::string>& list_gc);
	size_t list_KineticReactions(std::list<std::string>& list_kr);
	size_t list_SolidSolutions(std::list<std::string>& list_comps, std::list<std::string>& list_names);
	size_t list_Surfaces(std::list<std::string>& surftype, std::list<std::string>& surf);
	size_t list_Exchangers(std::list<std::string>& ex);
	PHRQ_io* Get_phrq_io(void) { return this->phrq_io; }
	void Set_run_cells_one_step(const bool tf) { this->run_cells_one_step = tf; }


	std::map<int, cxxSolution>& Get_Rxn_solution_map() { return this->Rxn_solution_map; }
	std::map<int, cxxExchange>& Get_Rxn_exchange_map() { return this->Rxn_exchange_map; }
	std::map<int, cxxGasPhase>& Get_Rxn_gas_phase_map() { return this->Rxn_gas_phase_map; }
	std::map<int, cxxKinetics>& Get_Rxn_kinetics_map() { return this->Rxn_kinetics_map; }
	std::map<int, cxxPPassemblage>& Get_Rxn_pp_assemblage_map() { return this->Rxn_pp_assemblage_map; }
	std::map<int, cxxSSassemblage>& Get_Rxn_ss_assemblage_map() { return this->Rxn_ss_assemblage_map; }
	std::map<int, cxxSurface>& Get_Rxn_surface_map() { return this->Rxn_surface_map; }
	std::map<int, cxxMix>& Get_Rxn_mix_map() { return this->Rxn_mix_map; }
	std::map<int, cxxReaction>& Get_Rxn_reaction_map() { return this->Rxn_reaction_map; }
	std::map<int, cxxTemperature>& Get_Rxn_temperature_map() { return this->Rxn_temperature_map; }
	std::map<int, cxxPressure>& Get_Rxn_pressure_map() { return this->Rxn_pressure_map; }

protected:
	void init(void);

	//
	//Data members
	//
protected:
	PHRQ_io* phrq_io;
	PHRQ_io ioInstance;
	int same_model;

	LDBLE current_tc;
	LDBLE current_pa;
	LDBLE current_mu;
	bool mu_terms_in_logk;

	/* ----------------------------------------------------------------------
	*   STRUCTURES
	* ---------------------------------------------------------------------- */

	Model last_model;
	//struct punch punch;
	bool high_precision;

	/* ----------------------------------------------------------------------
	*   Temperatures
	* ---------------------------------------------------------------------- */

	std::map<int, cxxTemperature> Rxn_temperature_map;

	/* ----------------------------------------------------------------------
	*   Pressures
	* ---------------------------------------------------------------------- */
	std::map<int, cxxPressure> Rxn_pressure_map;

	/* ----------------------------------------------------------------------
	*   Surface
	* --------------------------------------------------------------------- */

	int g_iterations;
	LDBLE G_TOL;
	std::map <int, cxxSurface> Rxn_surface_map;
	std::map <LDBLE, LDBLE> charge_group_map;
	int change_surf_count;
	struct Change_Surf* change_surf;

	/* ----------------------------------------------------------------------
	*   Exchange
	* ---------------------------------------------------------------------- */
	std::map<int, cxxExchange> Rxn_exchange_map;

	/* ----------------------------------------------------------------------
	*   Kinetics
	* ---------------------------------------------------------------------- */
	std::map<int, cxxKinetics> Rxn_kinetics_map;
	bool use_kinetics_limiter;

	/*----------------------------------------------------------------------
	*   Save
	*---------------------------------------------------------------------- */
	std::map<std::string, double> save_values;
	std::map<std::string, std::string> save_strings;
	class save save;

	/*----------------------------------------------------------------------
	*   Use
	*---------------------------------------------------------------------- */
	cxxUse use;

	/*----------------------------------------------------------------------
	*   Copy
	*---------------------------------------------------------------------- */
	class copier copy_solution;
	class copier copy_pp_assemblage;
	class copier copy_exchange;
	class copier copy_surface;
	class copier copy_ss_assemblage;
	class copier copy_gas_phase;
	class copier copy_kinetics;
	class copier copy_mix;
	class copier copy_reaction;
	class copier copy_temperature;
	class copier copy_pressure;

	/*----------------------------------------------------------------------
	*   Inverse
	*---------------------------------------------------------------------- */
	std::vector<class inverse> inverse;
	int count_inverse;
	/*----------------------------------------------------------------------
	*   Rates
	*---------------------------------------------------------------------- */
	std::map<std::string, std::vector<double> > rate_parameters_pk;
	std::map<std::string, std::vector<double> > rate_parameters_svd;
	std::map<std::string, std::vector<double> > rate_parameters_hermanska;
	/*----------------------------------------------------------------------
	*   Mean gammas
	*---------------------------------------------------------------------- */
	std::map<std::string, cxxNameDouble> mean_gammas;
	/*----------------------------------------------------------------------
	*   Mix
	*---------------------------------------------------------------------- */
	std::map<int, cxxMix> Rxn_mix_map;
	std::map<int, cxxMix> Dispersion_mix_map;
	std::map<int, cxxMix> Rxn_solution_mix_map;
	std::map<int, cxxMix> Rxn_exchange_mix_map;
	std::map<int, cxxMix> Rxn_gas_phase_mix_map;
	std::map<int, cxxMix> Rxn_kinetics_mix_map;
	std::map<int, cxxMix> Rxn_pp_assemblage_mix_map;
	std::map<int, cxxMix> Rxn_ss_assemblage_mix_map;
	std::map<int, cxxMix> Rxn_surface_mix_map;
	/*
	* List new definitions
	*/
	std::set<int> Rxn_new_exchange;
	std::set<int> Rxn_new_gas_phase;
	std::set<int> Rxn_new_kinetics;     // not used
	std::set<int> Rxn_new_mix;          // not used
	std::set<int> Rxn_new_pp_assemblage;
	std::set<int> Rxn_new_pressure;     // not used
	std::set<int> Rxn_new_reaction;     // not used
	std::set<int> Rxn_new_solution;
	std::set<int> Rxn_new_ss_assemblage;
	std::set<int> Rxn_new_surface;
	std::set<int> Rxn_new_temperature;  // not used
	/*----------------------------------------------------------------------
	*   Irreversible reaction
	*---------------------------------------------------------------------- */
	std::map<int, cxxReaction> Rxn_reaction_map;

	/*----------------------------------------------------------------------
	*   Gas phase
	*---------------------------------------------------------------------- */
	std::map<int, cxxGasPhase> Rxn_gas_phase_map;

	/*----------------------------------------------------------------------
	*   Solid solution
	*---------------------------------------------------------------------- */
	std::map<int, cxxSSassemblage> Rxn_ss_assemblage_map;

	/*----------------------------------------------------------------------
	*   Pure-phase assemblage
	*---------------------------------------------------------------------- */
	std::map<int, cxxPPassemblage> Rxn_pp_assemblage_map;

	/*----------------------------------------------------------------------
	*   Species_list
	*---------------------------------------------------------------------- */
	std::vector<class species_list> species_list;

	/*----------------------------------------------------------------------
	*   Jacobian and Mass balance lists
	*---------------------------------------------------------------------- */
	std::vector<class list0> sum_jacob0;	/* array of pointers to targets and coefficients for array */

	std::vector<class list1> sum_mb1; /* array of pointers to sources and targets for mass
										balance summations with coef = 1.0 */
	std::vector<class list1> sum_jacob1;	/* array of pointers to sources and targets for array
											equations with coef = 1.0 */
	std::vector<class list2> sum_mb2; /* array of coefficients and pointers to sources and
									   targets for mass balance summations with coef != 1.0 */
	std::vector<class list2> sum_jacob2; /* array of coefficients and pointers to sources and
										  targets, coef != 1.0 */
	std::vector<class list2> sum_delta; /* array of pointers to sources, targets and coefficients for
										 summing deltas for mass balance equations */
										 /*----------------------------------------------------------------------
										 *   Solution
										 *---------------------------------------------------------------------- */
	std::map<int, cxxSolution> Rxn_solution_map;
	std::vector<cxxSolution> unnumbered_solutions;
	bool save_species;

	/*----------------------------------------------------------------------
	*   Global solution
	*---------------------------------------------------------------------- */
	std::string title_x;
	std::string last_title_x;
	int new_x;
	std::string description_x;
	LDBLE tc_x;
	LDBLE tk_x;
	LDBLE patm_x;
	LDBLE last_patm_x;
	LDBLE potV_x;
	bool numerical_fixed_volume;
	bool force_numerical_fixed_volume;
	//bool switch_numerical;
	LDBLE ph_x;
	LDBLE solution_pe_x;
	LDBLE mu_x;
	LDBLE ah2o_x;
	LDBLE total_h_x;
	LDBLE total_o_x;
	LDBLE cb_x;
	LDBLE total_ions_x;
	LDBLE mass_water_aq_x;
	LDBLE mass_water_surfaces_x;
	LDBLE mass_water_bulk_x;
	std::string units_x;
	std::map < std::string, CReaction > pe_x;
	std::map<std::string, cxxSolutionIsotope> isotopes_x;
	std::string default_pe_x;
	cxxSurface::DIFFUSE_LAYER_TYPE dl_type_x;
	LDBLE total_carbon;
	LDBLE total_co2;
	LDBLE total_alkalinity;
	LDBLE gfw_water;
	LDBLE step_x;
	LDBLE kin_time_x;

	/*----------------------------------------------------------------------
	*   Transport data
	*---------------------------------------------------------------------- */
	int count_cells;
	int count_shifts;
	int ishift;
	int bcon_first;
	int bcon_last;
	int correct_disp;
	LDBLE tempr;
	LDBLE timest;
	int simul_tr;
	LDBLE diffc;
	LDBLE heat_diffc;
	int cell;
	LDBLE mcd_substeps;
	class stag_data stag_data;
	int print_modulus;
	int punch_modulus;
	int dump_in;
	int dump_modulus;
	int transport_warnings;
	std::vector<class cell_data> cell_data;
	int old_cells, max_cells, all_cells;
	int multi_Dflag;		/* signals calc'n of multicomponent diffusion */
	int interlayer_Dflag;	/* multicomponent diffusion and diffusion through interlayer porosity */
	int implicit;	    /* implicit calculation of diffusion */
	LDBLE max_mixf;     /* the maximum value of the implicit mixfactor = De * Dt / (Dx^2) */
	LDBLE min_dif_LM;    /* the minimal log10(molality) for including a species in multicomponent diffusion */
	LDBLE default_Dw;		/* default species diffusion coefficient in water at 25oC, m2/s */
	int correct_Dw;         /* if true, Dw is adapted in calc_SC */
	LDBLE multi_Dpor;		/* uniform porosity of free porewater in solid medium */
	LDBLE interlayer_Dpor;	/* uniform porosity of interlayer space of montmorillonite in solid medium */
	LDBLE multi_Dpor_lim;	/* limiting free porewater porosity where transport stops */
	LDBLE interlayer_Dpor_lim;	/* limiting interlayer porosity where transport stops */
	LDBLE multi_Dn;		/* exponent to calculate pore water diffusion coefficient,
						Dp = Dw * (multi_Dpor)^multi_Dn */
	LDBLE interlayer_tortf;	/* tortuosity_factor in interlayer porosity,
							Dpil = Dw / interlayer_tortf */

	int cell_no, mixrun;
	/*----------------------------------------------------------------------
	*   Advection data
	*---------------------------------------------------------------------- */
	int count_ad_cells;
	int count_ad_shifts;
	int print_ad_modulus;
	int punch_ad_modulus;
	std::vector<int> advection_print, advection_punch;
	LDBLE advection_kin_time;
	LDBLE advection_kin_time_defined;
	int advection_warnings;

	/*----------------------------------------------------------------------
	*   Tidy data
	*---------------------------------------------------------------------- */
	int new_model, new_exchange, new_pp_assemblage, new_surface,
		new_reaction, new_temperature, new_mix, new_solution, new_gas_phase,
		new_inverse, new_punch, new_ss_assemblage, new_kinetics, new_copy,
		new_pitzer;

	/*----------------------------------------------------------------------
	*   Elements
	*---------------------------------------------------------------------- */
	std::vector<class element*> elements;
	class element* element_h_one;

	/*----------------------------------------------------------------------
	*   Element List
	*---------------------------------------------------------------------- */
	std::vector<class elt_list> elt_list;
	size_t count_elts;		/* number of elements in elt_list = position of next */
	/*----------------------------------------------------------------------
	*   Reaction
	*---------------------------------------------------------------------- */
	bool run_cells_one_step;
	/*----------------------------------------------------------------------
	*   Species
	*---------------------------------------------------------------------- */
	std::vector<class logk*> logk;

	std::string moles_per_kilogram_string;

	std::vector<class species*> s;
	std::vector< std::map < std::string, cxxSpeciesDL > > s_diff_layer;
	std::vector<class species*> s_x;

	class species* s_h2o;
	class species* s_hplus;
	class species* s_h3oplus;
	class species* s_eminus;
	class species* s_co3;
	class species* s_h2;
	class species* s_o2;

	/*----------------------------------------------------------------------
	*   Phases
	*---------------------------------------------------------------------- */
	std::vector<class phase*> phases;

	/*----------------------------------------------------------------------
	*   Master species
	*---------------------------------------------------------------------- */
	std::vector<class master*> master;

	/*----------------------------------------------------------------------
	*   Unknowns
	*---------------------------------------------------------------------- */
	std::vector<class unknown*> x;
	size_t count_unknowns;
	size_t sit_aqueous_unknowns;
	size_t max_unknowns;

	class unknown* ah2o_unknown;
	class unknown* alkalinity_unknown;
	class unknown* carbon_unknown;
	class unknown* charge_balance_unknown;
	class unknown* exchange_unknown;
	class unknown* mass_hydrogen_unknown;
	class unknown* mass_oxygen_unknown;
	class unknown* mb_unknown;
	class unknown* mu_unknown;
	class unknown* pe_unknown;
	class unknown* ph_unknown;
	class unknown* pure_phase_unknown;
	class unknown* solution_phase_boundary_unknown;
	class unknown* surface_unknown;
	class unknown* gas_unknown;
	class unknown* ss_unknown;
	std::vector<class unknown*> gas_unknowns;

	/*----------------------------------------------------------------------
	*   Reaction work space
	*---------------------------------------------------------------------- */
	class reaction_temp trxn;	/* structure array of working space while reading equations
								species names are in "temp_strings" */
	size_t count_trxn;		        /* number of reactants in trxn = position of next */

	std::vector<class unknown_list> mb_unknowns;

	/* ----------------------------------------------------------------------
	*   Print
	* ---------------------------------------------------------------------- */
	class prints pr;
	bool status_on;
	clock_t status_interval;
	clock_t status_timer;
	std::string status_string;
	int count_warnings;

	/* ----------------------------------------------------------------------
	*   RATES
	* ---------------------------------------------------------------------- */
	std::vector<class rate> rates;
	LDBLE rate_m, rate_m0, rate_time, rate_kin_time, rate_sim_time_start,
		rate_sim_time_end, rate_sim_time, rate_moles, initial_total_time;
	std::vector<LDBLE> rate_p;
	int count_rate_p;

	/* ----------------------------------------------------------------------
	*   USER PRINT COMMANDS
	* ---------------------------------------------------------------------- */
	class rate* user_print;
	int n_user_punch_index;

	int fpunchf_user_s_warning;
	char fpunchf_user_buffer[80];

#if defined MULTICHART
	ChartHandler chart_handler;
public:
	ChartHandler& Get_chart_handler(void)
	{
		return chart_handler;
	}
	const ChartHandler& Get_chart_handler(void)const
	{
		return chart_handler;
	}
protected:
#endif

	/* ----------------------------------------------------------------------
	*   GLOBAL DECLARATIONS
	* ---------------------------------------------------------------------- */
	const char* error_string;
	int simulation;
	int state;
	int reaction_step;
	int transport_step;
	int transport_start;
	int advection_step;
	int stop_program;
	int incremental_reactions;

	double MIN_LM;
	double LOG_ZERO_MOLALITY;
	double MIN_TOTAL;
	double MIN_TOTAL_SS;
	double MIN_RELATED_SURFACE;
	double MIN_RELATED_LOG_ACTIVITY;

	int count_strings;
	int max_strings;

	std::vector<double> my_array, delta, residual;

	int input_error;

	Keywords::KEYWORDS next_keyword;
	int parse_error;
	int paren_count;
	int iterations;
	int gamma_iterations;
	size_t density_iterations;
	LDBLE kgw_kgs;
	int run_reactions_iterations;
	int overall_iterations;

	int max_line;
	char* line;
	char* line_save;

	LDBLE LOG_10;

	int debug_model;
	int debug_prep;
	int debug_mass_action;
	int debug_mass_balance;
	int debug_set;
	int debug_diffuse_layer;
	int debug_inverse;

	LDBLE inv_tol_default;
	int itmax;
	int max_tries;
	LDBLE ineq_tol;
	LDBLE convergence_tolerance;
	LDBLE step_size;
	LDBLE pe_step_size;
	LDBLE step_size_now;
	LDBLE pe_step_size_now;
	LDBLE pp_scale;
	LDBLE pp_column_scale;
	int diagonal_scale;	/* 0 not used, 1 used */
	int mass_water_switch;
	int delay_mass_water;
	int equi_delay;
	bool dampen_ah2o;
	LDBLE censor;
	int aqueous_only;
	int negative_concentrations;
	int calculating_deriv;
	int numerical_deriv;

	int count_total_steps;
	int phast;
	bool output_newline;
	inline void Set_output_newline(bool tf) { this->output_newline = tf; }
	inline bool Get_output_newline() { return this->output_newline; }
	double a_llnl, b_llnl, bdot_llnl;
	std::vector<double> llnl_temp, llnl_adh, llnl_bdh, llnl_bdot, llnl_co2_coefs;

	//char *selected_output_file_name;
	std::map<int, SelectedOutput> SelectedOutput_map;
	SelectedOutput* current_selected_output;

	std::map <int, UserPunch> UserPunch_map;
	UserPunch* current_user_punch;

	char* dump_file_name;
	int remove_unstable_phases;
	std::string screen_string;
#ifdef PHREEQCI_GUI
	class spread_sheet g_spread_sheet;
#endif
	int spread_length;

	/* ---------------------------------------------------------------------- */
	/*
	*   Map definitions
	*/

	std::map<std::string, std::string*> strings_map;
	std::map<std::string, class element*> elements_map;
	std::map<std::string, class species*> species_map;
	std::map<std::string, class phase*> phases_map;
	std::map<std::string, class logk*> logk_map;
	std::map<std::string, class master_isotope*> master_isotope_map;

#if defined(PHREEQCI_GUI)
#include "../../phreeqci_gui.h"
#endif /* defined(PHREEQCI_GUI) */
	/* ----------------------------------------------------------------------
	*   ISOTOPES
	* ---------------------------------------------------------------------- */
	std::vector<class master_isotope*> master_isotope;
	int initial_solution_isotopes;
	std::vector<class calculate_value*> calculate_value;
	std::map<std::string, class calculate_value*> calculate_value_map;
public:
	std::map<std::string, class calculate_value*>&  GetCalculateValueMap() 
	{
		return this->calculate_value_map;
	}
protected:
	std::vector<class isotope_ratio*> isotope_ratio;
	std::map<std::string, class isotope_ratio*> isotope_ratio_map;
	std::vector<class isotope_alpha*> isotope_alpha;
	std::map<std::string, class isotope_alpha*> isotope_alpha_map;
	int phreeqc_mpi_myself;
	int first_read_input;
	std::string user_database;

	//int have_punch_name;
	/* VP: Density Start */
	int print_density;
	/* VP: Density End */

	int print_viscosity;
	LDBLE viscos, viscos_0, viscos_0_25; // viscosity of the solution, of pure water, of pure water at 25 C
	LDBLE density_x;
	LDBLE solution_volume_x;
	LDBLE solution_mass_x;
	LDBLE cell_pore_volume;
	LDBLE cell_porosity;
	LDBLE cell_volume;
	LDBLE cell_saturation;
	std::vector<class system_species> sys;
	LDBLE sys_tot;

	LDBLE V_solutes, rho_0, rho_0_sat, kappa_0, p_sat/*, ah2o_x0*/;
	LDBLE SC; // specific conductance mS/cm
	LDBLE eps_r; // relative dielectric permittivity
	LDBLE DH_A, DH_B, DH_Av; // Debye-Hueckel A, B and Av
	LDBLE QBrn; // Born function d(ln(eps_r))/dP / eps_r * 41.84004, for supcrt calc'n of molal volume
	LDBLE ZBrn; // Born function (-1/eps_r + 1) * 41.84004, for supcrt calc'n of molal volume
	LDBLE dgdP; // dg / dP, pressure derivative of g-function, for supcrt calc'n of molal volume

	int need_temp_msg;
	LDBLE solution_mass, solution_volume;

	/* phqalloc.cpp ------------------------------- */
	PHRQMemHeader* s_pTail;

	/* Basic */
	PBasic* basic_interpreter;

	double (*basic_callback_ptr) (double x1, double x2, const char* str, void* cookie);
	void* basic_callback_cookie;
#ifdef IPHREEQC_NO_FORTRAN_MODULE
	double (*basic_fortran_callback_ptr) (double* x1, double* x2, const char* str, size_t l);
#else
	double (*basic_fortran_callback_ptr) (double* x1, double* x2, const char* str, int l);
#endif
#if defined(SWIG) || defined(SWIG_IPHREEQC)
	class BasicCallback* basicCallback;
	void SetCallback(BasicCallback* cb) { basicCallback = cb; }
#endif

	/* cl1.cpp ------------------------------- */
	std::vector<double> x_arg, res_arg, scratch;
	/* gases.cpp ------------------------------- */
	LDBLE a_aa_sum, b2, b_sum, R_TK;
	std::map < std::pair<std::string, std::string>, double > gas_binary_parameters;

	/* input.cpp ------------------------------- */
	int check_line_return;
	int reading_db;

	/* integrate.cpp ------------------------------- */
	LDBLE midpoint_sv;
	LDBLE z_global, xd_global, alpha_global;

	/* inverse.cpp ------------------------------- */
	size_t max_row_count, max_column_count;
	int carbon;
	std::vector<const char*> col_name, row_name;
	size_t count_rows, count_optimize;
	size_t col_phases, col_redox, col_epsilon, col_ph, col_water,
		col_isotopes, col_phase_isotopes;
	size_t row_mb, row_fract, row_charge, row_carbon, row_isotopes,
		row_epsilon, row_isotope_epsilon, row_water;
	std::vector<double> inv_zero, array1, inv_res, inv_delta1, delta2, 
		delta3, inv_cu, delta_save;
	std::vector<double> min_delta, max_delta;
	std::vector<int> inv_iu, inv_is;
	size_t klmd, nklmd, n2d;
	int kode, iter;
	LDBLE toler, error, max_pct, scaled_error;
	class master* master_alk;
	std::vector<int> row_back, col_back;
	std::vector<unsigned long> good, bad, minimal;
	size_t max_good, max_bad, max_minimal;
	int count_good, count_bad, count_minimal, count_calls;
	unsigned long soln_bits, phase_bits, current_bits, temp_bits;
	FILE* netpath_file;
	int count_inverse_models, count_pat_solutions;
	int min_position[32], max_position[32], now[32];
	std::vector <std::string> inverse_heading_names;

	/* kinetics.cpp ------------------------------- */
public:
	int count_pp, count_pg, count_ss;
	void* cvode_kinetics_ptr;
	int cvode_test;
	int cvode_error;
	int cvode_n_user;
	int cvode_n_reactions;
	realtype cvode_step_fraction;
	realtype cvode_rate_sim_time;
	realtype cvode_rate_sim_time_start;
	realtype cvode_last_good_time;
	realtype cvode_prev_good_time;
	N_Vector cvode_last_good_y;
	N_Vector cvode_prev_good_y;
	M_Env kinetics_machEnv;
	N_Vector kinetics_y, kinetics_abstol;
	void* kinetics_cvode_mem;
	cxxSSassemblage* cvode_ss_assemblage_save;
	cxxPPassemblage* cvode_pp_assemblage_save;
protected:
	std::vector<double> m_temp, m_original, rk_moles, x0_moles;
	int set_and_run_attempt;

	/* model.cpp ------------------------------- */
	int gas_in;
	LDBLE min_value;
	std::vector<double> normal, ineq_array, res, cu, zero, delta1;
	std::vector<int> iu, is, back_eq;

	/* phrq_io_output.cpp ------------------------------- */
	int forward_output_to_log;

	/* phreeqc_files.cpp ------------------------------- */
	std::string default_data_base;
	/* Pitzer  */
	int pitzer_model, sit_model, pitzer_pe;
	int full_pitzer, always_full_pitzer, ICON, IC;
	LDBLE COSMOT;
	LDBLE AW;
	LDBLE VP, DW0;
	std::vector<class pitz_param*> pitz_params;
	std::map< std::string, size_t > pitz_param_map;
	std::vector<class theta_param*> theta_params;
	int use_etheta;
	LDBLE OTEMP, OPRESS;
	LDBLE A0;
	class pitz_param* aphi/* = NULL*/;
	std::vector<class species*> spec;
	class species** cations, ** anions, ** neutrals; // pointers to spec
	int count_cations, count_anions, count_neutrals;
	int MAXCATIONS, FIRSTANION, MAXNEUTRAL;
	class pitz_param* mcb0, * mcb1, * mcc0;
	std::vector<int> IPRSNT;
	std::vector<double> M, LGAMMA;
	LDBLE BK[23], DK[23];

	LDBLE dummy;

	/* print.cpp ------------------------------- */

	/* read.cpp */
	const char* prev_next_char;
#if defined PHREEQ98 
	int shifts_as_points;
#endif

	/* read_class.cxx */
	dumper dump_info;
	StorageBinList delete_info;
	runner run_info;
	char* sformatf_buffer;
	size_t sformatf_buffer_size;

	/* readtr.cpp */
	std::string dump_file_name_cpp;

	/* sit.cpp ------------------------------- */
	std::vector<class pitz_param*> sit_params;
	std::map< std::string, size_t > sit_param_map;
	LDBLE sit_A0;
	int sit_count_cations, sit_count_anions, sit_count_neutrals;
	int sit_MAXCATIONS, sit_FIRSTANION, sit_MAXNEUTRAL;
	std::vector<int> sit_IPRSNT;
	std::vector<double> sit_M, sit_LGAMMA;
	std::vector<int> s_list, cation_list, neutral_list, anion_list, ion_list, param_list;

	/* tidy.cpp ------------------------------- */
	LDBLE a0, a1, kc, kb;

	/* tally.cpp ------------------------------- */
	class tally_buffer* t_buffer;
	size_t tally_count_component;
	//class tally* tally_table;
	std::vector<class tally> tally_table;
	size_t count_tally_table_columns;
	size_t count_tally_table_rows;

	/* transport.cpp ------------------------------- */
	class sol_D* sol_D;
	class sol_D* sol_D_dbg;
	class J_ij* J_ij, * J_ij_il;
	int J_ij_count_spec;

	class M_S* m_s;
	int count_m_s;
	LDBLE tot1_h, tot1_o, tot2_h, tot2_o;
	LDBLE diffc_max, diffc_tr, J_ij_sum;
	int transp_surf;
	LDBLE* heat_mix_array;
	LDBLE* temp1, * temp2;
	int nmix, heat_nmix;
	LDBLE heat_mix_f_imm, heat_mix_f_m;
	int warn_MCD_X, warn_fixed_Surf;
	LDBLE current_x, current_A, fix_current; // current: coulomb / s, Ampere, fixed current (Ampere)

	/* utilities.cpp ------------------------------- */
	int spinner;
	std::map<std::string, double> gfw_map;
	std::map<const char*, int> rates_map;

	/* new after release of Version 3 */
	std::map<std::string, std::vector < std::string> > sum_species_map;
	std::map<std::string, std::vector < std::string> > sum_species_map_db;

	friend class PBasic;
	friend class ChartObject;
	friend class IPhreeqc;
	friend class TestIPhreeqc;
	friend class TestSelectedOutput;
	friend class IPhreeqcMMS;
	friend class IPhreeqcPhast;
	friend class PhreeqcRM;

	std::vector<int> keycount;  // used to mark keywords that have been read 

public:
	static const class const_iso iso_defaults[];
	static const int count_iso_defaults;
};
#endif /* _INC_PHREEQC_H */

#ifndef _INC_ISFINITE_H
#define _INC_ISFINITE_H
/*********************************
isfinite handling
(Note: Should NOT be guarded)
**********************************/

#if defined (PHREEQ98) || defined (_MSC_VER)
#  define HAVE_FINITE
#  define finite _finite
#else  /*defined (PHREEQ98) || defined (_MSC_VER)*/
#  if defined(DJGPP)
#    define HAVE_FINITE
#  endif
#endif /*defined (PHREEQ98) || defined (_MSC_VER)*/

#if defined(HAVE_ISFINITE)
#  if __GNUC__ && (__cplusplus >= 201103L)
#    define PHR_ISFINITE(x) std::isfinite(x)
#  else
#  define PHR_ISFINITE(x) std::isfinite(x) /* changed when <math.h> was changed to <cmath> */
#  endif
#elif defined(HAVE_FINITE)
#  define PHR_ISFINITE(x) finite(x)
#elif defined(HAVE_ISNAN)
#  define PHR_ISFINITE(x) ( ((x) == 0.0) || ((!std::isnan(x)) && ((x) != (2.0 * (x)))) )
#else
#  define PHR_ISFINITE(x) ( ((x) == 0.0) || (((x) == (x)) && ((x) != (2.0 * (x)))) )
#endif
#endif // _INC_ISFINITE_H

#ifndef _INC_UTILITIES_NAMESPACE_H
#define _INC_UTILITIES_NAMESPACE_H
namespace Utilities
{
	LDBLE get_nan(void);

	// operations on maps of entities (Solution, Exchange, ...)
	template < typename T >
	void Rxn_dump_raw(const T& b, std::ostream& s_oss, unsigned int indent)
	{
		typename T::const_iterator it;
		for (it = b.begin(); it != b.end(); ++it)
		{
			// Adding logic to dump only non-negative entities
			//if (it->second.Get_n_user() >= 0)
			if (it->first >= 0 && it->second.Get_n_user() >= 0)
			{
				it->second.dump_raw(s_oss, indent);
			}
		}
		return;
	}

	template < typename T >
	void Rxn_dump_raw_range(const T& b, std::ostream& s_oss, int start, int end, unsigned int indent)
	{
		typename T::const_iterator it;
		for (int i = start; i <= end; i++)
		{
			if (i < 0) continue;
			it = b.find(i);
			if (it != b.end())
			{
				it->second.dump_raw(s_oss, indent);
			}
		}
		return;
	}

	template < typename T >
	T* Rxn_find(std::map < int, T >& b, int i)
	{
		if (b.find(i) != b.end())
		{
			return (&(b.find(i)->second));
		}
		else
		{
			return (NULL);
		}
	}

	template < typename T >
	int Rxn_next_user_number(std::map < int, T >& b)
	{
		int ret = 0;
		if (b.size() != 0)
		{
			ret = b.rbegin()->first + 1;
		}
		return ret;
	}

	template < typename T >
	T* Rxn_copy(std::map < int, T >& b, int i, int j)
	{
		typename std::map < int, T >::iterator it;
		it = b.find(i);
		if (it != b.end())
		{
			b[j] = it->second;
			it = b.find(j);
			it->second.Set_n_user(j);
			it->second.Set_n_user_end(j);
			return &(it->second);
		}
		else
		{
			return (NULL);
		}
	}

	template < typename T >
	void Rxn_copies(std::map < int, T >& b, int n_user, int n_user_end)
	{
		if (n_user_end <= n_user) return;
		typename std::map < int, T >::iterator it;
		it = b.find(n_user);
		if (it != b.end())
		{
			for (int j = n_user + 1; j <= n_user_end; j++)
			{
				b[j] = it->second;
				it = b.find(j);
				it->second.Set_n_user(j);
				it->second.Set_n_user_end(j);
			}
		}
	}
	template < typename T >
	int Rxn_read_raw(std::map < int, T >& m, std::set < int >& s, Phreeqc* phreeqc_cookie)
	{
		typename std::map < int, T >::iterator it;
		assert(!phreeqc_cookie->reading_database());

		T entity(phreeqc_cookie->Get_phrq_io());

		CParser parser(phreeqc_cookie->Get_phrq_io());
		entity.read_raw(parser);

		// Store
		if (entity.Get_base_error_count() == 0)
		{
			m[entity.Get_n_user()] = entity;
		}

		// Make copies if necessary
		Utilities::Rxn_copies(m, entity.Get_n_user(), entity.Get_n_user_end());
		for (int i = entity.Get_n_user(); i <= entity.Get_n_user_end(); i++)
		{
			s.insert(i);
		}
		return phreeqc_cookie->cleanup_after_parser(parser);
	}

	template < typename T >
	int Rxn_read_modify(std::map < int, T >& m, std::set < int >& s, Phreeqc* phreeqc_cookie)
	{
		typename std::map < int, T >::iterator it;

		CParser parser(phreeqc_cookie->Get_phrq_io());

		std::string key_name;
		std::string::iterator b = parser.line().begin();
		std::string::iterator e = parser.line().end();
		CParser::copy_token(key_name, b, e);

		cxxNumKeyword nk;
		nk.read_number_description(parser);
		T* entity_ptr = Utilities::Rxn_find(m, nk.Get_n_user());
		if (!entity_ptr)
		{
			std::ostringstream errstr;
			errstr << "Could not find " << key_name << " " << nk.Get_n_user() << ", ignoring modify data.\n";
			phreeqc_cookie->warning_msg(errstr.str().c_str());
			//phreeqc_cookie->error_msg(errstr.str().c_str(), PHRQ_io::OT_STOP);

			// Don't throw, read data into dummy entity, then ignore
			T entity;
			entity_ptr = &entity;
			entity_ptr->read_raw(parser, false);
			return phreeqc_cookie->cleanup_after_parser(parser);
		}

		entity_ptr->read_raw(parser, false);
		entity_ptr->Set_n_user(nk.Get_n_user());
		entity_ptr->Set_n_user_end(nk.Get_n_user_end());
		entity_ptr->Set_description(nk.Get_description());
		s.insert(entity_ptr->Get_n_user());

		return phreeqc_cookie->cleanup_after_parser(parser);
	}

	template < typename T >
	int SB_read_modify(std::map < int, T >& m, CParser& parser)
	{
		typename std::map < int, T >::iterator it;

		std::string key_name;
		std::string::iterator b = parser.line().begin();
		std::string::iterator e = parser.line().end();
		CParser::copy_token(key_name, b, e);

		cxxNumKeyword nk;
		nk.read_number_description(parser);
		T* entity_ptr = Utilities::Rxn_find(m, nk.Get_n_user());
		if (!entity_ptr)
		{
			std::ostringstream errstr;
			errstr << "Could not find " << key_name << " " << nk.Get_n_user() << ", ignoring modify data.\n";
			//io->warning_msg(errstr.str().c_str());

			// Don't throw, read data into dummy entity, then ignore
			T entity;
			entity_ptr = &entity;
			entity_ptr->read_raw(parser, false);
			return FALSE;
		}

		entity_ptr->read_raw(parser, false);
		entity_ptr->Set_n_user(nk.Get_n_user());
		entity_ptr->Set_n_user_end(nk.Get_n_user_end());
		entity_ptr->Set_description(nk.Get_description());

		return TRUE;
	}

	template < typename T >
	void Rxn_mix(std::map <int, cxxMix>& mix_map, std::map < int, T >& entity_map, Phreeqc* phreeqc_cookie)
	{
		std::map<int, cxxMix>::iterator mix_it;
		for (mix_it = mix_map.begin(); mix_it != mix_map.end(); mix_it++)
		{
			T entity(entity_map, mix_it->second, mix_it->second.Get_n_user(), phreeqc_cookie->Get_phrq_io());
			entity_map[mix_it->second.Get_n_user()] = entity;
			Utilities::Rxn_copies(entity_map, mix_it->second.Get_n_user(), mix_it->second.Get_n_user_end());
		}
		mix_map.clear();
	}

} // namespace Utilities


#if defined(PHREEQCI_GUI)
void PhreeqcIWait(Phreeqc* phreeqc);
#endif

#if !defined(NDEBUG) && defined(WIN32_MEMORY_DEBUG)
#define   string_duplicate(s)             _string_duplicate(s, __FILE__, __LINE__)
#endif
#if defined(_DEBUG)
char* _string_duplicate(const char* token, const char* szFileName, int nLine);
#endif

#endif //_INC_UTILITIES_NAMESPACE_H

#ifndef _INC_MISSING_SNPRINTF_H
#define _INC_MISSING_SNPRINTF_H

// Section _INC_MISSING_SNPRINTF_H is based on
// https://stackoverflow.com/questions/2915672/snprintf-and-visual-studio-2010

#if defined(_MSC_VER) && (_MSC_VER < 1900)
#if (_MSC_VER <= 1700) // VS2012
namespace std {
	__inline bool isnan(double num) {
		return _isnan(num) != 0;
	}
}
#endif
#include <stdarg.h>

#define snprintf c99_snprintf
#define vsnprintf c99_vsnprintf

#pragma warning( push )
// warning C4793: 'vararg' : causes native code generation
#pragma warning( disable : 4793 )

__inline int c99_vsnprintf(char *outBuf, size_t size, const char *format, va_list ap)
{
    int count = -1;

    if (size != 0)
        count = _vsnprintf_s(outBuf, size, _TRUNCATE, format, ap);
    if (count == -1)
        count = _vscprintf(format, ap);

    return count;
}

__inline int c99_snprintf(char *outBuf, size_t size, const char *format, ...)
{
    int count;
    va_list ap;

    va_start(ap, format);
    count = c99_vsnprintf(outBuf, size, format, ap);
    va_end(ap);

    return count;
}
#pragma warning( pop )
#endif // defined(_MSC_VER) && (_MSC_VER < 1900)

#endif //_INC_MISSING_SNPRINTF_H
