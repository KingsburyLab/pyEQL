#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "IPhreeqc.h"
#include "Phreeqc.h" /* snprintf */
#include "CVar.hxx"

using ::testing::HasSubstr;

#if defined(_WIN32)
#include <windows.h>
#endif

#include <fstream>
#include <string>
#include <string.h> // strstr
#include <cmath>
#include <cfloat>
#include <stdlib.h>

#include "FileTest.h"


IPQ_RESULT SOLUTION(int n, double C, double Ca, double Na);
IPQ_RESULT EQUILIBRIUM_PHASES(int n, const char* phase, double si, double amount);
IPQ_RESULT USER_PUNCH(int n, const char* element, int max);
IPQ_RESULT SELECTED_OUTPUT(int n);
IPQ_RESULT DUMP(int n);

void TestFileOnOff(const char* FILENAME_FORMAT, int output_file_on, int error_file_on, int log_file_on, int selected_output_file_on, int dump_file_on);


const char ex15_dat[] =
"SOLUTION_MASTER_SPECIES\n"
"C        CO2            2.0     61.0173         12.0111\n"
"Cl       Cl-            0.0     Cl              35.453\n"
"Co       Co+2           0.0     58.93           58.93   \n"
"E        e-             0.0     0.0             0.0\n"
"H        H+             -1.     1.008           1.008\n"
"H(0)     H2             0.0     1.008\n"
"H(1)     H+             -1.     1.008\n"
"N        NH4+           0.0     14.0067         14.0067\n"
"Na       Na+            0.0     Na              22.9898\n"
"Nta      Nta-3          3.0     1.              1.\n"
"O        H2O            0.0     16.00           16.00\n"
"O(-2)    H2O            0.0     18.016\n"
"O(0)     O2             0.0     16.00\n"
"SOLUTION_SPECIES\n"
"2H2O = O2 + 4H+ + 4e- \n"
"        log_k   -86.08; -gamma  1e7   0.0\n"
"2 H+ + 2 e- = H2\n"
"        log_k   -3.15;  -gamma  1e7   0.0\n"
"H+ = H+\n"
"        log_k   0.0;    -gamma  1e7   0.0\n"
"e- = e-\n"
"        log_k   0.0;    -gamma  1e7   0.0\n"
"H2O = H2O\n"
"        log_k   0.0;    -gamma  1e7   0.0\n"
"CO2 = CO2\n"
"        log_k   0.0;    -gamma  1e7   0.0\n"
"Na+ = Na+\n"
"        log_k   0.0;    -gamma  1e7   0.0\n"
"Cl- = Cl-\n"
"        log_k   0.0;    -gamma  1e7   0.0\n"
"Co+2 = Co+2\n"
"        log_k   0.0;    -gamma  1e7   0.0\n"
"NH4+ = NH4+\n"
"        log_k   0.0;    -gamma  1e7   0.0\n"
"Nta-3 = Nta-3\n"
"        log_k   0.0;    -gamma  1e7   0.0\n"
"Nta-3 + 3H+ = H3Nta\n"
"        log_k   14.9;   -gamma  1e7   0.0\n"
"Nta-3 + 2H+ = H2Nta-\n"
"        log_k   13.3;   -gamma  1e7   0.0\n"
"Nta-3 + H+ = HNta-2\n"
"        log_k   10.3;   -gamma  1e7   0.0\n"
"Nta-3 + Co+2 = CoNta-\n"
"        log_k   11.7;   -gamma  1e7   0.0\n"
"2 Nta-3 + Co+2 = CoNta2-4\n"
"        log_k   14.5;   -gamma  1e7   0.0\n"
"Nta-3 + Co+2 + H2O = CoOHNta-2 + H+\n"
"        log_k   0.5;    -gamma  1e7   0.0\n"
"Co+2 + H2O = CoOH+ + H+\n"
"        log_k   -9.7;   -gamma  1e7   0.0\n"
"Co+2 + 2H2O = Co(OH)2 + 2H+\n"
"        log_k   -22.9;  -gamma  1e7   0.0\n"
"Co+2 + 3H2O = Co(OH)3- + 3H+\n"
"        log_k   -31.5;  -gamma  1e7   0.0\n"
"CO2 + H2O = HCO3- + H+\n"
"        log_k   -6.35;  -gamma  1e7   0.0\n"
"CO2 + H2O = CO3-2 + 2H+\n"
"        log_k   -16.68; -gamma  1e7   0.0\n"
"NH4+ = NH3 + H+\n"
"        log_k   -9.3;   -gamma  1e7   0.0\n"
"H2O = OH- +  H+\n"
"        log_k   -14.0;  -gamma  1e7   0.0\n"
"END\n";


TEST(TestIPhreeqcLib, TestCreateIPhreeqc)
{
	int n = ::CreateIPhreeqc();
	ASSERT_GE(n, 0);

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestDestroyIPhreeqc)
{
	for (int i = 0; i < 10; ++i)
	{
		ASSERT_EQ(IPQ_BADINSTANCE, ::DestroyIPhreeqc(i));
	}

	int n = ::CreateIPhreeqc();
	ASSERT_GE(n, 0);
	ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
}

TEST(TestIPhreeqcLib, TestLoadDatabase)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	std::string FILES[] = { "phreeqc.dat", "pitzer.dat", "wateq4f.dat",
		"Amm.dat", "frezchem.dat", "iso.dat",
		"llnl.dat", "minteq.dat", "minteq.v4.dat",
		"sit.dat","ColdChem.dat","core10.dat",
		"Tipping_Hurley.dat"
	};

	for (int j = 0; j < sizeof(FILES) / sizeof(std::string); ++j)
	{
		for (int i = 0; i < 10; ++i)
		{
			ASSERT_EQ(true, ::FileExists(FILES[j].c_str()));
			ASSERT_TRUE(::FileSize(FILES[j].c_str()) > 0);
			ASSERT_EQ(0, ::LoadDatabase(n, FILES[j].c_str()));
		}
	}
	ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
}

TEST(TestIPhreeqcLib, TestLoadDatabaseString)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	for (int i = 0; i < 10; ++i)
	{
		ASSERT_EQ(0, ::LoadDatabaseString(n, ex15_dat));
	}

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestLoadDatabaseMissingFile)
{
	ASSERT_EQ(false, ::FileExists("missing.file"));

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	for (int i = 0; i < 10; ++i)
	{
		ASSERT_EQ(1, ::LoadDatabase(n, "missing.file"));
	}

	const char expected[] =
		"ERROR: LoadDatabase: Unable to open:\"missing.file\".\n";

	const char* err = ::GetErrorString(n);

	ASSERT_EQ(std::string(expected), std::string(err));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestSetErrorOn)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(1, ::GetErrorOn(n));			// initial setting is true
	ASSERT_EQ(IPQ_OK, ::SetErrorOn(n, 0));
	ASSERT_EQ(0, ::GetErrorOn(n));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestLoadDatabaseWithErrors)
{
	ASSERT_EQ(true, ::FileExists("missing_e.dat"));

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	for (int i = 0; i < 10; ++i)
	{
		ASSERT_EQ(true, ::FileExists("missing_e.dat"));
		ASSERT_TRUE(::FileSize("missing_e.dat") > 0);
		ASSERT_EQ(6, ::LoadDatabase(n, "missing_e.dat"));

		static const char* expected =
			"ERROR: Could not reduce equation to primary master species, CH4.\n"
			"ERROR: Could not reduce equation to primary master species, Cu+.\n"
			"ERROR: Could not reduce equation to primary master species, Fe+3.\n"
			"ERROR: Could not reduce equation to primary master species, H2.\n"
			"ERROR: Could not reduce equation to primary master species, Mn+3.\n"
			"ERROR: Could not reduce equation to primary master species, NH4+.\n"
			"ERROR: Could not reduce equation to primary master species, N2.\n"
			"ERROR: Could not reduce equation to primary master species, NO2-.\n"
			"ERROR: Could not reduce equation to primary master species, O2.\n"
			"ERROR: Could not reduce equation to primary master species, HS-.\n"
			"ERROR: Could not reduce equation to secondary master species, e-.\n"
			"ERROR: Non-master species in secondary reaction, e-.\n"
			"ERROR: No master species for element e.\n"
			"ERROR: Could not find primary master species for e.\n"
			"ERROR: No master species for element e.\n"
			"ERROR: Could not reduce equation to secondary master species, Hausmannite.\n"
			"ERROR: Could not reduce equation to secondary master species, Manganite.\n"
			"ERROR: Could not reduce equation to secondary master species, Pyrite.\n"
			"ERROR: Could not reduce equation to secondary master species, Pyrolusite.\n"
			"ERROR: Could not reduce equation to secondary master species, Sulfur.\n"
			"ERROR: e-, primary master species for E-, not defined.\n"
			"ERROR: Calculations terminating due to input errors.\n";

		const char* err = ::GetErrorString(n);

		ASSERT_EQ(std::string(expected), std::string(err));
	}

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestRunAccumulated)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0,      ::LoadDatabase(n, "phreeqc.dat"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "solution 12"));
	ASSERT_EQ(0,      ::GetOutputFileOn(n));
	ASSERT_EQ(0,      ::GetErrorFileOn(n));
	ASSERT_EQ(0,      ::GetLogFileOn(n));
	ASSERT_EQ(0,      ::GetSelectedOutputFileOn(n));
	ASSERT_EQ(0,      ::GetDumpFileOn(n));
	ASSERT_EQ(0,      ::GetDumpStringOn(n));
	ASSERT_EQ(0,      ::RunAccumulated(n));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestRunWithErrors)
{
	const char dump_file[] = "error.inp";

	// remove dump file if it exists
	//
	FileTest dump(dump_file);
	ASSERT_TRUE(dump.RemoveExisting());

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0,      ::LoadDatabase(n, "phreeqc.dat"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "SOLUTION 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	pH	7"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Na	1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "PHASES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Fix_H+"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	H+ = H+"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	log_k	0"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "EQUILIBRIUM_PHASES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Fix_H+ -10 HCl	10"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "END"));
	ASSERT_EQ(0,      ::GetOutputFileOn(n));
	ASSERT_EQ(0,      ::GetErrorFileOn(n));
	ASSERT_EQ(0,      ::GetLogFileOn(n));
	ASSERT_EQ(0,      ::GetSelectedOutputFileOn(n));
	ASSERT_EQ(0,      ::GetDumpFileOn(n));
	ASSERT_EQ(0,      ::GetDumpStringOn(n));
	ASSERT_EQ(1,      ::RunAccumulated(n));
				      
	const char expected[] =
		"ERROR: Numerical method failed on all combinations of convergence parameters, cell/soln/mix 1\n";
	const char* err = ::GetErrorString(n);

	ASSERT_EQ(std::string(expected), std::string(err));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
	ASSERT_TRUE(dump.VerifyExists());
}

TEST(TestIPhreeqcLib, TestRunFile)
{
	static const char dump_file[] = "error.inp";

	// remove dump file if it exists
	//
	if (::FileExists(dump_file))
	{
		ASSERT_TRUE(::DeleteFile(dump_file));
	}
	ASSERT_EQ(false, ::FileExists(dump_file));

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(true, ::FileExists("phreeqc.dat"));
	ASSERT_TRUE(::FileSize("phreeqc.dat") > 0);
	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));
	ASSERT_EQ(0, ::GetOutputFileOn(n));
	ASSERT_EQ(0, ::GetErrorFileOn(n));
	ASSERT_EQ(0, ::GetLogFileOn(n));
	ASSERT_EQ(0, ::GetSelectedOutputFileOn(n));
	ASSERT_EQ(0, ::GetDumpFileOn(n));
	ASSERT_EQ(0, ::GetDumpStringOn(n));
	ASSERT_EQ(true, ::FileExists("conv_fail.in"));
	ASSERT_TRUE(::FileSize("conv_fail.in") > 0);
	ASSERT_EQ(1, ::RunFile(n, "conv_fail.in"));

	static const char expected[] =
		"ERROR: Numerical method failed on all combinations of convergence parameters, cell/soln/mix 1\n";
	const char* err = ::GetErrorString(n);

	ASSERT_EQ(std::string(expected), std::string(err));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}

	// Note: should this file exist since GetDumpFileOn is false?
	ASSERT_EQ(true, ::FileExists(dump_file));
	ASSERT_TRUE(::FileSize(dump_file) > 0);
	ASSERT_TRUE(::DeleteFile(dump_file));
}

TEST(TestIPhreeqcLib, TestRunString)
{
	const char input[] =
		"TITLE Example 1.--Add uranium and speciate seawater.\n"
		"SOLUTION 1  SEAWATER FROM NORDSTROM ET AL. (1979)\n"
		"        units   ppm\n"
		"        pH      8.22\n"
		"        pe      8.451\n"
		"        density 1.023\n"
		"        temp    25.0\n"
		"        redox   O(0)/O(-2)\n"
		"        Ca              412.3\n"
		"        Mg              1291.8\n"
		"        Na              10768.0\n"
		"        K               399.1\n"
		"        Fe              0.002\n"
		"        Mn              0.0002  pe\n"
		"        Si              4.28\n"
		"        Cl              19353.0\n"
		"        Alkalinity      141.682 as HCO3\n"
		"        S(6)            2712.0\n"
		"        N(5)            0.29    gfw   62.0\n"
		"        N(-3)           0.03    as    NH4\n"
		"        U               3.3     ppb   N(5)/N(-3)\n"
		"        O(0)            1.0     O2(g) -0.7\n"
		"SOLUTION_MASTER_SPECIES\n"
		"        U       U+4     0.0     238.0290     238.0290\n"
		"        U(4)    U+4     0.0     238.0290\n"
		"        U(5)    UO2+    0.0     238.0290\n"
		"        U(6)    UO2+2   0.0     238.0290\n"
		"SOLUTION_SPECIES\n"
		"        #primary master species for U\n"
		"        #is also secondary master species for U(4)\n"
		"        U+4 = U+4\n"
		"                log_k          0.0\n"
		"        U+4 + 4 H2O = U(OH)4 + 4 H+\n"
		"                log_k          -8.538\n"
		"                delta_h        24.760 kcal\n"
		"        U+4 + 5 H2O = U(OH)5- + 5 H+\n"
		"                log_k          -13.147\n"
		"                delta_h        27.580 kcal\n"
		"        #secondary master species for U(5)\n"
		"        U+4 + 2 H2O = UO2+ + 4 H+ + e-\n"
		"                log_k          -6.432\n"
		"                delta_h        31.130 kcal\n"
		"        #secondary master species for U(6)\n"
		"        U+4 + 2 H2O = UO2+2 + 4 H+ + 2 e-\n"
		"                log_k          -9.217\n"
		"                delta_h        34.430 kcal\n"
		"        UO2+2 + H2O = UO2OH+ + H+\n"
		"                log_k          -5.782\n"
		"                delta_h        11.015 kcal\n"
		"        2UO2+2 + 2H2O = (UO2)2(OH)2+2 + 2H+\n"
		"                log_k          -5.626\n"
		"                delta_h        -36.04 kcal\n"
		"        3UO2+2 + 5H2O = (UO2)3(OH)5+ + 5H+\n"
		"                log_k          -15.641\n"
		"                delta_h        -44.27 kcal\n"
		"        UO2+2 + CO3-2 = UO2CO3\n"
		"                log_k          10.064\n"
		"                delta_h        0.84 kcal\n"
		"        UO2+2 + 2CO3-2 = UO2(CO3)2-2\n"
		"                log_k          16.977\n"
		"                delta_h        3.48 kcal\n"
		"        UO2+2 + 3CO3-2 = UO2(CO3)3-4\n"
		"                log_k          21.397\n"
		"                delta_h        -8.78 kcal\n"
		"PHASES\n"
		"        Uraninite\n"
		"        UO2 + 4 H+ = U+4 + 2 H2O\n"
		"        log_k          -3.490\n"
		"        delta_h        -18.630 kcal\n"
		"END\n"
		"\n";

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	char OUTPUT_FILE[80];
	snprintf(OUTPUT_FILE, sizeof(OUTPUT_FILE), "phreeqc.%d.out", n);

	if (::FileExists(OUTPUT_FILE))
	{
		ASSERT_TRUE(::DeleteFile(OUTPUT_FILE));
	}
	ASSERT_EQ(false, ::FileExists(OUTPUT_FILE));

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));
	::SetOutputFileOn(n, 1);
	::SetErrorFileOn(n, 0);
	::SetLogFileOn(n, 0);
	::SetSelectedOutputFileOn(n, 0);
	::SetDumpFileOn(n, 0);
	ASSERT_EQ(false, ::FileExists(OUTPUT_FILE));
	ASSERT_EQ(0, ::RunString(n, input));
	ASSERT_EQ(true, ::FileExists(OUTPUT_FILE));
	ASSERT_TRUE(::FileSize(OUTPUT_FILE) > 0);
	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
	if (::FileExists(OUTPUT_FILE))
	{
		ASSERT_TRUE(::DeleteFile(OUTPUT_FILE));
	}
}

TEST(TestIPhreeqcLib, TestGetSelectedOutputRowCount)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(true, ::FileExists("llnl.dat"));
	ASSERT_TRUE(::FileSize("llnl.dat") > 0);
	ASSERT_EQ(0, ::LoadDatabase(n, "llnl.dat"));

	int max = 6;

	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));
	ASSERT_EQ(IPQ_OK, ::EQUILIBRIUM_PHASES(n, "calcite", 0.0, 0.010));
	ASSERT_EQ(IPQ_OK, ::USER_PUNCH(n, "Ca", max));
	ASSERT_EQ(0, ::GetOutputFileOn(n));
	ASSERT_EQ(0, ::GetErrorFileOn(n));
	ASSERT_EQ(0, ::GetLogFileOn(n));
	ASSERT_EQ(0, ::GetSelectedOutputFileOn(n));
	ASSERT_EQ(0, ::GetDumpFileOn(n));
	ASSERT_EQ(0, ::GetDumpStringOn(n));
	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(3, ::GetSelectedOutputRowCount(n)); // rows + header

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}
#define USE_VAR
TEST(TestIPhreeqcLib, TestGetSelectedOutputValue)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "llnl.dat"));

	int max = 6;

	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));
	ASSERT_EQ(IPQ_OK, ::EQUILIBRIUM_PHASES(n, "calcite", 0.0, 0.010));
	ASSERT_EQ(IPQ_OK, ::USER_PUNCH(n, "Ca", max));
	ASSERT_EQ(0, ::GetOutputFileOn(n));
	ASSERT_EQ(0, ::GetErrorFileOn(n));
	ASSERT_EQ(0, ::GetLogFileOn(n));
	ASSERT_EQ(0, ::GetSelectedOutputFileOn(n));
	ASSERT_EQ(0, ::GetDumpFileOn(n));
	ASSERT_EQ(0, ::GetDumpStringOn(n));
	ASSERT_EQ(0, ::RunAccumulated(n));

	/*
	EXPECTED selected.out:
			 sim	       state	        soln	      dist_x	        time	        step	          pH	          pe	           C	          Ca	          Na	     m_CO3-2	     m_CaOH+	    m_NaCO3-	    la_CO3-2	    la_CaOH+	   la_NaCO3-	     Calcite	   d_Calcite	   si_CO2(g)	 si_Siderite	    pressure	   total mol	      volume	    g_CO2(g)	     g_N2(g)	    k_Albite	   dk_Albite	    k_Pyrite	   dk_Pyrite	     s_CaSO4	     s_SrSO4	      1.name	      1.type	     1.moles	      2.name	      2.type	     2.moles	      3.name	      3.type	     3.moles	      4.name	      4.type	     4.moles	      5.name	      5.type	     5.moles	      6.name	      6.type	     6.moles
			   1	      i_soln	           1	         -99	         -99	         -99	           7	           4	 1.0000e-003	 1.0000e-003	 1.0000e-003	 4.2975e-007	 1.1819e-009	 1.1881e-009	-6.4686e+000	-8.9530e+000	-8.9507e+000	 0.0000e+000	 0.0000e+000	     -2.2870	   -999.9990	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	        Ca+2	          aq	 9.9178e-004	     CaHCO3+	          aq	 7.5980e-006	       CaCO3	          aq	 6.2155e-007	       CaOH+	          aq	 1.1819e-009
			   1	       react	           1	         -99	           0	           1	     7.86135	       10.18	 1.1556e-003	 1.1556e-003	 1.0000e-003	 4.2718e-006	 9.7385e-009	 1.1620e-008	-5.4781e+000	-8.0388e+000	-7.9621e+000	 9.8444e-003	-1.5555e-004	     -3.0192	   -999.9990	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	     calcite	        equi	 9.8444e-003	        Ca+2	          aq	 1.1371e-003	     CaHCO3+	          aq	 1.1598e-005	       CaCO3	          aq	 6.8668e-006	       CaOH+	          aq	 9.7385e-009
	*/

#if defined(USE_VAR)
	VAR v;
	::VarInit(&v);
#else
	CVar v;
#endif

	ASSERT_EQ(IPQ_INVALIDROW, ::GetSelectedOutputValue(n, -1, 0, &v));
	ASSERT_EQ(TT_ERROR, v.type);
	ASSERT_EQ(VR_INVALIDROW, v.vresult);

	ASSERT_EQ(IPQ_INVALIDROW, ::GetSelectedOutputValue(n, ::GetSelectedOutputRowCount(n), 0, &v));
	ASSERT_EQ(TT_ERROR, v.type);
	ASSERT_EQ(VR_INVALIDROW, v.vresult);

	ASSERT_EQ(IPQ_INVALIDCOL, ::GetSelectedOutputValue(n, 0, -1, &v));
	ASSERT_EQ(TT_ERROR, v.type);
	ASSERT_EQ(VR_INVALIDCOL, v.vresult);

	ASSERT_EQ(IPQ_INVALIDCOL, ::GetSelectedOutputValue(n, 0, ::GetSelectedOutputColumnCount(n), &v));
	ASSERT_EQ(TT_ERROR, v.type);
	ASSERT_EQ(VR_INVALIDCOL, v.vresult);


	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, 0, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("sim"), std::string(v.sVal));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, 1, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("state"), std::string(v.sVal));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, 2, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("soln"), std::string(v.sVal));

	//{{{{{{
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, 3, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("dist_x"), std::string(v.sVal));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, 4, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("time"), std::string(v.sVal));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, 5, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("step"), std::string(v.sVal));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, 6, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("pH"), std::string(v.sVal));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, 7, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("pe"), std::string(v.sVal));

	int col = 7;

	// -totals C Ca Na
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("C(mol/kgw)"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("Ca(mol/kgw)"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("Na(mol/kgw)"), std::string(v.sVal));

	// -molalities CO3-2  CaOH+  NaCO3-
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("m_CO3-2(mol/kgw)"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("m_CaOH+(mol/kgw)"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("m_NaCO3-(mol/kgw)"), std::string(v.sVal));

	// -activities CO3-2  CaOH+  NaCO3-
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("la_CO3-2"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("la_CaOH+"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("la_NaCO3-"), std::string(v.sVal));

	// -equilibrium_phases Calcite
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("Calcite"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("d_Calcite"), std::string(v.sVal));


	// -saturation_indices CO2(g) Siderite
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("si_CO2(g)"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("si_Siderite"), std::string(v.sVal));

	// -gases CO2(g) N2(g)
	//                      pressure "total mol" volume g_CO2(g) g_N2(g)
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("pressure"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("total mol"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("volume"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("g_CO2(g)"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("g_N2(g)"), std::string(v.sVal));

	// -kinetic_reactants Albite Pyrite
	//                               k_Albite dk_Albite k_Pyrite dk_Pyrite
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("k_Albite"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("dk_Albite"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("k_Pyrite"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("dk_Pyrite"), std::string(v.sVal));

	// -solid_solutions CaSO4 SrSO4
	//                              s_CaSO4 s_SrSO4
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("s_CaSO4"), std::string(v.sVal));
	++col;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("s_SrSO4"), std::string(v.sVal));

	for (int i = 0; i < max; ++i)
	{
		std::ostringstream oss1, oss2, oss3;

		// 1.name
		//
		ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col + 1 + (i * 3), &v));
		ASSERT_EQ(TT_STRING, v.type);
		oss1 << i + 1 << ".name";
		ASSERT_EQ(oss1.str(), std::string(v.sVal));

		// 1.type
		//
		ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col + 2 + (i * 3), &v));
		ASSERT_EQ(TT_STRING, v.type);
		oss2 << i + 1 << ".type";
		ASSERT_EQ(oss2.str(), std::string(v.sVal));

		// 1.moles
		//
		ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, col + 3 + (i * 3), &v));
		ASSERT_EQ(TT_STRING, v.type);
		oss3 << i + 1 << ".moles";
		ASSERT_EQ(oss3.str(), std::string(v.sVal));
	}

	// sim
	//
	col = 0;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_LONG, v.type);
	ASSERT_EQ(1L, v.lVal);
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_LONG, v.type);
	ASSERT_EQ(1L, v.lVal);

	// state
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("i_soln"), std::string(v.sVal));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("react"), std::string(v.sVal));

	// soln
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_LONG, v.type);
	ASSERT_EQ(1L, v.lVal);
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_LONG, v.type);
	ASSERT_EQ(1L, v.lVal);

	// dist_x -- sometimes as double sometimes as long (depends on state)
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_EQ(-99., v.dVal);
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_EQ(-99., v.dVal);


	// time -- sometimes as double sometimes as long (depends on state)
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_EQ(-99., v.dVal);
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	//CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v.dVal, ::pow(10., -DBL_DIG));
	ASSERT_NEAR(0.0, v.dVal, ::pow(10., -DBL_DIG));
	ASSERT_DOUBLE_EQ(0.0, v.dVal);

	// step
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_LONG, v.type);
	ASSERT_EQ(-99L, v.lVal);
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_LONG, v.type);
	ASSERT_EQ(1L, v.lVal);


	// pH
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(7.0, v.dVal, ::pow(10., -DBL_DIG));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(7.861354, v.dVal, ::pow(10., -6));

	// pe
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(4.0, v.dVal, ::pow(10., -DBL_DIG));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	// COMMENT: {8/8/2013 3:38:12 PM}	ASSERT_NEAR( 9.90852, v.dVal, ::pow(10., -1) );

		//
		// -totals C Ca Na
		//

		// C(mol/kgw)
		//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(1.0000e-003, v.dVal, ::pow(10., -DBL_DIG));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(1.1556e-003, v.dVal, ::pow(10., -7));


	// Ca(mol/kgw)
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(1.0000e-003, v.dVal, ::pow(10., -DBL_DIG));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(1.1556e-003, v.dVal, ::pow(10., -7));


	// Na(mol/kgw)
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(1.0000e-003, v.dVal, ::pow(10., -DBL_DIG));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(1.0000e-003, v.dVal, ::pow(10., -7));

	// -molalities CO3-2  CaOH+  NaCO3-
	col += 3;

	// -activities CO3-2  CaOH+  NaCO3-
	col += 3;

	// -equilibrium_phases Calcite
	col += 2;

	// -saturation_indices CO2(g) Siderite
	col += 2;

	// -gases CO2(g) N2(g)
	col += 5;

	// -kinetic_reactants Albite Pyrite
	//                               k_Albite dk_Albite k_Pyrite dk_Pyrite
	col += 4;

	// -solid_solutions CaSO4 SrSO4
	//                              s_CaSO4 s_SrSO4
	col += 2;


	// 1.name
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("Ca+2"), std::string(v.sVal));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("Calcite"), std::string(v.sVal));

	// 1.type
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("aq"), std::string(v.sVal));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("equi"), std::string(v.sVal));

	// 1.moles
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(9.9177923E-04, v.dVal, ::pow(10., -11));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(9.8444477E-03, v.dVal, ::pow(10., -10));

	// 2.name
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("CaHCO3+"), std::string(v.sVal));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("Ca+2"), std::string(v.sVal));

	// 2.type
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("aq"), std::string(v.sVal));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("aq"), std::string(v.sVal));

	// 2.moles
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(7.5980e-006, v.dVal, ::pow(10., -10));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(1.1371e-003, v.dVal, ::pow(10., -7));


	// 3.name
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("CaCO3"), std::string(v.sVal));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("CaHCO3+"), std::string(v.sVal));

	// 3.type
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("aq"), std::string(v.sVal));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("aq"), std::string(v.sVal));

	// 3.moles
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(6.2155e-007, v.dVal, ::pow(10., -11));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(1.1598e-005, v.dVal, ::pow(10., -9));



	// 4.name
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("CaOH+"), std::string(v.sVal));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("CaCO3"), std::string(v.sVal));

	// 4.type
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("aq"), std::string(v.sVal));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("aq"), std::string(v.sVal));

	// 4.moles
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(1.1819e-009, v.dVal, ::pow(10., -13));
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(6.8668e-006, v.dVal, ::pow(10., -10));


	// 5.name
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_EMPTY, v.type);
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("CaOH+"), std::string(v.sVal));

	// 5.type
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_EMPTY, v.type);
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("aq"), std::string(v.sVal));

	// 5.moles
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_EMPTY, v.type);
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(9.7385e-009, v.dVal, ::pow(10., -13));


	// 6.name
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_EMPTY, v.type);
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_EMPTY, v.type);

	// 6.type
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_EMPTY, v.type);
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_EMPTY, v.type);

	// 6.moles
	//
	++col;
	//   i_soln
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, col, &v));
	ASSERT_EQ(TT_EMPTY, v.type);
	//   react
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, col, &v));
	ASSERT_EQ(TT_EMPTY, v.type);

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}

#if defined(USE_VAR)
	::VarClear(&v);
#endif
}

IPQ_RESULT
SOLUTION(int n, double C, double Ca, double Na)
{
	std::ostringstream oss;

	oss << "SOLUTION 1\n";
	oss << "C " << C << "\n";
	oss << "Ca " << Ca << "\n";
	oss << "Na " << Na << "\n";

	return ::AccumulateLine(n, oss.str().c_str());
}

IPQ_RESULT
EQUILIBRIUM_PHASES(int n, const char* phase, double si, double amount)
{
	std::ostringstream oss;

	oss << "EQUILIBRIUM_PHASES\n";
	oss << phase << " " << si << " " << amount << "\n";
	return ::AccumulateLine(n, oss.str().c_str());
}

IPQ_RESULT
USER_PUNCH(int n, const char* element, int max)
{
	std::ostringstream oss;

	oss << "USER_PUNCH\n";

	oss << "-head ";
	for (int i = 1; i <= max; ++i)
	{
		oss << i << ".name " << i << ".type " << i << ".moles ";
	}
	oss << "\n";
	oss << "-start" << "\n";
	oss << "10 n = sys(\"" << element << "\"" << ", count, names$, types$, moles)" << "\n";
	oss << "20 n = " << max << "\n";
	oss << "30 if count < " << max << " then n = count" << "\n";
	oss << "40 for i = 1 to count" << "\n";
	oss << "50 PUNCH names$(i), types$(i), moles(i)" << "\n";
	oss << "60 next i" << "\n";
	oss << "70 list" << "\n";
	oss << "-end" << "\n";
	oss << "SELECTED_OUTPUT" << "\n";
	oss << "-totals C Ca Na" << "\n";
	oss << "-molalities CO3-2  CaOH+  NaCO3-" << "\n";
	oss << "-activities CO3-2  CaOH+  NaCO3-" << "\n";
	oss << "-equilibrium_phases Calcite" << "\n";
	oss << "-saturation_indices CO2(g) Siderite" << "\n";
	oss << "-gases CO2(g) N2(g)" << "\n";
	oss << "-kinetic_reactants Albite Pyrite" << "\n";
	oss << "-solid_solutions CaSO4 SrSO4" << "\n";

	return ::AccumulateLine(n, oss.str().c_str());
}

IPQ_RESULT
USER_PUNCH_NEH(int n)
{
	std::ostringstream oss;

	oss << "USER_PUNCH\n";

	oss << "-head head0 head1 head2\n";
	oss << "-start" << "\n";
	oss << "10 PUNCH \"have0\", \"have1\", \"have2\"" << "\n";
	oss << "20 PUNCH \"missing0\", \"missing1\", \"missing2\"" << "\n";
	oss << "-end" << "\n";
	oss << "SELECTED_OUTPUT" << "\n";
	oss << "-totals C Ca Na" << "\n";
	oss << "-molalities CO3-2  CaOH+  NaCO3-" << "\n";
	oss << "-activities CO3-2  CaOH+  NaCO3-" << "\n";
	oss << "-equilibrium_phases Calcite" << "\n";
	oss << "-saturation_indices CO2(g) Siderite" << "\n";
	oss << "-gases CO2(g) N2(g)" << "\n";
	oss << "-kinetic_reactants Albite Pyrite" << "\n";
	oss << "-solid_solutions CaSO4 SrSO4" << "\n";

	return ::AccumulateLine(n, oss.str().c_str());
}


IPQ_RESULT
SELECTED_OUTPUT(int n)
{
	std::ostringstream oss;

	oss << "SELECTED_OUTPUT" << "\n";
	oss << "-totals C Ca Na" << "\n";

	return ::AccumulateLine(n, oss.str().c_str());
}

IPQ_RESULT
DUMP(int n)
{
	std::ostringstream oss;
	oss << "DUMP" << "\n";
	oss << "-solution 1" << "\n";
	return ::AccumulateLine(n, oss.str().c_str());
}

TEST(TestIPhreeqcLib, TestGetSelectedOutputColumnCount)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "llnl.dat"));
	ASSERT_EQ(0, ::GetSelectedOutputColumnCount(n));
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));
	ASSERT_EQ(IPQ_OK, ::EQUILIBRIUM_PHASES(n, "calcite", 1.0, 1.0));
	ASSERT_EQ(IPQ_OK, ::USER_PUNCH(n, "Ca", 10));
	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(62, ::GetSelectedOutputColumnCount(n));
	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestAddError)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// make sure initialized to empty
	//
	const char* err = ::GetErrorString(n);
	ASSERT_EQ(std::string(""), std::string(err));

	// make sure initialized to empty
	//
	const char* expected = "TESTING AddError\n";
	ASSERT_EQ(1, ::AddError(n, expected));

	// check 1
	//
	err = ::GetErrorString(n);
	ASSERT_EQ(std::string(expected), std::string(err));

	// check increment
	//
	const char* expected2 = "XXXXXX\n";
	ASSERT_EQ(2, ::AddError(n, expected2));

	// check concatenation
	//
	err = ::GetErrorString(n);
	ASSERT_EQ(std::string(expected) + std::string(expected2), std::string(err));


	// clear errors
	//
	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// make sure back to empty
	//
	err = ::GetErrorString(n);
	ASSERT_EQ(std::string(""), std::string(err));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestAccumulateLine)
{
	// TODO
}

TEST(TestIPhreeqcLib, TestOutputErrorString)
{
	// TODO
}

TEST(TestIPhreeqcLib, TestRunWithCallback)
{
	// TODO
}

TEST(TestIPhreeqcLib, TestRunNoDatabaseLoaded)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(1, ::RunAccumulated(n));

	const char expected[] =
		"ERROR: RunAccumulated: No database is loaded\n";
	const char* err = ::GetErrorString(n);

	ASSERT_EQ(std::string(expected), std::string(err));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestCase1)
{
	// Case 1 (see do_run)
	// pr.punch == TRUE
	// punch.new_def == FALSE
	// output_isopen(OUTPUT_PUNCH) == FALSE
	// selected_output_on == TRUE

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	char SELECTED_OUT[80];
	snprintf(SELECTED_OUT, sizeof(SELECTED_OUT), "selected_1.%d.out", n);

	// remove punch file if it exists
	if (::FileExists(SELECTED_OUT))
	{
		ASSERT_TRUE(::DeleteFile(SELECTED_OUT));
	}
	ASSERT_EQ(false, ::FileExists(SELECTED_OUT));
	ASSERT_EQ(true, ::FileExists("phreeqc.dat"));
	ASSERT_TRUE(::FileSize("phreeqc.dat") > 0);
	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));
	ASSERT_EQ(IPQ_OK, ::USER_PUNCH(n, "Ca", 10));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 1));
	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(true, ::FileExists(SELECTED_OUT));
	ASSERT_TRUE(::FileSize(SELECTED_OUT) > 0);
	ASSERT_EQ(62, ::GetSelectedOutputColumnCount(n));
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 1));
	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(true, ::FileExists(SELECTED_OUT));
	ASSERT_TRUE(::FileSize(SELECTED_OUT) > 0);
	ASSERT_EQ(62, ::GetSelectedOutputColumnCount(n));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestCase2)
{
	// Case 2 (see do_run)
	// pr.punch == TRUE
	// punch.new_def == TRUE
	// output_isopen(OUTPUT_PUNCH) == FALSE
	// selected_output_on == TRUE

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	char SELECTED_OUT[80];
	snprintf(SELECTED_OUT, sizeof(SELECTED_OUT), "selected_1.%d.out", n);

	// remove punch files if they exists
	//
	if (::FileExists(SELECTED_OUT))
	{
		ASSERT_TRUE(::DeleteFile(SELECTED_OUT));
	}
	if (::FileExists("case2.punch"))
	{
		ASSERT_TRUE(::DeleteFile("case2.punch"));
	}
	ASSERT_EQ(false, ::FileExists(SELECTED_OUT));
	ASSERT_EQ(false, ::FileExists("case2.punch"));
	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));
	ASSERT_EQ(IPQ_OK, ::USER_PUNCH(n, "Ca", 10));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "-file case2.punch")); // force have_punch_name to TRUE (see read_selected_ouput)
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 1));
	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(false, ::FileExists(SELECTED_OUT));
	ASSERT_EQ(true, ::FileExists("case2.punch"));
	ASSERT_TRUE(::FileSize("case2.punch") > 0);
	ASSERT_EQ(62, ::GetSelectedOutputColumnCount(n));


	// remove punch files if they exists
	//
	if (::FileExists(SELECTED_OUT))
	{
		ASSERT_TRUE(::DeleteFile(SELECTED_OUT));
	}
	if (::FileExists("case2.punch"))
	{
		ASSERT_TRUE(::DeleteFile("case2.punch"));
	}
	ASSERT_EQ(false, ::FileExists(SELECTED_OUT));
	ASSERT_EQ(false, ::FileExists("case2.punch"));
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));
	ASSERT_EQ(IPQ_OK, ::USER_PUNCH(n, "Ca", 10));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 1));
	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(false, ::FileExists(SELECTED_OUT));
	ASSERT_EQ(true, ::FileExists("case2.punch"));
	ASSERT_TRUE(::FileSize("case2.punch") > 0);
	ASSERT_EQ(62, ::GetSelectedOutputColumnCount(n));

	if (::FileExists("case2.punch"))
	{
		ASSERT_TRUE(::DeleteFile("case2.punch"));
	}
	ASSERT_EQ(false, ::FileExists("case2.punch"));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestPrintSelectedOutputFalse)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	// remove punch files if they exists
	//
	if (::FileExists("selected.out"))
	{
		::DeleteFile("selected.out");
	}
	ASSERT_EQ(false, ::FileExists("selected.out"));

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add selected output block
	ASSERT_EQ(IPQ_OK, ::SELECTED_OUTPUT(n));

	// turn off selected output
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "PRINT; -selected_output false \n"));

	// run
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 1));
	ASSERT_EQ(0, ::RunAccumulated(n));

	ASSERT_EQ(0, ::GetSelectedOutputColumnCount(n));
	ASSERT_EQ(0, ::GetSelectedOutputRowCount(n));


	// reset pr.punch to TRUE
	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add selected output block
	ASSERT_EQ(IPQ_OK, ::SELECTED_OUTPUT(n));

	// run
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 1));
	ASSERT_EQ(0, ::RunAccumulated(n));

	ASSERT_EQ(11, ::GetSelectedOutputColumnCount(n));
	ASSERT_EQ(2, ::GetSelectedOutputRowCount(n));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestOutputFileOnOff)
{
	int onoff[5];
	onoff[0] = 1;  // output_file_on
	onoff[1] = 0;  // error_file_on
	onoff[2] = 0;  // log_file_on
	onoff[3] = 0;  // selected_output_file_on
	onoff[4] = 0;  // dump_file_on
	TestFileOnOff("phreeqc.%d.out", onoff[0], onoff[1], onoff[2], onoff[3], onoff[4]);
}

TEST(TestIPhreeqcLib, TestErrorFileOnOff)
{
	int onoff[5];
	onoff[0] = 0;  // output_file_on
	onoff[1] = 1;  // error_file_on
	onoff[2] = 0;  // log_file_on
	onoff[3] = 0;  // selected_output_file_on
	onoff[4] = 0;  // dump_file_on
	TestFileOnOff("phreeqc.%d.err", onoff[0], onoff[1], onoff[2], onoff[3], onoff[4]);
}

TEST(TestIPhreeqcLib, TestLogFileOnOff)
{
	int onoff[5];
	onoff[0] = 0;  // output_file_on
	onoff[1] = 0;  // error_file_on
	onoff[2] = 1;  // log_file_on
	onoff[3] = 0;  // selected_output_file_on
	onoff[4] = 0;  // dump_file_on
	TestFileOnOff("phreeqc.%d.log", onoff[0], onoff[1], onoff[2], onoff[3], onoff[4]);
}

TEST(TestIPhreeqcLib, TestDumpFileOnOff)
{
	int onoff[5];
	onoff[0] = 0;  // output_file_on
	onoff[1] = 0;  // error_file_on
	onoff[2] = 0;  // log_file_on
	onoff[3] = 0;  // selected_output_file_on
	onoff[4] = 1;  // dump_file_on
	TestFileOnOff("dump.%d.out", onoff[0], onoff[1], onoff[2], onoff[3], onoff[4]);
}

TEST(TestIPhreeqcLib, TestSelOutFileOnOff)
{
	int onoff[5];
	onoff[0] = 0;  // output_file_on
	onoff[1] = 0;  // error_file_on
	onoff[2] = 0;  // log_file_on
	onoff[3] = 1;  // selected_output_file_on
	onoff[4] = 0;  // dump_file_on
	TestFileOnOff("selected_1.%d.out", onoff[0], onoff[1], onoff[2], onoff[3], onoff[4]);
}

//void TestIPhreeqcLib::TestFileOnOff(const char* FILENAME_FORMAT, int output_file_on, int error_file_on, int log_file_on, int selected_output_file_on, int dump_file_on)
void TestFileOnOff(const char* FILENAME_FORMAT, int output_file_on, int error_file_on, int log_file_on, int selected_output_file_on, int dump_file_on)
{
	int dump_string_on = 0;

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	char FILENAME[80];
	snprintf(FILENAME, sizeof(FILENAME), FILENAME_FORMAT, n);

	// remove FILENAME if it exists
	//
	if (::FileExists(FILENAME))
	{
		::DeleteFile(FILENAME);
	}
	ASSERT_EQ(false, ::FileExists(FILENAME));

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add selected output block
	ASSERT_EQ(IPQ_OK, ::SELECTED_OUTPUT(n));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// run all off
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(false, ::FileExists(FILENAME));
	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add selected output block
	ASSERT_EQ(IPQ_OK, ::SELECTED_OUTPUT(n));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// run
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, dump_file_on));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, error_file_on));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, log_file_on));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, output_file_on));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, selected_output_file_on));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, dump_string_on));
	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(true, ::FileExists(FILENAME));
	ASSERT_TRUE(::DeleteFile(FILENAME));
	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add selected output block
	ASSERT_EQ(IPQ_OK, ::SELECTED_OUTPUT(n));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// run
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(false, ::FileExists(FILENAME));
	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add selected output block
	ASSERT_EQ(IPQ_OK, ::SELECTED_OUTPUT(n));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// run
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, dump_file_on));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, error_file_on));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, log_file_on));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, output_file_on));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, selected_output_file_on));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, dump_string_on));
	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(true, ::FileExists(FILENAME));
	ASSERT_TRUE(::DeleteFile(FILENAME));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}


TEST(TestIPhreeqcLib, TestLongHeadings)
{
	char long_header[] = "this_is_a_long_header_0123456789012345678901234567890123456789";
	char long_value[] = "this_is_a_long_value_01234567890123456789012345678901234567890";

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	std::ostringstream oss;
	oss << "SOLUTION" << "\n";

	oss << "SELECTED_OUTPUT" << "\n";
	oss << "-reset false" << "\n";

	oss << "USER_PUNCH" << "\n";
	oss << "-head " << long_header << "\n";
	oss << "-start" << "\n";
	oss << "10 PUNCH \"" << long_value << "\"\n";
	oss << "-end" << "\n";
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, oss.str().c_str()));

	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(0, ::RunAccumulated(n));

	ASSERT_EQ(2, ::GetSelectedOutputRowCount(n));
	ASSERT_EQ(1, ::GetSelectedOutputColumnCount(n));

	CVar v;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, 0, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string(long_header), std::string(v.sVal));

	CVar v1;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, 0, &v1));
	ASSERT_EQ(TT_STRING, v1.type);
	ASSERT_EQ(std::string(long_value), std::string(v1.sVal));

	ASSERT_EQ(IPQ_INVALIDCOL, ::GetSelectedOutputValue(n, 1, 1, &v1));
	ASSERT_EQ(IPQ_INVALIDROW, ::GetSelectedOutputValue(n, 2, 0, &v1));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestDatabaseKeyword)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 1));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(1, ::RunFile(n, "dump"));

	const char* exp_errs =
		"ERROR: Gas not found in PHASES database, Amm(g).\n"
		"ERROR: Calculations terminating due to input errors.\n";

	const char* err = ::GetErrorString(n);
	ASSERT_EQ(std::string(exp_errs), std::string(err));

	const char* exp_warn =
		"WARNING: DATABASE keyword is ignored by IPhreeqc.\n"
		"WARNING: Cell-lengths were read for 1 cells. Last value is used till cell 100.\n"
		"WARNING: No dispersivities were read; disp = 0 assumed.\n"
		"WARNING: Could not find element in database, Amm.\n"
		"	Concentration is set to zero.\n";

	const char* warn = ::GetWarningString(n);
	ASSERT_EQ(std::string(exp_warn), std::string(warn));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestDumpString)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// run
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 1));
	ASSERT_EQ(0, ::RunAccumulated(n));

	const char* dump_str = ::GetDumpString(n);

	EXPECT_THAT(dump_str, HasSubstr("SOLUTION_RAW"));
	EXPECT_THAT(dump_str, HasSubstr("-temp"));
	EXPECT_THAT(dump_str, HasSubstr("-pressure"));
	EXPECT_THAT(dump_str, HasSubstr("-potential"));
	EXPECT_THAT(dump_str, HasSubstr("-total_h"));
	EXPECT_THAT(dump_str, HasSubstr("-total_o"));
	EXPECT_THAT(dump_str, HasSubstr("-cb"));
	EXPECT_THAT(dump_str, HasSubstr("-density"));
	EXPECT_THAT(dump_str, HasSubstr("-viscosity"));
	EXPECT_THAT(dump_str, HasSubstr("-totals"));
	EXPECT_THAT(dump_str, HasSubstr(" C(4) "));
	EXPECT_THAT(dump_str, HasSubstr(" Ca "));
	EXPECT_THAT(dump_str, HasSubstr(" H(0) "));
	EXPECT_THAT(dump_str, HasSubstr(" Na "));
	EXPECT_THAT(dump_str, HasSubstr("-pH"));
	EXPECT_THAT(dump_str, HasSubstr("-pe"));
	EXPECT_THAT(dump_str, HasSubstr("-mu"));
	EXPECT_THAT(dump_str, HasSubstr("-ah2o"));
	EXPECT_THAT(dump_str, HasSubstr("-mass_water"));
	EXPECT_THAT(dump_str, HasSubstr("-total_alkalinity"));
	EXPECT_THAT(dump_str, HasSubstr("-activities"));
	EXPECT_THAT(dump_str, HasSubstr(" C(-4) "));
	EXPECT_THAT(dump_str, HasSubstr(" C(4) "));
	EXPECT_THAT(dump_str, HasSubstr(" Ca "));
	EXPECT_THAT(dump_str, HasSubstr(" E "));
	EXPECT_THAT(dump_str, HasSubstr(" H(0) "));
	EXPECT_THAT(dump_str, HasSubstr(" Na "));
	EXPECT_THAT(dump_str, HasSubstr(" O(0) "));
	EXPECT_THAT(dump_str, HasSubstr("-gammas"));
	EXPECT_THAT(dump_str, HasSubstr("USE mix none"));
	EXPECT_THAT(dump_str, HasSubstr("USE reaction none"));
	EXPECT_THAT(dump_str, HasSubstr("USE reaction_temperature none"));
	EXPECT_THAT(dump_str, HasSubstr("USE reaction_pressure none"));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestGetDumpStringLineCount)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// run
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 1));
	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(35, ::GetDumpStringLineCount(n));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestGetDumpStringLine)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// run
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 1));
	ASSERT_EQ(0,      ::RunAccumulated(n));
	ASSERT_EQ(35,     ::GetDumpStringLineCount(n));

	int line = 0;

	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("SOLUTION_RAW"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-temp"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-pressure"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-potential"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-total_h"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-total_o"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-cb"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-density"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-viscosity"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-viscos_0"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-totals"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr(" C(4) "));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr(" Ca "));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr(" H(0) "));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr(" Na "));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-pH"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-pe"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-mu"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-ah2o"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-mass_water"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-soln_vol"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-total_alkalinity"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-activities"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr(" C(-4) "));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr(" C(4) "));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr(" Ca "));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr(" E "));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr(" H(0) "));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr(" Na "));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr(" O(0) "));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("-gammas"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("USE mix none"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("USE reaction none"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("USE reaction_temperature none"));
	EXPECT_THAT(::GetDumpStringLine(n, line++), HasSubstr("USE reaction_pressure none"));

	// remaining lines should be empty
	ASSERT_EQ(std::string(""), std::string(::GetDumpStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetDumpStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetDumpStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetDumpStringLine(n, line++)));

	// negative lines should be empty
	ASSERT_EQ(std::string(""), std::string(::GetDumpStringLine(n, -1)));
	ASSERT_EQ(std::string(""), std::string(::GetDumpStringLine(n, -2)));
	ASSERT_EQ(std::string(""), std::string(::GetDumpStringLine(n, -3)));
	ASSERT_EQ(std::string(""), std::string(::GetDumpStringLine(n, -4)));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestGetComponentCount)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// run
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(3, ::GetComponentCount(n));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestGetComponent)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// run
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(3, ::GetComponentCount(n));

	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, -2)));
	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, -1)));

	ASSERT_EQ(std::string("C"), std::string(::GetComponent(n, 0)));
	ASSERT_EQ(std::string("Ca"), std::string(::GetComponent(n, 1)));
	ASSERT_EQ(std::string("Na"), std::string(::GetComponent(n, 2)));

	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, 3)));
	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, 4)));
	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, 5)));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// run
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(3, ::GetComponentCount(n));

	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, -2)));
	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, -1)));

	ASSERT_EQ(std::string("C"), std::string(::GetComponent(n, 0)));
	ASSERT_EQ(std::string("Ca"), std::string(::GetComponent(n, 1)));
	ASSERT_EQ(std::string("Na"), std::string(::GetComponent(n, 2)));

	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, 3)));
	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, 4)));
	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, 5)));

	// clear using LoadDatabase
	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));
	ASSERT_EQ(0, ::GetComponentCount(n));

	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, -2)));
	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, -1)));

	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, 0)));
	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, 1)));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));
	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(3, ::GetComponentCount(n));

	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, -2)));
	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, -1)));

	ASSERT_EQ(std::string("C"), std::string(::GetComponent(n, 0)));
	ASSERT_EQ(std::string("Ca"), std::string(::GetComponent(n, 1)));
	ASSERT_EQ(std::string("Na"), std::string(::GetComponent(n, 2)));

	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, 3)));
	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, 4)));
	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, 5)));


	// clear using LoadDatabaseString
	ASSERT_EQ(0, ::LoadDatabaseString(n, ex15_dat));
	ASSERT_EQ(0, ::GetComponentCount(n));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));
	ASSERT_EQ(0, ::RunAccumulated(n));

	// ex15.dat doesn't have Ca
	ASSERT_EQ(2, ::GetWarningStringLineCount(n));
	ASSERT_EQ(std::string("WARNING: Could not find element in database, Ca."), std::string(::GetWarningStringLine(n, 0)));
	ASSERT_EQ(std::string("\tConcentration is set to zero."), std::string(::GetWarningStringLine(n, 1)));

	ASSERT_EQ(2, ::GetComponentCount(n));

	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, -2)));
	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, -1)));

	ASSERT_EQ(std::string("C"), std::string(::GetComponent(n, 0)));
	ASSERT_EQ(std::string("Na"), std::string(::GetComponent(n, 1)));

	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, 2)));
	ASSERT_EQ(std::string(""), std::string(::GetComponent(n, 3)));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestGetErrorStringLine)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(1, ::RunFile(n, "dump"));
#if 0
	ASSERT_EQ(8, ::GetErrorStringLineCount(n));
#else
	ASSERT_EQ(2, ::GetErrorStringLineCount(n));
#endif

	int line = 0;
#if 0
	ASSERT_EQ(std::string("WARNING: DATABASE keyword is ignored by IPhreeqc."), std::string(::GetErrorStringLine(n, line++)));
	ASSERT_EQ(std::string("WARNING: Cell-lengths were read for 1 cells. Last value is used till cell 100."), std::string(::GetErrorStringLine(n, line++)));
	ASSERT_EQ(std::string("WARNING: No dispersivities were read; disp = 0 assumed."), std::string(::GetErrorStringLine(n, line++)));
#endif
	ASSERT_EQ(std::string("ERROR: Gas not found in PHASES database, Amm(g)."), std::string(::GetErrorStringLine(n, line++)));
#if 0
	ASSERT_EQ(std::string("WARNING: Could not find element in database, Amm."), std::string(::GetErrorStringLine(n, line++)));
	ASSERT_EQ(std::string("\tConcentration is set to zero."), std::string(::GetErrorStringLine(n, line++)));
#endif
	ASSERT_EQ(std::string("ERROR: Calculations terminating due to input errors."), std::string(::GetErrorStringLine(n, line++)));
#if 0
	ASSERT_EQ(std::string("Stopping."), std::string(::GetErrorStringLine(n, line++)));
#endif

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestErrorFileOn)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	char FILENAME[80];
	snprintf(FILENAME, sizeof(FILENAME), "phreeqc.%d.err", n);

	if (::FileExists(FILENAME))
	{
		::DeleteFile(FILENAME);
	}

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 1));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(1, ::RunFile(n, "dump"));
	ASSERT_EQ(true, ::FileExists(FILENAME));

	std::string lines[10];
	std::ifstream ifs(FILENAME);

	size_t i = 0;
	while (i < sizeof(lines) / sizeof(lines[0]) && std::getline(ifs, lines[i]))
	{
		++i;
	}

	ASSERT_EQ((size_t)8, i);

	ASSERT_EQ(std::string("WARNING: DATABASE keyword is ignored by IPhreeqc."), lines[0]);
	ASSERT_EQ(std::string("WARNING: Cell-lengths were read for 1 cells. Last value is used till cell 100."), lines[1]);
	ASSERT_EQ(std::string("WARNING: No dispersivities were read; disp = 0 assumed."), lines[2]);
	ASSERT_EQ(std::string("ERROR: Gas not found in PHASES database, Amm(g)."), lines[3]);
	ASSERT_EQ(std::string("WARNING: Could not find element in database, Amm."), lines[4]);
	ASSERT_EQ(std::string("\tConcentration is set to zero."), lines[5]);
	ASSERT_EQ(std::string("ERROR: Calculations terminating due to input errors."), lines[6]);
	ASSERT_EQ(std::string("Stopping."), lines[7]);
	ASSERT_EQ(std::string(""), lines[8]);
	ASSERT_EQ(std::string(""), lines[9]);

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestLogFileOn)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	char FILENAME[80];
	snprintf(FILENAME, sizeof(FILENAME), "phreeqc.%d.log", n);

	if (::FileExists(FILENAME))
	{
		::DeleteFile(FILENAME);
	}

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 1));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));

	ASSERT_EQ(1, ::RunFile(n, "dump"));
	ASSERT_EQ(true, ::FileExists(FILENAME));

	std::string lines[10];
	std::ifstream ifs(FILENAME);

	size_t i = 0;
	while (i < sizeof(lines) / sizeof(lines[0]) && std::getline(ifs, lines[i]))
	{
		++i;
	}

	ASSERT_EQ((size_t)6, i);

	ASSERT_EQ(std::string("WARNING: Cell-lengths were read for 1 cells. Last value is used till cell 100."), lines[0]);
	ASSERT_EQ(std::string("WARNING: No dispersivities were read; disp = 0 assumed."), lines[1]);
	ASSERT_EQ(std::string("ERROR: Gas not found in PHASES database, Amm(g)."), lines[2]);
	ASSERT_EQ(std::string("WARNING: Could not find element in database, Amm."), lines[3]);
	ASSERT_EQ(std::string("\tConcentration is set to zero."), lines[4]);
	ASSERT_EQ(std::string("ERROR: Calculations terminating due to input errors."), lines[5]);
	ASSERT_EQ(std::string(""), lines[8]);
	ASSERT_EQ(std::string(""), lines[9]);

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestGetWarningStringLine)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(1, ::RunFile(n, "dump"));

	ASSERT_EQ(5, ::GetWarningStringLineCount(n));
	ASSERT_EQ(std::string("WARNING: DATABASE keyword is ignored by IPhreeqc."), std::string(::GetWarningStringLine(n, 0)));
	ASSERT_EQ(std::string("WARNING: Cell-lengths were read for 1 cells. Last value is used till cell 100."), std::string(::GetWarningStringLine(n, 1)));
	ASSERT_EQ(std::string("WARNING: No dispersivities were read; disp = 0 assumed."), std::string(::GetWarningStringLine(n, 2)));
	ASSERT_EQ(std::string("WARNING: Could not find element in database, Amm."), std::string(::GetWarningStringLine(n, 3)));
	ASSERT_EQ(std::string("	Concentration is set to zero."), std::string(::GetWarningStringLine(n, 4)));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestPitzer)
{
	ASSERT_EQ(true, ::FileExists("pitzer.dat"));
	ASSERT_TRUE(::FileSize("pitzer.dat") > 0);

	int id = ::CreateIPhreeqc();
	ASSERT_TRUE(id >= 0);

	ASSERT_EQ(0, ::LoadDatabase(id, "pitzer.dat"));

	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "SOLUTION 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "units mol/kgw"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "pH    7"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "temp  25"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "Mg    0.1    charge"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "Cl    0"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "S(6)  1.0612244897959"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "EQUILIBRIUM_PHASES 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "Halite       0 10"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "Sylvite      0 10"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "Bischofite   0 0"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "Carnallite   0 0"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "Epsomite     0 0"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "Kieserite    0 0"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "END"));

	// this fails -r4375
	// fixed in -r4380
	ASSERT_EQ(0, ::RunAccumulated(id));

	if (id >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(id));
	}
}

TEST(TestIPhreeqcLib, TestClearAccumulatedLines)
{
	ASSERT_EQ(true, ::FileExists("wateq4f.dat"));
	ASSERT_TRUE(::FileSize("wateq4f.dat") > 0);

	int id = ::CreateIPhreeqc();
	ASSERT_TRUE(id >= 0);

	ASSERT_EQ(0, ::LoadDatabase(id, "wateq4f.dat"));

	// phreeqc can now handle pH of -2 (-r5885)
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "SOLUTION 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "	pH	7"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "	Na	1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "PHASES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "	Fix_H+"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "	H+ = H+"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "	log_k	0"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "EQUILIBRIUM_PHASES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "	Fix_H+ -10 HCl	10"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "END"));

	ASSERT_EQ(1, ::RunAccumulated(id));

	ASSERT_EQ(1, ::GetErrorStringLineCount(id));

	const char expected[] =
		"ERROR: Numerical method failed on all combinations of convergence parameters, cell/soln/mix 1\n";
	const char* err = ::GetErrorString(id);

	ASSERT_EQ(std::string(expected), std::string(err));

	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "SOLUTION 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "pH 2"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "END"));

	// RunAccumulated clears accumulated lines on next call to AccumulateLine
	//
	ASSERT_EQ(0, ::RunAccumulated(id));
	ASSERT_EQ(0, ::GetErrorStringLineCount(id));

	ASSERT_EQ(IPQ_OK, ::ClearAccumulatedLines(id));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "SOLUTION 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "pH 2"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "END"));

	ASSERT_EQ(0, ::RunAccumulated(id));

	if (id >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(id));
	}
}

TEST(TestIPhreeqcLib, TestSetDumpFileName)
{
	char DUMP_FILENAME[80];
	snprintf(DUMP_FILENAME, sizeof(DUMP_FILENAME), "dump.%06d.out", ::rand());
	if (::FileExists(DUMP_FILENAME))
	{
		::DeleteFile(DUMP_FILENAME);
	}

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// run
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 1));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileName(n, DUMP_FILENAME));

	ASSERT_EQ(std::string(DUMP_FILENAME), std::string(::GetDumpFileName(n)));

	ASSERT_EQ(0, ::RunAccumulated(n));

	ASSERT_EQ(std::string(DUMP_FILENAME), std::string(::GetDumpFileName(n)));

	ASSERT_EQ(true, ::FileExists(DUMP_FILENAME));

	std::string lines[35];
	std::ifstream ifs(DUMP_FILENAME);

	size_t i = 0;
	while (i < sizeof(lines) / sizeof(lines[0]) && std::getline(ifs, lines[i]))
	{
		++i;
	}
	ifs.close();

	int line = 0;
	EXPECT_THAT(lines[line++], HasSubstr("SOLUTION_RAW"));
	EXPECT_THAT(lines[line++], HasSubstr("-temp"));
	EXPECT_THAT(lines[line++], HasSubstr("-pressure"));
	EXPECT_THAT(lines[line++], HasSubstr("-potential"));
	EXPECT_THAT(lines[line++], HasSubstr("-total_h"));
	EXPECT_THAT(lines[line++], HasSubstr("-total_o"));
	EXPECT_THAT(lines[line++], HasSubstr("-cb"));
	EXPECT_THAT(lines[line++], HasSubstr("-density"));
	EXPECT_THAT(lines[line++], HasSubstr("-viscosity"));
	EXPECT_THAT(lines[line++], HasSubstr("-viscos_0"));
	EXPECT_THAT(lines[line++], HasSubstr("-totals"));
	EXPECT_THAT(lines[line++], HasSubstr(" C(4) "));
	EXPECT_THAT(lines[line++], HasSubstr(" Ca "));
	EXPECT_THAT(lines[line++], HasSubstr(" H(0) "));
	EXPECT_THAT(lines[line++], HasSubstr(" Na "));
	EXPECT_THAT(lines[line++], HasSubstr("-pH"));
	EXPECT_THAT(lines[line++], HasSubstr("-pe"));
	EXPECT_THAT(lines[line++], HasSubstr("-mu"));
	EXPECT_THAT(lines[line++], HasSubstr("-ah2o"));
	EXPECT_THAT(lines[line++], HasSubstr("-mass_water"));
	EXPECT_THAT(lines[line++], HasSubstr("-soln_vol"));
	EXPECT_THAT(lines[line++], HasSubstr("-total_alkalinity"));
	EXPECT_THAT(lines[line++], HasSubstr("-activities"));
	EXPECT_THAT(lines[line++], HasSubstr(" C(-4) "));
	EXPECT_THAT(lines[line++], HasSubstr(" C(4) "));
	EXPECT_THAT(lines[line++], HasSubstr(" Ca "));
	EXPECT_THAT(lines[line++], HasSubstr(" E "));
	EXPECT_THAT(lines[line++], HasSubstr(" H(0) "));
	EXPECT_THAT(lines[line++], HasSubstr(" Na "));
	EXPECT_THAT(lines[line++], HasSubstr(" O(0) "));
	EXPECT_THAT(lines[line++], HasSubstr("-gammas"));
	EXPECT_THAT(lines[line++], HasSubstr("USE mix none"));
	EXPECT_THAT(lines[line++], HasSubstr("USE reaction none"));
	EXPECT_THAT(lines[line++], HasSubstr("USE reaction_temperature none"));
	EXPECT_THAT(lines[line++], HasSubstr("USE reaction_pressure none"));

	if (::FileExists(DUMP_FILENAME))
	{
		::DeleteFile(DUMP_FILENAME);
	}

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestSetOutputFileName)
{
	char OUTPUT_FILENAME[80];
	snprintf(OUTPUT_FILENAME, sizeof(OUTPUT_FILENAME), "output.%06d.out", ::rand());
	if (::FileExists(OUTPUT_FILENAME))
	{
		::DeleteFile(OUTPUT_FILENAME);
	}

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat.old"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// run
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 1));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileName(n, OUTPUT_FILENAME));

	ASSERT_EQ(std::string(OUTPUT_FILENAME), std::string(::GetOutputFileName(n)));


	ASSERT_EQ(0, ::RunAccumulated(n));

	ASSERT_EQ(true, ::FileExists(OUTPUT_FILENAME));

	std::string lines[200];
	std::ifstream ifs(OUTPUT_FILENAME);

	size_t i = 0;
	while (i < sizeof(lines) / sizeof(lines[0]) && std::getline(ifs, lines[i]))
	{
		++i;
	}

#ifndef TESTING
	ASSERT_EQ((size_t)100, i);
#else
	ASSERT_EQ((size_t)96, i);
#endif

	int line = 0;

	EXPECT_THAT(lines[line++], HasSubstr("------------------------------------"));
	EXPECT_THAT(lines[line++], HasSubstr("Reading input data for simulation 1."));
	EXPECT_THAT(lines[line++], HasSubstr("------------------------------------"));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("	SOLUTION 1"));
	EXPECT_THAT(lines[line++], HasSubstr("	C 1"));
	EXPECT_THAT(lines[line++], HasSubstr("	Ca 1"));
	EXPECT_THAT(lines[line++], HasSubstr("	Na 1"));
	EXPECT_THAT(lines[line++], HasSubstr("	DUMP"));
	EXPECT_THAT(lines[line++], HasSubstr("	-solution 1"));
	EXPECT_THAT(lines[line++], HasSubstr("-------------------------------------------"));
	EXPECT_THAT(lines[line++], HasSubstr("Beginning of initial solution calculations."));
	EXPECT_THAT(lines[line++], HasSubstr("-------------------------------------------"));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("Initial solution 1.	"));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("-----------------------------Solution composition--------"));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("	Elements           Molality       Moles"));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("	C                "));
	EXPECT_THAT(lines[line++], HasSubstr("	Ca               "));
	EXPECT_THAT(lines[line++], HasSubstr("	Na               "));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("----------------------------Description of solution------"));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("                                       pH  =  "));
	EXPECT_THAT(lines[line++], HasSubstr("                                       pe  =  "));
	EXPECT_THAT(lines[line++], HasSubstr("                        Activity of water  =  "));
	EXPECT_THAT(lines[line++], HasSubstr("                 Ionic strength (mol/kgw)  =  "));
	EXPECT_THAT(lines[line++], HasSubstr("                       Mass of water (kg)  =  "));
	EXPECT_THAT(lines[line++], HasSubstr("                 Total alkalinity (eq/kg)  =  "));
	EXPECT_THAT(lines[line++], HasSubstr("                       Total CO2 (mol/kg)  =  "));
	EXPECT_THAT(lines[line++], HasSubstr("                         Temperature (°C)  =  "));
	EXPECT_THAT(lines[line++], HasSubstr("                  Electrical balance (eq)  =  "));
	EXPECT_THAT(lines[line++], HasSubstr(" Percent error, 100*(Cat-|An|)/(Cat+|An|)  =  "));
	EXPECT_THAT(lines[line++], HasSubstr("                               Iterations  =  "));
	EXPECT_THAT(lines[line++], HasSubstr("                                  Total H  = "));
	EXPECT_THAT(lines[line++], HasSubstr("                                  Total O  = "));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("----------------------------Distribution of species----------------------------"));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("                                               Log       Log       Log    mole V"));
	EXPECT_THAT(lines[line++], HasSubstr("   Species          Molality    Activity  Molality  Activity     Gamma   cm³/mol"));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("   OH- "));
	EXPECT_THAT(lines[line++], HasSubstr("   H+ "));
	EXPECT_THAT(lines[line++], HasSubstr("   H2O "));
	EXPECT_THAT(lines[line++], HasSubstr("C(-4) "));
	EXPECT_THAT(lines[line++], HasSubstr("   CH4 "));
	EXPECT_THAT(lines[line++], HasSubstr("C(4) "));
	EXPECT_THAT(lines[line++], HasSubstr("   HCO3- "));
	EXPECT_THAT(lines[line++], HasSubstr("   CO2 "));
	EXPECT_THAT(lines[line++], HasSubstr("   CaHCO3+ "));
	EXPECT_THAT(lines[line++], HasSubstr("   CaCO3 "));
	EXPECT_THAT(lines[line++], HasSubstr("   CO3-2 "));
	EXPECT_THAT(lines[line++], HasSubstr("   NaHCO3 "));
	EXPECT_THAT(lines[line++], HasSubstr("   NaCO3- "));
	EXPECT_THAT(lines[line++], HasSubstr("Ca "));
	EXPECT_THAT(lines[line++], HasSubstr("   Ca+2 "));
	EXPECT_THAT(lines[line++], HasSubstr("   CaHCO3+ "));
	EXPECT_THAT(lines[line++], HasSubstr("   CaCO3 "));
	EXPECT_THAT(lines[line++], HasSubstr("   CaOH+ "));
	EXPECT_THAT(lines[line++], HasSubstr("H(0) "));
	EXPECT_THAT(lines[line++], HasSubstr("   H2 "));
	EXPECT_THAT(lines[line++], HasSubstr("Na "));
	EXPECT_THAT(lines[line++], HasSubstr("   Na+ "));
	EXPECT_THAT(lines[line++], HasSubstr("   NaHCO3 "));
	EXPECT_THAT(lines[line++], HasSubstr("   NaCO3- "));
	EXPECT_THAT(lines[line++], HasSubstr("   NaOH "));
	EXPECT_THAT(lines[line++], HasSubstr("O(0) "));
	EXPECT_THAT(lines[line++], HasSubstr("   O2 "));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("------------------------------Saturation indices-------------------------------"));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("  Phase               SI** log IAP   log K(298 K,   1 atm)"));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("  Aragonite"));
	EXPECT_THAT(lines[line++], HasSubstr("  Calcite"));
	EXPECT_THAT(lines[line++], HasSubstr("  CH4(g)"));
	EXPECT_THAT(lines[line++], HasSubstr("  CO2(g)"));
	EXPECT_THAT(lines[line++], HasSubstr("  H2(g)"));
	EXPECT_THAT(lines[line++], HasSubstr("  H2O(g)"));
	EXPECT_THAT(lines[line++], HasSubstr("  O2(g)"));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("**For a gas, SI = log10(fugacity). Fugacity = pressure * phi / 1 atm."));
	EXPECT_THAT(lines[line++], HasSubstr("  For ideal gases, phi = 1."));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("------------------"));
	EXPECT_THAT(lines[line++], HasSubstr("End of simulation."));
	EXPECT_THAT(lines[line++], HasSubstr("------------------"));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr("------------------------------------"));
	EXPECT_THAT(lines[line++], HasSubstr("Reading input data for simulation 2."));
	EXPECT_THAT(lines[line++], HasSubstr("------------------------------------"));
	EXPECT_THAT(lines[line++], HasSubstr(""));
#ifndef TESTING
	EXPECT_THAT(lines[line++], HasSubstr("----------------"));
	EXPECT_THAT(lines[line++], HasSubstr("End of Run after "));
	EXPECT_THAT(lines[line++], HasSubstr("----------------"));
#endif
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr(""));
	EXPECT_THAT(lines[line++], HasSubstr(""));
///#endif

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}

	if (::FileExists(OUTPUT_FILENAME))
	{
		::DeleteFile(OUTPUT_FILENAME);
	}
}

TEST(TestIPhreeqcLib, TestOutputStringOnOff)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(false, ::GetOutputStringOn(n) != 0);

	ASSERT_EQ(IPQ_OK, ::SetOutputStringOn(n, 1));
	ASSERT_EQ(true, ::GetOutputStringOn(n) != 0);

	ASSERT_EQ(IPQ_OK, ::SetOutputStringOn(n, 0));
	ASSERT_EQ(false, ::GetOutputStringOn(n) != 0);

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestGetOutputString)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	char OUTPUT_FILENAME[80];
	snprintf(OUTPUT_FILENAME, sizeof(OUTPUT_FILENAME), "output.%06d.out", ::rand());
	if (::FileExists(OUTPUT_FILENAME))
	{
		::DeleteFile(OUTPUT_FILENAME);
	}
	ASSERT_EQ(false, ::FileExists(OUTPUT_FILENAME));

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// run
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 1));
	ASSERT_EQ(IPQ_OK, ::SetOutputStringOn(n, 1));

	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileName(n, OUTPUT_FILENAME));

	ASSERT_EQ(0, ::RunAccumulated(n));

	ASSERT_EQ(true, ::FileExists(OUTPUT_FILENAME));

	{
		std::ifstream ifs(OUTPUT_FILENAME);
		std::string fline((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

		std::string sline(::GetOutputString(n));
		ASSERT_TRUE(sline.size() > 0);

		ASSERT_EQ(fline, sline);
	}

	if (::FileExists(OUTPUT_FILENAME))
	{
		::DeleteFile(OUTPUT_FILENAME);
	}
	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestGetOutputStringLineCount)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::GetOutputStringLineCount(n));
	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat.old"));
	ASSERT_EQ(0, ::GetOutputStringLineCount(n));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// run
	ASSERT_EQ(IPQ_OK, ::SetOutputStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));

	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(0, ::GetOutputStringLineCount(n));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	ASSERT_EQ(IPQ_OK, ::SetOutputStringOn(n, 1));

	ASSERT_EQ(0, ::RunAccumulated(n));
#ifndef TESTING
	ASSERT_EQ(100, ::GetOutputStringLineCount(n));
#else
	ASSERT_EQ(96, ::GetOutputStringLineCount(n));
#endif

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	ASSERT_EQ(IPQ_OK, ::SetOutputStringOn(n, 0));

	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(0, ::GetOutputStringLineCount(n));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestGetOutputStringLine)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::GetOutputStringLineCount(n));
	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat.old"));
	ASSERT_EQ(0, ::GetOutputStringLineCount(n));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// run
	ASSERT_EQ(IPQ_OK, ::SetOutputStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));

	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(0, ::GetOutputStringLineCount(n));

	int line = 0;
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, line++)));

	// negative lines should be empty
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, -1)));
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, -2)));
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, -3)));
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, -4)));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	ASSERT_EQ(IPQ_OK, ::SetOutputStringOn(n, 1));

	ASSERT_EQ(0, ::RunAccumulated(n));
#ifndef TESTING
	ASSERT_EQ(100, ::GetOutputStringLineCount(n));
#else
	ASSERT_EQ(96, ::GetOutputStringLineCount(n));
#endif

	ASSERT_EQ(std::string("------------------------------------"), std::string(::GetOutputStringLine(n, 0)));
	ASSERT_EQ(std::string("Reading input data for simulation 1."), std::string(::GetOutputStringLine(n, 1)));
	ASSERT_EQ(std::string("------------------------------------"), std::string(::GetOutputStringLine(n, 2)));
	ASSERT_EQ(std::string("-----------------------------Solution composition------------------------------"), std::string(::GetOutputStringLine(n, 16)));
	ASSERT_EQ(std::string("----------------------------Description of solution----------------------------"), std::string(::GetOutputStringLine(n, 24)));
	ASSERT_EQ(std::string("----------------------------Distribution of species----------------------------"), std::string(::GetOutputStringLine(n, 40)));
	ASSERT_EQ(std::string("------------------------------Saturation indices-------------------------------"), std::string(::GetOutputStringLine(n, 73)));
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, 100)));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	ASSERT_EQ(IPQ_OK, ::SetOutputStringOn(n, 0));

	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(0, ::GetOutputStringLineCount(n));

	line = 0;
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, line++)));

	// negative lines should be empty
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, -1)));
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, -2)));
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, -3)));
	ASSERT_EQ(std::string(""), std::string(::GetOutputStringLine(n, -4)));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestSetLogFileName)
{
	char LOG_FILENAME[80];
	snprintf(LOG_FILENAME, sizeof(LOG_FILENAME), "log.%06d.out", ::rand());
	if (::FileExists(LOG_FILENAME))
	{
		::DeleteFile(LOG_FILENAME);
	}

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// add knobs
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "KNOBS"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "\t-logfile TRUE"));

	// run
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 1));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileName(n, LOG_FILENAME));


	ASSERT_EQ(0, ::RunAccumulated(n));

	ASSERT_EQ(true, ::FileExists(LOG_FILENAME));

	std::string lines[33];
	std::ifstream ifs(LOG_FILENAME);

	size_t i = 0;
	while (i < sizeof(lines) / sizeof(lines[0]) && std::getline(ifs, lines[i]))
	{
		++i;
	}

#ifndef TESTING
	ASSERT_EQ((size_t)25, i);
#else
	ASSERT_EQ((size_t)21, i);
#endif

	int line = 0;
	ASSERT_EQ(std::string("-------------------------------------------"), lines[line++]);
	ASSERT_EQ(std::string("Beginning of initial solution calculations."), lines[line++]);
	ASSERT_EQ(std::string("-------------------------------------------"), lines[line++]);
	ASSERT_EQ(std::string(""), lines[line++]);
	ASSERT_EQ(std::string("Initial solution 1.	"), lines[line++]);
	ASSERT_EQ(std::string(""), lines[line++]);
	ASSERT_EQ(std::string("Iterations in revise_guesses: 2"), lines[line++]);
	ASSERT_EQ(std::string(""), lines[line++]);
	ASSERT_EQ(std::string("Number of infeasible solutions: 0"), lines[line++]);
	ASSERT_EQ(std::string("Number of basis changes: 0"), lines[line++]);
	ASSERT_EQ(std::string(""), lines[line++]);
	ASSERT_EQ(std::string("Number of iterations: 8"), lines[line++]);
	ASSERT_EQ(std::string(""), lines[line++]);
	ASSERT_EQ(std::string("------------------"), lines[line++]);
	ASSERT_EQ(std::string("End of simulation."), lines[line++]);
	ASSERT_EQ(std::string("------------------"), lines[line++]);
	ASSERT_EQ(std::string(""), lines[line++]);
	ASSERT_EQ(std::string("------------------------------------"), lines[line++]);
	ASSERT_EQ(std::string("Reading input data for simulation 2."), lines[line++]);
	ASSERT_EQ(std::string("------------------------------------"), lines[line++]);
	ASSERT_EQ(std::string(""), lines[line++]);
	line++;
	line++;
	line++;
	ASSERT_EQ(std::string(""), lines[line++]);

	if (::FileExists(LOG_FILENAME))
	{
		::DeleteFile(LOG_FILENAME);
	}

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestLogStringOnOff)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(false, ::GetLogStringOn(n) != 0);

	ASSERT_EQ(IPQ_OK, ::SetLogStringOn(n, 1));
	ASSERT_EQ(true, ::GetLogStringOn(n) != 0);

	ASSERT_EQ(IPQ_OK, ::SetLogStringOn(n, 0));
	ASSERT_EQ(false, ::GetLogStringOn(n) != 0);

	ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
}

TEST(TestIPhreeqcLib, TestGetLogString)
{
	char LOG_FILENAME[80];
	snprintf(LOG_FILENAME, sizeof(LOG_FILENAME), "log.%06d.out", ::rand());
	if (::FileExists(LOG_FILENAME))
	{
		::DeleteFile(LOG_FILENAME);
	}
	ASSERT_EQ(false, ::FileExists(LOG_FILENAME));

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// add knobs
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "KNOBS"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "\t-logfile TRUE"));


	// run
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 1));
	ASSERT_EQ(IPQ_OK, ::SetLogStringOn(n, 1));

	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));

	ASSERT_EQ(IPQ_OK, ::SetLogFileName(n, LOG_FILENAME));
	ASSERT_EQ(std::string(LOG_FILENAME), std::string(::GetLogFileName(n)));

	ASSERT_EQ(0, ::RunAccumulated(n));

	ASSERT_EQ(std::string(LOG_FILENAME), std::string(::GetLogFileName(n)));

	ASSERT_EQ(true, ::FileExists(LOG_FILENAME));

	{
		std::ifstream ifs(LOG_FILENAME);
		std::string fline((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

		std::string sline(::GetLogString(n));
		ASSERT_TRUE(sline.size() > 0);

		ASSERT_EQ(fline, sline);
	}

	if (::FileExists(LOG_FILENAME))
	{
		::DeleteFile(LOG_FILENAME);
	}

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestGetLogStringLineCount)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::GetLogStringLineCount(n));

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	ASSERT_EQ(0, ::GetLogStringLineCount(n));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// add knobs
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "KNOBS"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "\t-logfile TRUE"));

	// run
	ASSERT_EQ(IPQ_OK, ::SetLogStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));

	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(0, ::GetLogStringLineCount(n));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// add knobs
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "KNOBS"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "\t-logfile TRUE"));

	ASSERT_EQ(IPQ_OK, ::SetLogStringOn(n, 1));

	ASSERT_EQ(0, ::RunAccumulated(n));
#ifndef TESTING
	ASSERT_EQ(29, ::GetLogStringLineCount(n));
#else
	ASSERT_EQ(25, ::GetLogStringLineCount(n));
#endif

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// add knobs
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "KNOBS"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "\t-logfile TRUE"));

	ASSERT_EQ(IPQ_OK, ::SetLogStringOn(n, 0));

	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(0, ::GetLogStringLineCount(n));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestGetLogStringLine)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);
	ASSERT_EQ(0, ::GetLogStringLineCount(n));

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	ASSERT_EQ(0, ::GetLogStringLineCount(n));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// add knobs
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "KNOBS"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "\t-logfile TRUE"));

	// run
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));

	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(0, ::GetLogStringLineCount(n));

	int line = 0;
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));

	// negative lines should be empty
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, -1)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, -2)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, -3)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, -4)));

	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// add knobs
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "KNOBS"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "\t-logfile TRUE"));

	ASSERT_EQ(IPQ_OK, ::SetLogStringOn(n, 1));

	ASSERT_EQ(0, ::RunAccumulated(n));
#ifndef TESTING
	ASSERT_EQ(29, ::GetLogStringLineCount(n));
#else
	ASSERT_EQ(25, ::GetLogStringLineCount(n));
#endif

	line = 0;
	ASSERT_EQ(std::string("------------------------------------"), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("Reading input data for simulation 1."), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("------------------------------------"), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("-------------------------------------------"), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("Beginning of initial solution calculations."), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("-------------------------------------------"), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("Initial solution 1.	"), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("Iterations in revise_guesses: 2"), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("Number of infeasible solutions: 0"), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("Number of basis changes: 0"), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("Number of iterations: 8"), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("------------------"), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("End of simulation."), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("------------------"), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("------------------------------------"), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("Reading input data for simulation 2."), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string("------------------------------------"), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));
	line++;
	line++;
	line++;
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));

	// add solution block
	// add solution block
	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));

	// add dump block
	ASSERT_EQ(IPQ_OK, ::DUMP(n));

	// add knobs
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "KNOBS"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "\t-logfile TRUE"));

	ASSERT_EQ(IPQ_OK, ::SetLogStringOn(n, 0));


	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(0, ::GetLogStringLineCount(n));

	line = 0;
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, line++)));

	// negative lines should be empty
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, -1)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, -2)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, -3)));
	ASSERT_EQ(std::string(""), std::string(::GetLogStringLine(n, -4)));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestSetErrorFileName)
{
	char ERR_FILENAME[80];
	snprintf(ERR_FILENAME, sizeof(ERR_FILENAME), "error.%s.out", "TestIPhreeqcLib-TestSetErrorFileName");
	if (::FileExists(ERR_FILENAME))
	{
		::DeleteFile(ERR_FILENAME);
	}

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "SOLUTION 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	pH	7"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Na	1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "PHASES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Fix_H+"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	H+ = H+"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	log_k	0"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "EQUILIBRIUM_PHASES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Fix_H+ -10 HCl	10"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "END"));

	// run
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 1));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 1));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileName(n, ERR_FILENAME));


	ASSERT_EQ(1, ::RunAccumulated(n));

	ASSERT_EQ(true, ::FileExists(ERR_FILENAME));

	std::string lines[100];
	{
		std::ifstream ifs(ERR_FILENAME);

		size_t i = 0;
		while (i < sizeof(lines) / sizeof(lines[0]) && std::getline(ifs, lines[i]))
		{
			++i;
		}

		ASSERT_EQ((size_t)90, i);
	}

	ASSERT_EQ(std::string("WARNING: Maximum iterations exceeded, 100"), lines[0]);
	ASSERT_EQ(std::string("WARNING: Numerical method failed with this set of convergence parameters."), lines[2]);
	ASSERT_EQ(std::string("ERROR: Numerical method failed on all combinations of convergence parameters, cell/soln/mix 1"), lines[88]);
	ASSERT_EQ(std::string("Stopping."), lines[89]);

	if (::FileExists(ERR_FILENAME))
	{
		::DeleteFile(ERR_FILENAME);
	}

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestErrorStringOnOff)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(true, ::GetErrorStringOn(n) != 0);

	ASSERT_EQ(IPQ_OK, ::SetErrorStringOn(n, 1));
	ASSERT_EQ(true, ::GetErrorStringOn(n) != 0);

	ASSERT_EQ(IPQ_OK, ::SetErrorStringOn(n, 0));
	ASSERT_EQ(false, ::GetErrorStringOn(n) != 0);

	ASSERT_EQ(IPQ_OK, ::SetErrorStringOn(n, 1));
	ASSERT_EQ(true, ::GetErrorStringOn(n) != 0);

	ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
}

TEST(TestIPhreeqcLib, TestGetErrorString)
{
	char ERR_FILENAME[80];
	snprintf(ERR_FILENAME, sizeof(ERR_FILENAME), "error.%06d.out", ::rand());
	if (::FileExists(ERR_FILENAME))
	{
		::DeleteFile(ERR_FILENAME);
	}
	ASSERT_EQ(false, ::FileExists(ERR_FILENAME));

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "SOLUTION 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	pH	7"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Na	1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "PHASES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Fix_H+"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	H+ = H+"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	log_k	0"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "EQUILIBRIUM_PHASES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Fix_H+ -10 HCl	10"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "END"));

	// run
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 1));
	ASSERT_EQ(IPQ_OK, ::SetErrorStringOn(n, 1));

	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));

	ASSERT_EQ(IPQ_OK, ::SetErrorFileName(n, ERR_FILENAME));
	ASSERT_EQ(std::string(ERR_FILENAME), std::string(::GetErrorFileName(n)));

	ASSERT_EQ(1, ::RunAccumulated(n));

	ASSERT_EQ(std::string(ERR_FILENAME), std::string(::GetErrorFileName(n)));

	ASSERT_EQ(true, ::FileExists(ERR_FILENAME));

	{
#if 0
		std::ifstream ifs(ERR_FILENAME);
		std::string fline((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
#else
		std::string fline("ERROR: Numerical method failed on all combinations of convergence parameters, cell/soln/mix 1\n");
#endif

		std::string sline(::GetErrorString(n));
		ASSERT_TRUE(sline.size() > 0);

		ASSERT_EQ(fline, sline);
	}

	if (::FileExists(ERR_FILENAME))
	{
		::DeleteFile(ERR_FILENAME);
	}

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestGetErrorStringLineCount)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::GetErrorStringLineCount(n));

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	ASSERT_EQ(0, ::GetErrorStringLineCount(n));

	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "SOLUTION 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	pH	7"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Na	1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "PHASES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Fix_H+"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	H+ = H+"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	log_k	0"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "EQUILIBRIUM_PHASES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Fix_H+ -10 HCl	10"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "END"));

	ASSERT_EQ(true, ::GetErrorStringOn(n) != 0);
	ASSERT_EQ(1, ::RunAccumulated(n));
	ASSERT_EQ(1, ::GetErrorStringLineCount(n));

	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "SOLUTION 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	pH	7"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Na	1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "PHASES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Fix_H+"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	H+ = H+"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	log_k	0"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "EQUILIBRIUM_PHASES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Fix_H+ -10 HCl	10"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "END"));

	ASSERT_EQ(IPQ_OK, ::SetErrorStringOn(n, 1));
	ASSERT_EQ(1, ::RunAccumulated(n));
	ASSERT_EQ(1, ::GetErrorStringLineCount(n));

	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "SOLUTION 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	pH	7"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Na	1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "PHASES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Fix_H+"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	H+ = H+"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	log_k	0"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "EQUILIBRIUM_PHASES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "	Fix_H+ -10 HCl	10"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "END"));

	ASSERT_EQ(IPQ_OK, ::SetErrorStringOn(n, 0));

	ASSERT_EQ(1, ::RunAccumulated(n));
	ASSERT_EQ(0, ::GetErrorStringLineCount(n));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestSetSelectedOutputFileName)
{
	char SELOUT_FILENAME[80];
	snprintf(SELOUT_FILENAME, sizeof(SELOUT_FILENAME), "selected_output.%06d.out", ::rand());
	if (::FileExists(SELOUT_FILENAME))
	{
		::DeleteFile(SELOUT_FILENAME);
	}

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "llnl.dat"));

	int max = 6;

	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));
	ASSERT_EQ(IPQ_OK, ::EQUILIBRIUM_PHASES(n, "calcite", 0.0, 0.010));
	ASSERT_EQ(IPQ_OK, ::USER_PUNCH(n, "Ca", max));

	ASSERT_EQ(0, ::GetOutputFileOn(n));
	ASSERT_EQ(0, ::GetErrorFileOn(n));
	ASSERT_EQ(0, ::GetLogFileOn(n));
	ASSERT_EQ(0, ::GetSelectedOutputFileOn(n));
	ASSERT_EQ(0, ::GetDumpFileOn(n));
	ASSERT_EQ(0, ::GetDumpStringOn(n));

	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 1));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileName(n, SELOUT_FILENAME));

	ASSERT_EQ(0, ::RunAccumulated(n));

	ASSERT_EQ(true, ::FileExists(SELOUT_FILENAME));

	/*
	EXPECTED selected.out:
			 sim	       state	        soln	      dist_x	        time	        step	          pH	          pe	           C	          Ca	          Na	     m_CO3-2	     m_CaOH+	    m_NaCO3-	    la_CO3-2	    la_CaOH+	   la_NaCO3-	     Calcite	   d_Calcite	   si_CO2(g)	 si_Siderite	    pressure	   total mol	      volume	    g_CO2(g)	     g_N2(g)	    k_Albite	   dk_Albite	    k_Pyrite	   dk_Pyrite	     s_CaSO4	     s_SrSO4	      1.name	      1.type	     1.moles	      2.name	      2.type	     2.moles	      3.name	      3.type	     3.moles	      4.name	      4.type	     4.moles	      5.name	      5.type	     5.moles	      6.name	      6.type	     6.moles
			   1	      i_soln	           1	         -99	         -99	         -99	           7	           4	 1.0000e-003	 1.0000e-003	 1.0000e-003	 4.2975e-007	 1.1819e-009	 1.1881e-009	-6.4686e+000	-8.9530e+000	-8.9507e+000	 0.0000e+000	 0.0000e+000	     -2.2870	   -999.9990	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	        Ca+2	          aq	 9.9178e-004	     CaHCO3+	          aq	 7.5980e-006	       CaCO3	          aq	 6.2155e-007	       CaOH+	          aq	 1.1819e-009
			   1	       react	           1	         -99	           0	           1	     7.86135	       10.18	 1.1556e-003	 1.1556e-003	 1.0000e-003	 4.2718e-006	 9.7385e-009	 1.1620e-008	-5.4781e+000	-8.0388e+000	-7.9621e+000	 9.8444e-003	-1.5555e-004	     -3.0192	   -999.9990	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	 0.0000e+000	     calcite	        equi	 9.8444e-003	        Ca+2	          aq	 1.1371e-003	     CaHCO3+	          aq	 1.1598e-005	       CaCO3	          aq	 6.8668e-006	       CaOH+	          aq	 9.7385e-009
	*/

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}

	if (::FileExists(SELOUT_FILENAME))
	{
		::DeleteFile(SELOUT_FILENAME);
	}
}

TEST(TestIPhreeqcLib, TestSelectedOutputStringOnOff)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(false, ::GetSelectedOutputFileOn(n) != 0);

	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 1));
	ASSERT_EQ(true, ::GetSelectedOutputFileOn(n) != 0);

	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(false, ::GetSelectedOutputFileOn(n) != 0);

	ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
}

TEST(TestIPhreeqcLib, TestGetSelectedOutputString)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "llnl.dat"));

	int max = 6;

	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));
	ASSERT_EQ(IPQ_OK, ::EQUILIBRIUM_PHASES(n, "calcite", 0.0, 0.010));
	ASSERT_EQ(IPQ_OK, ::USER_PUNCH(n, "Ca", max));

	ASSERT_EQ(0, ::GetOutputFileOn(n));
	ASSERT_EQ(0, ::GetErrorFileOn(n));
	ASSERT_EQ(0, ::GetLogFileOn(n));
	ASSERT_EQ(0, ::GetSelectedOutputFileOn(n));
	ASSERT_EQ(0, ::GetDumpFileOn(n));
	ASSERT_EQ(0, ::GetDumpStringOn(n));

	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputStringOn(n, 1));

	ASSERT_EQ(0, ::RunAccumulated(n));

	std::string sline(::GetSelectedOutputString(n));

	ASSERT_TRUE(::strstr(sline.c_str(), "sim\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "state\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "soln\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "dist_x\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "time\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "step\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "pH\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "pe\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "C\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "Ca\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "Na\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "m_CO3-2\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "m_CaOH+\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "m_NaCO3-\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "la_CO3-2\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "la_CaOH+\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "la_NaCO3-\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "Calcite\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "d_Calcite\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "si_CO2(g)\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "si_Siderite\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "pressure\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "total mol\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "volume\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "g_CO2(g)\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "g_N2(g)\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "k_Albite\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "dk_Albite\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "k_Pyrite\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "dk_Pyrite\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "s_CaSO4\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "s_SrSO4\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "1.name\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "1.type\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "1.moles\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "2.name\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "2.type\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "2.moles\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "3.name\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "3.type\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "3.moles\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "4.name\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "4.type\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "4.moles\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "5.name\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "5.type\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "5.moles\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "6.name\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "6.type\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "6.moles\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "\n") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "i_soln\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "react\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "Ca+2\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "aq\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "CaHCO3+\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "CaCO3\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "CaOH+\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "Calcite\t") != NULL);
	ASSERT_TRUE(::strstr(sline.c_str(), "equi\t") != NULL);

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestGetSelectedOutputStringLineCount)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::GetSelectedOutputFileOn(n));

	ASSERT_EQ(0, ::LoadDatabase(n, "llnl.dat"));

	ASSERT_EQ(0, ::GetSelectedOutputFileOn(n));

	int max = 6;

	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));
	ASSERT_EQ(IPQ_OK, ::EQUILIBRIUM_PHASES(n, "calcite", 0.0, 0.010));
	ASSERT_EQ(IPQ_OK, ::USER_PUNCH(n, "Ca", max));

	ASSERT_EQ(0, ::GetOutputFileOn(n));
	ASSERT_EQ(0, ::GetErrorFileOn(n));
	ASSERT_EQ(0, ::GetLogFileOn(n));
	ASSERT_EQ(0, ::GetSelectedOutputFileOn(n));
	ASSERT_EQ(0, ::GetDumpFileOn(n));
	ASSERT_EQ(0, ::GetDumpStringOn(n));

	ASSERT_EQ(false, ::GetSelectedOutputStringOn(n) != 0);
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputStringOn(n, 1));

	ASSERT_EQ(0, ::RunAccumulated(n));

	ASSERT_EQ(3, ::GetSelectedOutputStringLineCount(n));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestGetSelectedOutputStringLine)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::GetSelectedOutputFileOn(n));

	ASSERT_EQ(0, ::LoadDatabase(n, "llnl.dat"));

	ASSERT_EQ(0, ::GetSelectedOutputFileOn(n));

	int max = 6;

	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));
	ASSERT_EQ(IPQ_OK, ::EQUILIBRIUM_PHASES(n, "calcite", 0.0, 0.010));
	ASSERT_EQ(IPQ_OK, ::USER_PUNCH(n, "Ca", max));

	ASSERT_EQ(0, ::GetOutputFileOn(n));
	ASSERT_EQ(0, ::GetErrorFileOn(n));
	ASSERT_EQ(0, ::GetLogFileOn(n));
	ASSERT_EQ(0, ::GetSelectedOutputFileOn(n));
	ASSERT_EQ(0, ::GetDumpFileOn(n));
	ASSERT_EQ(0, ::GetDumpStringOn(n));

	ASSERT_EQ(false, ::GetSelectedOutputStringOn(n) != 0);

	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(0, ::GetSelectedOutputStringLineCount(n));

	int line = 0;
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, line++)));
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, line++)));

	// negative lines should be empty
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, -1)));
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, -2)));
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, -3)));

	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));
	ASSERT_EQ(IPQ_OK, ::EQUILIBRIUM_PHASES(n, "calcite", 0.0, 0.010));
	ASSERT_EQ(IPQ_OK, ::USER_PUNCH(n, "Ca", max));

	ASSERT_EQ(false, ::GetSelectedOutputStringOn(n) != 0);
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputStringOn(n, 1));

	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(3, ::GetSelectedOutputStringLineCount(n));

	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "sim\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "state\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "soln\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "dist_x\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "time\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "step\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "pH\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "pe\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "C\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "Ca\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "Na\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "m_CO3-2\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "m_CaOH+\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "m_NaCO3-\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "la_CO3-2\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "la_CaOH+\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "la_NaCO3-\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "Calcite\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "d_Calcite\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "si_CO2(g)\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "si_Siderite\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "pressure\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "total mol\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "volume\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "g_CO2(g)\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "g_N2(g)\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "k_Albite\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "dk_Albite\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "k_Pyrite\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "dk_Pyrite\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "s_CaSO4\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "s_SrSO4\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "1.name\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "1.type\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "1.moles\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "2.name\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "2.type\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "2.moles\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "3.name\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "3.type\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "3.moles\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "4.name\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "4.type\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "4.moles\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "5.name\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "5.type\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "5.moles\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "6.name\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "6.type\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "6.moles\t") != NULL);

	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 1), "i_soln\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 1), "Ca+2\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 1), "aq\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 1), "CaHCO3+\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 1), "CaCO3\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 1), "CaOH+\t") != NULL);

	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 2), "react\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 2), "Calcite\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 2), "equi\t") != NULL);

	// after obj.GetSelectedOutputStringLineCount() should be empty
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, ::GetSelectedOutputStringLineCount(n))));
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, ::GetSelectedOutputStringLineCount(n) + 1)));
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, ::GetSelectedOutputStringLineCount(n) + 2)));

	// negative lines should be empty
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, -1)));
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, -2)));
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, -3)));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestGetSelectedOutputStringLineNotEnoughHeadings)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::GetSelectedOutputFileOn(n));

	ASSERT_EQ(0, ::LoadDatabase(n, "llnl.dat"));

	ASSERT_EQ(0, ::GetSelectedOutputFileOn(n));

	ASSERT_EQ(IPQ_OK, ::SOLUTION(n, 1.0, 1.0, 1.0));
	ASSERT_EQ(IPQ_OK, ::EQUILIBRIUM_PHASES(n, "calcite", 0.0, 0.010));
	ASSERT_EQ(IPQ_OK, ::USER_PUNCH_NEH(n));

	ASSERT_EQ(0, ::GetOutputFileOn(n));
	ASSERT_EQ(0, ::GetErrorFileOn(n));
	ASSERT_EQ(0, ::GetLogFileOn(n));
	ASSERT_EQ(0, ::GetSelectedOutputFileOn(n));
	ASSERT_EQ(0, ::GetDumpFileOn(n));
	ASSERT_EQ(0, ::GetDumpStringOn(n));

	ASSERT_EQ(false, ::GetSelectedOutputStringOn(n) != 0);
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputStringOn(n, 1));

	ASSERT_EQ(0, ::RunAccumulated(n));
	ASSERT_EQ(3, ::GetSelectedOutputStringLineCount(n));

	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "sim\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "state\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "soln\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "dist_x\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "time\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "step\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "pH\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "pe\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "C\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "Ca\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "Na\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "m_CO3-2\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "m_CaOH+\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "m_NaCO3-\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "la_CO3-2\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "la_CaOH+\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "la_NaCO3-\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "Calcite\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "d_Calcite\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "si_CO2(g)\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "si_Siderite\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "pressure\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "total mol\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "volume\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "g_CO2(g)\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "g_N2(g)\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "k_Albite\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "dk_Albite\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "k_Pyrite\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "dk_Pyrite\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "s_CaSO4\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "s_SrSO4\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "head0\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "head1\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 0), "head2\t") != NULL);

	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 1), "i_soln\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 1), "have0\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 1), "have1\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 1), "have2\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 1), "missing0\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 1), "missing2\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 1), "missing2\t") != NULL);

	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 2), "react\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 2), "have0\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 2), "have1\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 2), "have2\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 2), "missing0\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 2), "missing2\t") != NULL);
	ASSERT_TRUE(::strstr(::GetSelectedOutputStringLine(n, 2), "missing2\t") != NULL);

	// after obj.GetSelectedOutputStringLineCount() should be empty
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, ::GetSelectedOutputStringLineCount(n))));
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, ::GetSelectedOutputStringLineCount(n) + 1)));
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, ::GetSelectedOutputStringLineCount(n) + 2)));

	// negative lines should be empty
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, -1)));
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, -2)));
	ASSERT_EQ(std::string(""), std::string(::GetSelectedOutputStringLine(n, -3)));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestLongUser_Punch)
{
	const char input[] =
		"PRINT\n"
		" -selected_output t\n"
		"SOLUTION\n"
		"SELECTED_OUTPUT\n"
		" -reset false\n"
		"USER_PUNCH\n"
		"1 REM 255 CHARACTER STRING\n"
		"10 temp$ = \"XXXXX123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345\"\n"
		"20 PUNCH temp$\n";

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputStringOn(n, 1));
	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));
	ASSERT_EQ(0, ::RunString(n, input));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestBasicSURF)
{
	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));

	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "SURFACE_MASTER_SPECIES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "		 Surfa Surfa"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "		 Surfb Surfb"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "SURFACE_SPECIES"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "		 Surfa = Surfa"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "		 log_k 0"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "		 Surfb = Surfb"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "		 log_k 0"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "		 Surfa + Zn+2 = SurfaZn+2"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "		 log_k  5."));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "		 Surfb + Zn+2 = SurfbZn+2"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "		 log_k  6."));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "		 Surfa + Cu+2 = SurfaCu+2"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "		 log_k  4.5"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "		 Surfb + Cu+2 = SurfbCu+2"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "		 log_k  6.5"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "SOLUTION 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "   pH        8"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "   units     mol/kgw"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "   Fe(3)     1e-2"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "   Zn        1e-4"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "   Cu        1e-5"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "   Na        1e-1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "   Cl        1e-1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "EQUILIBRIUM_PHASES 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "   Fe(OH)3(a) 0 0"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "SELECTED_OUTPUT"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "USER_PUNCH"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "    -headings Hfo-Zn Surfa-Zn Surfb-Zn Surfa-Cu Surfb-Cu"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "-start"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "10 PUNCH SURF(\"Zn\",\"Hfo\")"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "20 PUNCH SURF(\"Zn\",\"Surfa\")"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "30 PUNCH SURF(\"Zn\",\"Surfb\")"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "40 PUNCH SURF(\"Cu\",\"Surfa\")"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "50 PUNCH SURF(\"Cu\",\"Surfb\")"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "-end"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "SURFACE 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "    Hfo_sOH Fe(OH)3(a)      equilibrium_phase 0.005  53300"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "    Hfo_wOH Fe(OH)3(a)      equilibrium_phase 0.2"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "    Surfa  0.2 100. 2"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "    Surfb  0.1 100. 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(n, "END"));


	ASSERT_EQ(IPQ_OK, ::SetOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetOutputStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetErrorFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetLogFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpStringOn(n, 0));
	ASSERT_EQ(IPQ_OK, ::SetDumpFileOn(n, 0));


	ASSERT_EQ(0, ::RunAccumulated(n));

	ASSERT_EQ(13, ::GetSelectedOutputColumnCount(n));
	ASSERT_EQ(3, ::GetSelectedOutputRowCount(n));

	VAR v;
	VarInit(&v);

	const int offset = 8;

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, offset + 0, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("Hfo-Zn"), std::string(v.sVal));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, offset + 1, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("Surfa-Zn"), std::string(v.sVal));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, offset + 2, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("Surfb-Zn"), std::string(v.sVal));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, offset + 3, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("Surfa-Cu"), std::string(v.sVal));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 0, offset + 4, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("Surfb-Cu"), std::string(v.sVal));


	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, offset + 0, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(0.0, v.dVal, ::pow(10., -FLT_DIG));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, offset + 1, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(0.0, v.dVal, ::pow(10., -FLT_DIG));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, offset + 2, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(0.0, v.dVal, ::pow(10., -FLT_DIG));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, offset + 3, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(0.0, v.dVal, ::pow(10., -FLT_DIG));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 1, offset + 4, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(0.0, v.dVal, ::pow(10., -FLT_DIG));


	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, offset + 0, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(6.3861e-005, v.dVal, ::pow(10., -FLT_DIG));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, offset + 1, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(1.7868e-005, v.dVal, ::pow(10., -FLT_DIG));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, offset + 2, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(1.8248e-005, v.dVal, ::pow(10., -FLT_DIG));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, offset + 3, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(4.6216e-009, v.dVal, ::pow(10., -FLT_DIG));

	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(n, 2, offset + 4, &v));
	ASSERT_EQ(TT_DOUBLE, v.type);
	ASSERT_NEAR(4.7201e-008, v.dVal, ::pow(10., -FLT_DIG));

	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
}

TEST(TestIPhreeqcLib, TestIEEE)
{
	// COMMENT: {1/18/2013 6:34:57 PM}	int n = ::CreateIPhreeqc();
	// COMMENT: {1/18/2013 6:34:57 PM}	ASSERT_TRUE(n >= 0);
	// COMMENT: {1/18/2013 6:34:57 PM}
	// COMMENT: {1/18/2013 6:34:57 PM}	ASSERT_EQ( 0,      ::LoadDatabase(n, "wateq4f.dat") );
	// COMMENT: {1/18/2013 6:34:57 PM}	//ASSERT_EQ( 0,      ::LoadDatabase(n, "cdmusic_hiemstra.dat") );
	// COMMENT: {1/18/2013 6:34:57 PM}	//::LoadDatabase(n, "cdmusic_hiemstra.dat");
	// COMMENT: {1/18/2013 6:34:57 PM}	//const char* errdb = ::GetErrorString(n);
	// COMMENT: {1/18/2013 6:34:57 PM}	//fprintf(stderr, "%s\n", errdb);
	// COMMENT: {1/18/2013 6:34:57 PM}
	// COMMENT: {1/18/2013 6:34:57 PM}	ASSERT_EQ( IPQ_OK, ::SetOutputFileOn(n, 0) );
	// COMMENT: {1/18/2013 6:34:57 PM}	ASSERT_EQ( IPQ_OK, ::SetOutputStringOn(n, 0) );
	// COMMENT: {1/18/2013 6:34:57 PM}	ASSERT_EQ( IPQ_OK, ::SetErrorFileOn(n, 0) );
	// COMMENT: {1/18/2013 6:34:57 PM}	ASSERT_EQ( IPQ_OK, ::SetLogFileOn(n, 0) );
	// COMMENT: {1/18/2013 6:34:57 PM}	ASSERT_EQ( IPQ_OK, ::SetSelectedOutputFileOn(n, 0) );
	// COMMENT: {1/18/2013 6:34:57 PM}	ASSERT_EQ( IPQ_OK, ::SetDumpStringOn(n, 0) );
	// COMMENT: {1/18/2013 6:34:57 PM}	ASSERT_EQ( IPQ_OK, ::SetDumpFileOn(n, 0) );
	// COMMENT: {1/18/2013 6:34:57 PM}
	// COMMENT: {1/18/2013 6:34:57 PM}	::RunFile(n, "IEEE");
	// COMMENT: {1/18/2013 6:34:57 PM}	const char* err = ::GetErrorString(n);
	// COMMENT: {1/18/2013 6:34:57 PM}	fprintf(stderr, "%s\n", err);
	// COMMENT: {1/18/2013 6:34:57 PM}
	// COMMENT: {1/18/2013 6:34:57 PM}	if (n >= 0)
	// COMMENT: {1/18/2013 6:34:57 PM}	{
	// COMMENT: {1/18/2013 6:34:57 PM}		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	// COMMENT: {1/18/2013 6:34:57 PM}	}
}

TEST(TestIPhreeqcLib, TestDelete)
{
	const char input[] =
		"SOLUTION 1 # definition of initial condition 1\n"
		"COPY cell 1 7405 # copy cell 1 to placeholder cell with index larger than the number of cells in the model domain\n"
		"END\n"
		"DELETE # delete initial condition 1 to allow for a redefinition of all reactions\n"
		"-cell 1\n"
		"END\n"
		"# define other initial conditions and copy to another placeholder cell\n"
		"\n"
		"COPY cell 7405    1 # copy back from placeholder cell to domain cell 1\n"
		"END\n"
		"MIX    1 # mix according to initial moisture content\n"
		"   1 0.25\n"
		"END\n"
		"RUN_CELLS\n"
		"-cells 1\n"
		"-start_time 0\n"
		"-time_step 0\n"
		"DELETE # remove mix reaction in subsequent runs\n"
		"-mix 1\n"
		"END\n"
		"RUN_CELLS\n"
		"-cells 1\n";

	int n = ::CreateIPhreeqc();
	ASSERT_TRUE(n >= 0);

	char OUTPUT_FILE[80];
	snprintf(OUTPUT_FILE, sizeof(OUTPUT_FILE), "phreeqc.%d.out", n);

	if (::FileExists(OUTPUT_FILE))
	{
		ASSERT_TRUE(::DeleteFile(OUTPUT_FILE));
	}
	ASSERT_EQ(false, ::FileExists(OUTPUT_FILE));

	ASSERT_EQ(0, ::LoadDatabase(n, "phreeqc.dat"));
	::SetOutputFileOn(n, 0);
	::SetErrorFileOn(n, 0);
	::SetLogFileOn(n, 0);
	::SetSelectedOutputFileOn(n, 0);
	::SetDumpFileOn(n, 0);
	ASSERT_EQ(false, ::FileExists(OUTPUT_FILE));
	ASSERT_EQ(0, ::RunString(n, input));
	ASSERT_EQ(false, ::FileExists(OUTPUT_FILE));
	if (n >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(n));
	}
	if (::FileExists(OUTPUT_FILE))
	{
		ASSERT_TRUE(::DeleteFile(OUTPUT_FILE));
	}
}

TEST(TestIPhreeqcLib, TestMultiPunchCSelectedOutput)
{
	VAR var;
	::VarInit(&var);

	int id = ::CreateIPhreeqc();
	ASSERT_TRUE(id >= 0);

	ASSERT_EQ(0, ::LoadDatabase(id, "phreeqc.dat.90a6449"));
	ASSERT_EQ(0, ::RunFile(id, "multi_punch"));

	ASSERT_EQ(6, ::GetSelectedOutputRowCount(id));
	ASSERT_EQ(35, ::GetSelectedOutputColumnCount(id));

	// headings
	int ncol = 0;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("sim"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("state"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("soln"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("dist_x"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("time"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("step"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("pH"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("pe"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("reaction"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("temp(C)"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Alk(eq/kgw)"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("mu"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("mass_H2O"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("charge(eq)"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("pct_err"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Na(mol/kgw)"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Ca(mol/kgw)"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("m_Na+(mol/kgw)"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("m_HCO3-(mol/kgw)"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("la_Ca+2"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("la_CO3-2"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("CO2(g)"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("d_CO2(g)"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("dolomite"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("d_dolomite"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("si_Halite"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("pressure"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("total mol"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("volume"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("g_N2(g)"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("k_Calcite"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("dk_Calcite"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("s_Anhydrite"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("s_Barite"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("V_TOTAL_C"), std::string(var.sVal));

	// sim
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 0, &var));   ASSERT_EQ((long)8, var.lVal);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 0, &var));   ASSERT_EQ((long)10, var.lVal);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 0, &var));   ASSERT_EQ((long)11, var.lVal);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 0, &var));   ASSERT_EQ((long)12, var.lVal);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 0, &var));   ASSERT_EQ((long)14, var.lVal);

	// state
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 1, &var));   ASSERT_EQ(std::string("react"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 1, &var));   ASSERT_EQ(std::string("react"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 1, &var));   ASSERT_EQ(std::string("i_soln"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 1, &var));   ASSERT_EQ(std::string("i_soln"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 1, &var));   ASSERT_EQ(std::string("react"), std::string(var.sVal));

	// pH
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 6, &var));   ASSERT_NEAR(7.30475, var.dVal, ::pow(10., -2));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 6, &var));   ASSERT_NEAR(7.29765, var.dVal, ::pow(10., -2));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 6, &var));   ASSERT_NEAR(6.99738, var.dVal, ::pow(10., -2));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 6, &var));   ASSERT_NEAR(6.99698, var.dVal, ::pow(10., -2));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 6, &var));   ASSERT_NEAR(7.2942, var.dVal, ::pow(10., -2));

	// V_TOTAL_C
#ifdef SKIP_TEST
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 34, &var));   ASSERT_NEAR(4.3729e-003, var.dVal, ::pow(10., -6));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 34, &var));   ASSERT_NEAR(4.3090e-003, var.dVal, ::pow(10., -6));
#endif
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 34, &var));   ASSERT_NEAR(0.0000e+000, var.dVal, ::pow(10., -6));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 34, &var));   ASSERT_NEAR(0.0000e+000, var.dVal, ::pow(10., -6));
#ifdef SKIP_TEST
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 34, &var));   ASSERT_NEAR(4.2784e-003, var.dVal, ::pow(10., -6));
#endif

	// edge cases
	int r = ::GetSelectedOutputRowCount(id);
	int c = ::GetSelectedOutputColumnCount(id);
	ASSERT_EQ(IPQ_INVALIDROW, ::GetSelectedOutputValue(id, -1, 0, &var));
	ASSERT_EQ(TT_ERROR, var.type);
	ASSERT_EQ(VR_INVALIDROW, var.vresult);

	ASSERT_EQ(IPQ_INVALIDROW, ::GetSelectedOutputValue(id, r, 0, &var));  ASSERT_EQ(TT_ERROR, var.type);  ASSERT_EQ(VR_INVALIDROW, var.vresult);
	ASSERT_EQ(IPQ_INVALIDCOL, ::GetSelectedOutputValue(id, 0, -1, &var));  ASSERT_EQ(TT_ERROR, var.type);  ASSERT_EQ(VR_INVALIDCOL, var.vresult);
	ASSERT_EQ(IPQ_INVALIDCOL, ::GetSelectedOutputValue(id, 0, c, &var));  ASSERT_EQ(TT_ERROR, var.type);  ASSERT_EQ(VR_INVALIDCOL, var.vresult);

	::SetCurrentSelectedOutputUserNumber(id, 2);
	ASSERT_EQ(7, ::GetSelectedOutputRowCount(id));
	ASSERT_EQ(16, ::GetSelectedOutputColumnCount(id));

	// headings
	ncol = 0;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("si_Halite"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("si_Calcite"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("DUMMY_1"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("DUMMY_2"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Sum_resid"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Sum_Delta/U"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("MaxFracErr"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Soln_2"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Soln_2_min"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Soln_2_max"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Soln_3"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Soln_3_min"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Soln_3_max"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Halite"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Halite_min"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Halite_max"), std::string(var.sVal));

	// si_Halite
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 0, &var));   ASSERT_NEAR(-7.70857, var.dVal, ::pow(10., -2));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 0, &var));   ASSERT_NEAR(-7.67087, var.dVal, ::pow(10., -2));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 0, &var));   ASSERT_NEAR(-7.6362, var.dVal, ::pow(10., -2));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 0, &var));   ASSERT_NEAR(-999.999, var.dVal, ::pow(10., -2));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 0, &var));   ASSERT_NEAR(-7.60092, var.dVal, ::pow(10., -2));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 0, &var));   ASSERT_NEAR(-7.60411, var.dVal, ::pow(10., -2));

	// si_Calcite
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 1, &var));   ASSERT_NEAR(0.702316, var.dVal, ::pow(10., -2));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 1, &var));   ASSERT_NEAR(0.695856, var.dVal, ::pow(10., -2));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 1, &var));   ASSERT_NEAR(0.689518, var.dVal, ::pow(10., -2));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 1, &var));   ASSERT_NEAR(-999.999, var.dVal, ::pow(10., -2));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 1, &var));   ASSERT_NEAR(-999.999, var.dVal, ::pow(10., -2));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 1, &var));   ASSERT_NEAR(0.683300, var.dVal, ::pow(10., -2));

	// DUMMY_1
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 2, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 2, &var));   ASSERT_EQ(std::string("Dummy1"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 2, &var));   ASSERT_EQ(std::string("Dummy1"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 2, &var));   ASSERT_EQ(std::string("Dummy1"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 2, &var));   ASSERT_EQ(std::string("Dummy1"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 2, &var));   ASSERT_EQ(std::string("Dummy1"), std::string(var.sVal));

	// DUMMY_2
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 3, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 3, &var));   ASSERT_EQ(std::string("Dummy2"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 3, &var));   ASSERT_EQ(std::string("Dummy2"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 3, &var));   ASSERT_EQ(std::string("Dummy2"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 3, &var));   ASSERT_EQ(std::string("Dummy2"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 3, &var));   ASSERT_EQ(std::string("Dummy2"), std::string(var.sVal));

	// Sum_resid
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 4, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 4, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 4, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 4, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 4, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 4, &var));   ASSERT_NEAR(4.12e-13, var.dVal, ::pow(10., log10(4.12e-13) - 2));

	// Sum_Delta/U
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 5, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 5, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 5, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 5, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 5, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 5, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));

	// MaxFracErr
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 6, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 6, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 6, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 6, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 6, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 6, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));

	// Soln_2
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 7, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 7, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 7, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 7, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 7, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 7, &var));   ASSERT_NEAR(1, var.dVal, ::pow(10., -3));

	// Soln_2_min
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 8, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 8, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 8, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 8, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 8, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 8, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));

	// Soln_2_max
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 9, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 9, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 9, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 9, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 9, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 9, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));

	// Soln_3
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 10, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 10, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 10, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 10, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 10, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 10, &var));   ASSERT_NEAR(1, var.dVal, ::pow(10., -3));

	// Soln_3_min
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 11, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 11, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 11, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 11, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 11, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 11, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));

	// Soln_3_max
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 12, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 12, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 12, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 12, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 12, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 12, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));

	// Halite
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 13, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 13, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 13, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 13, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 13, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 13, &var));   ASSERT_NEAR(0.001, var.dVal, ::pow(10., -3));

	// Halite_min
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 14, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 14, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 14, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 14, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 14, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 14, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));

	// Halite_max
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 15, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 2, 15, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 3, 15, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 4, 15, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 5, 15, &var));   ASSERT_EQ(TT_EMPTY, var.type);
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 6, 15, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));

	// edge cases
	r = ::GetSelectedOutputRowCount(id);
	c = ::GetSelectedOutputColumnCount(id);
	ASSERT_EQ(IPQ_INVALIDROW, ::GetSelectedOutputValue(id, -1, 0, &var));  ASSERT_EQ(TT_ERROR, var.type);  ASSERT_EQ(VR_INVALIDROW, var.vresult);
	ASSERT_EQ(IPQ_INVALIDROW, ::GetSelectedOutputValue(id, r, 0, &var));  ASSERT_EQ(TT_ERROR, var.type);  ASSERT_EQ(VR_INVALIDROW, var.vresult);
	ASSERT_EQ(IPQ_INVALIDCOL, ::GetSelectedOutputValue(id, 0, -1, &var));  ASSERT_EQ(TT_ERROR, var.type);  ASSERT_EQ(VR_INVALIDCOL, var.vresult);
	ASSERT_EQ(IPQ_INVALIDCOL, ::GetSelectedOutputValue(id, 0, c, &var));  ASSERT_EQ(TT_ERROR, var.type);  ASSERT_EQ(VR_INVALIDCOL, var.vresult);

	::SetCurrentSelectedOutputUserNumber(id, 3);
	ASSERT_EQ(2, ::GetSelectedOutputRowCount(id));
	ASSERT_EQ(12, ::GetSelectedOutputColumnCount(id));

	// headings
	ncol = 0;
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Sum_resid"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Sum_Delta/U"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("MaxFracErr"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Soln_2"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Soln_2_min"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Soln_2_max"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Soln_3"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Soln_3_min"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Soln_3_max"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Halite"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Halite_min"), std::string(var.sVal));
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 0, ncol++, &var));   ASSERT_EQ(std::string("Halite_max"), std::string(var.sVal));

	// Sum_resid
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 0, &var));   ASSERT_NEAR(4.12E-13, var.dVal, ::pow(10., log10(4.12E-13) - 2));

	// Sum_Delta/U
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 1, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));

	// MaxFracErr
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 2, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));

	// Soln_2
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 3, &var));   ASSERT_NEAR(1, var.dVal, ::pow(10., -3));

	// Soln_2_min
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 4, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));

	// Soln_2_max
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 5, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));

	// Soln_3
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 6, &var));   ASSERT_NEAR(1, var.dVal, ::pow(10., -3));

	// Soln_3_min
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 7, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));

	// Soln_3_max
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 8, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));

	// Halite
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 9, &var));   ASSERT_NEAR(0.001, var.dVal, ::pow(10., -3));

	// Halite_min
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 10, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));

	// Halite_max
	ASSERT_EQ(IPQ_OK, ::GetSelectedOutputValue(id, 1, 11, &var));   ASSERT_NEAR(0, var.dVal, ::pow(10., -3));


	// edge cases
	r = ::GetSelectedOutputRowCount(id);
	c = ::GetSelectedOutputColumnCount(id);
	ASSERT_EQ(IPQ_INVALIDROW, ::GetSelectedOutputValue(id, -1, 0, &var));  ASSERT_EQ(TT_ERROR, var.type);  ASSERT_EQ(VR_INVALIDROW, var.vresult);
	ASSERT_EQ(IPQ_INVALIDROW, ::GetSelectedOutputValue(id, r, 0, &var));  ASSERT_EQ(TT_ERROR, var.type);  ASSERT_EQ(VR_INVALIDROW, var.vresult);
	ASSERT_EQ(IPQ_INVALIDCOL, ::GetSelectedOutputValue(id, 0, -1, &var));  ASSERT_EQ(TT_ERROR, var.type);  ASSERT_EQ(VR_INVALIDCOL, var.vresult);
	ASSERT_EQ(IPQ_INVALIDCOL, ::GetSelectedOutputValue(id, 0, c, &var));  ASSERT_EQ(TT_ERROR, var.type);  ASSERT_EQ(VR_INVALIDCOL, var.vresult);

	ASSERT_EQ(IPQ_INVALIDARG, ::SetCurrentSelectedOutputUserNumber(id, -1));
	ASSERT_EQ(IPQ_OK, ::SetCurrentSelectedOutputUserNumber(id, 0));
	ASSERT_EQ(0, ::GetCurrentSelectedOutputUserNumber(id));

	if (id >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(id));
	}
}

TEST(TestIPhreeqcLib, TestGetSelectedOutputCount)
{
	int id = ::CreateIPhreeqc();
	ASSERT_TRUE(id >= 0);

	ASSERT_EQ(0, ::GetSelectedOutputCount(id));
	ASSERT_EQ(0, ::LoadDatabase(id, "phreeqc.dat"));
	ASSERT_EQ(0, ::GetSelectedOutputCount(id));
	ASSERT_EQ(0, ::RunFile(id, "multi_punch"));
	ASSERT_EQ(3, ::GetSelectedOutputCount(id));

	if (id >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(id));
	}
}

TEST(TestIPhreeqcLib, TestGetNthSelectedOutputUserNumber)
{
	int id = ::CreateIPhreeqc();
	ASSERT_TRUE(id >= 0);

	ASSERT_EQ(0, ::LoadDatabase(id, "phreeqc.dat"));
	ASSERT_EQ(0, ::RunFile(id, "multi_punch"));

	ASSERT_EQ(3, ::GetSelectedOutputCount(id));

	ASSERT_EQ(1, ::GetNthSelectedOutputUserNumber(id, 0));
	ASSERT_EQ(2, ::GetNthSelectedOutputUserNumber(id, 1));
	ASSERT_EQ(3, ::GetNthSelectedOutputUserNumber(id, 2));

	// edge cases
	ASSERT_EQ((int)IPQ_INVALIDARG, ::GetNthSelectedOutputUserNumber(id, -1));
	ASSERT_EQ((int)IPQ_INVALIDARG, ::GetNthSelectedOutputUserNumber(id, 4));

	if (id >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(id));
	}
}

TEST(TestIPhreeqcLib, TestGetCurrentSelectedOutputUserNumber)
{
	int id = ::CreateIPhreeqc();
	ASSERT_TRUE(id >= 0);
	ASSERT_EQ(1, ::GetCurrentSelectedOutputUserNumber(id));

	ASSERT_EQ(0, ::LoadDatabase(id, "phreeqc.dat"));
	ASSERT_EQ(1, ::GetCurrentSelectedOutputUserNumber(id));
	ASSERT_EQ(0, ::RunFile(id, "multi_punch"));

	ASSERT_EQ(1, ::GetCurrentSelectedOutputUserNumber(id));
	ASSERT_EQ(6, ::GetSelectedOutputRowCount(id));
	ASSERT_EQ(35, ::GetSelectedOutputColumnCount(id));

	::SetCurrentSelectedOutputUserNumber(id, 2);
	ASSERT_EQ(2, ::GetCurrentSelectedOutputUserNumber(id));
	ASSERT_EQ(7, ::GetSelectedOutputRowCount(id));
	ASSERT_EQ(16, ::GetSelectedOutputColumnCount(id));

	::SetCurrentSelectedOutputUserNumber(id, 3);
	ASSERT_EQ(3, ::GetCurrentSelectedOutputUserNumber(id));
	ASSERT_EQ(2, ::GetSelectedOutputRowCount(id));
	ASSERT_EQ(12, ::GetSelectedOutputColumnCount(id));

	// edge cases
	ASSERT_EQ(IPQ_INVALIDARG, ::SetCurrentSelectedOutputUserNumber(id, -1));
	ASSERT_EQ(3, ::GetCurrentSelectedOutputUserNumber(id));
	ASSERT_EQ(IPQ_OK, ::SetCurrentSelectedOutputUserNumber(id, 0));
	ASSERT_EQ(0, ::GetCurrentSelectedOutputUserNumber(id));

	// unload database
	ASSERT_EQ(0, ::LoadDatabase(id, "phreeqc.dat"));
	ASSERT_EQ(1, ::GetCurrentSelectedOutputUserNumber(id));

	if (id >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(id));
	}
}

TEST(TestIPhreeqcLib, TestMultiSetSelectedOutputFileName)
{
	FileTest set1("state.sel");
	ASSERT_TRUE(set1.RemoveExisting());

	FileTest set2("si.sel");
	ASSERT_TRUE(set2.RemoveExisting());

	int id = ::CreateIPhreeqc();
	ASSERT_TRUE(id >= 0);

	ASSERT_EQ(0, ::LoadDatabase(id, "phreeqc.dat"));

	ASSERT_EQ(IPQ_OK, ::SetCurrentSelectedOutputUserNumber(id, 1));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputStringOn(id, true));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(id, true));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileName(id, set1.GetName().c_str()));

	ASSERT_EQ(IPQ_OK, ::SetCurrentSelectedOutputUserNumber(id, 2));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputStringOn(id, true));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileOn(id, true));
	ASSERT_EQ(IPQ_OK, ::SetSelectedOutputFileName(id, set2.GetName().c_str()));

	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "TITLE Temperature dependence of solubility"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "      of gypsum and anhydrite             "));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "SOLUTION 1 Pure water                     "));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "        pH      7.0                       "));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "        temp    25.0                      "));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "EQUILIBRIUM_PHASES 1                      "));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "        Gypsum          0.0     1.0       "));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "        Anhydrite       0.0     1.0       "));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "REACTION_TEMPERATURE 1                    "));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "        25.0 75.0 in 51 steps             "));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "SELECTED_OUTPUT 1                         "));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "        -temperature                      "));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "SELECTED_OUTPUT 2                         "));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "        -si     anhydrite  gypsum         "));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "END                                       "));

	ASSERT_EQ(0, ::RunAccumulated(id));

	ASSERT_TRUE(set1.VerifyExists());
	ASSERT_TRUE(set2.VerifyExists());

	if (id >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(id));
	}
}

TEST(TestIPhreeqcLib, TestWissmeier20131203)
{
	int id = ::CreateIPhreeqc();
	ASSERT_TRUE(id >= 0);

	ASSERT_EQ(0, ::LoadDatabase(id, "phreeqc.dat"));

	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "selected_output 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "-totals O"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, ""));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "END"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "selected_output 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "-totals H"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, ""));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "solution"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "END"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "selected_output"));

	// original asserts here
	ASSERT_EQ(0, ::RunAccumulated(id));

	ASSERT_EQ(1, ::GetSelectedOutputCount(id));

	ASSERT_EQ(9, ::GetSelectedOutputColumnCount(id));
	ASSERT_EQ(2, ::GetSelectedOutputRowCount(id));

	if (id >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(id));
	}
}

TEST(TestIPhreeqcLib, TestWissmeier20131203_2)
{
	int id = ::CreateIPhreeqc();
	ASSERT_TRUE(id >= 0);

	ASSERT_EQ(0, ::LoadDatabase(id, "phreeqc.dat"));

	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "selected_output 22"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "-totals O"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, ""));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "END"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "selected_output 22"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "-totals H"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, ""));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "solution"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "END"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "selected_output 22"));

	// original asserts here
	ASSERT_EQ(0, ::RunAccumulated(id));

	ASSERT_EQ(1, ::GetSelectedOutputCount(id));

	ASSERT_EQ(IPQ_OK, ::SetCurrentSelectedOutputUserNumber(id, 22));
	ASSERT_EQ(1, ::GetSelectedOutputColumnCount(id));
	ASSERT_EQ(2, ::GetSelectedOutputRowCount(id));

	if (id >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(id));
	}
}

TEST(TestIPhreeqcLib, TestWissmeier20131203_3)
{
	int id = ::CreateIPhreeqc();
	ASSERT_TRUE(id >= 0);

	ASSERT_EQ(0, ::LoadDatabase(id, "phreeqc.dat"));

	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "selected_output 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "-reset false"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "-totals O"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, ""));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "END"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "selected_output 1"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "-reset false"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "-totals H"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, ""));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "solution"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "END"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "selected_output"));
	ASSERT_EQ(IPQ_OK, ::AccumulateLine(id, "-reset false"));

	// original asserts here
	ASSERT_EQ(0, ::RunAccumulated(id));

	ASSERT_EQ(1, ::GetSelectedOutputCount(id));

	ASSERT_EQ(1, ::GetSelectedOutputColumnCount(id));
	ASSERT_EQ(2, ::GetSelectedOutputRowCount(id));

	if (id >= 0)
	{
		ASSERT_EQ(IPQ_OK, ::DestroyIPhreeqc(id));
	}
}

TEST(TestIPhreeqcLib, TestIsZeroInitialized)
{
	std::map<void*, void*> test_map;
	ASSERT_TRUE(test_map.find((void*)0) == test_map.end());
	ASSERT_TRUE(test_map[(void*)0] == (void*)0);

	int i = int();
	long l = long();
	float f = float();
	double d = double();

	// these would fail
	/**
	int i;
	long l;
	float f;
	double d;
	**/

	ASSERT_EQ((int)0, i);
	ASSERT_EQ((long)0, l);
	ASSERT_EQ((float)0, f);
	ASSERT_EQ((double)0, d);
}
