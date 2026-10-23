#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <IPhreeqc.hpp>
#include <sstream>
#include <iomanip>

class MyData
{
public:
	void Mydata() 
	{
		this->IPhreeqc_ptr = NULL;
		this->year         = 0.0;
	}
	~MyData() {};

	// data members
	IPhreeqc *          IPhreeqc_ptr;
	std::vector < VAR > results;
	double              year;

	// methods
	void ExtractWrite(int cell)
	{
		results.clear();
		results.resize(IPhreeqc_ptr->GetSelectedOutputColumnCount());
		for (int j = 0; j < IPhreeqc_ptr->GetSelectedOutputColumnCount(); ++j) 
		{
			VarInit(&results[j]);
			IPhreeqc_ptr->GetSelectedOutputValue(1, j, &results[j]);
		}
		std::cerr << "Cell "     << cell << "  " << (int) results[7].dVal 
			<< std::setprecision(2) << std::fixed
			<< "\n\tpH: "        << results[5].dVal 
			<< "\tSR(calcite): " << results[6].dVal << "\n";
	}

	void EHandler(void)
	{
		IPhreeqc_ptr->OutputErrorString();
		exit(EXIT_FAILURE);	
	}

	static double MyCallback(double x1, double x2, const char * str1, void *my_ptr)
	{
		/*
		Use of a callback is optional.

		The callback provides a way to obtain data from a Basic program
		through the variables x1, x2, and str1, and send data to a 
		Basic program through the return value of the callback.

		The void pointer mydata can be used to obtain data from the
		calling program; in this example, it points to a structure.

		The callback function is called whenever CALLBACK(x1, x2, str$)  
		is used in a Basic program (usually USER_PUNCH). See file "ic".
		*/
		if (strcmp(str1, "Year") == 0)
		{
			fprintf(stderr, "\nCallback for cell %d: pH %8.2f\n", (int) x1, x2);
			return ((MyData *) my_ptr)->year;
		}
		return -1;
	}

	void Exec(void)
	{
		// Create module
		this->IPhreeqc_ptr = new IPhreeqc;
		// Load database
		if (this->IPhreeqc_ptr->LoadDatabase("phreeqc.dat") != 0) this->EHandler();
		// Set callback
		this->IPhreeqc_ptr->SetBasicCallback(this->MyCallback, (void *) this);
		// Define initial conditions and selected output 
		this->year = 2014.0;
		if (this->IPhreeqc_ptr->RunFile("ic") != 0) this->EHandler();
		// Run cell 1
		if (this->IPhreeqc_ptr->RunString("RUN_CELLS; -cells; 1; END") != 0) this->EHandler();
		// Extract/write results 
		this->ExtractWrite(1);

		// Advect cell 1 solution to cell 2
		this->year += 1.0;
		// Define new solution composition for cell 2 from cell 1 selected-output
		std::ostringstream oss;
		oss << "SOLUTION_MODIFY 2" << "\n";
		oss << "   -cb      "      << this->results[0].dVal << "\n"; 
		oss << "   -total_h "      << this->results[1].dVal << "\n"; 
		oss << "   -total_o "      << this->results[2].dVal << "\n"; 
		oss << "   -totals  "      << "\n"; 
		oss << "      C     "      << this->results[3].dVal << "\n"; 
		oss << "      Ca    "      << this->results[4].dVal << "\n"; 
		// run cell 2
		oss << "RUN_CELLS; -cells; 2; END\n";
		if (this->IPhreeqc_ptr->RunString(oss.str().c_str()) != 0) this->EHandler();
		// Extract/write results
		this->ExtractWrite(2);

		// Destroy module 
		delete this->IPhreeqc_ptr;
	}
};

int main(void)
{
	MyData mydata;
	mydata.Exec();
}