#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <IPhreeqc.h>

int id;
int vt[8];
double dv[8];
char sv[8][100];
char buffer[100];
struct MyData
{
	double year;
};
void ExtractWrite(int cell)
{
	VAR v;
	int j;
	VarInit(&v);
	for (j = 0; j < 8; ++j) {
		GetSelectedOutputValue(id, 1, j, &v);
		vt[j] = v.type;
		switch (vt[j]) {
		case TT_DOUBLE:
			dv[j] = v.dVal;
			snprintf(sv[j], sizeof(sv[j]), "%23.15e", v.dVal);
			break;
		case TT_STRING:
			strcpy(sv[j], v.sVal);
			break;
		}
		VarClear(&v);
	}
	printf("Cell %d %d\n\tpH: %4.2f\tSR(calcite): %4.2f\n", cell, (int) dv[7], dv[5], dv[6]);
}

void EHandler(void)
{
	OutputErrorString(id);
	exit(EXIT_FAILURE);	
}

const char *ConCat(const char *str1, const char *str2)
{
	strcpy(buffer, str1);
	return strcat(buffer, str2);
}

double MyCallback(double x1, double x2, const char * str1, void *mydata)
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
		return ((struct MyData *) mydata)->year;
	}
	return -1;
}
int main(void)
{
	struct MyData mydata;
	mydata.year = 2012.0;
	/* Create module, load database, define initial conditions and selected output */
	id = CreateIPhreeqc();
	if (LoadDatabase(id, "phreeqc.dat") != 0) EHandler();
	if (SetBasicCallback(id, MyCallback, &mydata) != 0) EHandler();
	if (RunFile(id, "ic") != 0) EHandler();

	/* Run cell 1, extract/write result */
	if (RunString(id, "RUN_CELLS; -cells; 1; END") != 0) EHandler();
	ExtractWrite(1);

	/* Advect cell 1 solution to cell 2, run cell 2, extract/write results */
	AccumulateLine(id, ConCat("SOLUTION_MODIFY 2",         ""   ));
	AccumulateLine(id, ConCat("   -cb      ",              sv[0]));
	AccumulateLine(id, ConCat("   -total_h ",              sv[1]));
	AccumulateLine(id, ConCat("   -total_o ",              sv[2]));
	AccumulateLine(id, ConCat("   -totals  ",              ""   ));
	AccumulateLine(id, ConCat("      C     ",              sv[3]));
	AccumulateLine(id, ConCat("      Ca    ",              sv[4]));
	mydata.year += 1.0;
	AccumulateLine(id, ConCat("RUN_CELLS; -cells; 2; END", ""   ));
	if (RunAccumulated(id) != 0) EHandler();
	ExtractWrite(2);

	/* Destroy module */
	if (DestroyIPhreeqc(id) != IPQ_OK) EHandler();
	exit(EXIT_SUCCESS);
}