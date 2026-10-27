#include <IPhreeqc.hpp>
#include <iostream>

int main (void)
{
    IPhreeqc iphreeqc;
    if (iphreeqc.LoadDatabase("phreeqc.dat") != 0)
    {
        std::cout << iphreeqc.GetErrorString();
        return EXIT_FAILURE;
    }
    iphreeqc.SetOutputFileOn(true);
    iphreeqc.SetOutputFileName("ex2.out");
    if (iphreeqc.RunFile("ex2") != 0)
    {
        std::cout << iphreeqc.GetErrorString();
        return EXIT_FAILURE;
    }
    return 0;
}
