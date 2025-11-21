#include <cstdlib>
#include <iostream>
#include <IPhreeqc.hpp>

template <class TClass> class TTestGetSet
{
private:
  bool (TClass::*_get)(void)const;
  void (TClass::*_set)(bool);
  TClass* _p;
public:
  TTestGetSet(TClass* p, bool(TClass::*get)(void)const, void(TClass::*set)(bool))
    {
      _p = p;
      _get = get;
      _set = set;
    }
  int Test(void)
    {
      if ((*_p.*_get)())
      {
        return EXIT_FAILURE;
      }

      (*_p.*_set)(true);

      if (!(*_p.*_get)())
      {
        return EXIT_FAILURE;
      }

      (*_p.*_set)(false);

      if ((*_p.*_get)())
      {
        return EXIT_FAILURE;
      }

      return EXIT_SUCCESS;
    }

  int TestInitOn(void)
  {
      if (!(*_p.*_get)())
      {
          return EXIT_FAILURE;
      }

      (*_p.*_set)(false);

      if ((*_p.*_get)())
      {
          return EXIT_FAILURE;
      }

      (*_p.*_set)(true);

      if (!(*_p.*_get)())
      {
          return EXIT_FAILURE;
      }

      return EXIT_SUCCESS;
  }
};

int
main(int argc, const char* argv[])
{
  IPhreeqc iphreeqc;

  // Dump
  TTestGetSet<IPhreeqc> testDump(&iphreeqc, &IPhreeqc::GetDumpFileOn, &IPhreeqc::SetDumpFileOn);
  if (testDump.Test() != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }

  // Dump string
  TTestGetSet<IPhreeqc> testDumpString(&iphreeqc, &IPhreeqc::GetDumpStringOn, &IPhreeqc::SetDumpStringOn);
  if (testDumpString.Test() != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }

  // Error file
  TTestGetSet<IPhreeqc> testErrorFile(&iphreeqc, &IPhreeqc::GetErrorFileOn, &IPhreeqc::SetErrorFileOn);
  if (testErrorFile.Test() != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }

  // Error
  TTestGetSet<IPhreeqc> testError(&iphreeqc, &IPhreeqc::GetErrorOn, &IPhreeqc::SetErrorOn);
  if (testError.TestInitOn() != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }


  // Log
  TTestGetSet<IPhreeqc> testLog(&iphreeqc, &IPhreeqc::GetLogFileOn, &IPhreeqc::SetLogFileOn);
  if (testLog.Test() != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }

  // Output
  TTestGetSet<IPhreeqc> testOutput(&iphreeqc, &IPhreeqc::GetOutputFileOn, &IPhreeqc::SetOutputFileOn);
  if (testOutput.Test() != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }

  // Selected output
  TTestGetSet<IPhreeqc> testSelectedOutput(&iphreeqc, &IPhreeqc::GetSelectedOutputFileOn, &IPhreeqc::SetSelectedOutputFileOn);
  if (testSelectedOutput.Test() != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }

  if (iphreeqc.LoadDatabase("phreeqc.dat") != 0)
  {
    std::cout << iphreeqc.GetErrorString();
    return EXIT_FAILURE;
  }

  if (iphreeqc.RunFile("ex2") != 0)
  {
    std::cout << iphreeqc.GetErrorString();
    return EXIT_FAILURE;
  }

  VAR v;
  ::VarInit(&v);
  for (int r = 0; r < iphreeqc.GetSelectedOutputRowCount(); ++r)
  {
    for (int c = 0; c < iphreeqc.GetSelectedOutputColumnCount(); ++c)
    {
      if (iphreeqc.GetSelectedOutputValue(r, c, &v) != VR_OK)
      {
        return EXIT_FAILURE;
      }
	  if (::VarClear(&v) != VR_OK)
      {
        return EXIT_FAILURE;
      }
    }
  }

  return EXIT_SUCCESS;
}
