#include <stdlib.h>
#include <IPhreeqc.h>

typedef int (*getFunc)(int);
typedef IPQ_RESULT (*setFunc)(int, int);
int TestGetSet(int, getFunc, setFunc);
int TestGetSetInitOn(int, getFunc, setFunc);

int
main(int argc, const char* argv[])
{
  int id;
  int r, c;
  VAR v;

  id = CreateIPhreeqc();
  if (id < 0)
  {
    return EXIT_FAILURE;
  }

  /* Dump */
  if (TestGetSet(id, GetDumpFileOn, SetDumpFileOn))
  {
    return EXIT_FAILURE;
  }

  /* Dump string */
  if (TestGetSet(id, GetDumpStringOn, SetDumpStringOn))
  {
    return EXIT_FAILURE;
  }

  /* Error file */
  if (TestGetSet(id, GetErrorFileOn, SetErrorFileOn))
  {
    return EXIT_FAILURE;
  }

  /* Error */
  if (TestGetSetInitOn(id, GetErrorOn, SetErrorOn))
  {
    return EXIT_FAILURE;
  }

  /* Log */
  if (TestGetSet(id, GetLogFileOn, SetLogFileOn))
  {
    return EXIT_FAILURE;
  }

  /* Output */
  if (TestGetSet(id, GetOutputFileOn, SetOutputFileOn))
  {
    return EXIT_FAILURE;
  }

  /* Selected output */
  if (TestGetSet(id, GetSelectedOutputFileOn, SetSelectedOutputFileOn))
  {
    return EXIT_FAILURE;
  }

  if (LoadDatabase(id, "phreeqc.dat") != 0)
  {
    OutputErrorString(id);
    return EXIT_FAILURE;
  }

  if (RunFile(id, "ex2") != 0)
  {
    OutputErrorString(id);
    return EXIT_FAILURE;
  }

  VarInit(&v);
  for (r = 0; r < GetSelectedOutputRowCount(id); ++r)
  {
    for (c = 0; c < GetSelectedOutputColumnCount(id); ++c)
    {
      if (GetSelectedOutputValue(id, r, c, &v) != IPQ_OK)
      {
        return EXIT_FAILURE;
      }
      VarClear(&v);
    }
  }
  
  if (DestroyIPhreeqc(id) != IPQ_OK)
  {
    OutputErrorString(id);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int 
TestGetSet(int id, getFunc gf, setFunc sf)
{
  if (gf(id))
  {  
    return EXIT_FAILURE;
  }
  
  if (sf(id, 1) != IPQ_OK)
  {
    return EXIT_FAILURE;
  }
  
  if (!gf(id))
  {
    return EXIT_FAILURE;
  }
  
  if (sf(id,0) != IPQ_OK)
  {
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}

int
TestGetSetInitOn(int id, getFunc gf, setFunc sf)
{
    if (!gf(id))
    {
        return EXIT_FAILURE;
    }

    if (sf(id, 0) != IPQ_OK)
    {
        return EXIT_FAILURE;
    }

    if (gf(id))
    {
        return EXIT_FAILURE;
    }

    if (sf(id, 1) != IPQ_OK)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
