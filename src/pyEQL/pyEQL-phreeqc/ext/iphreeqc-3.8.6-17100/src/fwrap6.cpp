#ifdef IPHREEQC_NO_FORTRAN_MODULE
#include "IPhreeqc.h"
#include "fwrap.h"

#if defined(_WIN32) && !defined(_M_AMD64)

#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:stdref /names:uppercase /assume:underscore
//
#define IPQ_DECL __stdcall
#define IPQ_CASE_UND(name, NAME, name_, NAME_) NAME_

#include "fimpl.h"

#if defined(__cplusplus)
}
#endif

#endif // _WIN32
#endif

