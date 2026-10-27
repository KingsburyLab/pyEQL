// ISolution.cxx: implementation of the cxxSolutionxx class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <cassert>				// assert
#include <algorithm>			// std::sort
#include <sstream>

#include "Utils.h" // define first
#include "Phreeqc.h"
#include "ISolution.h"
#include "phqalloc.h"

#if defined(PHREEQCI_GUI)
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxISolution::cxxISolution(PHRQ_io *io)
:
units("mMol/kgw")
{
	default_pe = "pe";
	CReaction temp_pe_reactions;
	pe_reactions[default_pe] = temp_pe_reactions;
	this->calc_density = false;

}
cxxISolution::~cxxISolution()
{
}
