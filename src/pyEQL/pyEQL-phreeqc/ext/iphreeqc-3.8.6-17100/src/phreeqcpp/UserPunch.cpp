#include "UserPunch.h"
#include "Phreeqc.h"

#if defined(PHREEQCI_GUI)
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#endif

UserPunch::UserPunch(int n, PHRQ_io *io)
:	cxxNumKeyword(io)
{
	this->PhreeqcPtr    = NULL;
	this->rate          = NULL;
}


UserPunch::~UserPunch(void)
{
	if (this->rate != NULL)
	{
		if (this->PhreeqcPtr != NULL)
		{
			this->PhreeqcPtr->rate_free(this->rate);
			delete this->rate;
		}
	}
	this->PhreeqcPtr = NULL;
	this->rate = NULL;
}
