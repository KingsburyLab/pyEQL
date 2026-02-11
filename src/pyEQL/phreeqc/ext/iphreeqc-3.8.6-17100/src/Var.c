#include "Var.h"
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// The VarInit function initializes the VAR
// by setting the type field to VT_EMPTY
//
void VarInit(VAR* pvar)
{
	pvar->type    = TT_EMPTY;
	pvar->sVal    = 0;
	pvar->vresult = VR_OK;
}

VRESULT VarClear(VAR* pvar)
{
	switch (pvar->type)
	{
	case TT_EMPTY:
		break;
	case TT_LONG:
		break;
	case TT_DOUBLE:
		break;
	case TT_STRING:
		VarFreeString(pvar->sVal);
		break;
	case TT_ERROR:
		break;
	default:
		assert(0);
		return VR_BADVARTYPE;
	}
	VarInit(pvar);
	return VR_OK;
}

VRESULT VarCopy(VAR* pvarDest, const VAR* pvarSrc)
{
	VarClear(pvarDest);

	pvarDest->type = pvarSrc->type;
	switch (pvarSrc->type)
	{
	case TT_EMPTY:
		break;
	case TT_LONG:
		pvarDest->lVal = pvarSrc->lVal;
		break;
	case TT_DOUBLE:
		pvarDest->dVal = pvarSrc->dVal;
		break;
	case TT_STRING:
		pvarDest->sVal = VarAllocString(pvarSrc->sVal);
		if (pvarDest->sVal == NULL && pvarSrc->sVal != NULL) {
			pvarDest->type = TT_ERROR;
			pvarDest->vresult = VR_OUTOFMEMORY;
			return pvarDest->vresult;
		}
		break;
	case TT_ERROR:
		pvarDest->vresult = pvarSrc->vresult;
		break;
	default:
		assert(0);
		return VR_BADVARTYPE;
	}
	return VR_OK;
}

char* VarAllocString(const char* pSrc)
{
	char* psz;
	if (!pSrc) return NULL;
	psz = (char*) malloc(strlen(pSrc) + 1);
	if(psz != NULL) strcpy(psz, pSrc);
	return psz;
}

void VarFreeString(char* pSrc)
{
	if (pSrc) free(pSrc);
}
