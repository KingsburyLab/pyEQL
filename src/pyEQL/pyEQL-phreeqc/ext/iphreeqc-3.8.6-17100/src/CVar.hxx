#ifndef INC_CVAR_HXX
#define INC_CVAR_HXX

#include "Var.h"
#include <ostream>          // std::ostream

class CVar : public VAR
{
// Constructors
public:
	CVar(void) throw()
	{
		::VarInit(this);
	}
	~CVar(void) //throw()
	{
		Clear();
	}
	CVar(const VAR& varSrc)
	{
		this->type = TT_EMPTY;
		InternalCopy(&varSrc);
	}
	CVar(const CVar& varSrc)
	{
		this->type = TT_EMPTY;
		InternalCopy(&varSrc);
	}
	CVar(double dblSrc) throw()
	{
		this->type = TT_DOUBLE;
		this->dVal = dblSrc;
	}
	CVar(const char* pszSrc)
	{
		this->type = TT_EMPTY;
		*this = pszSrc;
	}
	CVar(long nSrc) throw()
	{
		this->type = TT_LONG;
		this->lVal = nSrc;
	}

// Assignment Operators
public:
	CVar& operator=(const CVar& varSrc)
	{
		InternalCopy(&varSrc);
		return *this;
	}
	CVar& operator=(const VAR& varSrc)
	{
		InternalCopy(&varSrc);
		return *this;
	}
	CVar& operator=(double dblSrc) throw()
	{
		if (type != TT_DOUBLE)
		{
			Clear();
			type = TT_DOUBLE;
		}
		this->dVal = dblSrc;
		return *this;
	}
	CVar& operator=(const char* pszSrc)
	{
		Clear();
		this->type = TT_STRING;
		this->sVal = ::VarAllocString(pszSrc);

		if (this->sVal == NULL && pszSrc != NULL)
		{
			this->type = TT_ERROR;
			this->vresult = VR_OUTOFMEMORY;
		}
		return *this;
	}

	friend std::ostream& operator<< (std::ostream &os, const CVar &a);


// Operations
public:
	VRESULT Clear(void) { return ::VarClear(this); }
	VRESULT Copy(const VAR* pSrc) { return ::VarCopy(this, const_cast<VAR*>(pSrc)); }

// Implementation
public:
	VRESULT InternalClear()
	{
		VRESULT vr = Clear();
		if (vr != VR_OK)
		{
			this->type = TT_ERROR;
			this->vresult = vr;
		}
		return vr;
	}

	void InternalCopy(const VAR* pSrc)
	{
		VRESULT vr = Copy(pSrc);
		if (vr != VR_OK)
		{
			this->type = TT_ERROR;
			this->vresult = vr;
		}
	}

};

inline std::ostream& operator<< (std::ostream &os, const CVar &a)
{
	switch (a.type) {
	case TT_EMPTY:
		os << "(TT_EMPTY)";
		break;
	case TT_ERROR:
		switch (a.vresult) {
		case VR_OK:
			os << "VR_OK";
			break;
		case VR_OUTOFMEMORY:
			os << "VR_OUTOFMEMORY";
			break;
		case VR_BADVARTYPE:
			os << "VR_BADVARTYPE";
			break;
		case VR_INVALIDARG:
			os << "VR_INVALIDARG";
			break;
		case VR_INVALIDROW:
			os << "VR_INVALIDROW";
			break;
		case VR_INVALIDCOL:
			os << "VR_INVALIDCOL";
			break;
		}
		os << "(TT_ERROR)";
		break;
	case TT_LONG:
		os << a.lVal;
		os << "(TT_LONG)";
		break;
	case TT_DOUBLE:
		os << a.dVal;
		os << "(TT_DOUBLE)";
		break;
	case TT_STRING:
		os << "\"" << a.sVal << "\"";
		os << "(TT_STRING)";
		break;
	default:
		os << "(BAD)";
		break;
	}
	return os;
}

#endif // INC_CVAR_HXX
