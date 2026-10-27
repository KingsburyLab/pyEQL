#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <stdlib.h>				// ::tolower
#include <ctype.h>				// ::tolower
#include <algorithm>			// std::transform
#include <iostream>				// std::cout std::cerr
#include <string.h>

#include "Utils.h"
#include "Parser.h"
#include "float.h"
#include <cmath>

#if defined(_MSC_VER) && (_MSC_VER <= 1400) // VS2005
#  define nullptr NULL
#endif

#if __cplusplus < 201103L // Check if C++ standard is pre-C++11
#  ifndef nullptr
#    define nullptr NULL
#  endif
#endif

#if defined(PHREEQCI_GUI)
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#endif

////////////////////////////////////////////////////////////////////////////
int
Utilities::strcmp_nocase_arg1(const char *str1, const char *str2)
////////////////////////////////////////////////////////////////////////////
{
	//
	// Compare two strings disregarding case
	//
	int c1, c2;
	while ((c1 =::tolower(*str1++)) == (c2 = *str2++))
	{
		if (c1 == '\0')
			return (0);
	}
	if (c1 < c2)
		return (-1);
	return (1);
}

////////////////////////////////////////////////////////////////////////////
int
Utilities::strcmp_nocase(const char *str1, const char *str2)
////////////////////////////////////////////////////////////////////////////
{
	//
	// Compare two strings disregarding case
	//
	int c1, c2;
	while ((c1 =::tolower(*str1++)) == (c2 =::tolower(*str2++)))
	{
		if (c1 == '\0')
			return (0);
	}
	if (c1 < c2)
		return (-1);
	return (1);
}

////////////////////////////////////////////////////////////////////////////
void
Utilities::str_tolower(std::string & str)
////////////////////////////////////////////////////////////////////////////
{
	std::transform(str.begin(), str.end(), str.begin(), tolower);
}

////////////////////////////////////////////////////////////////////////////
void
Utilities::str_toupper(std::string & str)
////////////////////////////////////////////////////////////////////////////
{
	std::transform(str.begin(), str.end(), str.begin(), toupper);
}
////////////////////////////////////////////////////////////////////////////
std::string 
Utilities::pad_right(const std::string & str, size_t l)
////////////////////////////////////////////////////////////////////////////
{
	std::string new_str(str);
	size_t length = new_str.size();
	if (length < l)
	{
		new_str = new_str.insert(length, l - length, ' ');
	}
	return new_str;
}


////////////////////////////////////////////////////////////////////////////
bool
Utilities::replace(const char *str1, const char *str2, std::string & str)
////////////////////////////////////////////////////////////////////////////
{
	std::string::size_type n = str.find(str1, 0);
	if (n == std::string::npos)
		return false;

	str.replace(n, ::strlen(str1), str2);
	return true;
}

////////////////////////////////////////////////////////////////////////////
void
Utilities::squeeze_white(std::string & s_l)
////////////////////////////////////////////////////////////////////////////
{
	std::string str;
	std::string::iterator beg = s_l.begin();
	std::string::iterator end = s_l.end();
	//CParser::copy_token(str, beg, end);
	std::string::iterator pos;
	for (pos = beg; pos != end; pos++)
	{
		int c = *pos;
		if (!::isspace(c))
		{
			str += c;
		}
	}
	s_l = str;
}
////////////////////////////////////////////////////////////////////////////
double 
Utilities::convert_time(double t, std::string in, std::string out)
////////////////////////////////////////////////////////////////////////////
{
	Utilities::str_tolower(in);

	// convert t to seconds
	if (in.substr(0,1) == "m")
	{
		t = t * 60.;
	}
	if (in.substr(0,1) == "h")
	{
		t = t * 3600.;
	}
	if (in.substr(0,1) == "d")
	{
		t = t * 3600. * 24.;
	}
	if (in.substr(0,1) == "y")
	{
		t = t * 3600. * 24. * 365.25;
	}
	// convert to output units
	if (out.substr(0,1) == "m")
	{
		t = t / 60.;
	}
	if (out.substr(0,1) == "h")
	{
		t = t / 3600.;
	}
	if (out.substr(0,1) == "d")
	{
		t = t / ( 3600. * 24.);
	}
	if (out.substr(0,1) == "y")
	{
		t = t / (3600. * 24. * 365.25);
	}
	return t;

}
LDBLE 
Utilities::safe_exp(LDBLE t)
////////////////////////////////////////////////////////////////////////////
{
	LDBLE f = 1.442695*t; // convert to exp for 2.0

	if (f > DBL_MAX_EXP - 50.0)
	{
		return pow(2, DBL_MAX_EXP - 50.0);
	}
	else if (f < DBL_MIN_EXP + 50.0)
	{
		return pow(2, DBL_MIN_EXP + 50.0);
	}
	return exp(t);
}
size_t Utilities::
strcpy_safe(char* dest, size_t max, const char* src)
{
	size_t lsrc = 0;
	try
	{
		if (dest == nullptr || src == nullptr)
		{
#if !defined(R_SO)
			std::cerr << "nullptr in Utilities::strcpy_safe." << std::endl;
#endif
			throw;
		}
		lsrc = strlen(src);
		if (lsrc + 1 > max)
		{
#if !defined(R_SO)
			std::cerr << "Buffer overrun in Utilities::strcpy_safe." << std::endl;
#endif
			throw;
		}
		memcpy(dest, src, (lsrc + 1) * sizeof(char));
	}
	catch (...)
	{
		throw;
	}
	return lsrc;
}
size_t Utilities::
strcat_safe(char* dest, size_t max, const char* src)
{
	size_t ldest = 0, lsrc = 0;
	try
	{
		if (dest == nullptr || src == nullptr)
		{
#if !defined(R_SO)
			std::cerr << "nullptr in Utilities::strcat_safe." << std::endl;
#endif
			throw;
		}
		lsrc = strlen(src);
		ldest = strlen(dest);
		if (ldest + lsrc + 1 > max)
		{
#if !defined(R_SO)
			std::cerr << "Buffer overrun in Utilities::strcat_safe." << std::endl;
#endif
			throw;
		}
		memcpy(&dest[ldest], src, (lsrc + 1) * sizeof(char));
	}
	catch (...)
	{
		throw;
	}
	return ldest + lsrc;
}
//+NAN LDBLE: 7ff8000000000000
//-NAN LDBLE: fff8000000000000
/*
LDBLE Utilities::get_nan(void)
{
	unsigned long long raw = 0x7ff0000000000000;
	LDBLE d = *( LDBLE* )&raw;
	return(d);

}
*/
