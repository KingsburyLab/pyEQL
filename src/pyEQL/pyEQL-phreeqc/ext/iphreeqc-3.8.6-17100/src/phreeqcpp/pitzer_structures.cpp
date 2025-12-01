#include "Phreeqc.h"
#include "phqalloc.h"

#include <list>
#include <string>

#if defined(PHREEQCI_GUI)
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#endif

/* **********************************************************************
 *
 *   Routines related to structure "pitz_param"
 *
 * ********************************************************************** */

/* ---------------------------------------------------------------------- */
class pitz_param * Phreeqc::
pitz_param_read(char *string, int n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read pitzer parameter info from string
 *   n is number of species (character values)
 *          
 */
	int l, i, j, k;
	const char* cptr;
	char token[2 * MAX_LENGTH];
	class pitz_param pzp, *pzp_ptr;

	if (n != 2 && n != 3 && n != 0)
		return (NULL);
	if (string == NULL)
		return (NULL);

	cptr = string;
	if (copy_token(token, &cptr, &l) == EMPTY)
		return (NULL);
	cptr = string;
	for (i = 0; i < n; i++)
	{
		int j = copy_token(token, &cptr, &l);
		if (j == EMPTY)
			return (NULL);
		if (j != UPPER && token[0] != '(')
		{
			input_error++;
			std::ostringstream err;
			err << "Wrong number of species for a Pitzer parameter.\n"  << line;
			error_msg(err.str().c_str(), CONTINUE);
		}
		pzp.species[i] = string_hsave(token);
	}
	k = 0;
	for (i = 0; i < 6; i++)
	{
		if (copy_token(token, &cptr, &l) == EMPTY)
			break;
		j = sscanf(token, SCANFORMAT, &pzp.a[i]);
		if (j <= 0)
			break;
		k++;
	}
	if (k <= 0)
		return (NULL);
	pzp_ptr = new class pitz_param;
	*pzp_ptr = pzp;
	return (pzp_ptr);
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
pitz_param_store(class pitz_param *pzp_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *  Does linear search of pitz_params for same type and species
 *  Returns -1 if not found, index number in pitz_params if found
 */
	int i;
	if (pzp_ptr == NULL)
		return;
	if (pzp_ptr->type == TYPE_Other)
		return;
	std::set< std::string > header;
	for (i = 0; i < 3; i++)
	{
		if (pzp_ptr->species[i] != NULL) header.insert(pzp_ptr->species[i]);
	}

	std::ostringstream key_str;
	key_str << pzp_ptr->type << " ";
	std::set< std::string >::iterator it = header.begin();
	for(; it != header.end(); ++it)
	{
		key_str << *it << " ";
	}
	std::string key = key_str.str().c_str();
	std::map< std::string, size_t>::iterator jit = pitz_param_map.find(key);
	if (jit != pitz_param_map.end())
	{
		if (pzp_ptr->species[2] != NULL)
		{
			error_string = sformatf( "Redefinition of parameter, %s %s %s\n", 
			pzp_ptr->species[0], pzp_ptr->species[1], pzp_ptr->species[2]);
		}
		else
		{
			error_string = sformatf( "Redefinition of parameter, %s %s\n", 
			pzp_ptr->species[0], pzp_ptr->species[1]);
		}
	    warning_msg(error_string);
		delete pitz_params[(*jit).second]; 
		pitz_params[(*jit).second] = pzp_ptr;
	}
	else
	{
		size_t count_pitz_param = pitz_params.size();
		pitz_params.resize(count_pitz_param + 1);
		pitz_params[count_pitz_param] = pzp_ptr;
		pitz_param_map[key] = count_pitz_param;
	}
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
sit_param_store(class pitz_param *pzp_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *  Does linear search of pitz_params for same type and species
 *  Returns -1 if not found, index number in pitz_params if found
 */
	int i;
	if (pzp_ptr == NULL)
		return;
	if (pzp_ptr->type == TYPE_Other)
		return;

	std::set< std::string > header;
	for (i = 0; i < 3; i++)
	{
		if (pzp_ptr->species[i] != NULL) header.insert(pzp_ptr->species[i]);
	}

	std::ostringstream key_str;
	key_str << pzp_ptr->type << " ";
	std::set< std::string >::iterator it = header.begin();
	for(; it != header.end(); ++it)
	{
		key_str << *it << " ";
	}
	std::string key = key_str.str().c_str();

	std::map< std::string, size_t>::iterator jit = sit_param_map.find(key);
	if (jit != sit_param_map.end())
	{
		if (pzp_ptr->species[2] != NULL)
		{
			error_string = sformatf( "Redefinition of parameter, %s %s %s\n", 
			pzp_ptr->species[0], pzp_ptr->species[1], pzp_ptr->species[2]);
		}
		else
		{
			error_string = sformatf( "Redefinition of parameter, %s %s\n", 
			pzp_ptr->species[0], pzp_ptr->species[1]);
		}
	    warning_msg(error_string);
		delete sit_params[(*jit).second]; 
		sit_params[(*jit).second] = pzp_ptr;
	}
	else
	{
		size_t count_sit_param = sit_params.size();
		sit_params.resize(count_sit_param + 1);
		sit_params[count_sit_param] = pzp_ptr;
		sit_param_map[key] = count_sit_param;
	}
}
class pitz_param* Phreeqc::
pitz_param_copy(const class pitz_param* src)
{
	if (src == NULL) return NULL;
	class pitz_param* dest = new class pitz_param;
	*dest = *src;
	for (size_t i = 0; i < 3; i++)
	{
		if (src->species[i] != NULL)
		{
			dest->species[i] = string_hsave(src->species[i]);
		}
	}
	dest->thetas = NULL;
	return dest;
}

/* **********************************************************************
 *
 *   Routines related to structure "theta_parm"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
class theta_param * Phreeqc::
theta_param_search(LDBLE zj, LDBLE zk)
/* ---------------------------------------------------------------------- */
{
/*
 *  Does linear search of theta_params for same charge
 *  Returns NULL if not found, index number in theta_params if found
 */
	int i;
	for (i = 0; i < (int)theta_params.size(); i++)
	{
		if ((theta_params[i]->zj == zj && theta_params[i]->zk == zk) ||
			(theta_params[i]->zj == zk && theta_params[i]->zk == zj))
		{
			return theta_params[i];
		}
	}
	return NULL;
}
