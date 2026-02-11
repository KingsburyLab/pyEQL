#if !(defined(WIN32) && !defined(__GNUC__))
#include <assert.h>
#define _ASSERTE assert
#endif
#include <stdlib.h>
#include "Phreeqc.h"
#include "PBasic.h"
#include "phqalloc.h"
#include "NameDouble.h"
#include "Utils.h"
#include "Solution.h"

/* Run-time library for PhreeqcPtr->use with "p2c", the Pascal to C translator */

/* "p2c"  Copyright (C) 1989, 1990, 1991 Free Software Foundation.
 * By Dave Gillespie, daveg@csvax.cs.caltech.edu.  Version 1.20.
 * This file may be copied, modified, etc. in any way.  It is not restricted
 * by the licence agreement accompanying p2c itself.
 */

#define STOP 1
#define CONTINUE 0
#define Isspace(c)  isspace(c)	/* or "((c) == ' ')" if preferred */
#define toklength       20
typedef long chset[9];

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

/* Output from p2c, the Pascal-to-C translator */
/* From input file "basic.p" */

PBasic::PBasic(Phreeqc * ptr, PHRQ_io *phrq_io)
	: PHRQ_base(phrq_io)
{
	if (ptr == NULL)
	{
		error_msg("No Phreeqc instance in PBasic constructor\n", 1);
	}
	PhreeqcPtr = ptr;
	inbuf = NULL;
	linebase = NULL;
	varbase = NULL;
	loopbase = NULL;
	curline = 0;
	stmtline = NULL;
	dataline = NULL;
	stmttok = NULL;
	datatok = NULL;
	buf = NULL;
	exitflag = false;
	EXCP_LINE = 0;
	P_escapecode = 0;
	P_ioresult = 0;
	parse_all = false;
	phreeqci_gui = false;
	parse_whole_program = true;
#if defined(PHREEQCI_GUI)
	hInfiniteLoop = 0;
	nIDErrPrompt = 0;
#else
	nIDErrPrompt = (PBasic::IDErr)0;
#endif
	nErrLineNumber = 0;
	punch_tab = true;
	skip_punch = false;
	// Basic commands initialized at bottom of file
}
PBasic::~PBasic(void)
{
}

int PBasic::
basic_compile(const char *commands, void **lnbase, void **vbase, void **lpbase)
{								/*main */
	int l;
	const char *ptr;

	P_escapecode = 0;
	P_ioresult = 0;
	inbuf = (char *) PhreeqcPtr->PHRQ_calloc(PhreeqcPtr->max_line, sizeof(char));
	if (inbuf == NULL)
		PhreeqcPtr->malloc_error();
	linebase = NULL;
	varbase = NULL;
	loopbase = NULL;
	exitflag = false;
	ptr = commands;
	do
	{
		try
		{
			ptr = commands;
			do
			{
				if (sget_logical_line(&ptr, &l, inbuf) == EOF)
				{
					strcpy(inbuf, "bye");
				}
				parseinput(&buf);
				if (curline == 0)
				{
					stmtline = NULL;
					stmttok = buf;
					if (stmttok != NULL)
						exec();
					disposetokens(&buf);
				}
			}
			while (!(exitflag || P_eof()));
		}
		catch (const PBasicStop&)
		{
			if (P_escapecode != -20)
			{
				if (phreeqci_gui)
				{
					_ASSERTE(false);
				}
				else
				{
					char * error_string = PhreeqcPtr->sformatf("%d/%d", (int) P_escapecode,
						(int) P_ioresult);
					PhreeqcPtr->warning_msg(error_string);
				}
			}
			else
			{
				if (phreeqci_gui)
				{
					_ASSERTE(false);
				}
				else
				{
					/*    putchar('\n');*/
					output_msg("\n");
				}
			}
		}
		catch (const PhreeqcStop&)
		{
			// clean up memory
			disposetokens(&buf);
			PhreeqcPtr->PHRQ_free(inbuf);
			*lnbase = (void *) linebase;
			*vbase = (void *) varbase;
			*lpbase = (void *) loopbase;
			throw;  // rethrow
		}
	}
	while (!(exitflag || P_eof()));
	/*  exit(EXIT_SUCCESS); */
	PhreeqcPtr->PHRQ_free(inbuf);
	*lnbase = (void *) linebase;
	*vbase = (void *) varbase;
	*lpbase = (void *) loopbase;
	return (P_escapecode);
}

int PBasic::
basic_renumber(char *commands, void **lnbase, void **vbase, void **lpbase)
{								/*main */
	int l, i;
	const char *ptr;

	P_escapecode = 0;
	P_ioresult = 0;
	inbuf = (char *) PhreeqcPtr->PHRQ_calloc(PhreeqcPtr->max_line, sizeof(char));
	if (inbuf == NULL)
		PhreeqcPtr->malloc_error();
	linebase = NULL;
	varbase = NULL;
	loopbase = NULL;
	exitflag = false;
	ptr = commands;
	do
	{
		try
		{
			i = 0;
			ptr = commands;
			do
			{
				if (sget_logical_line(&ptr, &l, inbuf) == EOF)
				{
					i++;
					if (i == 1)
					{
						strcpy(inbuf, "renum");
					}
					else if (i == 2)
					{
						strcpy(inbuf, "list");
					}
					else if (i == 3)
					{
						strcpy(inbuf, "new");
					}
					else if (i == 4)
					{
						strcpy(inbuf, "bye");
					}
				}
				parseinput(&buf);
				if (curline == 0)
				{
					stmtline = NULL;
					stmttok = buf;
					if (stmttok != NULL)
						exec();
					disposetokens(&buf);
				}
			}
			while (!(exitflag || P_eof()));
		}
		catch (const PBasicStop&)
		{
			if (P_escapecode != -20)
			{
				char * error_string = PhreeqcPtr->sformatf( "%d/%d", (int) P_escapecode,
					(int) P_ioresult);
				PhreeqcPtr->warning_msg(error_string);
			}
			else
			{
				/*    putchar('\n');*/
				output_msg("\n");
			}
		}
	}
	while (!(exitflag || P_eof()));
	/*  exit(EXIT_SUCCESS); */
	PhreeqcPtr->PHRQ_free(inbuf);
	*lnbase = (void *) linebase;
	*vbase = (void *) varbase;
	*lpbase = (void *) loopbase;

	return (P_escapecode);
}

int PBasic::
basic_run(char *commands, void *lnbase, void *vbase, void *lpbase)
{								/*main */
	int l;
	const char *ptr;
	P_escapecode = 0;
	P_ioresult = 0;
	inbuf = (char *) PhreeqcPtr->PHRQ_calloc(PhreeqcPtr->max_line, sizeof(char));
	if (inbuf == NULL)
		PhreeqcPtr->malloc_error();
	linebase = NULL;
	varbase = NULL;
	loopbase = NULL;
	exitflag = false;
	ptr = commands;
	linebase = (linerec *) lnbase;
	varbase = (varrec *) vbase;
	loopbase = (looprec *) lpbase;
	do
	{
		try
		{
			do
			{
				if (sget_logical_line(&ptr, &l, inbuf) == EOF)
				{
					strcpy(inbuf, "bye");
				}
				parseinput(&buf);
				if (curline == 0)
				{
					stmtline = NULL;
					stmttok = buf;
					if (stmttok != NULL)
						exec();
					disposetokens(&buf);
				}
			}
			while (!(exitflag || P_eof()));
		}
		catch (const PBasicStop&)
		{
			if (P_escapecode != -20)
			{
				if (phreeqci_gui)
				{
					_ASSERTE(FALSE);
				}
				else
				{
					char * error_string = PhreeqcPtr->sformatf( "%d/%d", (int) P_escapecode,
						(int) P_ioresult);
					PhreeqcPtr->warning_msg(error_string);
				}
			}
			else
			{
				/*    putchar('\n');*/
				output_msg("\n");
			}
		}
	}
	while (!(exitflag || P_eof()));

	/*  exit(EXIT_SUCCESS); */
	PhreeqcPtr->PHRQ_free(inbuf);

	// Cleanup after run
	clearvars();
	clearloops();
	restoredata();

	return (P_escapecode);
}

int PBasic::
basic_main(const char *commands)
{								/*main */
	int l;
	const char *ptr;

	P_escapecode = 0;
	P_ioresult = 0;
	inbuf = (char *) PhreeqcPtr->PHRQ_calloc(PhreeqcPtr->max_line, sizeof(char));
	if (inbuf == NULL)
		PhreeqcPtr->malloc_error();
	linebase = NULL;
	varbase = NULL;
	loopbase = NULL;
#ifdef SKIP
	printf("Chipmunk BASIC 1.0\n\n");
#endif
	exitflag = false;
	ptr = commands;
	do
	{
		try
		{
			do
			{
				if (sget_logical_line(&ptr, &l, inbuf) == EOF)
				{
					strcpy(inbuf, "bye");
				}
				parseinput(&buf);
				if (curline == 0)
				{
					stmtline = NULL;
					stmttok = buf;
					if (stmttok != NULL)
						exec();
					disposetokens(&buf);
				}
			}
			while (!(exitflag || P_eof()));
		}
		catch (const PBasicStop&)
		{
			if (P_escapecode != -20)
			{
				char * error_string = PhreeqcPtr->sformatf( "%d/%d", (int) P_escapecode,
					(int) P_ioresult);
				PhreeqcPtr->warning_msg(error_string);
			}
			else
			{
				/*    putchar('\n');*/
				output_msg("\n");
			}
		}
	}
	while (!(exitflag || P_eof()));
	return 1;
/*  exit(EXIT_SUCCESS); */
}

/* ---------------------------------------------------------------------- */
int PBasic::
sget_logical_line(const char **ptr, int *l, char *return_line)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads file fp until end of line, ";", or eof
 *   stores characters in line_save
 *   reallocs line_save and line if more space is needed
 *
 *   returns:
 *	   EOF on empty line on end of file or
 *	   OK otherwise
 *	   *l returns length of line
 */
	int i;
	char c;
	i = 0;
	if (**ptr == '\0')
		return (EOF);
	for (;;)
	{
		c = **ptr;
		if (c == '\0')
			break;
		(*ptr)++;
		if (c == ';' || c == '\n')
			break;
		return_line[i++] = c;
	}
	return_line[i] = '\0';
	*l = i;
	return (1);
}
void PBasic::
restoredata(void)
{
	dataline = NULL;
	datatok = NULL;
}

void PBasic::
clearloops(void)
{
	looprec *l;

	while (loopbase != NULL)
	{
		l = loopbase->next;
		PhreeqcPtr->PHRQ_free(loopbase);
		loopbase = l;
	}
}

void PBasic::
clearvar(varrec * v)
{
	if (v->numdims != 0)
	{
		if (v->stringvar == 0)
		{
			PhreeqcPtr->PHRQ_free(v->UU.U0.arr);
			v->UU.U0.arr = NULL;
		}
		else
		{
			free_dim_stringvar(v);
		}
	}
	else if (v->stringvar && v->UU.U1.sv != NULL)
	{
		PhreeqcPtr->PHRQ_free(v->UU.U1.sv);
	}
	v->numdims = 0;
	if (v->stringvar)
	{
		v->UU.U1.sv = NULL;
		v->UU.U1.sval = &v->UU.U1.sv;
	}
	else
	{
		v->UU.U0.rv = 0.0;
		v->UU.U0.val = &v->UU.U0.rv;
	}
}

void PBasic::
clearvars(void)
{
	varrec *v;

	v = varbase;
	while (v != NULL)
	{
		clearvar(v);
		v = v->next;
	}
}

char * PBasic::
numtostr(char * Result, LDBLE n)
{
	char *l_s;
	long i;

	l_s = (char *) PhreeqcPtr->PHRQ_calloc(PhreeqcPtr->max_line, sizeof(char));
	if (l_s == NULL)
	{
		PhreeqcPtr->malloc_error();
#if !defined(R_SO)
		exit(4);
#endif
	}
	l_s[PhreeqcPtr->max_line - 1] = '\0';
/*  if ((n != 0 && fabs(n) < 1e-2) || fabs(n) >= 1e12) { */
	if (ceil(n) == floor(n))
	{
		//if (PhreeqcPtr->current_selected_output != NULL &&
		//	!PhreeqcPtr->current_selected_output->Get_high_precision())
		//{
		//	snprintf(l_s, PhreeqcPtr->max_line, "%12.0f", (double) n);
		//}
		//else
		//{
		//	snprintf(l_s, PhreeqcPtr->max_line, "%20.0f", (double) n);
		//}
		bool temp_high_precision = (PhreeqcPtr->current_selected_output != NULL) ? 
			PhreeqcPtr->current_selected_output->Get_high_precision() : 
			PhreeqcPtr->high_precision;
		if (!temp_high_precision)
		{
			snprintf(l_s, PhreeqcPtr->max_line, "%12.0f", (double) n);
		}
		else
		{
			snprintf(l_s, PhreeqcPtr->max_line, "%20.0f", (double) n);
		}
	}
	else
	{
		bool temp_high_precision = (PhreeqcPtr->current_selected_output != NULL) ? 
			PhreeqcPtr->current_selected_output->Get_high_precision() : 
			PhreeqcPtr->high_precision;
		if (!temp_high_precision)
		{
			snprintf(l_s, PhreeqcPtr->max_line, "%12.4e", (double) n);
		}
		else
		{
			snprintf(l_s, PhreeqcPtr->max_line, "%20.12e", (double) n);
		}
	}
	i = (int) strlen(l_s) + 1;
	l_s[i - 1] = '\0';
/* p2c: basic.p, line 237:
 * Note: Modification of string length may translate incorrectly [146] */
	strcpy(Result, l_s);
	PhreeqcPtr->free_check_null(l_s);
	return (Result);
/*  } else {
    if (PhreeqcPtr->punch.high_precision == FALSE) snprintf(l_s, PhreeqcPtr->max_line, "%30.10f", n);
      else snprintf(l_s, PhreeqcPtr->max_line, "%30.12f", n);
    i = strlen(l_s) + 1;
    do {
      i--;
    } while (l_s[i - 1] == '0');
    if (l_s[i - 1] == '.')
      i--;
    l_s[i] = '\0';
 * p2c: basic.p, line 248:
 * Note: Modification of string length may translate incorrectly [146] *
     Utilities::strcpy_safe(Result, MAX_LENGTH, strltrim(l_s));
	 return Result;
  } */
}

void PBasic::
parse(char * l_inbuf, tokenrec ** l_buf)
{
	long i, j, begin, len, m, lp, q;
	char token[toklength + 1] = {0};
	tokenrec *t, *tptr;
	varrec *v;
	char ch;
	char *ptr;

	tptr = NULL;
	*l_buf = NULL;
	i = 1;
	lp = q = 0;
	do
	{
		ch = ' ';
		while (i <= (int) strlen(l_inbuf) && (ch == ' ' || ch == '\t'))
		{
			ch = l_inbuf[i - 1];
			i++;
		}
		if (ch != ' ')
		{
			t = (tokenrec *) PhreeqcPtr->PHRQ_calloc(1, sizeof(tokenrec));
			if (t == NULL)
				PhreeqcPtr->malloc_error();
			if (tptr == NULL)
				*l_buf = t;
			else
				tptr->next = t;
			if (phreeqci_gui)
			{
				t->n_sz = 0;
				t->sz_num = 0;
			}
			tptr = t;
			t->next = NULL;
			switch (ch)
			{

			case '"':
			case '\'':
				q += 1;
				t->kind = tokstr;
				j = 0;
				len = (int) strlen(l_inbuf);
				begin = i;
				while (i <= len && l_inbuf[i - 1] != ch)
				{
					++j;
					++i;
				}
				if (l_inbuf[i - 1] == ch) q -= 1;
				m = 256;
				if (j + 1 > m)
					m = j + 1;
				t->UU.sp = (char *) PhreeqcPtr->PHRQ_calloc(m, sizeof(char));
				t->sp_sz = m;
				if (t->UU.sp == NULL)
				{
					PhreeqcPtr->malloc_error();
#if !defined(R_SO)
					exit(4);
#endif
				}
				strncpy(t->UU.sp, l_inbuf + begin - 1, j);
				t->UU.sp[j] = '\0';
/* p2c: basic.p, line 415:
 * Note: Modification of string length may translate incorrectly [146] */
				i++;
				break;

			case '+':
				t->kind = tokplus;
				break;

			case '-':
				t->kind = tokminus;
				break;

			case '*':
				t->kind = toktimes;
				break;

			case '/':
				t->kind = tokdiv;
				break;

			case '^':
				t->kind = tokup;
				break;

			case '(':
			case '[':
				t->kind = toklp;
				lp += 1;
				break;

			case ')':
			case ']':
				t->kind = tokrp;
				lp -= 1;
				break;

			case ',':
				t->kind = tokcomma;
				break;

			case ';':
				t->kind = toksemi;
				break;

			case ':':
				t->kind = tokcolon;
				break;

			case '?':
				t->kind = tokprint;
				break;

			case '=':
				t->kind = tokeq;
				break;

			case '<':
				if (i <= (int) strlen(l_inbuf) && l_inbuf[i - 1] == '=')
				{
					t->kind = tokle;
					i++;
				}
				else if (i <= (int) strlen(l_inbuf) && l_inbuf[i - 1] == '>')
				{
					t->kind = tokne;
					i++;
				}
				else
					t->kind = toklt;
				break;

			case '>':
				if (i <= (int) strlen(l_inbuf) && l_inbuf[i - 1] == '=')
				{
					t->kind = tokge;
					i++;
				}
				else
					t->kind = tokgt;
				break;

			default:
				if (isalpha((int) ch))
				{
					i--;
					j = 0;
					token[toklength] = '\0';
					while (i <= (int) strlen(l_inbuf) &&
						   (l_inbuf[i - 1] == '$' || l_inbuf[i - 1] == '_' ||
							isalnum((int) l_inbuf[i - 1])))
					{
						if (j < toklength)
						{
							j++;
							token[j - 1] = l_inbuf[i - 1];
						}
						i++;
					}
					token[j] = '\0';
/* p2c: basic.p, line 309:
 * Note: Modification of string length may translate incorrectly [146] */

/*
 *   Search list
 */
					PhreeqcPtr->str_tolower(token);
					std::map<const std::string, BASIC_TOKEN>::const_iterator item;
					item = command_tokens.find(token);
					if (item != command_tokens.end())
					{
						t->kind = item->second;
						if (t->kind == tokrem)
						{
							m = (int) strlen(l_inbuf) + 1;
							if (m < 256)
								m = 256;
							t->UU.sp = (char *) PhreeqcPtr->PHRQ_calloc(m, sizeof(char));
							t->sp_sz = m;
							if (t->UU.sp == NULL)
							{
								PhreeqcPtr->malloc_error();
#if !defined(R_SO)
								exit(4);
#endif
							}
							snprintf(t->UU.sp, t->sp_sz, "%.*s",
									(int) (strlen(l_inbuf) - i + 1),
									l_inbuf + i - 1);
							i = (int) strlen(l_inbuf) + 1;
						}
					}
					else
					{
						t->kind = tokvar;
						v = varbase;
						while (v != NULL && strcmp(v->name, token))
							v = v->next;
						if (v == NULL)
						{
							v = (varrec *) PhreeqcPtr->PHRQ_calloc(1, sizeof(varrec));
							if (v == NULL)
							{
								PhreeqcPtr->malloc_error();
#if !defined(R_SO)
								exit(4);
#endif
							}
							v->UU.U0.arr = NULL;
							v->next = varbase;
							varbase = v;
							strcpy(v->name, token);
							v->numdims = 0;
							if (token[strlen(token) - 1] == '$')
							{
								v->stringvar = true;
								v->UU.U1.sv = NULL;
								v->UU.U1.sval = &v->UU.U1.sv;
							}
							else
							{
								v->stringvar = false;
								v->UU.U0.rv = 0.0;
								v->UU.U0.val = &v->UU.U0.rv;
							}
						}
						t->UU.vp = v;
					}
				}
				else if (isdigit((int) ch) || ch == '.')
				{
					t->kind = toknum;
					i--;
					t->UU.num = strtod(&l_inbuf[i - 1], &ptr);
					if (&l_inbuf[i - 1] == ptr)
					{
						/*
						   Note: the following causes an infinite loop:
						   X = ..9
						 */
						t->kind = toksnerr;
						t->UU.snch = ch;
						i++;
						break;
					}
#if defined(PHREEQCI_GUI)
					if (phreeqci_gui)
					{
						_ASSERTE(t->n_sz == 0);
						_ASSERTE(t->sz_num == NULL);
						t->n_sz = max(23, ptr - &inbuf[i - 1]);
						t->sz_num =
							(char *) PhreeqcPtr->PHRQ_calloc((t->n_sz + 1), sizeof(char));
						if (t->sz_num == NULL)
							PhreeqcPtr->malloc_error();
						if (ptr > &inbuf[i - 1])
						{
							strncpy(t->sz_num, &inbuf[i - 1],
								(ptr - &inbuf[i - 1]));
							t->sz_num[ptr - &inbuf[i - 1]] = '\0';
						}
						else
						{
							t->sz_num[0] = '\0';
						}
					}
#endif
					i += (int) (ptr - &l_inbuf[i - 1]);
				}
				else
				{
					t->kind = toksnerr;
					t->UU.snch = ch;
				}
				break;
			}
		}
	}
	while (i <= (int) strlen(l_inbuf));
	if (q) {
		if (phreeqci_gui)
		{
			_ASSERTE(nIDErrPrompt == 0);
			_ASSERTE(P_escapecode == 0);
			nIDErrPrompt = IDS_ERR_MISSING_Q;
			P_escapecode = -20;
			return;
		}
		else
		{
			char * error_string = PhreeqcPtr->sformatf( " missing \" or \' in BASIC line\n %ld %s", curline, l_inbuf);
			error_msg(error_string, STOP);
		}
	}
	if (lp > 0) {
		if (phreeqci_gui)
		{
			_ASSERTE(nIDErrPrompt == 0);
			_ASSERTE(P_escapecode == 0);
			nIDErrPrompt = IDS_ERR_MISSING_RP;
			P_escapecode = -20;
			return;
		}
		else
		{
			char * error_string = PhreeqcPtr->sformatf( " missing ) or ] in BASIC line\n %ld %s", curline, l_inbuf);
			error_msg(error_string, STOP);
		}
	}
	else if (lp < 0) {
		if (phreeqci_gui)
		{
			_ASSERTE(nIDErrPrompt == 0);
			_ASSERTE(P_escapecode == 0);
			nIDErrPrompt = IDS_ERR_MISSING_RP;
			P_escapecode = -20;
			return;
		}
		else
		{
			char * error_string = PhreeqcPtr->sformatf( " missing ( or [ in BASIC line\n %ld %s", curline, l_inbuf);
			error_msg(error_string, STOP);
		}
	}
}

#undef toklength

void PBasic::
listtokens(FILE * f, tokenrec * l_buf)
{
	bool ltr;
	char STR1[256] = {0};
	char *string;
	ltr = false;
	while (l_buf != NULL)
	{
		if ((l_buf->kind >= (long) toknot && l_buf->kind <= (long) tokrenum) ||
			l_buf->kind == (long) toknum || l_buf->kind == (long) tokvar ||
			l_buf->kind >= (long) toktc)
		{
			if (ltr)
				/*putc(' ', f); */
				output_msg(" ");
			ltr = (bool) (l_buf->kind != toknot);
		}
		else
			ltr = false;
		switch (l_buf->kind)
		{

		case tokvar:
			/*fputs(l_buf->UU.vp->name, f); */
			PhreeqcPtr->output_msg(PhreeqcPtr->sformatf("%s", l_buf->UU.vp->name));
			break;

		case toknum:
			/*fputs(numtostr(STR1, l_buf->UU.num), f); */
			string = numtostr(STR1, l_buf->UU.num);
			PhreeqcPtr->string_trim(string);
			output_msg(PhreeqcPtr->sformatf("%s", string));
			break;

		case tokstr:
			if (strchr(l_buf->UU.sp, '\"'))
			{
				output_msg(PhreeqcPtr->sformatf("\'%s\'", l_buf->UU.sp));
			}
			else
			{
				output_msg(PhreeqcPtr->sformatf("\"%s\"", l_buf->UU.sp));
			}
			break;

		case toksnerr:
			output_msg(PhreeqcPtr->sformatf("{%c}", l_buf->UU.snch));
			break;

		case tokplus:
			/*putc('+', f); */
			output_msg("+");
			break;

		case tokminus:
			/*putc('-', f); */
			output_msg("-");
			break;

		case toktimes:
			/*putc('*', f); */
			output_msg("*");
			break;

		case tokdiv:
			/*putc('/', f); */
			output_msg("/");
			break;

		case tokup:
			/*putc('^', f); */
			output_msg("^");
			break;

		case toklp:
			/*putc('(', f); */
			output_msg("(");
			break;

		case tokrp:
			/*putc(')', f); */
			output_msg(")");
			break;

		case tokcomma:
			/*putc(',', f); */
			output_msg(",");
			break;

		case toksemi:
			/*putc(';', f); */
			output_msg(";");
			break;

		case tokcolon:
			output_msg(" : ");
			break;

		case tokeq:
			output_msg(" = ");
			break;

		case toklt:
			output_msg(" < ");
			break;

		case tokgt:
			output_msg(" > ");
			break;

		case tokle:
			output_msg(" <= ");
			break;

		case tokge:
			output_msg(" >= ");
			break;

		case tokne:
			output_msg(" <> ");
			break;

		case tokand:
			output_msg(" AND ");
			break;

		case tokor:
			output_msg(" OR ");
			break;

		case tokxor:
			output_msg(" XOR ");
			break;

		case tokmod:
			output_msg(" MOD ");
			break;

		case toknot:
			output_msg("NOT ");
			break;

		case toksqr:
			output_msg("SQR");
			break;

		case toksqrt:
			output_msg("SQRT");
			break;

		case toksin:
			output_msg("SIN");
			break;

		case tokcos:
			output_msg("COS");
			break;

		case toktan:
			output_msg("TAN");
			break;

		case tokarctan:
			output_msg("ARCTAN");
			break;

		case toklog:
			output_msg("LOG");
			break;

		case tokexp:
			output_msg("EXP");
			break;

		case tokabs:
			output_msg("ABS");
			break;

		case toksgn:
			output_msg("SGN");
			break;

		case tokstr_:
			output_msg("STR$");
			break;

		case tokval:
			output_msg("VAL");
			break;

		case tokchr_:
			output_msg("CHR$");
			break;

		case tokasc:
			output_msg("ASC");
			break;

		case toklen:
			output_msg("LEN");
			break;

		case tokmid_:
			output_msg("MID$");
			break;

		case tokpeek:
			output_msg("PEEK");
			break;

		case tokrem:
			output_msg(PhreeqcPtr->sformatf("REM%s", l_buf->UU.sp));
			break;

		case toklet:
			output_msg("LET");
			break;

		case tokinput:
			output_msg("INPUT");
			break;

		case tokgoto:
			output_msg("GOTO");
			break;

		case tokif:
			output_msg("IF");
			break;

		case tokend:
			output_msg("END");
			break;

		case tokstop:
			output_msg("STOP");
			break;

		case tokfor:
			output_msg("FOR");
			break;

		case toknext:
			output_msg("NEXT");
			break;

		case tokwhile:
			output_msg("WHILE");
			break;

		case tokwend:
			output_msg("WEND");
			break;

		case tokgosub:
			output_msg("GOSUB");
			break;

		case tokreturn:
			output_msg("RETURN");
			break;

		case tokread:
			output_msg("READ");
			break;

		case tokdata:
			output_msg("DATA");
			break;

		case tokrestore:
			output_msg("RESTORE");
			break;

		case tokgotoxy:
			output_msg("GOTOXY");
			break;

		case tokon:
			output_msg("ON");
			break;

		case tokdim:
			output_msg("DIM");
			break;

		case tokpoke:
			output_msg("POKE");
			break;

		case toklist:
			output_msg("LIST");
			break;

		case tokrun:
			output_msg("RUN");
			break;

		case toknew:
			output_msg("NEW");
			break;

		case tokload:
			output_msg("LOAD");
			break;

		case tokmerge:
			output_msg("MERGE");
			break;

		case tokbye:
			output_msg("BYE");
			break;

		case tokdel:
			output_msg("DEL");
			break;

		case tokrenum:
			output_msg("RENUM");
			break;

		case tokthen:
			output_msg(" THEN ");
			break;

		case tokelse:
			output_msg(" ELSE ");
			break;

		case tokto:
			output_msg(" TO ");
			break;
		case tokstep:
			output_msg(" STEP ");
			break;
/*
* PHREEQC functions
*/
		case tokact:
			output_msg("ACT");
			break;
		case tokadd_heading:
			output_msg("ADD_HEADING");
			break;
		case tokalk:
			output_msg("ALK");
			break;
		case tokaphi:
			output_msg("APHI"); // mole volume of a phase 
			break;
		case tokcalc_value:
			output_msg("CALC_VALUE");
			break;
		case tokceil:
			output_msg("CEIL");
			break;
		case tokcell_no:
			output_msg("CELL_NO");
			break;
		case tokchange_por:
			output_msg("CHANGE_POR");
			break;
		case tokchange_surf:
			output_msg("CHANGE_SURF");
			break;
		case tokcharge_balance:
			output_msg("CHARGE_BALANCE");
			break;
		case tokcurrent_a:
			output_msg("CURRENT_A");
			break;
		case tokdebye_length:
			output_msg("DEBYE_LENGTH"); // Debye-Hueckel length
			break;
		case tokdelta_h_phase:
			output_msg("DELTA_H_PHASE");
			break;
		case tokdelta_h_species:
			output_msg("DELTA_H_SPECIES");
			break;
		case tokdescription:
			output_msg("DESCRIPTION");
			break;
		case tokdh_a:
			output_msg("DH_A"); // Debye-Hueckel A
			break;
		case tokdh_a0:
			output_msg("DH_A0");
			break;
		case tokdh_av:
			output_msg("DH_Av"); // Debye-Hueckel Av
			break;
		case tokdh_b:
			output_msg("DH_B"); // Debye-Hueckel B
			break;
		case tokdh_bdot:
			output_msg("DH_BDOT");
			break;
		case tokdiff_c:
			output_msg("DIFF_C");
			break;
		case tokdist:
			output_msg("DIST");
			break;
		case tokedl:
			output_msg("EDL");
			break;
		case tokedl_species:
			output_msg("EDL_SPECIES");
			break;
		case tokeol_:
			output_msg("EOL$");
			break;
		case tokeol_notab_:
			output_msg("EOL_NOTAB$");
			break;
		case tokeps_r:
			output_msg("EPS_R"); // dielectric constant
			break;
		case tokeq_frac:
		case tokequiv_frac:
			output_msg("EQ_FRAC");
			break;
		case tokequi:
			output_msg("EQUI");
			break;
		case tokequi_delta:
			output_msg("EQUI_DELTA");
			break;
		case tokerase:
			output_msg("ERASE");
			break;
		case tokexists:
			output_msg("EXISTS");
			break;
		case tokfloor:
			output_msg("FLOOR");
			break;
		case tokgamma:
			output_msg("GAMMA");
			break;
		case tokgas:
			output_msg("GAS");
			break;
		case tokgas_p:
			output_msg("GAS_P");
			break;
		case tokgas_vm:
			output_msg("GAS_VM");
			break;
		case tokget:
			output_msg("GET");
			break;
		case tokget_:
			output_msg("GET$");
			break;
		case tokget_por:
			output_msg("GET_POR");
			break;
		case tokgfw:
			output_msg("GFW"); // gram formula weight of a formula
			break;
#if defined MULTICHART
		case tokgraph_x:
			output_msg("GRAPH_X");
			break;
		case tokgraph_y:
			output_msg("GRAPH_Y");
			break;
		case tokgraph_sy:
			output_msg("GRAPH_SY");
			break;
#endif
		case tokinstr:
			output_msg("INSTR");
			break;
		case tokiso:
			output_msg("ISO");
			break;
		case tokiso_unit:
			output_msg("ISO_UNIT");
			break;
		case tokiterations:
			output_msg("ITERATIONS");
			break;
		case tokkappa:
			output_msg("KAPPA"); // compressibility of pure water, d(rho)/d(P) / rho
			break;
		case tokkin:
			output_msg("KIN");
			break;
		case tokkin_delta:
			output_msg("KIN_DELTA");
			break;
		case tokkin_time:
			output_msg("KIN_TIME");
			break;
		case tokkinetics_formula:
		case tokkinetics_formula_:
			output_msg("KINETICS_FORMULA$");
			break;
		case tokla:
			output_msg("LA");
			break;
		case toklg:
			output_msg("LG");
			break;
		case toklist_s_s:
			output_msg("LIST_S_S");
			break;
		case toklk_named:
			output_msg("LK_NAMED");
			break;
		case toklk_phase:
			output_msg("LK_PHASE");
			break;
		case toklk_species:
			output_msg("LK_SPECIES");
			break;
		case toklm:
			output_msg("LM");
			break;
		case tokltrim:
			output_msg("LTRIM");
			break;
		case tokm:
			output_msg("M");
			break;
		case tokm0:
			output_msg("M0");
			break;
		case tokmcd_jtot:
			output_msg("MCD_JTOT");
			break;
		case tokmcd_jconc:
			output_msg("MCD_JCONC");
			break;
		case tokmisc1:
			output_msg("MISC1");
			break;
		case tokmisc2:
			output_msg("MISC2");
			break;
		case tokmol:
			output_msg("MOL");
			break;
		case tokmu:
			output_msg("MU");
			break;
		case tokno_newline_:
			output_msg("NO_NEWLINE$");
			break;
		case tokosmotic:
			output_msg("OSMOTIC");
			break;
		case tokpad_:
		case tokpad:
			output_msg("PAD");
			break;
		case tokparm:
			output_msg("PARM");
			break;
		case tokrate_pk:
			output_msg("RATE_PK");
			break;
		case tokrate_svd:
			output_msg("RATE_SVD");
			break;
		case tokrate_hermanska:
			output_msg("RATE_HERMANSKA");
			break;
		case tokmeang:
			output_msg("MEANG");
			break;
		case tokpercent_error:
			output_msg("PERCENT_ERROR");
			break;
		case tokphase_formula:
		case tokphase_formula_:
			output_msg("PHASE_FORMULA$");
			break;
		case tokphase_vm:
			output_msg("PHASE_VM"); // mole volume of a phase 
			break;
#if defined MULTICHART
		case tokplot_xy:
			output_msg("PLOT_XY");
			break;
#endif
		case tokpot_v:
			output_msg("POT_V");
			break;
		case tokpr_p:
			output_msg("PR_P");
			break;
		case tokpr_phi:
			output_msg("PR_PHI");
			break;
		case tokpressure:
			output_msg("PRESSURE");
			break;
		case tokprint:
			output_msg("PRINT");
			break;
		case tokpunch:
			output_msg("PUNCH");
			break;
		case tokput:
			output_msg("PUT");
			break;
		case tokput_:
			output_msg("PUT$");
			break;
		case tokqbrn:
			output_msg("QBrn"); // Q_Born, d(eps_r)/d(P)/(eps_r^2)
			break;
		case tokrho:
			output_msg("RHO");
			break;
		case tokrho_0:
			output_msg("RHO_0");
			break;
		case tokrtrim:
			output_msg("RTRIM");
			break;
		case tokrxn:
			output_msg("RXN");
			break;
		case toks_s:
			output_msg("S_S");
			break;
		case toksave:
			output_msg("SAVE");
			break;
		case toksc:
			output_msg("SC");
			break;
		case toksetdiff_c:
			output_msg("SETDIFF_C");
			break;
		case toksi:
			output_msg("SI");
		case toksim_no:
			output_msg("SIM_NO");
			break;
		case toksim_time:
			output_msg("SIM_TIME");
			break;
		case toksoln_vol:
			output_msg("SOLN_VOL"); // volume of solution
			break;
		case tokspecies_formula:
		case tokspecies_formula_:
			output_msg("SPECIES_FORMULA$");
			break;
		case tokphase_equation:
		case tokphase_equation_:
			output_msg("PHASE_EQUATION$");
			break;
		case tokspecies_equation:
		case tokspecies_equation_:
			output_msg("SPECIES_EQUATION$");
			break;
		case toksr:
			output_msg("SR");
			break;
		case tokstep_no:
			output_msg("STEP_NO");
			break;
		case tokstr_e_:
			output_msg("STR_E$");
			break;
		case tokstr_f_:
			output_msg("STR_F$");
			break;
		case toksum_gas:
			output_msg("SUM_GAS");
			break;
		case toksum_s_s:
			output_msg("SUM_s_s");
			break;
		case toksum_species:
			output_msg("SUM_SPECIES");
		case toksurf:
			output_msg("SURF");
			break;
		case toksys:
			output_msg("SYS");
			break;
		case tokt_sc:
			output_msg("T_SC");
			break;
		case tokf_visc:
			output_msg("F_VISC");
			break;
		case toktc:
			output_msg("TC");
			break;
		case toktime:
			output_msg("TIME");
			break;
		case toktitle:
			output_msg("TITLE");
			break;
		case toktk:
			output_msg("TK");
			break;
		case toktot:
			output_msg("TOT");
			break;
		case toktotal_time:
			output_msg("TOTAL_TIME");
			break;
		case toktotmole:
		case toktotmol:
		case toktotmoles:
			output_msg("TOTMOLE");
			break;
		case toktrim:
			output_msg("TRIM");
			break;
		case tokviscos:
			output_msg("VISCOS");
			break;
		case tokviscos_0:
			output_msg("VISCOS_0");
			break;
		case tokvm:
			output_msg("VM"); // mole volume of aqueous solute
			break;
/*
* End PHREEQC functions
*/	
		case toksa_declercq: // Undocumented function
			output_msg("SA_DECLERCQ");
			break;
		case tokcallback: // PHAST function
			output_msg("CALLBACK");
			break;
		case tokcell_pore_volume: // PHAST function
		case tokporevolume:
			output_msg("POREVOLUME");
			break;
		case tokcell_porosity: // PHAST function
			output_msg("CELL_POROSITY");
			break;
		case tokcell_saturation: // PHAST function
			output_msg("CELL_SATURATION");
			break;
		case tokcell_volume: // PHAST function
			output_msg("CELL_VOLUME");
			break;
		case toktransport_cell_no: // PHAST function
			output_msg("TRANSPORT_CELL_NO");
			break;
		case tokvelocity_x: // PHAST function
			output_msg("VELOCITY_X");
			break;
		case tokvelocity_y: // PHAST function
			output_msg("VELOCITY_Y");
			break;
		case tokvelocity_z: // PHAST function
			output_msg("VELOCITY_Z");
			break;
		case toklog10:
			output_msg("LOG10");
			break;;
		}
		l_buf = l_buf->next;
	}
}

void PBasic::
disposetokens(tokenrec ** tok)
{
	tokenrec *tok1;

	while (*tok != NULL)
	{
		tok1 = (*tok)->next;
		if (phreeqci_gui)
		{
			if ((*tok)->kind == (long) toknum)
			{
				PhreeqcPtr->PHRQ_free((*tok)->sz_num);
			}
#ifdef _DEBUG
			else
			{
				_ASSERTE((*tok)->sz_num == NULL);
			}
#endif /* _DEBUG */
		}
		if ((*tok)->kind == (long) tokrem || (*tok)->kind == (long) tokstr)
		{
			(*tok)->UU.sp = (char *) PhreeqcPtr->free_check_null((*tok)->UU.sp);
		}
		*tok = (tokenrec *) PhreeqcPtr->free_check_null(*tok);
		*tok = tok1;
	}
}

void PBasic::
parseinput(tokenrec ** l_buf)
{
	linerec *l, *l0, *l1;

	while (PhreeqcPtr->replace("\t", " ", inbuf));
	while (PhreeqcPtr->replace("\r", " ", inbuf));
	PhreeqcPtr->string_trim(inbuf);
	curline = 0;
	while (*inbuf != '\0' && isdigit((int) inbuf[0]))
	{
		curline = curline * 10 + inbuf[0] - 48;
		memmove(inbuf, inbuf + 1, strlen(inbuf));
	}
	parse(inbuf, l_buf);
	if (curline == 0)
		return;
	l = linebase;
	l0 = NULL;
	while (l != NULL && l->num < curline)
	{
		l0 = l;
		l = l->next;
	}
	if (l != NULL && l->num == curline)
	{
		l1 = l;
		l = l->next;
		if (l0 == NULL)
			linebase = l;
		else
			l0->next = l;
		disposetokens(&l1->txt);
		PhreeqcPtr->PHRQ_free(l1);
	}
	if (*l_buf != NULL)
	{
		l1 = (linerec *) PhreeqcPtr->PHRQ_calloc(1, sizeof(linerec));
		if (l1 == NULL)
		{
			PhreeqcPtr->malloc_error();
#if !defined(R_SO)
			exit(4);
#endif
		}
		l1->next = l;
		if (l0 == NULL)
			linebase = l1;
		else
			l0->next = l1;
		l1->num = curline;
		l1->txt = *l_buf;
		strncpy(l1->inbuf, inbuf, MAX_LINE);
		l1->inbuf[MAX_LINE-1] = '\0';
	}
	clearloops();
	restoredata();
}

void PBasic::
errormsg(const char * l_s)
{
	if (phreeqci_gui)
	{
		/* set nIDErrPrompt before calling errormsg see snerr */
		_ASSERTE(nIDErrPrompt != 0);
	}
	else
	{
		error_msg(l_s, CONTINUE);
	}
	_Escape(42);
}

void PBasic::
	snerr(const char * l_s)
{
	char str[MAX_LENGTH] = {0};
	Utilities::strcpy_safe(str, MAX_LENGTH, "Syntax_error ");
	if (phreeqci_gui)
	{
		_ASSERTE(nIDErrPrompt == 0);
		nIDErrPrompt = IDS_ERR_SYNTAX;
	}
	Utilities::strcat_safe(str, MAX_LENGTH, l_s);
	Utilities::strcat_safe(str, MAX_LENGTH, " in line: ");
	if (strcmp(inbuf, "run"))
		Utilities::strcat_safe(str, MAX_LENGTH, inbuf);
	errormsg(str);
}

void PBasic::
	tmerr(const char * l_s)
{
	char str[MAX_LENGTH] = {0};
	Utilities::strcpy_safe(str, MAX_LENGTH, "Type mismatch error");
	if (phreeqci_gui)
	{
		_ASSERTE(nIDErrPrompt == 0);
		nIDErrPrompt = IDS_ERR_MISMATCH;
	}
	Utilities::strcat_safe(str, MAX_LENGTH, l_s);
	Utilities::strcat_safe(str, MAX_LENGTH, " in line: ");
	if (strcmp(inbuf, "run"))
		Utilities::strcat_safe(str, MAX_LENGTH, inbuf);
	errormsg(str);
}

void PBasic::
	badsubscr(void)
{
	if (phreeqci_gui)
	{
		_ASSERTE(nIDErrPrompt == 0);
		nIDErrPrompt = IDS_ERR_BAD_SUBSCRIPT;
	}
	errormsg("Bad subscript");
}

LDBLE PBasic::
realfactor(struct LOC_exec *LINK)
{
	valrec n;
	n = factor(LINK);
	if (n.stringval)
		tmerr(": found characters, not a number");
	return (n.UU.val);
}

char * PBasic::
strfactor(struct LOC_exec * LINK)
{
	valrec n;
	n = factor(LINK);
	if (!n.stringval)
		//tmerr(": chemical name is not enclosed in \"  \"" );
		tmerr(": Expected quoted string or character variable." );
	return (n.UU.sval);
}

char * PBasic::
stringfactor(char * Result, struct LOC_exec * LINK)
{
	valrec n;

	n = factor(LINK);
	if (!n.stringval)
		//tmerr(": chemical name is not enclosed in \"  \"" );
		tmerr(": Expected quoted string or character variable." );
	strcpy(Result, n.UU.sval);
	PhreeqcPtr->PHRQ_free(n.UU.sval);
	return Result;
}

const char * PBasic::
stringfactor(std::string & Result, struct LOC_exec * LINK)
{
	valrec n;

	n = factor(LINK);
	if (!n.stringval)
		//tmerr(": chemical name is not enclosed in \"  \"" );
		tmerr(": Expected quoted string or character variable." );
	Result = n.UU.sval;
	PhreeqcPtr->PHRQ_free(n.UU.sval);
	return Result.c_str();
}

long PBasic::
intfactor(struct LOC_exec *LINK)
{
	return ((long) floor(realfactor(LINK) + 0.5));
}

LDBLE PBasic::
realexpr(struct LOC_exec *LINK)
{
	valrec n;

	n = expr(LINK);
	if (n.stringval)
		tmerr(": found characters, not a number");
	return (n.UU.val);
}

char * PBasic::
strexpr(struct LOC_exec * LINK)
{
	valrec n;

	n = expr(LINK);
	if (!n.stringval)
		//tmerr(": chemical name is not enclosed in \"  \"" );
		tmerr(": Expected quoted string or character variable." );
	return (n.UU.sval);
}

char * PBasic::
stringexpr(char * Result, struct LOC_exec * LINK)
{
	valrec n;

	n = expr(LINK);
	if (!n.stringval)
		//tmerr(": chemical name is not enclosed in \"  \"" );
		tmerr(": Expected quoted string or character variable." );
	strcpy(Result, n.UU.sval);
	PhreeqcPtr->PHRQ_free(n.UU.sval);
	return Result;
}

long PBasic::
intexpr(struct LOC_exec *LINK)
{
	return ((long) floor(realexpr(LINK) + 0.5));
}

void PBasic::
require(int k, struct LOC_exec *LINK)
{
	char str[MAX_LENGTH] = {0};
	if (LINK->t == NULL || LINK->t->kind != k)
	{
		std::map<const std::string, BASIC_TOKEN>::const_iterator item;
		for (item = command_tokens.begin(); item != command_tokens.end(); item++)
		{
			if (item->second == k)
				break;
		}

		if (item == command_tokens.end())
			snerr(": missing unknown command");
		else {
			Utilities::strcpy_safe(str, MAX_LENGTH, ": missing ");
			Utilities::strcat_safe(str, MAX_LENGTH, item->first.c_str());
			snerr(str);
		}
#if !defined(R_SO)
		exit(4);
#endif
	}
	LINK->t = LINK->t->next;
}


void PBasic::
skipparen(struct LOC_exec *LINK)
{
	do
	{
		if (LINK->t == NULL)
		{
			snerr(": parenthesis missing");
#if !defined(R_SO)
			exit(4);
#endif
		}
		if (LINK->t->kind == tokrp || LINK->t->kind == tokcomma)
			goto _L1;
		if (LINK->t->kind == toklp)
		{
			LINK->t = LINK->t->next;
			skipparen(LINK);
		}
		LINK->t = LINK->t->next;
	}
	while (true);
  _L1:;
}

varrec * PBasic::
findvar(struct LOC_exec *LINK)
{
	varrec *v;
	long i, j, k;
	tokenrec *tok;
	long FORLIM;

	if (LINK->t == NULL || LINK->t->kind != tokvar)
	{
		snerr(": can`t find variable");
#if !defined(R_SO)
		exit(4);
#endif
	}
	v = LINK->t->UU.vp;
	LINK->t = LINK->t->next;
	if (LINK->t == NULL || LINK->t->kind != toklp)
	{
		if (v->numdims != 0)
			badsubscr();
		return v;
	}
	if (v->numdims == 0)
	{
		tok = LINK->t;
		i = 0;
		j = 1;
		do
		{
			if (i >= maxdims)
				badsubscr();
			LINK->t = LINK->t->next;
			skipparen(LINK);
			j *= 11;
			i++;
			v->dims[i - 1] = 11;
		}
		while (LINK->t->kind != tokrp);
		v->numdims = (char) i;
		if (v->stringvar)
		{
			v->UU.U1.sarr = (char **) PhreeqcPtr->PHRQ_malloc(j * sizeof(char *));
			if (v->UU.U1.sarr == NULL)
				PhreeqcPtr->malloc_error();
			for (k = 0; k < j; k++)
				v->UU.U1.sarr[k] = NULL;
		}
		else
		{
			v->UU.U0.arr = (LDBLE *) PhreeqcPtr->PHRQ_malloc(j * sizeof(LDBLE));
			if (v->UU.U0.arr == NULL)
				PhreeqcPtr->malloc_error();
			for (k = 0; k < j; k++)
				v->UU.U0.arr[k] = 0.0;
		}
		LINK->t = tok;
	}
	k = 0;
	LINK->t = LINK->t->next;
	FORLIM = v->numdims;
	for (i = 1; i <= FORLIM; i++)
	{
		j = intexpr(LINK);
		if ((unsigned long) j >= (unsigned long) v->dims[i - 1])
			badsubscr();
		k = k * v->dims[i - 1] + j;
		if (i < v->numdims)
			require(tokcomma, LINK);
	}
	require(tokrp, LINK);
	if (v->stringvar)
		v->UU.U1.sval = &v->UU.U1.sarr[k];
	else
		v->UU.U0.val = &v->UU.U0.arr[k];
	return v;
}

valrec PBasic::
factor(struct LOC_exec * LINK)
{
	char string[MAX_LENGTH] = {0};
	cxxSolution *soln_ptr;
	varrec *v;
	tokenrec *facttok;
	valrec n;

	long i, j, m;
	tokenrec *tok, *tok1;
	char *l_s;
	LDBLE l_dummy;
	int i_rate;
	union
	{
		long i;
		char *c;
	} trick;
	LDBLE TEMP;
	std::string STR1, STR2;
	const char *elt_name, *surface_name, *mytemplate, *name;
	varrec *count_varrec = NULL, *names_varrec = NULL, *types_varrec =
		NULL, *moles_varrec = NULL;
	char **names_arg, **types_arg;
	LDBLE *moles_arg;
	int arg_num;
	LDBLE count_species;
	const char *string1, *string2;

	if (LINK->t == NULL)
		snerr(": missing variable or command");
	facttok = LINK->t;
	LINK->t = LINK->t->next;
	n.stringval = false;
	switch (facttok->kind)
	{

	case toknum:
		n.UU.val = facttok->UU.num;
		break;

	case tokstr:
		n.stringval = true;
		m = (int)strlen(facttok->UU.sp) + 1;
		if (m < 256)
			m = 256;
		n.UU.sval = (char*)PhreeqcPtr->PHRQ_calloc(m, sizeof(char));
		if (n.UU.sval == NULL)
			PhreeqcPtr->malloc_error();
		strcpy(n.UU.sval, facttok->UU.sp);
		break;

	case tokvar:
		LINK->t = facttok;
		v = findvar(LINK);
		n.stringval = v->stringvar;
		if (n.stringval)
		{
			if (*v->UU.U1.sval != NULL)
			{
				m = (int)strlen(*v->UU.U1.sval) + 1;
				if (m < 256)
					m = 256;
			}
			else
			{
				m = 256;
			}
			n.UU.sval = (char*)PhreeqcPtr->PHRQ_calloc(m, sizeof(char));
			if (n.UU.sval == NULL)
				PhreeqcPtr->malloc_error();
			if (*v->UU.U1.sval != NULL)
			{
				strcpy(n.UU.sval, *v->UU.U1.sval);
			}

		}
		else
			n.UU.val = *v->UU.U0.val;
		break;

	case toklp:
		n = expr(LINK);
		require(tokrp, LINK);
		break;

	case tokminus:
		n.UU.val = -realfactor(LINK);
		break;

	case tokplus:
	{
		n.UU.val = realfactor(LINK);
	}
	break;

	case toknot:
	{
		n.UU.val = ~intfactor(LINK);
	}
	break;

	case toksqr:
	{
		TEMP = realfactor(LINK);
		n.UU.val = TEMP * TEMP;
	}
	break;

	case toksqrt:
	{
		n.UU.val = sqrt(realfactor(LINK));
	}
	break;
	/*
	* PHREEQC functions=============================================
	*/
	case tokact:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->activity(str);
	}
	break;

	case tokadd_heading:
	{
		require(toklp, LINK);
		name = strexpr(LINK);
		require(tokrp, LINK);
		if (PhreeqcPtr->current_user_punch != NULL)
		{
			PhreeqcPtr->current_user_punch->Get_headings().push_back(name);
			n.UU.val = (parse_all) ? 1 : (double)PhreeqcPtr->current_user_punch->Get_headings().size();
		}
		else {
			n.UU.val = 0;
		}
		name = (char*)PhreeqcPtr->free_check_null((void*) name);
	}
	break;

	case tokalk:
	{
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->total_alkalinity / PhreeqcPtr->mass_water_aq_x;
	}
	break;

	case tokaphi:
	{
		n.UU.val = PhreeqcPtr->A0;
	}
	break;

	case tokcalc_value:
	{
		require(toklp, LINK);
		name = stringfactor(STR1, LINK);
		require(tokrp, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->get_calculate_value(name);
	}
	break;

	case tokceil:
	{
		n.UU.val = ceil(realfactor(LINK));
	}
	break;

	case tokcell_no:
	{
		if (parse_all)
		{
			n.UU.val = 1;
			break;
		}
		n.UU.val = PhreeqcPtr->solution_number();
	}
	break;

	case tokcharge_balance:
	{
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->cb_x;
	}
	break;

	case tokcurrent_a:
	{
		//n.UU.val = (parse_all) ? 1 : PhreeqcPtr->current_x;
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->current_A;
	}
	break;

	case tokdebye_length:
	{
		double debye_length = (PhreeqcPtr->eps_r * EPSILON_ZERO * R_KJ_DEG_MOL * 1000.0 * PhreeqcPtr->tk_x)
			/ (2. * F_C_MOL * F_C_MOL * PhreeqcPtr->mu_x * 1000.);
		n.UU.val = sqrt(debye_length);
	}
	break;

	case tokdelta_h_phase:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->calc_deltah_p(str);
	}
	break;

	case tokdelta_h_species:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->calc_deltah_s(str);
	}
	break;

	case tokdescription:
	{
		n.stringval = true;
		if (PhreeqcPtr->state == REACTION)
		{
			if (PhreeqcPtr->use.Get_mix_in())
			{
				snprintf(string, sizeof(string), "Mix %d", PhreeqcPtr->use.Get_n_mix_user());
				n.UU.sval = PhreeqcPtr->string_duplicate(string);
			}
			else
			{
				soln_ptr = Utilities::Rxn_find(PhreeqcPtr->Rxn_solution_map,
					PhreeqcPtr->use.Get_n_solution_user());
				if (soln_ptr != NULL)
				{
					n.UU.sval = PhreeqcPtr->string_duplicate(soln_ptr->Get_description().c_str());
				}
				else
				{
					n.UU.sval = PhreeqcPtr->string_duplicate("Unknown");
				}
			}
		}
		else if (PhreeqcPtr->state == ADVECTION || PhreeqcPtr->state == TRANSPORT || PhreeqcPtr->state == PHAST)
		{
			snprintf(string, sizeof(string), "Cell %d", PhreeqcPtr->cell_no);
			n.UU.sval = PhreeqcPtr->string_duplicate(string);
		}
		else
		{
			if (PhreeqcPtr->use.Get_solution_ptr() != NULL)
			{
				n.UU.sval = PhreeqcPtr->string_duplicate(PhreeqcPtr->use.Get_solution_ptr()->Get_description().c_str());
			}
			else
			{
				n.UU.sval = PhreeqcPtr->string_duplicate("Unknown");
			}
		}
		while (PhreeqcPtr->replace("\t", " ", n.UU.sval));
	}
	break;

	case tokdh_a:
	{
		if (PhreeqcPtr->llnl_temp.size() > 0)
		{
			n.UU.val = PhreeqcPtr->a_llnl;
		}
		else
		{
			n.UU.val = PhreeqcPtr->DH_A;
		}
	}
	break;

	case tokdh_a0:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->dh_a0(str);
	}
	break;

	case tokdh_av:
	{
		n.UU.val = PhreeqcPtr->DH_Av;
	}
	break;

	case tokdh_b:
	{
		if (PhreeqcPtr->llnl_temp.size() > 0)
		{
			n.UU.val = PhreeqcPtr->b_llnl;
		}
		else
		{
			n.UU.val = PhreeqcPtr->DH_B;
		}
	}
	break;
	
	case tokdh_bdot:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->dh_bdot(str);
	}
	break;

	case tokdiff_c:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->diff_c(str);
	}
	break;

	case tokdist:
	{
		if (PhreeqcPtr->state == PHAST)
		{
			n.UU.val = 0;
		}
		else if (PhreeqcPtr->state == TRANSPORT)
		{
			n.UU.val = (parse_all) ? 1 : PhreeqcPtr->cell_data[PhreeqcPtr->cell].mid_cell_x;
		}
		else if (PhreeqcPtr->state == ADVECTION)
		{
			n.UU.val = (parse_all) ? 1 : (LDBLE)PhreeqcPtr->use.Get_n_solution_user();
		}
		else
		{
			n.UU.val = 0;
		}
	}
	break;

	case tokedl:
	{
		require(toklp, LINK);
		elt_name = stringfactor(STR1, LINK);
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			surface_name = stringfactor(STR2, LINK);
		}
		else
		{
			surface_name = NULL;
		}
		require(tokrp, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->diff_layer_total(elt_name, surface_name);
	}
	break;

	case tokedl_species:
	{
		double area = 0.0, thickness = 0.0;
		require(toklp, LINK);
		const char* surf_name = stringfactor(STR1, LINK);
		require(tokcomma, LINK);
		// variable for number of species
		count_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || !count_varrec || count_varrec->stringvar != 0)
		{
			snerr(": Missing or wrong type count variable.");
#if !defined(R_SO)
			exit(4);
#endif
		}
		// variable for species names
		LINK->t = LINK->t->next;
		require(tokcomma, LINK);
		names_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || !names_varrec || names_varrec->stringvar != 1)
		{
			snerr(": Missing or wrong type name variable.");
#if !defined(R_SO)
			exit(4);
#endif
		}
		// variable for species concentrations
		LINK->t = LINK->t->next;
		require(tokcomma, LINK);
		moles_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || moles_varrec->stringvar != 0)
			snerr(": Missing or wrong type moles variable.");
		// variable for area
		LINK->t = LINK->t->next;
		require(tokcomma, LINK);
		varrec* area_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || area_varrec->stringvar != 0)
			snerr(": Missing or wrong type area variable.");
		// variable for thickness
		LINK->t = LINK->t->next;
		require(tokcomma, LINK);
		varrec* thickness_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || thickness_varrec->stringvar != 0)
			snerr(": Missing or wrong type thickness variable.");
		LINK->t = LINK->t->next;
		require(tokrp, LINK);

		free_dim_stringvar(names_varrec);
		PhreeqcPtr->free_check_null(moles_varrec->UU.U0.arr);
		moles_varrec->UU.U0.arr = NULL;

		// Call subroutine
		if (parse_all)
		{
			PhreeqcPtr->sys_tot = 0;
			//PhreeqcPtr->count_sys = 1000;
			//int count_sys = PhreeqcPtr->count_sys;
			size_t count_sys = 1000;
			names_arg = (char**)PhreeqcPtr->PHRQ_calloc((count_sys + 1), sizeof(char*));
			if (names_arg == NULL)
			{
				PhreeqcPtr->malloc_error();
#if !defined(R_SO)
				exit(4);
#endif
			}
			moles_arg = (LDBLE*)PhreeqcPtr->PHRQ_calloc((count_sys + 1), sizeof(LDBLE));
			if (moles_arg == NULL)
			{
				PhreeqcPtr->malloc_error();
#if !defined(R_SO)
				exit(4);
#endif
			}
			names_arg[0] = NULL;
			moles_arg[0] = 0;
			count_species = (LDBLE)count_sys;
			n.UU.val = 0;
		}
		else
		{
			//n.UU.val = PhreeqcPtr->system_total(elt_name, &count_species, &(names_arg),
			//	&(types_arg), &(moles_arg));
			n.UU.val = PhreeqcPtr->edl_species(surf_name, &count_species, &(names_arg), &(moles_arg), &area, &thickness);
		}
		/*
		*  fill in varrec structures
		*/
		*count_varrec->UU.U0.val = count_species;
		names_varrec->UU.U1.sarr = names_arg;
		moles_varrec->UU.U0.arr = moles_arg;
		*area_varrec->UU.U0.val = area;
		*thickness_varrec->UU.U0.val = thickness;

		for (i = 0; i < maxdims; i++)
		{
			names_varrec->dims[i] = 0;
			moles_varrec->dims[i] = 0;
		}
		names_varrec->dims[0] = (long)(*count_varrec->UU.U0.val) + 1;
		moles_varrec->dims[0] = (long)(*count_varrec->UU.U0.val) + 1;
		names_varrec->numdims = 1;
		moles_varrec->numdims = 1;
	}
	break;

	case tokeol_:
	{
		n.stringval = true;
		n.UU.sval = (char*)PhreeqcPtr->PHRQ_calloc(256, sizeof(char));
		if (n.UU.sval == NULL)
			PhreeqcPtr->malloc_error();
		strcpy(n.UU.sval, "\n");
	}
	break;

	case tokeol_notab_:
	{
		n.stringval = true;
		n.UU.sval = (char*)PhreeqcPtr->PHRQ_calloc(256, sizeof(char));
		if (n.UU.sval == NULL)
			PhreeqcPtr->malloc_error();
		strcpy(n.UU.sval, "\n");
		punch_tab = false;
	}
	break;

	case tokeps_r:
	{
		n.UU.val = PhreeqcPtr->eps_r;
	}
	break;

	case tokeq_frac:
	case tokequiv_frac:
	{
		// left parenthesis
		require(toklp, LINK);

		// species name
		std::string species_name(stringfactor(STR1, LINK));

		require(tokcomma, LINK);

		// equivalents
		count_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || count_varrec->stringvar != 0)
			snerr(": Cannot find equivalents variable");

		LINK->t = LINK->t->next;
		require(tokcomma, LINK);

		// exchange or surface element 
		varrec* elt_varrec = NULL;
		elt_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || elt_varrec->stringvar != 1)
			snerr(": Cannot find element string variable");
		free_dim_stringvar(elt_varrec);
		*elt_varrec->UU.U1.sval = (char*)PhreeqcPtr->free_check_null(*elt_varrec->UU.U1.sval);

		// right parenthesis
		LINK->t = LINK->t->next;
		require(tokrp, LINK);

		// set function value
		LDBLE eq;
		std::string elt_name;

		// return equivalent fraction
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->equivalent_fraction(species_name.c_str(), &eq, elt_name);

		// set equivalents
		*count_varrec->UU.U0.val = (parse_all) ? 1 : eq;

		// set element name
		size_t l = elt_name.size();
		l = l < 256 ? 256 : l + 1;
		char* token = (char*)PhreeqcPtr->PHRQ_malloc(l * sizeof(char));
		 Utilities::strcpy_safe(token, l, elt_name.c_str());
		*elt_varrec->UU.U1.sval = token;
	}
	break;

	case tokequi:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->equi_phase(str);
	}
	break;

	case tokequi_delta:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->equi_phase_delta(str);
	}
	break;

	case tokexists:
	{
		std::ostringstream oss;
		require(toklp, LINK);

		/* get first subscript */
		if (LINK->t != NULL && LINK->t->kind != tokrp)
		{
			i = intexpr(LINK);
			oss << i << ",";
		}

		/* get other subscripts */
		for (;;)
		{
			if (LINK->t != NULL && LINK->t->kind == tokcomma)
			{
				LINK->t = LINK->t->next;
				j = intexpr(LINK);
				oss << j << ",";
			}
			else
			{
				/* get right parentheses */
				require(tokrp, LINK);
				break;
			}
		}
		if (parse_all)
		{
			n.UU.val = 1;
		}
		else
		{
			std::map<std::string, double>::iterator it = PhreeqcPtr->save_values.find(oss.str());
			n.UU.val = (it == PhreeqcPtr->save_values.end()) ? 0 : 1;
		}
	}
	break;

	case tokfloor:
	{
		n.UU.val = floor(realfactor(LINK));
	}
	break;

	case tokmcd_jtot:
	{
		double f = 0.0;
		const char* str = stringfactor(STR1, LINK);
		if (PhreeqcPtr->state == TRANSPORT && PhreeqcPtr->multi_Dflag)
		{
			f = PhreeqcPtr->flux_mcd(str, 1);
		}
		n.UU.val = (parse_all) ? 1 : f;
	}
	break;
	case tokmcd_jconc:
	{	
		double f = 0.0;
		const char* str = stringfactor(STR1, LINK);
		if (PhreeqcPtr->state == TRANSPORT && PhreeqcPtr->multi_Dflag)
		{
			f = PhreeqcPtr->flux_mcd(str, 2);
		}
		n.UU.val = (parse_all) ? 1 : f;
	}
	break;

	case tokgamma:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->activity_coefficient(str);
	}
	break;

	case tokgas:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->find_gas_comp(str);
	}
	break;

	case tokgas_p:
	{
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->find_gas_p();
	}
	break;

	case tokgas_vm:
	{
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->find_gas_vm();
	}
	break;

	case tokget_:
	{
		std::ostringstream oss;
		require(toklp, LINK);

		/* get first subscript */
		if (LINK->t != NULL && LINK->t->kind != tokrp)
		{
			i = intexpr(LINK);
			oss << i << ",";
		}

		/* get other subscripts */
		for (;;)
		{
			if (LINK->t != NULL && LINK->t->kind == tokcomma)
			{
				LINK->t = LINK->t->next;
				j = intexpr(LINK);
				oss << j << ",";
			}
			else
			{
				/* get right parentheses */
				require(tokrp, LINK);
				break;
			}
		}
		if (parse_all)
		{
			n.UU.val = 1;
		}
		else
		{
			n.stringval = true;
			n.UU.sval = (char*)PhreeqcPtr->PHRQ_calloc(256, sizeof(char));
			if (n.UU.sval == NULL)
				PhreeqcPtr->malloc_error();
			std::map<std::string, std::string>::iterator it = PhreeqcPtr->save_strings.find(oss.str());
			n.UU.sval = (it == PhreeqcPtr->save_strings.end()) ? strcpy(n.UU.sval, "unknown") :
				strcpy(n.UU.sval, it->second.c_str());
		}
		break;
	}

	case tokget:
	{
		std::ostringstream oss;
		require(toklp, LINK);

		/* get first subscript */
		if (LINK->t != NULL && LINK->t->kind != tokrp)
		{
			i = intexpr(LINK);
			oss << i << ",";
		}

		/* get other subscripts */
		for (;;)
		{
			if (LINK->t != NULL && LINK->t->kind == tokcomma)
			{
				LINK->t = LINK->t->next;
				j = intexpr(LINK);
				oss << j << ",";
			}
			else
			{
				/* get right parentheses */
				require(tokrp, LINK);
				break;
			}
		}
		if (parse_all)
		{
			n.UU.val = 1;
		}
		else
		{
			std::map<std::string, double>::iterator it = PhreeqcPtr->save_values.find(oss.str());
			n.UU.val = (it == PhreeqcPtr->save_values.end()) ? 0 : it->second;
		}
		break;
	}
	case tokget_por:
	{
		i = intfactor(LINK);
		if (parse_all)
		{
			n.UU.val = 1;
		}
		else
		{
			if (PhreeqcPtr->phast != TRUE)
			{
				if (i <= 0 || i > PhreeqcPtr->count_cells * (1 + PhreeqcPtr->stag_data.count_stag) + 1
					|| i == PhreeqcPtr->count_cells + 1)
				{
					/*		warning_msg("Note... no porosity for boundary solutions."); */
					n.UU.val = 0;
					break;
				}
				else
					n.UU.val = PhreeqcPtr->cell_data[i].por;
				break;
			}
			else
			{
				n.UU.val = PhreeqcPtr->cell_porosity;
				break;
			}
		}
	}
	break;

	case tokgfw:
	{
		const char* str = stringfactor(STR1, LINK);
		LDBLE gfw;
		PhreeqcPtr->compute_gfw(str, &gfw);
		n.UU.val = (parse_all) ? 1 : gfw;
	}
	break;

	case tokinstr:
	{
		require(toklp, LINK);
		string1 = stringfactor(STR1, LINK);
		require(tokcomma, LINK);
		string2 = stringfactor(STR2, LINK);
		require(tokrp, LINK);
		{
			const char* cptr = strstr(string1, string2);
			if (cptr == NULL)
			{
				n.UU.val = 0;
			}
			else
			{
				n.UU.val = ((LDBLE)(cptr - string1)) + 1;
			}
		}
	}
	break;

	case tokiso:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->iso_value(str);
	}
	break;

	case tokiso_unit:
	{
		n.stringval = true;
		require(toklp, LINK);
		string1 = stringfactor(STR1, LINK);
		require(tokrp, LINK);
		trim(STR1);
		n.UU.sval = (parse_all) ? PhreeqcPtr->string_duplicate("unknown") : PhreeqcPtr->iso_unit(STR1.c_str());
	}
	break;

	case tokiterations:
	{
		n.UU.val = (parse_all) ? 0 : PhreeqcPtr->overall_iterations;
	}
	break;

	case tokkappa:
	{
		n.UU.val = PhreeqcPtr->kappa_0;
	}
	break;

	case tokkin:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->kinetics_moles(str);
	}
	break;

	case tokkin_delta:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->kinetics_moles_delta(str);
	}
	break;

	case tokkin_time:
	{
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->rate_kin_time;
	}
	break;

	case tokkinetics_formula:
	case tokkinetics_formula_:
	{
		require(toklp, LINK);
		std::string kinetics_name(stringfactor(STR1, LINK));
		varrec* elts_varrec = NULL, * coef_varrec = NULL;
		cxxNameDouble stoichiometry;
		/*
		*  Parse arguments
		*/
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			/* kinetics_formula("calcite", count, elt, coef) */
			/* return formula */
			/*int c; */
			/*  struct varrec *count_varrec, *names_varrec, *types_varrec, *moles_varrec; */
			/*  struct varrec *count_varrec, *elt_varrec, *coef_varrec; */
			/* return number of species */
			LINK->t = LINK->t->next;
			count_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || count_varrec->stringvar != 0)
				snerr(": Cannot find count variable");

			/* return number of names of elements */
			LINK->t = LINK->t->next;
			require(tokcomma, LINK);
			elts_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || elts_varrec->stringvar != 1)
				snerr(": Cannot find element string variable");

			/* return coefficients of species */
			LINK->t = LINK->t->next;
			require(tokcomma, LINK);
			coef_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || coef_varrec->stringvar != 0)
				snerr(": Cannot find coefficient variable");
			LINK->t = LINK->t->next;
			arg_num = 4;
		}
		else
		{
			arg_num = 1;
		}
		require(tokrp, LINK);

		if (arg_num > 1)
		{
			free_dim_stringvar(elts_varrec);
			PhreeqcPtr->free_check_null(coef_varrec->UU.U0.arr);
			coef_varrec->UU.U0.arr = NULL;
		}
		/*
		*  Call subroutine
		*/
		std::string form = PhreeqcPtr->kinetics_formula(kinetics_name, stoichiometry);

		// put formula as return value
		n.stringval = true;
		n.UU.sval = PhreeqcPtr->string_duplicate(form.c_str());

		/*
		*  fill in varrec structure
		*/

		if (arg_num > 1)
		{
			size_t count = stoichiometry.size();
			*count_varrec->UU.U0.val = (LDBLE)count;
			/*
			* malloc space
			*/
			elts_varrec->UU.U1.sarr = (char**)PhreeqcPtr->PHRQ_malloc((count + 1) * sizeof(char*));
			if (elts_varrec->UU.U1.sarr == NULL)
				PhreeqcPtr->malloc_error();
			coef_varrec->UU.U0.arr = (LDBLE*)PhreeqcPtr->PHRQ_malloc((count + 1) * sizeof(LDBLE));
			if (coef_varrec->UU.U0.arr == NULL)
				PhreeqcPtr->malloc_error();

			// first position not used
			elts_varrec->UU.U1.sarr[0] = NULL;
			coef_varrec->UU.U0.arr[0] = 0;

			// set dims for Basic array
			for (i = 0; i < maxdims; i++)
			{
				elts_varrec->dims[i] = 0;
				coef_varrec->dims[i] = 0;
			}
			// set dims for first dimension and number of dims
			elts_varrec->dims[0] = (long)(count + 1);
			coef_varrec->dims[0] = (long)(count + 1);
			elts_varrec->numdims = 1;
			coef_varrec->numdims = 1;

			// fill in arrays
			i = 1;
			for (cxxNameDouble::iterator it = stoichiometry.begin(); it != stoichiometry.end(); it++)
			{
				elts_varrec->UU.U1.sarr[i] = PhreeqcPtr->string_duplicate((it->first).c_str());
				coef_varrec->UU.U0.arr[i] = it->second;
				i++;
			}
		}
	}
	break;

	case tokla:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->log_activity(str);
	}
	break;

	case toklg:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->log_activity_coefficient(str);
	}
	break;

	case toklist_s_s:
	{
		/* list_s_s("calcite", count, name$, moles) */
		/* return total moles */
		require(toklp, LINK);
		std::string s_s_name(stringfactor(STR1, LINK));
		cxxNameDouble composition;
		/*
		*  Parse arguments
		*/
		arg_num = -1;
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			count_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || count_varrec->stringvar != 0)
				snerr(": Cannot find count variable");

			/* return number of names of components */
			LINK->t = LINK->t->next;
			require(tokcomma, LINK);
			names_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || names_varrec->stringvar != 1)
				snerr(": Cannot find component string variable");

			/* return number of moles  of components */
			LINK->t = LINK->t->next;
			require(tokcomma, LINK);
			moles_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || moles_varrec->stringvar != 0)
				snerr(": Cannot find moles of component variable");
			LINK->t = LINK->t->next;
			arg_num = 4;
		}
		else
		{
			snerr(": Expected 4 arguments for list_s_s");
#if !defined(R_SO)
			exit(4);
#endif
		}
		require(tokrp, LINK);

		if (arg_num > 1)
		{
			free_dim_stringvar(names_varrec);
			if (moles_varrec)
			{
				PhreeqcPtr->free_check_null(moles_varrec->UU.U0.arr);
				moles_varrec->UU.U0.arr = NULL;
			}
		}
		/*
		*  Call subroutine
		*/
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->list_ss(s_s_name, composition);

		/*
		*  fill in varrec structure
		*/

		if (arg_num > 1)
		{
			size_t count = composition.size();
			*count_varrec->UU.U0.val = (LDBLE)count;
			/*
			* malloc space
			*/
			names_varrec->UU.U1.sarr = (char**)PhreeqcPtr->PHRQ_malloc((count + 1) * sizeof(char*));
			if (names_varrec->UU.U1.sarr == NULL)
			{
				PhreeqcPtr->malloc_error();
#if !defined(R_SO)
				exit(4);
#endif
			}
			moles_varrec->UU.U0.arr = (LDBLE*)PhreeqcPtr->PHRQ_malloc((count + 1) * sizeof(LDBLE));
			if (moles_varrec->UU.U0.arr == NULL)
			{
				PhreeqcPtr->malloc_error();
#if !defined(R_SO)
				exit(4);
#endif
			}

			// first position not used
			names_varrec->UU.U1.sarr[0] = NULL;
			moles_varrec->UU.U0.arr[0] = 0;

			// set dims for Basic array
			for (i = 0; i < maxdims; i++)
			{
				names_varrec->dims[i] = 0;
				moles_varrec->dims[i] = 0;
			}
			// set dims for first dimension and number of dims
			names_varrec->dims[0] = (long)(count + 1);
			moles_varrec->dims[0] = (long)(count + 1);
			names_varrec->numdims = 1;
			moles_varrec->numdims = 1;

			// fill in arrays
			i = 1;
			std::vector< std::pair<std::string, LDBLE> > sort_comp = composition.sort_second();
			size_t j;
			for (j = 0; j != sort_comp.size(); j++)
			{
				names_varrec->UU.U1.sarr[i] = PhreeqcPtr->string_duplicate(sort_comp[j].first.c_str());
				moles_varrec->UU.U0.arr[i] = sort_comp[j].second;
				i++;
			}

		}
	}
	break;

	case toklk_named:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->calc_logk_n(str);
	}
	break;

	case toklk_phase:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->calc_logk_p(str);
	}
	break;

	case toklk_species:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->calc_logk_s(str);
	}
	break;

	case toklm:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->log_molality(str);
	}
	break;

	case tokltrim:
	{
		n.stringval = true;
		require(toklp, LINK);
		string1 = stringfactor(STR1, LINK);
		require(tokrp, LINK);
		trim_left(STR1);
		n.UU.sval = PhreeqcPtr->string_duplicate(STR1.c_str());
	}
	break;

	case tokm:
	{
		n.UU.val = PhreeqcPtr->rate_m;
	}
	break;

	case tokm0:
	{
		n.UU.val = PhreeqcPtr->rate_m0;
	}
	break;

	case tokmisc1:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->find_misc1(str);
	}
	break;

	case tokmisc2:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->find_misc2(str);
	}
	break;

	case tokmol:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->molality(str);
	}
	break;

	case tokmu:
	{
		n.UU.val = PhreeqcPtr->mu_x;
	}
	break;

	case tokno_newline_:
	{
		n.stringval = true;
		PhreeqcPtr->Set_output_newline(false);
		this->skip_punch = true;
	}
	break;

	case tokosmotic:
	{
		if (PhreeqcPtr->pitzer_model == TRUE || PhreeqcPtr->sit_model == TRUE)
		{
			n.UU.val = PhreeqcPtr->COSMOT;
		}
		else
		{
			n.UU.val = 0.0;
		}
	}
	break;

	case tokpad_:
	case tokpad:
	{
		char* str;
		n.stringval = true;
		require(toklp, LINK);
		str = strexpr(LINK);
		require(tokcomma, LINK);
		i = intexpr(LINK);
		require(tokrp, LINK);
		n.UU.sval = PhreeqcPtr->string_pad(str, i);
		PhreeqcPtr->PHRQ_free(str);
	}
	break;

	case tokparm:
	{
		i_rate = intfactor(LINK);
		if (parse_all)
		{
			n.UU.val = 1;
		}
		else
		{
			if (i_rate > PhreeqcPtr->count_rate_p || i_rate == 0)
			{
				errormsg("Parameter subscript out of range.");
			}
			n.UU.val = PhreeqcPtr->rate_p[(size_t)i_rate - 1];
		}
	}
	break;

	case tokrate_pk:
	{
		require(toklp, LINK);
		char* min_name = strexpr(LINK);
		require(tokrp, LINK);
		if (parse_all) {
			PhreeqcPtr->PHRQ_free(min_name);
			n.UU.val = 1;
			break;
		}
		std::string min_string = min_name;
		PhreeqcPtr->PHRQ_free(min_name);
		Utilities::str_tolower(min_string);
		std::map<std::string, std::vector<double> >::const_iterator it = PhreeqcPtr->rate_parameters_pk.find(min_string);
		if (it == PhreeqcPtr->rate_parameters_pk.end())
		{
			std::ostringstream oss;
			oss << "PK rate parameters not found for " << min_name << "\n";
			snerr(oss.str().c_str());
		}
		//if (it->second.size() != 8)
		//{
		//	std::ostringstream oss;
		//	oss << "RATE_PK requires 8 rate parameters, " << it->second.size() << " were found for " << min_name << "\n";
		//	snerr(oss.str().c_str());
		//}
		// temperature factor, gas constant
		double dif_temp = 1.0 / PhreeqcPtr->tk_x - 1.0 / 298.15;
		double dT_R = dif_temp / (2.303 * 8.314e-3);
		int Table = 0;
		double rate_H = 0.0, rate_H2O = 0.0, rate_OH = 0.0;
		double lgk_H = -30.0, lgk_H2O = -30.0, lgk_OH = -30.0;
		if (it->second.size() > 8)
			Table = (int) it->second.back();

		switch (Table)
		{
		case 0:
			if (it->second.size() != 8)
			{
				std::ostringstream oss;
				oss << "Expected 8 rate parameters, " << it->second.size() << " were found for " << min_name << "\n";
				snerr(oss.str().c_str());
			}
			break;
		case 33:
			if (it->second.size() != 9)
			{
				std::ostringstream oss;
				oss << "Expected 8 rate parameters for table 33 mineral. " << it->second.size() - 1 << " were found for " << min_name << ".\n";
				snerr(oss.str().c_str());
			}
			break;
		case 35:
			if (it->second.size() != 11)
			{
				std::ostringstream oss;
				oss << "Expected 10 rate parameters for table 35 mineral. " << it->second.size() - 1 << " were found for " << min_name << ".\n";
				snerr(oss.str().c_str());
			}
			break;
		default:
		{
			std::ostringstream oss;
			oss << "Unknown table value " << Table << " for " << min_name << ".";
			snerr(oss.str().c_str());
		}
		break;
		}
		switch (Table)
		{
		case 0:
			// rate by H+
			if ((lgk_H = it->second[0]) > -30)
			{
				double e_H = it->second[1];
				double nH = it->second[2];
				rate_H = pow(10.0, lgk_H - e_H * dT_R) * pow(PhreeqcPtr->activity("H+"), nH);
			}
			// rate by hydrolysis
			if ((lgk_H2O = it->second[3]) > -30)
			{
				double e_H2O = it->second[4];
				rate_H2O = pow(10.0, lgk_H2O - e_H2O * dT_R);
			}
			//	 rate by OH-
			if ((lgk_OH = it->second[5]) > -30)
			{
				double e_OH = it->second[6];
				double n_OH = it->second[7];
				rate_OH = pow(10.0, lgk_OH - e_OH * dT_R) * pow(PhreeqcPtr->activity("H+"), n_OH);
			}
			break;
		case 33:
			// rate by H+
			if ((lgk_H = it->second[0]) > -30)
			{
				double e_H = it->second[1];
				double nH = it->second[2];
				rate_H = pow(10.0, lgk_H - e_H * dT_R) * pow(PhreeqcPtr->activity("H+"), nH);
			}
			// rate by hydrolysis
			if ((lgk_H2O = it->second[3]) > -30)
			{
				double e_H2O = it->second[4];
				rate_H2O = pow(10.0, lgk_H2O - e_H2O * dT_R);
			}
			//	 rate by P_CO2
			if ((lgk_OH = it->second[5]) > -30)
			{
				double e_OH = it->second[6];
				double n_PCO2 = it->second[7];
				rate_OH = pow(10.0, lgk_OH - e_OH * dT_R) * pow(PhreeqcPtr->saturation_ratio("CO2(g)"), n_PCO2);
			}
			break;
		case 35:
			// rate by H+ and Fe+3
			if ((lgk_H = it->second[0]) > -30)
			{
				double e_H = it->second[1];
				double nH = it->second[2];
				double nFe = it->second[3];
				rate_H = pow(10.0, lgk_H - e_H * dT_R) * pow(PhreeqcPtr->activity("H+"), nH) * pow(PhreeqcPtr->activity("Fe+3"), nFe);
			}
			// rate by hydrolysis and O2
			if ((lgk_H2O = it->second[4]) > -30)
			{
				double e_H2O = it->second[5];
				double n_O2 = it->second[6];
				rate_H2O = pow(10.0, lgk_H2O - e_H2O * dT_R) * pow(PhreeqcPtr->activity("O2"), n_O2);
			}
			//	 rate by OH-
			if ((lgk_OH = it->second[7]) > -30)
			{
				double e_OH = it->second[8];
				double n_OH = it->second[9];
				rate_OH = pow(10.0, lgk_OH - e_OH * dT_R) * pow(PhreeqcPtr->activity("H+"), n_OH);
			}
			break;
		}
		// sum rates
		double rate = rate_H + rate_H2O + rate_OH;
		n.UU.val = rate;
		//	# affinity_factor m ^ 2 / mol roughness, lgkH  e_H  nH, lgkH2O e_H2O, lgkOH e_OH nOH
		//		# parm number  1       2          3, 4    5   6, 7     8, 9     10  11
		//		10  affinity = get(-99, 1) # retrieve number from memory
		//		20
		//		30  REM # specific area m2 / mol, surface roughness
		//		40  sp_area = get(-99, 2) : roughness = get(-99, 3)
		//		50
		//		60  REM # temperature factor, gas constant
		//		70  dif_temp = 1 / TK - 1 / 298 : R = 2.303 * 8.314e-3 : dT_R = dif_temp / R
		//		80
		//		90  REM # rate by H +
		//		100 lgk_H = get(-99, 4) : e_H = get(-99, 5) : nH = get(-99, 6)
		//		110 rate_H = 10 ^ (lgk_H - e_H * dT_R) * ACT("H+") ^ nH
		//		120
		//		130 REM # rate by hydrolysis
		//		140 lgk_H2O = get(-99, 7) : e_H2O = get(-99, 8)
		//		150 rate_H2O = 10 ^ (lgk_H2O - e_H2O * dT_R)
		//		160
		//		170 REM # rate by OH -
		//		180 lgk_OH = get(-99, 9) : e_OH = get(-99, 10) : nOH = get(-99, 11)
		//		190 rate_OH = 10 ^ (lgk_OH - e_OH * dT_R) * ACT("H+") ^ nOH
		//		200
		//		210 rate = rate_H + rate_H2O + rate_OH
		//		220 area = sp_area * M0 * (M / M0) ^ 0.67
		//		230
		//		240 rate = area * roughness * rate * affinity
		//		250 SAVE rate * TIME
		//		-end
	}
	break;
	case tokrate_svd:
	{
		require(toklp, LINK);
		char* min_name = strexpr(LINK);
		require(tokrp, LINK);
		if (parse_all) {
			PhreeqcPtr->PHRQ_free(min_name);
			n.UU.val = 1;
			break;
		}
		std::string min_string = min_name;
		PhreeqcPtr->PHRQ_free(min_name);
		Utilities::str_tolower(min_string);
		std::map<std::string, std::vector<double> >::const_iterator it = PhreeqcPtr->rate_parameters_svd.find(min_string);
		if (it == PhreeqcPtr->rate_parameters_svd.end())
		{
			std::ostringstream oss;
			oss << "SVD rate parameters not found for " << min_name << "\n";
			snerr(oss.str().c_str());
		}
		if (it->second.size() != 31)
		{
			std::ostringstream oss;
			oss << "RATE_SVD requires 31 rate parameters, " << it->second.size() << " were found for " << min_name << "\n";
			snerr(oss.str().c_str());
		}

		// temperature factor, gas constant
		double dif_temp = 1.0 / PhreeqcPtr->tk_x - 1.0 / 281.0;
		double e_H = it->second[0];
		double e_H2O = it->second[1];
		double e_CO2 = it->second[2];
		double e_OA = it->second[3];
		double e_OH = it->second[4];

		double BC = PhreeqcPtr->activity("Na+") + PhreeqcPtr->activity("K+") +
			PhreeqcPtr->activity("Mg+2") + PhreeqcPtr->activity("Ca+2");
		double aAl = PhreeqcPtr->activity("Al+3");
		double aSi = PhreeqcPtr->activity("H4SiO4") + PhreeqcPtr->activity("SiO2");
		double R = PhreeqcPtr->total("Organicmatter");
		//	rate by H +
		double pkH = it->second[5];
		double nH = it->second[6];
		double yAl = it->second[7];
		double CAl = it->second[8];
		double xBC = it->second[9];
		double CBC = it->second[10];
		double pk_H = pkH - 3.0 + e_H * dif_temp;
		CAl *= 1e-6;
		CBC *= 1e-6;
		double rate_H = pow(10.0, -pk_H) * pow(PhreeqcPtr->activity("H+"), nH) / 
			(pow(1.0 + aAl / CAl, yAl) * pow(1.0 + BC / CBC, xBC));
		//	rate by hydrolysis
		double pkH2O = it->second[11];
		yAl = it->second[12];
		CAl = it->second[13];
		xBC = it->second[14];
		CBC = it->second[15];
		double zSi = it->second[16];
		double CSi = it->second[17];
		CAl *= 1e-6;
		CBC *= 1e-6;
		CSi *= 1e-6;
		double pk_H2O = pkH2O - 3.0 + e_H2O * dif_temp;
		double rate_H2O = pow(10.0, -pk_H2O) / (pow(1.0 + aAl / CAl, yAl) * pow(1.0 + BC / CBC, xBC) * pow(1.0 + aSi / CSi, zSi));
		//	rate by CO2
		double pKCO2 = it->second[18];
		double nCO2 = it->second[19];
		double pk_CO2 = pKCO2 - 3.0 + e_CO2 * dif_temp;
		double rate_CO2 = pow(10.0, -pk_CO2) * pow(PhreeqcPtr->saturation_ratio("CO2(g)"), nCO2);
		//	rate by Organic Acids
		double pkOrg = it->second[20];
		double nOrg = it->second[21];
		double COrg = it->second[22];
		COrg *= 1e-6;
		double pk_Org = pkOrg - 3.0 + e_OA * dif_temp;
		double rate_Org = pow(10.0, -pkOrg) * pow(R / (1 + R / COrg), nOrg);
		//	rate by OH-
		double pkOH = it->second[23];
		double wOH = it->second[24];
		yAl = it->second[25];
		CAl = it->second[26];
		xBC = it->second[27];
		CBC = it->second[28];
		zSi = it->second[29];
		CSi = it->second[30];
		CAl *= 1e-6;
		CBC *= 1e-6;
		CSi *= 1e-6;
		double pk_OH = pkOH - 3.0 + e_OH * dif_temp;
		double rate_OH = pow(10.0, -pk_OH) * pow(PhreeqcPtr->activity("OH-"), wOH) /
			(pow(1.0 + aAl / CAl, yAl) * pow(1.0 + BC / CBC, xBC) * pow(1.0 + aSi / CSi, zSi));
		// sum rates
		double rate = rate_H + rate_H2O + rate_CO2 + rate_Org + rate_OH;
		n.UU.val = rate;
		//	Sverdrup_rate
		//		# in KINETICS, define 34 parms:
		//	#          affinity m ^ 2 / mol roughness, temperature_factors(TABLE 4) : e_H e_H2O e_CO2 e_OA e_OH, \
		//# (TABLE 3): pkH  nH yAl CAl xBC CBC,   pKH2O yAl CAl xBC CBC zSi CSi,  pKCO2 nCO2  pkOrg nOrg COrg, pkOH wOH yAl CAl xBC CBC zSi CSi
		//		10  affinity = get(-99, 1)
		//		20
		//		30  REM # specific area m2 / mol, surface roughness
		//		40  sp_area = get(-99, 2) : roughness = get(-99, 3)
		//		50
		//		60  REM # temperature factors
		//		70  dif_temp = 1 / TK - 1 / 281
		//		80  e_H = get(-99, 4) : e_H2O = get(-99, 5) : e_CO2 = get(-99, 6) : e_OA = get(-99, 7) : e_OH = get(-99, 8)
		//		90
		//		100 BC = ACT("Na+") + ACT("K+") + ACT("Mg+2") + ACT("Ca+2")
		//		110 aAl = act("Al+3")
		//		120 aSi = act("H4SiO4")
		//		130 R = tot("OrganicMatter")
		//		140
		//		150 REM # rate by H +
		//		160 pkH = get(-99, 9) : nH = get(-99, 10) : yAl = get(-99, 11) : CAl = get(-99, 12) : xBC = get(-99, 13) : CBC = get(-99, 14)
		//		170 pk_H = pkH - 3 + e_H * dif_temp
		//		180 CAl = CAl * 1e-6
		//		190 CBC = CBC * 1e-6
		//		200 rate_H = 10 ^ -pk_H * ACT("H+") ^ nH / ((1 + aAl / CAl) ^ yAl * (1 + BC / CBC) ^ xBC)
		//		210
		//		220 REM # rate by hydrolysis
		//		230 pkH2O = get(-99, 15) : yAl = get(-99, 16) : CAl = get(-99, 17) : xBC = get(-99, 18) : CBC = get(-99, 19) : zSi = get(-99, 20) : CSi = get(-99, 21)
		//		240 CAl = CAl * 1e-6
		//		250 CBC = CBC * 1e-6
		//		260 CSi = CSi * 1e-6
		//		270 pk_H2O = pkH2O - 3 + e_H2O * dif_temp
		//		280 rate_H2O = 10 ^ -pk_H2O / ((1 + aAl / CAl) ^ yAl * (1 + BC / CBC) ^ xBC * (1 + aSi / CSi) ^ zSi)
		//		290
		//		300 REM # rate by CO2
		//		310 pKCO2 = get(-99, 22) : nCO2 = get(-99, 23)
		//		320 pk_CO2 = pkCO2 - 3 + e_CO2 * dif_temp
		//		330 rate_CO2 = 10 ^ -pk_CO2 * SR("CO2(g)") ^ nCO2
		//		340
		//		350 REM # rate by Organic Acids
		//		360 pkOrg = get(-99, 24) : nOrg = get(-99, 25) : COrg = get(-99, 26)
		//		370 COrg = COrg * 1e-6
		//		380 pk_Org = pkOrg - 3 + e_OA * dif_temp
		//		390 rate_Org = 10 ^ -pk_Org * (R / (1 + R / COrg)) ^ nOrg
		//		400
		//		410 REM # rate by OH -
		//		420 pkOH = get(-99, 27) : wOH = get(-99, 28) : yAl = get(-99, 29) : CAl = get(-99, 30) : xBC = get(-99, 31) : CBC = get(-99, 32) : zSi = get(-99, 33) : CSi = get(-99, 34)
		//		430 CAl = CAl * 1e-6
		//		440 CBC = CBC * 1e-6
		//		450 CSi = CSi * 1e-6
		//		460 pk_OH = pkOH - 3 + e_OH * dif_temp
		//		470 rate_OH = 10 ^ -pk_OH * ACT("OH-") ^ wOH / ((1 + aAl / CAl) ^ yAl * (1 + BC / CBC) ^ xBC * (1 + aSi / CSi) ^ zSi)#  : print rate_OH
		//		480
		//		490 rate = rate_H + rate_H2O + rate_CO2 + rate_Org + rate_OH
		//		500 area = sp_area * M0 * (M / M0) ^ 0.67
		//		510
		//		520 rate = roughness * area * rate * affinity
		//		530 SAVE rate * TIME
		//		- end
	}
	break;

	case tokrate_hermanska:
	{
		require(toklp, LINK);
		char* min_name = strexpr(LINK);
		require(tokrp, LINK);
		if (parse_all) {
			PhreeqcPtr->PHRQ_free(min_name);
			n.UU.val = 1;
			break;
		}
		std::string min_string = min_name;
		PhreeqcPtr->PHRQ_free(min_name);
		Utilities::str_tolower(min_string);
		std::map<std::string, std::vector<double> >::const_iterator it = PhreeqcPtr->rate_parameters_hermanska.find(min_string);
		if (it == PhreeqcPtr->rate_parameters_hermanska.end())
		{
			std::ostringstream oss;
			oss << "Hermanska rate parameters not found for " << min_name << "\n";
			snerr(oss.str().c_str());
		}
		if (it->second.size() != 11)
		{
			std::ostringstream oss;
			oss << "RATE_HERMANSKA requires 11 rate parameters, " << it->second.size() << " were found for " << min_name << "\n";
			snerr(oss.str().c_str());
		}
		//	gas constant * Tk, act("H+")
		double RT = 8.314e-3 * PhreeqcPtr->tk_x;
		double aH = PhreeqcPtr->activity("H+");

		//	rate by H+
		double lgk_H = it->second[0];
		double Aa = it->second[1];
		double e_H = it->second[2];
		double nH = it->second[3];
		double rate_H = Aa * exp(-e_H / RT) * pow(aH, nH);

		//	rate by hydrolysis
		double rate_H2O = 0.0, lgk_H2O = it->second[4];
		if (lgk_H2O)
		{
			double Ab = it->second[5];
			double e_H2O = it->second[6];
			rate_H2O = Ab * exp(-e_H2O / RT);
		}

		//	rate by OH-
		//	180 lgk_OH = get(-99, 11) : Ac = get(-99, 12) : e_OH = get(-99, 13) : nOH = get(-99, 14)
		//	190 rate_OH = Ac * exp(-e_OH / RT) * aH ^ nOH
		double  rate_OH = 0.0, lgk_OH = it->second[7];
		if (lgk_OH)
		{
			double Ac = it->second[8];
			double e_OH = it->second[9];
			double nOH = it->second[10];
			rate_OH = Ac * exp(-e_OH / RT) * pow(aH, nOH);
		}
		// sum rates
		double rate = rate_H + rate_H2O + rate_OH;
		n.UU.val = rate;

//		Hermanska_rate
//			# in KINETICS, define 14 parms:
//		# parms affinity m ^ 2 / mol roughness, (TABLE 2) : (acid)logk25 Aa Ea na(neutral)logk25 Ab Eb(basic)logk25 Ac Ec nc
//# (Note that logk25 values are not used, they were transformed to A's.)
//	10  affinity = get(-99, 1) # retrieve number from memory
//	20
//	30  REM # specific area m2 / mol, surface roughness
//	40  sp_area = get(-99, 2) : roughness = get(-99, 3)
//	50
//	60  REM # gas constant * Tk, act("H+")
//	70  RT = 8.314e-3 * TK     : aH = act("H+")
//	80
//	90  REM # rate by H +
//	100 lgk_H = get(-99, 4) : Aa = get(-99, 5) : e_H = get(-99, 6) : nH = get(-99, 7)
//	110 rate_H = Aa * exp(-e_H / RT) * aH ^ nH
//	120
//	130 REM # rate by hydrolysis
//	140 lgk_H2O = get(-99, 8) : Ab = get(-99, 9) : e_H2O = get(-99, 10)
//	150 rate_H2O = Ab * exp(-e_H2O / RT)
//	160
//	170 REM # rate by OH -
//	180 lgk_OH = get(-99, 11) : Ac = get(-99, 12) : e_OH = get(-99, 13) : nOH = get(-99, 14)
//	190 rate_OH = Ac * exp(-e_OH / RT) * aH ^ nOH
//	200
//	210 rate = rate_H + rate_H2O + rate_OH
//	220 area = sp_area * M0 * (M / M0) ^ 0.67
//	230
//	240 rate = area * roughness * rate * affinity
//	250 SAVE rate * TIME
//	- end

	}
	break;

	case tokmeang:
	{
		require(toklp, LINK);
		char* min_name = strexpr(LINK);
		require(tokrp, LINK);
		if (parse_all) {
			PhreeqcPtr->PHRQ_free(min_name);
			n.UU.val = 1;
			break;
		}
		std::string min_string = min_name;
		PhreeqcPtr->PHRQ_free(min_name);
		Utilities::str_tolower(min_string);
		std::map<std::string, cxxNameDouble>::const_iterator it = PhreeqcPtr->mean_gammas.find(min_string);
		if (it == PhreeqcPtr->mean_gammas.end() || it->second.size() == 0)
		{
			std::ostringstream oss;
			oss << "No definition in MEAN_GAMMAS found for " << min_name << "\n";
			snerr(oss.str().c_str());
		}

		double mg = 1.0;
		double sum = 0.0;
		cxxNameDouble::const_iterator it_nd = it->second.begin();
		for (; it_nd != it->second.end(); it_nd++)
		{
			double g = PhreeqcPtr->activity_coefficient(it_nd->first.c_str());
			mg *= pow(g, it_nd->second);
			sum += it_nd->second;
		}
		mg = pow(mg, 1.0 / sum);
		n.UU.val = mg;
	}
	break;
	case tokpercent_error:
	{
		n.UU.val = (parse_all) ? 1 : 100 * PhreeqcPtr->cb_x / PhreeqcPtr->total_ions_x;
	}
	break;

	case tokphase_formula:
	case tokphase_formula_:
	{
		require(toklp, LINK);
		std::string phase_name(stringfactor(STR1, LINK));
		varrec* elts_varrec = NULL, * coef_varrec = NULL;
		cxxNameDouble stoichiometry;
		/*
		*  Parse arguments
		*/
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			/* phase_formula("calcite", count, elt, coef) */
			/* return formula */
			/*int c; */
			/*  struct varrec *count_varrec, *names_varrec, *types_varrec, *moles_varrec; */
			/*  struct varrec *count_varrec, *elt_varrec, *coef_varrec; */
			/* return number of species */
			LINK->t = LINK->t->next;
			count_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || count_varrec->stringvar != 0)
				snerr(": Cannot find count variable");

			/* return number of names of species */
			LINK->t = LINK->t->next;
			require(tokcomma, LINK);
			elts_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || elts_varrec->stringvar != 1)
				snerr(": Cannot find element string variable");

			/* return coefficients of species */
			LINK->t = LINK->t->next;
			require(tokcomma, LINK);
			coef_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || coef_varrec->stringvar != 0)
				snerr(": Cannot find coefficient variable");
			LINK->t = LINK->t->next;
			arg_num = 4;
		}
		else
		{
			arg_num = 1;
		}
		require(tokrp, LINK);

		if (arg_num > 1)
		{
			free_dim_stringvar(elts_varrec);
			PhreeqcPtr->free_check_null(coef_varrec->UU.U0.arr);
			coef_varrec->UU.U0.arr = NULL;
		}
		/*
		*  Call subroutine
		*/
		std::string form = PhreeqcPtr->phase_formula(phase_name, stoichiometry);

		// put formula as return value
		n.stringval = true;
		n.UU.sval = PhreeqcPtr->string_duplicate(form.c_str());

		/*
		*  fill in varrec structure
		*/

		if (arg_num > 1)
		{
			size_t count = stoichiometry.size();
			*count_varrec->UU.U0.val = (LDBLE)count;
			/*
			* malloc space
			*/
			elts_varrec->UU.U1.sarr = (char**)PhreeqcPtr->PHRQ_malloc((count + 1) * sizeof(char*));
			if (elts_varrec->UU.U1.sarr == NULL)
				PhreeqcPtr->malloc_error();
			coef_varrec->UU.U0.arr = (LDBLE*)PhreeqcPtr->PHRQ_malloc((count + 1) * sizeof(LDBLE));
			if (coef_varrec->UU.U0.arr == NULL)
				PhreeqcPtr->malloc_error();

			// first position not used
			elts_varrec->UU.U1.sarr[0] = NULL;
			coef_varrec->UU.U0.arr[0] = 0;

			// set dims for Basic array
			for (i = 0; i < maxdims; i++)
			{
				elts_varrec->dims[i] = 0;
				coef_varrec->dims[i] = 0;
			}
			// set dims for first dimension and number of dims
			elts_varrec->dims[0] = (long)(count + 1);
			coef_varrec->dims[0] = (long)(count + 1);
			elts_varrec->numdims = 1;
			coef_varrec->numdims = 1;

			// fill in arrays
			i = 1;
			for (cxxNameDouble::iterator it = stoichiometry.begin(); it != stoichiometry.end(); it++)
			{
				elts_varrec->UU.U1.sarr[i] = PhreeqcPtr->string_duplicate((it->first).c_str());
				coef_varrec->UU.U0.arr[i] = it->second;
				i++;
			}

		}
	}
	break;

	case tokphase_vm:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->phase_vm(str);
	}
	break;

	case tokpot_v:
	{
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->use.Get_solution_ptr()->Get_potV();
	}
	break;

	case tokpr_p:
	{
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->pr_pressure(stringfactor(STR1, LINK));
	}
	break;

	case tokpr_phi:
	{
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->pr_phi(stringfactor(STR1, LINK));
	}
	break;

	case tokpressure:
	{
		n.UU.val = PhreeqcPtr->pressure();
	}
	break;

	case tokqbrn:
	{
		n.UU.val = PhreeqcPtr->QBrn;
	}
	break;

	case tokrho:
	{
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->calc_dens();
	}
	break;

	case tokrho_0:
	{
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->rho_0;
	}
	break;

	case tokrtrim:
	{
		n.stringval = true;
		require(toklp, LINK);
		string1 = stringfactor(STR1, LINK);
		require(tokrp, LINK);
		trim_right(STR1);
		n.UU.sval = PhreeqcPtr->string_duplicate(STR1.c_str());
	}
	break;

	case tokrxn:
	{
		if (PhreeqcPtr->state == REACTION ||
			PhreeqcPtr->state == ADVECTION ||
			PhreeqcPtr->state == TRANSPORT)
		{
			n.UU.val = PhreeqcPtr->step_x;
		}
		else
		{
			n.UU.val = 0.0;
		}
	}
	break;

	case toks_s:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->find_ss_comp(str);
	}
	break;

	case toksc:
	{
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->calc_SC();
	}
	break;

	case toksetdiff_c:
	{
		double d, d_v_d = 0;

		require(toklp, LINK);

		const char* str = stringfactor(STR1, LINK);
		require(tokcomma, LINK);

		d = realexpr(LINK);
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			d_v_d = realexpr(LINK);
		}
		require(tokrp, LINK);

		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->setdiff_c(str, d, d_v_d);
	}
	break;

	case toksi:
	{
		const char* str = stringfactor(STR1, LINK);
		if (parse_all)
		{
			n.UU.val = 1;
		}
		else
		{
			PhreeqcPtr->saturation_index(str, &l_dummy, &n.UU.val);
		}
	}
	break;

	case toksim_no:
	{
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->simulation;
	}
	break;

	case toksim_time:
	{
		if (!PhreeqcPtr->use.Get_kinetics_in())
		{
			if (PhreeqcPtr->state == PHAST)
			{
				n.UU.val = PhreeqcPtr->rate_sim_time;
			}
			else if (PhreeqcPtr->state == TRANSPORT)
			{
				n.UU.val = PhreeqcPtr->transport_step * PhreeqcPtr->timest;
			}
			else if (PhreeqcPtr->state == ADVECTION)
			{
				if (PhreeqcPtr->advection_kin_time_defined == TRUE)
				{
					n.UU.val = PhreeqcPtr->advection_step * PhreeqcPtr->advection_kin_time;
				}
				else
				{
					n.UU.val = PhreeqcPtr->advection_step;
				}
			}
			else
			{
				n.UU.val = 0;
			}
		}
		else
		{
			n.UU.val = PhreeqcPtr->rate_sim_time;
		}
	}
	break;

	case toksoln_vol:
	{
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->calc_solution_volume();
	}
	break;

	case tokspecies_formula:
	case tokspecies_formula_:
	{
		require(toklp, LINK);
		std::string species_name(stringfactor(STR1, LINK));
		varrec* elts_varrec = NULL, * coef_varrec = NULL;
		cxxNameDouble stoichiometry;
		/*
		*  Parse arguments
		*/
		require(tokcomma, LINK);

		count_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || count_varrec->stringvar != 0)
			snerr(": Cannot find count variable");

		/* return number of names of species */
		LINK->t = LINK->t->next;
		require(tokcomma, LINK);
		elts_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || elts_varrec->stringvar != 1)
			snerr(": Cannot find element string variable");

		/* return coefficients of species */
		LINK->t = LINK->t->next;
		require(tokcomma, LINK);
		coef_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || coef_varrec->stringvar != 0)
			snerr(": Cannot find coefficient variable");
		LINK->t = LINK->t->next;

		require(tokrp, LINK);

		free_dim_stringvar(elts_varrec);
		PhreeqcPtr->free_check_null(coef_varrec->UU.U0.arr);
		coef_varrec->UU.U0.arr = NULL;
		/*
		*  Call subroutine
		*/
		std::string type = PhreeqcPtr->species_formula(species_name, stoichiometry);

		// put type as return value
		n.stringval = true;
		n.UU.sval = PhreeqcPtr->string_duplicate(type.c_str());

		/*
		*  fill in varrec structure
		*/

		size_t count = stoichiometry.size();
		*count_varrec->UU.U0.val = (LDBLE)count;
		/*
		* malloc space
		*/
		elts_varrec->UU.U1.sarr = (char**)PhreeqcPtr->PHRQ_malloc((count + 1) * sizeof(char*));
		if (elts_varrec->UU.U1.sarr == NULL)
		{
			PhreeqcPtr->malloc_error();
#if !defined(R_SO)
			exit(4);
#endif
		}
		coef_varrec->UU.U0.arr = (LDBLE*)PhreeqcPtr->PHRQ_malloc((count + 1) * sizeof(LDBLE));
		if (coef_varrec->UU.U0.arr == NULL)
		{
			PhreeqcPtr->malloc_error();
#if !defined(R_SO)
			exit(4);
#endif
		}

		// first position not used
		elts_varrec->UU.U1.sarr[0] = NULL;
		coef_varrec->UU.U0.arr[0] = 0;

		// set dims for Basic array
		for (i = 0; i < maxdims; i++)
		{
			elts_varrec->dims[i] = 0;
			coef_varrec->dims[i] = 0;
		}
		// set dims for first dimension and number of dims
		elts_varrec->dims[0] = (long)(count + 1);
		coef_varrec->dims[0] = (long)(count + 1);
		elts_varrec->numdims = 1;
		coef_varrec->numdims = 1;

		// fill in arrays
		i = 1;
		for (cxxNameDouble::iterator it = stoichiometry.begin(); it != stoichiometry.end(); it++)
		{
			elts_varrec->UU.U1.sarr[i] = PhreeqcPtr->string_duplicate((it->first).c_str());
			coef_varrec->UU.U0.arr[i] = it->second;
			i++;
		}
	}
	break;

	case tokphase_equation:
	case tokphase_equation_:
	{
		require(toklp, LINK);
		std::string phase_name(stringfactor(STR1, LINK));
		varrec* elts_varrec = NULL, * coef_varrec = NULL;
		std::vector<std::pair<std::string, double> > stoichiometry;
		/*
		*  Parse arguments
		*/
		require(tokcomma, LINK);

		count_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || count_varrec->stringvar != 0)
			snerr(": Cannot find count variable");

		/* return number of names of species */
		LINK->t = LINK->t->next;
		require(tokcomma, LINK);
		elts_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || elts_varrec->stringvar != 1)
			snerr(": Cannot find species string variable");

		/* return coefficients of species */
		LINK->t = LINK->t->next;
		require(tokcomma, LINK);
		coef_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || coef_varrec->stringvar != 0)
			snerr(": Cannot find coefficient variable");
		LINK->t = LINK->t->next;

		require(tokrp, LINK);

		free_dim_stringvar(elts_varrec);
		PhreeqcPtr->free_check_null(coef_varrec->UU.U0.arr);
		coef_varrec->UU.U0.arr = NULL;
		/*
		*  Call subroutine
		*/
		std::string eq = PhreeqcPtr->phase_equation(phase_name, stoichiometry);

		// put type as return value
		n.stringval = true;
		n.UU.sval = PhreeqcPtr->string_duplicate(eq.c_str());

		/*
		*  fill in varrec structure
		*/

		size_t count = stoichiometry.size();
		*count_varrec->UU.U0.val = (LDBLE)count;
		/*
		* malloc space
		*/
		elts_varrec->UU.U1.sarr = (char**)PhreeqcPtr->PHRQ_malloc((count + 1) * sizeof(char*));
		if (elts_varrec->UU.U1.sarr == NULL)
		{
			PhreeqcPtr->malloc_error();
#if !defined(R_SO)
			exit(4);
#endif
		}
		coef_varrec->UU.U0.arr = (LDBLE*)PhreeqcPtr->PHRQ_malloc((count + 1) * sizeof(LDBLE));
		if (coef_varrec->UU.U0.arr == NULL)
		{
			PhreeqcPtr->malloc_error();
#if !defined(R_SO)
			exit(4);
#endif
		}

		// first position not used
		elts_varrec->UU.U1.sarr[0] = NULL;
		coef_varrec->UU.U0.arr[0] = 0;

		// set dims for Basic array
		for (i = 0; i < maxdims; i++)
		{
			elts_varrec->dims[i] = 0;
			coef_varrec->dims[i] = 0;
		}
		// set dims for first dimension and number of dims
		elts_varrec->dims[0] = (long)(count + 1);
		coef_varrec->dims[0] = (long)(count + 1);
		elts_varrec->numdims = 1;
		coef_varrec->numdims = 1;

		// fill in arrays
		i = 1;
		for (std::vector<std::pair<std::string, double > >::iterator it = stoichiometry.begin(); it != stoichiometry.end(); it++)
		{
			elts_varrec->UU.U1.sarr[i] = PhreeqcPtr->string_duplicate((it->first).c_str());
			coef_varrec->UU.U0.arr[i] = it->second;
			i++;
		}
	}
	break;

	case tokspecies_equation:
	case tokspecies_equation_:
	{
		require(toklp, LINK);
		std::string species_name(stringfactor(STR1, LINK));
		varrec* elts_varrec = NULL, * coef_varrec = NULL;
		std::vector<std::pair<std::string, double> > stoichiometry;
		/*
		*  Parse arguments
		*/
		require(tokcomma, LINK);

		count_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || count_varrec->stringvar != 0)
			snerr(": Cannot find count variable");

		/* return number of names of species */
		LINK->t = LINK->t->next;
		require(tokcomma, LINK);
		elts_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || elts_varrec->stringvar != 1)
			snerr(": Cannot find species string variable");

		/* return coefficients of species */
		LINK->t = LINK->t->next;
		require(tokcomma, LINK);
		coef_varrec = LINK->t->UU.vp;
		if (LINK->t->kind != tokvar || coef_varrec->stringvar != 0)
			snerr(": Cannot find coefficient variable");
		LINK->t = LINK->t->next;

		require(tokrp, LINK);

		free_dim_stringvar(elts_varrec);
		PhreeqcPtr->free_check_null(coef_varrec->UU.U0.arr);
		coef_varrec->UU.U0.arr = NULL;
		/*
		*  Call subroutine
		*/
		std::string eq = PhreeqcPtr->species_equation(species_name, stoichiometry);

		// put type as return value
		n.stringval = true;
		n.UU.sval = PhreeqcPtr->string_duplicate(eq.c_str());

		/*
		*  fill in varrec structure
		*/

		size_t count = stoichiometry.size();
		*count_varrec->UU.U0.val = (LDBLE)count;
		/*
		* malloc space
		*/
		elts_varrec->UU.U1.sarr = (char**)PhreeqcPtr->PHRQ_malloc((count + 1) * sizeof(char*));
		if (elts_varrec->UU.U1.sarr == NULL)
		{
			PhreeqcPtr->malloc_error();
#if !defined(R_SO)
			exit(4);
#endif
		}
		coef_varrec->UU.U0.arr = (LDBLE*)PhreeqcPtr->PHRQ_malloc((count + 1) * sizeof(LDBLE));
		if (coef_varrec->UU.U0.arr == NULL)
		{
			PhreeqcPtr->malloc_error();
#if !defined(R_SO)
			exit(4);
#endif
		}

		// first position not used
		elts_varrec->UU.U1.sarr[0] = NULL;
		coef_varrec->UU.U0.arr[0] = 0;

		// set dims for Basic array
		for (i = 0; i < maxdims; i++)
		{
			elts_varrec->dims[i] = 0;
			coef_varrec->dims[i] = 0;
		}
		// set dims for first dimension and number of dims
		elts_varrec->dims[0] = (long)(count + 1);
		coef_varrec->dims[0] = (long)(count + 1);
		elts_varrec->numdims = 1;
		coef_varrec->numdims = 1;

		// fill in arrays
		i = 1;
		for (std::vector<std::pair<std::string, double > >::iterator it = stoichiometry.begin(); it != stoichiometry.end(); it++)
		{
			elts_varrec->UU.U1.sarr[i] = PhreeqcPtr->string_duplicate((it->first).c_str());
			coef_varrec->UU.U0.arr[i] = it->second;
			i++;
		}
	}
	break;

	case toksr:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->saturation_ratio(str);
	}
	break;

	case tokstep_no:
	{
		if (PhreeqcPtr->state == PHAST)
		{
			n.UU.val = 0;
		}
		else if (PhreeqcPtr->state == TRANSPORT)
		{
			n.UU.val = PhreeqcPtr->transport_step;
		}
		else if (PhreeqcPtr->state == ADVECTION)
		{
			n.UU.val = PhreeqcPtr->advection_step;
		}
		else if (PhreeqcPtr->state == REACTION)
		{
			n.UU.val = PhreeqcPtr->reaction_step;
		}
		else
		{
			n.UU.val = 0;
		}
	}
	break;

	case tokstr_e_:
	{
		// left parenthesis
		require(toklp, LINK);

		// set number
		LDBLE nmbr;
		nmbr = realexpr(LINK);

		// set total length
		require(tokcomma, LINK);
		int length = (int)realexpr(LINK);

		// set total length
		require(tokcomma, LINK);
		int width = (int)realexpr(LINK);

		// right parenthesis
		require(tokrp, LINK);

		// Make work space
		int max_length = length < 256 ? 256 : length;
		char* token = (char*)PhreeqcPtr->PHRQ_calloc((max_length + 1), sizeof(char));
		if (token == NULL) PhreeqcPtr->malloc_error();

		std::string std_num;
		{
			snprintf(token, max_length, "%*.*e", length, width, nmbr);
			std_num = token;
		}

		// set function value
		n.UU.sval = (char*)PhreeqcPtr->PHRQ_calloc(std_num.size() + 2, sizeof(char));
		if (n.UU.sval == NULL)
			PhreeqcPtr->malloc_error();
		strcpy(n.UU.sval, std_num.c_str());
		n.stringval = true;

		// free work space
		PhreeqcPtr->free_check_null(token);
	}
	break;

	case tokstr_f_:
	{
		// left parenthesis
		require(toklp, LINK);

		// set number
		LDBLE nmbr;
		nmbr = realexpr(LINK);

		// set total length
		require(tokcomma, LINK);
		int length = (int)realexpr(LINK);

		// set total length
		require(tokcomma, LINK);
		int width = (int)realexpr(LINK);

		// right parenthesis
		require(tokrp, LINK);

		// Make work space
		int max_length = length < 256 ? 256 : length;
		char* token = (char*)PhreeqcPtr->PHRQ_calloc((max_length + 1), sizeof(char));
		if (token == NULL) PhreeqcPtr->malloc_error();

		std::string std_num;
		{
			snprintf(token, max_length, "%*.*f", length, width, nmbr);
			std_num = token;
		}

		// set function value
		n.UU.sval = (char*)PhreeqcPtr->PHRQ_calloc(std_num.size() + 2, sizeof(char));
		if (n.UU.sval == NULL)
			PhreeqcPtr->malloc_error();
		strcpy(n.UU.sval, std_num.c_str());
		n.stringval = true;

		// free work space
		PhreeqcPtr->free_check_null(token);
	}
	break;

	case toksum_gas:
	{
		require(toklp, LINK);
		mytemplate = stringfactor(STR1, LINK);
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			elt_name = stringfactor(STR2, LINK);
		}
		else
		{
			elt_name = NULL;
		}
		require(tokrp, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->sum_match_gases(mytemplate, elt_name);
	}
	break;

	case toksum_s_s:
	{
		require(toklp, LINK);
		mytemplate = stringfactor(STR1, LINK);
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			elt_name = stringfactor(STR2, LINK);
		}
		else
		{
			elt_name = NULL;
		}
		require(tokrp, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->sum_match_ss(mytemplate, elt_name);
	}
	break;

	case toksum_species:
	{
		require(toklp, LINK);
		mytemplate = stringfactor(STR1, LINK);
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			elt_name = stringfactor(STR2, LINK);
		}
		else
		{
			elt_name = NULL;
		}
		require(tokrp, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->sum_match_species(mytemplate, elt_name);
	}
	break;

	case toksurf:
	{
		require(toklp, LINK);
		elt_name = stringfactor(STR1, LINK);
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			surface_name = stringfactor(STR2, LINK);
		}
		else
		{
			surface_name = NULL;
		}
		require(tokrp, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->surf_total(elt_name, surface_name);
	}
	break;

	case toksys:
	{
		int isort = 0;
		require(toklp, LINK);
		elt_name = stringfactor(STR1, LINK);
		/*
		 *  Parse arguments
		 */
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			/* return number of species */
			LINK->t = LINK->t->next;
			count_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || !count_varrec || count_varrec->stringvar != 0)
			{
				snerr(": can`t find variable");
#if !defined(R_SO)
				exit(4);
#endif
			}

			/* return number of names of species */
			LINK->t = LINK->t->next;
			require(tokcomma, LINK);
			names_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || !names_varrec || names_varrec->stringvar != 1)
			{
				snerr(": can`t find name of species");
#if !defined(R_SO)
				exit(4);
#endif
			}

			/* return number of types of species */
			LINK->t = LINK->t->next;
			require(tokcomma, LINK);
			types_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || types_varrec->stringvar != 1)
				snerr(": can`t find type of species");

			/* return number of moles  of species */
			LINK->t = LINK->t->next;
			require(tokcomma, LINK);
			moles_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || moles_varrec->stringvar != 0)
				snerr(": can`t find moles of species");
			LINK->t = LINK->t->next;
			if (LINK->t != NULL && LINK->t->kind == tokcomma)
			{
				LINK->t = LINK->t->next;
				isort = intexpr(LINK);
				arg_num = 5;
			}
			else
			{
				arg_num = 4;
			}
		}
		else
		{
			arg_num = 1;
		}
		require(tokrp, LINK);

		if (arg_num > 1)
		{
			free_dim_stringvar(names_varrec);
			free_dim_stringvar(types_varrec);
			PhreeqcPtr->free_check_null(moles_varrec->UU.U0.arr);
			moles_varrec->UU.U0.arr = NULL;
		}
		/*
		 *  Call subroutine
		 */
		 /*
			n.UU.val = system_total(elt_name, count_varrec->UU.U0.val, &(names_varrec->UU.U1.sarr), &(types_varrec->UU.U1.sarr), &(moles_varrec->UU.U0.arr));
		  */
		if (parse_all)
		{
			PhreeqcPtr->sys_tot = 0;
			//PhreeqcPtr->count_sys = 1000;
			//int count_sys = PhreeqcPtr->count_sys;
			size_t count_sys = 1000;
			names_arg = (char**)PhreeqcPtr->PHRQ_calloc((size_t)(count_sys + 1), sizeof(char*));
			if (names_arg == NULL)
			{
				PhreeqcPtr->malloc_error();
#if !defined(R_SO)
				exit(4);
#endif
			}
			types_arg = (char**)PhreeqcPtr->PHRQ_calloc((size_t)(count_sys + 1), sizeof(char*));
			if (types_arg == NULL)
			{
				PhreeqcPtr->malloc_error();
#if !defined(R_SO)
				exit(4);
#endif
			}
			moles_arg = (LDBLE*)PhreeqcPtr->PHRQ_calloc((size_t)(count_sys + 1), sizeof(LDBLE));
			if (moles_arg == NULL)
			{
				PhreeqcPtr->malloc_error();
#if !defined(R_SO)
				exit(4);
#endif
			}
			names_arg[0] = NULL;
			types_arg[0] = NULL;
			moles_arg[0] = 0;
			count_species = (LDBLE)count_sys;
			n.UU.val = 0;
		}
		else
		{
			n.UU.val = PhreeqcPtr->system_total(elt_name, &count_species, &(names_arg),
				&(types_arg), &(moles_arg), isort);
		}

		/*
		 *  fill in varrec structure
		 */
		if (arg_num > 1)
		{
			*count_varrec->UU.U0.val = count_species;
			names_varrec->UU.U1.sarr = names_arg;
			types_varrec->UU.U1.sarr = types_arg;
			moles_varrec->UU.U0.arr = moles_arg;

			for (i = 0; i < maxdims; i++)
			{
				names_varrec->dims[i] = 0;
				types_varrec->dims[i] = 0;
				moles_varrec->dims[i] = 0;
			}
			names_varrec->dims[0] = (long)(*count_varrec->UU.U0.val) + 1;
			types_varrec->dims[0] = (long)(*count_varrec->UU.U0.val) + 1;
			moles_varrec->dims[0] = (long)(*count_varrec->UU.U0.val) + 1;
			names_varrec->numdims = 1;
			types_varrec->numdims = 1;
			moles_varrec->numdims = 1;
		}
		else
		{
			for (i = 0; i < count_species + 1; i++)
			{
				PhreeqcPtr->free_check_null(names_arg[i]);
				PhreeqcPtr->free_check_null(types_arg[i]);
			}
			PhreeqcPtr->free_check_null(names_arg);
			PhreeqcPtr->free_check_null(types_arg);
			PhreeqcPtr->free_check_null(moles_arg);
		}
	}
	break;

	case tokt_sc:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->calc_t_sc(str);
	}
	break;

	case tokf_visc:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->calc_f_visc(str);
	}
	break;

	case toktc:
	{
		n.UU.val = PhreeqcPtr->tc_x;
	}
	break;

	case toktime:
	{
		n.UU.val = PhreeqcPtr->rate_time;
	}
	break;

	case toktitle:
	{
		n.stringval = true;
		if (strlen(PhreeqcPtr->last_title_x.c_str()) == 0)
		{
			PhreeqcPtr->last_title_x = " ";
		}
		n.UU.sval = PhreeqcPtr->string_duplicate(PhreeqcPtr->last_title_x.c_str());
		while (PhreeqcPtr->replace("\t", " ", n.UU.sval));
	}
	break;

	case toktk:
	{
		n.UU.val = PhreeqcPtr->tc_x + 273.15;
	}
	break;

	case toktot:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->total(str);
	}
	break;

	case toktotal_time:
	{
		if (!PhreeqcPtr->use.Get_kinetics_in())
		{
			if (PhreeqcPtr->state == PHAST)
			{
				n.UU.val = PhreeqcPtr->rate_sim_time_end;
			}
			else if (PhreeqcPtr->state == TRANSPORT)
			{
				n.UU.val = PhreeqcPtr->initial_total_time + PhreeqcPtr->transport_step * PhreeqcPtr->timest;
			}
			else if (PhreeqcPtr->state == ADVECTION)
			{
				n.UU.val =
					PhreeqcPtr->initial_total_time + PhreeqcPtr->advection_step * PhreeqcPtr->advection_kin_time;
			}
			else
			{
				n.UU.val = 0;
			}
		}
		else
		{
			n.UU.val = PhreeqcPtr->initial_total_time + PhreeqcPtr->rate_sim_time;
		}
	}
	break;

	case toktotmole:
	case toktotmol:
	case toktotmoles:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->total_mole(str);
	}
	break;

	case toktrim:
	{
		n.stringval = true;
		require(toklp, LINK);
		string1 = stringfactor(STR1, LINK);
		require(tokrp, LINK);
		STR1 = trim(STR1);
		n.UU.sval = PhreeqcPtr->string_duplicate(STR1.c_str());
	}
	break;

	case tokviscos:
	{
		if (PhreeqcPtr->print_viscosity)
			PhreeqcPtr->viscosity(nullptr);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->viscos;
	}
	break;

	case tokviscos_0:
	{
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->viscos_0;
	}
	break;

	case tokvm:
	{
		const char* str = stringfactor(STR1, LINK);
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->aqueous_vm(str);
	}
	break;
	/*
	* End of PHREEQC functions
	*/
	case toksa_declercq: // Undocumented function
	{
		double type, sa, d, m, m0, gfw;

		// left parenthesis
		require(toklp, LINK);

		// first double arugument, type
		type = realfactor(LINK);
		require(tokcomma, LINK);

		// second double arugument, Sa
		sa = realfactor(LINK);
		require(tokcomma, LINK);

		// third double arugument, Sa
		d = realfactor(LINK);
		require(tokcomma, LINK);

		// fourth double arugument, m
		m = realfactor(LINK);
		require(tokcomma, LINK);

		// fifth double arugument, m0
		m0 = realfactor(LINK);
		require(tokcomma, LINK);

		// sixth double arugument, gfw
		gfw = realfactor(LINK);
		require(tokrp, LINK);

		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->sa_declercq(type, sa, d, m, m0, gfw);
	}
	break;

	case tokcallback: // PHAST function
	{
		double x1, x2;
		char* str;

		// left parenthesis
		require(toklp, LINK);

		// first double arugument
		x1 = realfactor(LINK);
		require(tokcomma, LINK);

		// second double arugument
		x2 = realfactor(LINK);
		require(tokcomma, LINK);

		// string arugument
		str = strexpr(LINK);

		require(tokrp, LINK);

		// call callback Basic function

		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->basic_callback(x1, x2, str);
		PhreeqcPtr->PHRQ_free(str);
	}
	break;

	case tokcell_pore_volume: // PHAST function
	case tokporevolume:
	{
		double x1 = (double)PhreeqcPtr->solution_number();
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->basic_callback(x1, x1, "cell_pore_volume");
	}
	break;

	case tokcell_porosity: // PHAST function
	{
		double x1 = (double)PhreeqcPtr->solution_number();
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->basic_callback(x1, x1, "cell_porosity");
	}
	break;

	case tokcell_saturation: // PHAST function
	{
		double x1 = (double)PhreeqcPtr->solution_number();
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->basic_callback(x1, x1, "cell_saturation");
	}
	break;

	case tokcell_volume: // PHAST function
	{
		double x1 = (double)PhreeqcPtr->solution_number();
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->basic_callback(x1, x1, "cell_volume");
	}
	break;

	case toktransport_cell_no:  // PHAST function
	{
		double x1 = (double)PhreeqcPtr->solution_number();
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->basic_callback(x1, x1, "transport_cell_no");
	}
	break;

	case tokvelocity_x: // PHAST function
	{
		double x1 = (double)PhreeqcPtr->solution_number();
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->basic_callback(x1, x1, "velocity_x");
	}
	break;

	case tokvelocity_y: // PHAST function
	{
		double x1 = (double)PhreeqcPtr->solution_number();
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->basic_callback(x1, x1, "velocity_y");
	}
	break;

	case tokvelocity_z: // PHAST function
	{
		double x1 = (double)PhreeqcPtr->solution_number();
		n.UU.val = (parse_all) ? 1 : PhreeqcPtr->basic_callback(x1, x1, "velocity_z");
	}
	break;

	case toklog10:
	{
		LDBLE t = realfactor(LINK);
		{
			n.UU.val = log10(t);
		}
	}
	break;

	case toksin:
	{
		n.UU.val = sin(realfactor(LINK));
	}
	break;

	case tokcos:
	{
		n.UU.val = cos(realfactor(LINK));
	}
	break;

	case toktan:
	{
		n.UU.val = realfactor(LINK);
		n.UU.val = sin(n.UU.val) / cos(n.UU.val);
	}
	break;

	case tokarctan:
	{
		n.UU.val = atan(realfactor(LINK));
	}
	break;

	case toklog:
	{
		n.UU.val = log(realfactor(LINK));
	}
	break;

	case tokexp:
	{
		n.UU.val = exp(realfactor(LINK));
	}
	break;

	case tokabs:
	{
		n.UU.val = fabs(realfactor(LINK));
	}
	break;

	case toksgn:
	{
		n.UU.val = realfactor(LINK);
		n.UU.val = (double)(n.UU.val > 0) - (double)(n.UU.val < 0);
	}
	break;

	case tokstr_:
		n.stringval = true;
		n.UU.sval = (char *) PhreeqcPtr->PHRQ_calloc(256, sizeof(char));
		if (n.UU.sval == NULL)
			PhreeqcPtr->malloc_error();
		numtostr(n.UU.sval, realfactor(LINK));
		break;
	case tokval:
		l_s = strfactor(LINK);
		tok1 = LINK->t;
		parse(l_s, &LINK->t);
		tok = LINK->t;
		if (tok == NULL)
			n.UU.val = 0.0;
		else
			n = expr(LINK);
		disposetokens(&tok);
		LINK->t = tok1;
		PhreeqcPtr->PHRQ_free(l_s);
		break;

	case tokchr_:
		n.stringval = true;
		n.UU.sval = (char *) PhreeqcPtr->PHRQ_calloc(256, sizeof(char));
		if (n.UU.sval == NULL)
			PhreeqcPtr->malloc_error();
		strcpy(n.UU.sval, " ");
		n.UU.sval[0] = (char) intfactor(LINK);
		break;

	case tokasc:
		l_s = strfactor(LINK);
		if (*l_s == '\0')
			n.UU.val = 0.0;
		else
			n.UU.val = l_s[0];
		PhreeqcPtr->PHRQ_free(l_s);
		break;

	case tokmid_:
		n.stringval = true;
		require(toklp, LINK);
		n.UU.sval = strexpr(LINK);
		require(tokcomma, LINK);
		i = intexpr(LINK);
		if (i < 1)
			i = 1;
		j = (int) strlen(n.UU.sval);
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			j = intexpr(LINK);
		}
		{
			std::string str = n.UU.sval;
			str = str.substr((size_t)i - 1, (size_t)j);
			strcpy(n.UU.sval, str.c_str());
		}
		require(tokrp, LINK);
		break;
	case toklen:
		l_s = strfactor(LINK);
		n.UU.val = (double) strlen(l_s);
		PhreeqcPtr->PHRQ_free(l_s);
		break;

	case tokpeek:
/* p2c: basic.p, line 1029: Note: Range checking is OFF [216] */
		if (parse_all)
		{
			intfactor(LINK);
			n.UU.val = 1.0;
		}
		else
		{
			trick.i = intfactor(LINK);
			n.UU.val = *trick.c;
		}
/* p2c: basic.p, line 1032: Note: Range checking is ON [216] */
		break;

	default:
		snerr(": missing \" or (");
		break;
	}
	return n;
}

valrec PBasic::
upexpr(struct LOC_exec * LINK)
{
	valrec n, n2;

	n = factor(LINK);
	while (LINK->t != NULL && LINK->t->kind == tokup)
	{
		if (n.stringval)
			tmerr(": not a number before ^");
		LINK->t = LINK->t->next;
		n2 = upexpr(LINK);
		if (n2.stringval)
			tmerr(": not a number after ^");
		if (n.UU.val >= 0)
		{
			if (n.UU.val > 0)
			{
				n.UU.val = exp(n2.UU.val * log(n.UU.val));
			}
			continue;
		}
		if (n2.UU.val != (long) n2.UU.val)
		{
			tmerr(": negative number cannot be raised to a fractional power.");
		} 
		else
		{
			n.UU.val = exp(n2.UU.val * log(-n.UU.val));
			if (((long) n2.UU.val) & 1)
				n.UU.val = -n.UU.val;
		}
	}
	return n;
}

valrec PBasic::
term(struct LOC_exec * LINK)
{
	valrec n, n2;

	int k;

	n = upexpr(LINK);
	while (LINK->t != NULL && (unsigned long) LINK->t->kind < 32 &&
		   ((1L << ((long) LINK->t->kind)) & ((1L << ((long) toktimes)) |
											  (1L << ((long) tokdiv)) | (1L <<
																		 ((long) tokmod)))) != 0)
	{
		k = LINK->t->kind;
		LINK->t = LINK->t->next;
		n2 = upexpr(LINK);
		if (n.stringval || n2.stringval)
			tmerr(": found char, but need a number for * or /");
		if (k == tokmod)
		{
			/*      n.UU.val = (long)floor(n.UU.val + 0.5) % (long)floor(n2.UU.val + 0.5); */
			if (n.UU.val != 0)
			{
				n.UU.val =
					fabs(n.UU.val) / n.UU.val * fmod(fabs(n.UU.val) +
													 1e-14, n2.UU.val);
			}
			else
			{
				n.UU.val = 0;
			}
/* p2c: basic.p, line 1078:
 * Note: Using % for possibly-negative arguments [317] */
		}
		else if (k == toktimes)
			n.UU.val *= n2.UU.val;
		else if (n2.UU.val != 0)
		{
			n.UU.val /= n2.UU.val;
		}
		else
		{
			if (!parse_all)
			{
				char * error_string = PhreeqcPtr->sformatf( "Zero divide in BASIC line\n %ld %s.\nValue set to zero.", stmtline->num, stmtline->inbuf);
				PhreeqcPtr->warning_msg(error_string);
			}
			n.UU.val = 0;
		}
	}
	return n;
}

valrec PBasic::
sexpr(struct LOC_exec * LINK)
{
	valrec n, n2;

	int k, m;

	n = term(LINK);
	while (LINK->t != NULL && (unsigned long) LINK->t->kind < 32 &&
		   ((1L << ((long) LINK->t->kind)) &
			((1L << ((long) tokplus)) | (1L << ((long) tokminus)))) != 0)
	{
		k = LINK->t->kind;
		LINK->t = LINK->t->next;
		n2 = term(LINK);
		if (n.stringval != n2.stringval)
			tmerr(": found char, but need a number for + or - ");
		if (k == tokplus)
		{
			if (n.stringval)
			{
				m = 1;
				if (n.UU.sval)
				{
					m += (int) strlen(n.UU.sval);
				}
				if (n2.UU.sval)
				{
					m += (int) strlen(n2.UU.sval);
				}
				//m = (int) strlen(n.UU.sval) + (int) strlen(n2.UU.sval) + 1;
				if (m < 256)
					m = 256;

				n.UU.sval = (char *) PhreeqcPtr->PHRQ_realloc(n.UU.sval, (size_t) m * sizeof(char));
				if (n.UU.sval == NULL)
				{
					PhreeqcPtr->malloc_error();
				}
				else
				{
					if (n2.UU.sval)
					{
						strcat(n.UU.sval, n2.UU.sval);
						PhreeqcPtr->PHRQ_free(n2.UU.sval);
					}
				}
			}
			else
				n.UU.val += n2.UU.val;
		}
		else
		{
			if (n.stringval)
				tmerr(": found char, but need a number for - ");
			else
				n.UU.val -= n2.UU.val;
		}
	}
	return n;
}

valrec PBasic::
relexpr(struct LOC_exec * LINK)
{
	valrec n, n2;

	bool f;
	int k;

	n = sexpr(LINK);
	while (LINK->t != NULL && (unsigned long) LINK->t->kind < 32 &&
		   ((1L << ((long) LINK->t->kind)) &
			((1L << ((long) tokne + 1)) - (1L << ((long) tokeq)))) != 0)
	{
		k = LINK->t->kind;
		LINK->t = LINK->t->next;
		n2 = sexpr(LINK);
		if (n.stringval != n2.stringval)
			tmerr("");
		if (n.stringval)
		{
			f = (bool) ((!strcmp(n.UU.sval, n2.UU.sval)
							&& (unsigned long) k < 32
							&& ((1L << ((long) k)) &
								((1L << ((long) tokeq)) |
								 (1L << ((long) tokge)) | (1L <<
														   ((long) tokle))))
							!= 0) || (strcmp(n.UU.sval, n2.UU.sval) < 0
									  && (unsigned long) k < 32
									  && ((1L << ((long) k)) &
										  ((1L << ((long) toklt)) |
										   (1L << ((long) tokle)) | (1L <<
																	 ((long)
																	  tokne))))
									  != 0)
						   || (strcmp(n.UU.sval, n2.UU.sval) > 0
							   && (unsigned long) k < 32
							   && ((1L << ((long) k)) &
								   ((1L << ((long) tokgt)) |
									(1L << ((long) tokge)) | (1L <<
															  ((long)
															   tokne)))) !=
							   0));
			/* p2c: basic.p, line 2175: Note:
			 * Line breaker spent 0.0+1.00 seconds, 5000 tries on line 1518 [251] */
			PhreeqcPtr->PHRQ_free(n.UU.sval);
			PhreeqcPtr->PHRQ_free(n2.UU.sval);
		}
		else
			f = (bool) ((n.UU.val == n2.UU.val && (unsigned long) k < 32 &&
							((1L << ((long) k)) & ((1L << ((long) tokeq)) |
												   (1L << ((long) tokge)) |
												   (1L << ((long) tokle)))) !=
							0) || (n.UU.val < n2.UU.val
								   && (unsigned long) k < 32
								   && ((1L << ((long) k)) &
									   ((1L << ((long) toklt)) |
										(1L << ((long) tokle)) | (1L <<
																  ((long)
																   tokne))))
								   != 0) || (n.UU.val > n2.UU.val
											 && (unsigned long) k < 32
											 && ((1L << ((long) k)) &
												 ((1L << ((long) tokgt)) |
												  (1L << ((long) tokge)) | (1L
																			<<
																			((long) tokne)))) != 0));
		/* p2c: basic.p, line 2175: Note:
		 * Line breaker spent 0.0+2.00 seconds, 5000 tries on line 1532 [251] */
		n.stringval = false;
		n.UU.val = f;
	}
	return n;
}

valrec PBasic::
andexpr(struct LOC_exec * LINK)
{
	valrec n, n2;

	n = relexpr(LINK);
	while (LINK->t != NULL && LINK->t->kind == tokand)
	{
		LINK->t = LINK->t->next;
		n2 = relexpr(LINK);
		if (n.stringval || n2.stringval)
			tmerr("");
		n.UU.val = ((long) n.UU.val) & ((long) n2.UU.val);
	}
	return n;
}

valrec PBasic::
expr(struct LOC_exec * LINK)
{
	valrec n, n2;

	int k;

	n = andexpr(LINK);
	while (LINK->t != NULL && (unsigned long) LINK->t->kind < 32 &&
		   ((1L << ((long) LINK->t->kind)) &
			((1L << ((long) tokor)) | (1L << ((long) tokxor)))) != 0)
	{
		k = LINK->t->kind;
		LINK->t = LINK->t->next;
		n2 = andexpr(LINK);
		if (n.stringval || n2.stringval)
			tmerr("");
		if (k == tokor)
			n.UU.val = ((long) n.UU.val) | ((long) n2.UU.val);
		else
			n.UU.val = ((long) n.UU.val) ^ ((long) n2.UU.val);
	}
	return n;
}

void PBasic::
checkextra(struct LOC_exec *LINK)
{
	if (LINK->t != NULL)
	{
		if (phreeqci_gui)
		{
			_ASSERTE(nIDErrPrompt == 0);
			nIDErrPrompt = IDS_ERR_EXTRA;
		}
		errormsg("Extra information on line");
	}
}

bool PBasic::
iseos(struct LOC_exec *LINK)
{
	return ((bool) (LINK->t == NULL || LINK->t->kind == (long) tokelse ||
					   LINK->t->kind == (long) tokcolon));
}

void PBasic::
skiptoeos(struct LOC_exec *LINK)
{
	while (!iseos(LINK))
		LINK->t = LINK->t->next;
}

linerec * PBasic::
findline(long n)
{
	linerec *l;

	l = linebase;
	while (l != NULL && l->num != n)
		l = l->next;
	return l;
}
linerec * PBasic::
mustfindline(long n)
{
	linerec *l;

	l = findline(n);
	if (phreeqci_gui)
	{
		if (parse_whole_program)
		{
			if (l == NULL) 
			{
				_ASSERTE(nIDErrPrompt == 0);
				nIDErrPrompt = IDS_ERR_UNDEF_LINE;
				char * error_string = PhreeqcPtr->sformatf( "Undefined line %ld", n);
				errormsg(error_string);
			}
		}
	}
	else
	{
		if (l == NULL) 
		{
			char * error_string = PhreeqcPtr->sformatf( "Undefined line %ld", n);
			errormsg(error_string);
		}
	}
	return l;
}
void PBasic::
cmdend(struct LOC_exec *LINK)
{
	stmtline = NULL;
	LINK->t = NULL;
}

void PBasic::
cmdnew(struct LOC_exec *LINK)
{
	void *p;
	int i, k;

	cmdend(LINK);
	clearloops();
	restoredata();
	while (linebase != NULL)
	{
		p = linebase->next;
		disposetokens(&linebase->txt);
		PhreeqcPtr->PHRQ_free(linebase);
		linebase = (linerec *) p;
	}
	while (varbase != NULL)
	{
		p = varbase->next;
		if (varbase->stringvar)
		{
			if (varbase->numdims > 0)
			{
				k = 1;
				for (i = 0; i < varbase->numdims; i++)
				{
					k = k * (varbase->dims[i]);
				}
				for (i = 0; i < k; i++)
				{
					PhreeqcPtr->free_check_null(varbase->UU.U1.sarr[i]);
				}
				PhreeqcPtr->free_check_null(varbase->UU.U1.sarr);
			}
			else if (*varbase->UU.U1.sval != NULL)
			{
				*varbase->UU.U1.sval =
					(char *) PhreeqcPtr->free_check_null(*varbase->UU.U1.sval);
			}

		}
		else
		{
			PhreeqcPtr->free_check_null(varbase->UU.U0.arr);
			varbase->UU.U0.arr = NULL;
		}
		PhreeqcPtr->PHRQ_free(varbase);
		varbase = (varrec *) p;
	}
}

void PBasic::
cmdlist(struct LOC_exec *LINK)
{
	linerec *l;
	long n1, n2;

	do
	{
		n1 = 0;
		n2 = LONG_MAX;
		if (LINK->t != NULL && LINK->t->kind == toknum)
		{
			n1 = (long) LINK->t->UU.num;
			LINK->t = LINK->t->next;
			if (LINK->t == NULL || LINK->t->kind != tokminus)
				n2 = n1;
		}
		if (LINK->t != NULL && LINK->t->kind == tokminus)
		{
			LINK->t = LINK->t->next;
			if (LINK->t != NULL && LINK->t->kind == toknum)
			{
				n2 = (long) LINK->t->UU.num;
				LINK->t = LINK->t->next;
			}
			else
				n2 = LONG_MAX;
		}
		l = linebase;
		while (l != NULL && l->num <= n2)
		{
			if (l->num >= n1)
			{
				/* printf("%ld ", l->num); */
				/*	listtokens(stdout, l->txt); */
				/* putchar('\n'); */
				output_msg(PhreeqcPtr->sformatf("%ld ", l->num));
				listtokens(NULL, l->txt);
				output_msg("\n");
			}
			l = l->next;
		}
		if (!iseos(LINK))
			require(tokcomma, LINK);
	}
	while (!iseos(LINK));
}

void PBasic::
cmdload(bool merging, char * name, struct LOC_exec *LINK)
{
	FILE *f;
	tokenrec *l_buf;
	char STR1[256] = {0};
	char *TEMP;

	f = NULL;
	if (!merging)
		cmdnew(LINK);
	if (f != NULL)
	{
		snprintf(STR1, sizeof(STR1), "%s.TEXT", name);
		f = freopen(STR1, "r", f);
	}
	else
	{
		snprintf(STR1, sizeof(STR1), "%s.TEXT", name);
		f = fopen(STR1, "r");
	}
	if (f == NULL)
	{
		_EscIO(FileNotFound);
		return;
	}
	while (fgets(inbuf, 256, f) != NULL)
	{
		TEMP = strchr(inbuf, '\n');
		if (TEMP != NULL)
			*TEMP = 0;
		parseinput(&l_buf);
		if (curline == 0)
		{
			/*      printf("Bad line in file\n"); */
			output_msg("Bad line in file\n");
			disposetokens(&l_buf);
		}
	}
	if (f != NULL)
		fclose(f);
	f = NULL;
	if (f != NULL)
		fclose(f);
}

void PBasic::
cmdrun(struct LOC_exec *LINK)
{
	linerec *l;
	long i;
	char *l_s;

	l_s = (char *) PhreeqcPtr->PHRQ_calloc(PhreeqcPtr->max_line, sizeof(char));
	if (l_s == NULL)
		PhreeqcPtr->malloc_error();

	l = linebase;
	if (!iseos(LINK))
	{
		if (LINK->t->kind == toknum)
			/*l = mustfindline(intexpr(LINK), LINK); */
			l = mustfindline(intexpr(LINK));
		else
		{
			stringexpr(l_s, LINK);
			i = 0;
			if (!iseos(LINK))
			{
				require(tokcomma, LINK);
				i = intexpr(LINK);
			}
			checkextra(LINK);
			cmdload(false, l_s, LINK);
			if (i == 0)
				l = linebase;
			else
			l = mustfindline(i);
		}
	}
	stmtline = l;
	LINK->gotoflag = true;
	clearvars();
	clearloops();
	restoredata();
	PhreeqcPtr->free_check_null(l_s);
	return;
}

/* PhreeqcPtr->replace basic save command with transport of rate back to calc_kinetic_rate */
void PBasic::
cmdsave(struct LOC_exec *LINK)
{
	valrec n;

	while (!iseos(LINK))
	{
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);
		if (n.stringval)
		{
			snerr(": in SAVE command");
		}
		else
		{
			PhreeqcPtr->rate_moles = n.UU.val;
		}
	}
}

void PBasic::
cmdput(struct LOC_exec *LINK)
{
	int j;
	std::ostringstream oss;

	/* get parentheses */
	require(toklp, LINK);

	/* get first argument */
	double value = realexpr(LINK);

	for (;;)
	{
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			j = intexpr(LINK);
			oss << j << ",";
		}
		else
		{
			/* get right parentheses */
			require(tokrp, LINK);
			break;
		}
	}
	if (!parse_all)
	{
		PhreeqcPtr->save_values[oss.str()] = value;
	}
}

void PBasic::
cmdput_(struct LOC_exec* LINK)
{
	int j;
	std::ostringstream oss;

	/* get parentheses */
	require(toklp, LINK);

	/* get first argument */
	char* str = strexpr(LINK);
	std::string s_value = str;
	PhreeqcPtr->PHRQ_free(str);

	for (;;)
	{
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			j = intexpr(LINK);
			oss << j << ",";
		}
		else
		{
			/* get right parentheses */
			require(tokrp, LINK);
			break;
		}
	}
	if (!parse_all)
	{
		PhreeqcPtr->save_strings[oss.str()] = s_value;
	}
}
void PBasic::
cmdchange_por(struct LOC_exec *LINK)
{
	int j;
	LDBLE TEMP;
	require(toklp, LINK);
	/* get new porosity */
	TEMP = realexpr(LINK);
	require(tokcomma, LINK);
	/* get cell_no */
	j = intexpr(LINK);
	require(tokrp, LINK);
	if (j > 0 && j <= PhreeqcPtr->count_cells * (1 + PhreeqcPtr->stag_data.count_stag) + 1
		&& j != PhreeqcPtr->count_cells + 1)
		PhreeqcPtr->cell_data[j].por = TEMP;
}

void PBasic::
cmdchange_surf(struct LOC_exec *LINK)
{
/*
    change_surf("Hfo",    0.3,      "Sfo",      0,      5      )
	       (old_name, fraction, new_name, new_Dw, cell_no)
 */
	char *c1;
	int count;

	PhreeqcPtr->change_surf_count += 1;
	count = PhreeqcPtr->change_surf_count;
	if (PhreeqcPtr->change_surf[count - 1].next == FALSE)
		PhreeqcPtr->change_surf = PhreeqcPtr->change_surf_alloc(count + 1);

	require(toklp, LINK);
	/* get surface component name (change affects all comps of the same charge structure) */
	c1 = strexpr(LINK);
	PhreeqcPtr->change_surf[count - 1].comp_name = PhreeqcPtr->string_hsave(c1);
	PhreeqcPtr->PHRQ_free(c1);
	require(tokcomma, LINK);
	/* get fraction of comp to change */
	PhreeqcPtr->change_surf[count - 1].fraction = realexpr(LINK);
	require(tokcomma, LINK);
	/* get new surface component name */
	c1 = strexpr(LINK);
	PhreeqcPtr->change_surf[count - 1].new_comp_name = PhreeqcPtr->string_hsave(c1);
	PhreeqcPtr->PHRQ_free(c1);
	require(tokcomma, LINK);
	/* get new Dw (no transport if 0) */
	PhreeqcPtr->change_surf[count - 1].new_Dw = realexpr(LINK);
	require(tokcomma, LINK);
	/* get cell_no */
	PhreeqcPtr->change_surf[count - 1].cell_no = intexpr(LINK);
	require(tokrp, LINK);

	if (PhreeqcPtr->change_surf->cell_no == 0 || PhreeqcPtr->change_surf->cell_no == PhreeqcPtr->count_cells + 1)
		PhreeqcPtr->change_surf[count - 1].cell_no = -99;
}

void PBasic::
cmdbye(void)
{
	exitflag = true;
}

void PBasic::
cmddel(struct LOC_exec *LINK)
{
	linerec *l, *l0, *l1;
	long n1, n2;

	do
	{
		if (iseos(LINK))
			snerr(": no variable name after del");
		n1 = 0;
		n2 = LONG_MAX;
		if (LINK->t != NULL && LINK->t->kind == toknum)
		{
			n1 = (long) LINK->t->UU.num;
			LINK->t = LINK->t->next;
			if (LINK->t == NULL || LINK->t->kind != tokminus)
				n2 = n1;
		}
		if (LINK->t != NULL && LINK->t->kind == tokminus)
		{
			LINK->t = LINK->t->next;
			if (LINK->t != NULL && LINK->t->kind == toknum)
			{
				n2 = (long) LINK->t->UU.num;
				LINK->t = LINK->t->next;
			}
			else
				n2 = LONG_MAX;
		}
		l = linebase;
		l0 = NULL;
		while (l != NULL && l->num <= n2)
		{
			l1 = l->next;
			if (l->num >= n1)
			{
				if (l == stmtline)
				{
					cmdend(LINK);
					clearloops();
					restoredata();
				}
				if (l0 == NULL)
					linebase = l->next;
				else
					l0->next = l->next;
				disposetokens(&l->txt);
				PhreeqcPtr->PHRQ_free(l);
			}
			else
				l0 = l;
			l = l1;
		}
		if (!iseos(LINK))
			require(tokcomma, LINK);
	}
	while (!iseos(LINK));
}

void PBasic::
cmdrenum(struct LOC_exec *LINK)
{
	linerec *l, *l1;
	tokenrec *tok;
	long lnum, step;

	lnum = 10;
	step = 10;
	if (!iseos(LINK))
	{
		lnum = intexpr(LINK);
		if (!iseos(LINK))
		{
			require(tokcomma, LINK);
			step = intexpr(LINK);
		}
	}
	l = linebase;
	if (l == NULL)
		return;
	while (l != NULL)
	{
		l->num2 = lnum;
		lnum += step;
		l = l->next;
	}
	l = linebase;
	do
	{
		tok = l->txt;
		do
		{
			if (tok->kind == (long) tokdel || tok->kind == (long) tokrestore
				|| tok->kind == (long) toklist || tok->kind == (long) tokrun
				|| tok->kind == (long) tokelse || tok->kind == (long) tokthen
				|| tok->kind == (long) tokgosub
				|| tok->kind == (long) tokgoto)
			{
				while (tok->next != NULL && tok->next->kind == toknum)
				{
					tok = tok->next;
					lnum = (long) floor(tok->UU.num + 0.5);
					l1 = linebase;
					while (l1 != NULL && l1->num != lnum)
						l1 = l1->next;
					if (l1 == NULL)
					{
						/*      printf("Undefined line %ld in line %ld\n", lnum,
							   l->num2); */
						output_msg(PhreeqcPtr->sformatf("Undefined line %ld in line %ld\n", lnum, l->num2));
					}
					else
					{
#if defined(PHREEQCI_GUI)
						if (phreeqci_gui)
						{
							_snprintf(tok->sz_num, tok->n_sz, "%ld", l1->num2);
						}
#endif
						tok->UU.num = l1->num2;
					}
					if (tok->next != NULL && tok->next->kind == tokcomma)
						tok = tok->next;
				}
			}
			tok = tok->next;
		}
		while (tok != NULL);
		l = l->next;
	}
	while (l != NULL);
	l = linebase;
	while (l != NULL)
	{
		l->num = l->num2;
		l = l->next;
	}
}

void PBasic::
cmdprint(struct LOC_exec *LINK)
{
	bool semiflag;
	valrec n;

	char STR1[256] = {0};

	semiflag = false;
	while (!iseos(LINK))
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);
		if (n.stringval)
		{
			if (!skip_punch) {
				/*      fputs(n.UU.sval, stdout); */
				output_msg(PhreeqcPtr->sformatf("%s ", n.UU.sval));
			}
			n.UU.sval = (char*)PhreeqcPtr->free_check_null(n.UU.sval);
		}
		else
/*      printf("%s ", numtostr(STR1, n.UU.val)); */
			output_msg(PhreeqcPtr->sformatf("%s ", numtostr(STR1, n.UU.val)));
	}
	if (!semiflag && PhreeqcPtr->Get_output_newline())
/*    putchar('\n');*/
		output_msg("\n");
	PhreeqcPtr->Set_output_newline(true);
	skip_punch = false;
}

void PBasic::
cmdpunch(struct LOC_exec *LINK)
{
	valrec n;

	/*  char STR1[256]; */

	while (!iseos(LINK))
	{
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);
		bool temp_high_precision = (PhreeqcPtr->current_selected_output != NULL) ?
			PhreeqcPtr->current_selected_output->Get_high_precision() :
			PhreeqcPtr->high_precision;
		if (!this->skip_punch)
		{
			if (n.stringval)
			{
				/*      fputs(n.UU.sval, stdout); */
				{
					if (!temp_high_precision)
					{
						if (strlen(n.UU.sval) <= 12)
						{
							if (punch_tab)
							{
								PhreeqcPtr->fpunchf_user(PhreeqcPtr->n_user_punch_index, "%12.12s\t", n.UU.sval);
							}
							else {
								PhreeqcPtr->fpunchf_user(PhreeqcPtr->n_user_punch_index, "%12.12s", n.UU.sval);
							}
						}
						else
						{
							if (punch_tab)
							{
								PhreeqcPtr->fpunchf_user(PhreeqcPtr->n_user_punch_index, "%s\t", n.UU.sval);
							}
							else {
								PhreeqcPtr->fpunchf_user(PhreeqcPtr->n_user_punch_index, "%s", n.UU.sval);
							}
						}
					}
					else
					{
						if (strlen(n.UU.sval) <= 20)
						{
							if (punch_tab) {
								PhreeqcPtr->fpunchf_user(PhreeqcPtr->n_user_punch_index, "%20.20s\t", n.UU.sval);
							}
							else {
								PhreeqcPtr->fpunchf_user(PhreeqcPtr->n_user_punch_index, "%20.20s", n.UU.sval);
							}
						}
						else
						{
							if (punch_tab) {
								PhreeqcPtr->fpunchf_user(PhreeqcPtr->n_user_punch_index, "%s\t", n.UU.sval);
							}
							else {
								PhreeqcPtr->fpunchf_user(PhreeqcPtr->n_user_punch_index, "%s", n.UU.sval);
							}
						}
					}
				}
				n.UU.sval = (char*)PhreeqcPtr->free_check_null(n.UU.sval);
			}
			else if (!temp_high_precision)
			{
				PhreeqcPtr->fpunchf_user(PhreeqcPtr->n_user_punch_index, "%12.4e\t", (double)n.UU.val);
			}
			else
			{
				PhreeqcPtr->fpunchf_user(PhreeqcPtr->n_user_punch_index, "%20.12e\t", (double)n.UU.val);
			}
			punch_tab = true;
			++PhreeqcPtr->n_user_punch_index;
		}
		else
		{
			n.UU.sval = (char*)PhreeqcPtr->free_check_null(n.UU.sval);
		}
		this->skip_punch = false;
	}
}

#if defined PHREEQ98 
void PBasic::
cmdgraph_x(struct LOC_exec *LINK)
{
	bool semiflag;
	valrec n;

	char STR1[256];
	semiflag = false;
	while (!iseos(LINK))
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);
		if (colnr == 0)
		{
			rownr++;
		}
		if (n.stringval)
		{
/*      fputs(n.UU.sval, stdout); */
			GridChar(n.UU.sval, "x");
			PhreeqcPtr->PHRQ_free(n.UU.sval);
		}
		else
			GridChar(numtostr(STR1, n.UU.val), "x");
		colnr++;
	}
}

void PBasic::
cmdgraph_y(struct LOC_exec *LINK)
{
	bool semiflag;
	valrec n;

	char STR1[256];
	semiflag = false;
	while (!iseos(LINK))
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);
		if (colnr == 0)
		{
			rownr++;
		}
		if (n.stringval)
		{
/*      fputs(n.UU.sval, stdout); */
			GridChar(n.UU.sval, "y");
			PhreeqcPtr->PHRQ_free(n.UU.sval);
		}
		else
			GridChar(numtostr(STR1, n.UU.val), "y");
		colnr++;
	}
}

void PBasic::
cmdgraph_sy(struct LOC_exec *LINK)
{
	bool semiflag;
	valrec n;

	char STR1[256];
	semiflag = false;
	while (!iseos(LINK))
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);
		if (colnr == 0)
		{
			rownr++;
		}
		if (n.stringval)
		{
/*      fputs(n.UU.sval, stdout); */
			GridChar(n.UU.sval, "s");
			PhreeqcPtr->PHRQ_free(n.UU.sval);
		}
		else
			GridChar(numtostr(STR1, n.UU.val), "s");
		colnr++;
	}
}
#endif

void PBasic::
cmdlet(bool implied, struct LOC_exec *LINK)
{
	varrec *v;
	char *old, *mynew;
	LDBLE d_value;
	LDBLE *target;
	char **starget;
	target = NULL;
	starget = NULL;
	if (implied)
		LINK->t = stmttok;
	v = findvar(LINK);
	if (v->stringvar)
	{
		starget = v->UU.U1.sval;
	}
	else
	{
		target = v->UU.U0.val;
	}
	require(tokeq, LINK);
	if (!v->stringvar)
	{
		/* in case array is used on right hand side */
		d_value = realexpr(LINK);
		v->UU.U0.val = target;
		*v->UU.U0.val = d_value;
		/*  *v->UU.U0.val = realexpr(LINK); */
		return;
	}
	mynew = strexpr(LINK);
	v->UU.U1.sval = starget;
	old = *v->UU.U1.sval;
	*v->UU.U1.sval = mynew;
	if (old != NULL)
		PhreeqcPtr->PHRQ_free(old);
}

void PBasic::
cmdgoto(struct LOC_exec *LINK)
{
	stmtline = mustfindline(intexpr(LINK));
	LINK->t = NULL;
	LINK->gotoflag = true;
}

void PBasic::
cmdif(struct LOC_exec *LINK)
{
	LDBLE n;
	long i;

	n = realexpr(LINK);
	require(tokthen, LINK);
	if (n == 0)
	{
		i = 0;
		do
		{
			if (LINK->t != NULL)
			{
				if (LINK->t->kind == tokif)
					i++;
				if (LINK->t->kind == tokelse)
					i--;
				LINK->t = LINK->t->next;
			}
		}
		while (LINK->t != NULL && i >= 0);
	}
	if (LINK->t != NULL && LINK->t->kind == toknum)
		cmdgoto(LINK);
	else
		LINK->elseflag = true;
}

void PBasic::
cmdelse(struct LOC_exec *LINK)
{
	LINK->t = NULL;
}

bool PBasic::
skiploop(int up, int dn, struct LOC_exec *LINK)
{
	bool Result;
	long i;
	linerec *saveline;

	saveline = stmtline;
	i = 0;
	do
	{
		while (LINK->t == NULL)
		{
			if (stmtline == NULL || stmtline->next == NULL)
			{
				Result = false;
				stmtline = saveline;
				goto _L1;
			}
			stmtline = stmtline->next;
			LINK->t = stmtline->txt;
		}
		if (LINK->t->kind == up)
			i++;
		if (LINK->t->kind == dn)
			i--;
		LINK->t = LINK->t->next;
	}
	while (i >= 0);
	Result = true;
  _L1:
	return Result;
}

void PBasic::
cmdfor(struct LOC_exec *LINK)
{
	looprec *l, lr;
	linerec *saveline;
	long i, j;

	lr.UU.U0.vp = findvar(LINK);
	if (lr.UU.U0.vp->stringvar)
		snerr(": error in FOR command");
	require(tokeq, LINK);
	*lr.UU.U0.vp->UU.U0.val = realexpr(LINK);
	require(tokto, LINK);
	lr.UU.U0.max = realexpr(LINK);
	if (LINK->t != NULL && LINK->t->kind == tokstep)
	{
		LINK->t = LINK->t->next;
		lr.UU.U0.step = realexpr(LINK);
	}
	else
		lr.UU.U0.step = 1.0;
	lr.homeline = stmtline;
	lr.hometok = LINK->t;
	lr.kind = forloop;
	lr.next = loopbase;
	if ((lr.UU.U0.step >= 0 && *lr.UU.U0.vp->UU.U0.val > lr.UU.U0.max) ||
		(lr.UU.U0.step <= 0 && *lr.UU.U0.vp->UU.U0.val < lr.UU.U0.max))
	{
		saveline = stmtline;
		i = 0;
		j = 0;
		do
		{
			while (LINK->t == NULL)
			{
				if (stmtline == NULL || stmtline->next == NULL)
				{
					stmtline = saveline;
					if (phreeqci_gui)
					{
						_ASSERTE(nIDErrPrompt == 0);
						nIDErrPrompt = IDS_ERR_FOR_WO_NEXT;
					}
					errormsg("FOR without NEXT");
				}
				stmtline = stmtline->next;
				LINK->t = stmtline->txt;
			}
			if (LINK->t->kind == tokfor)
			{
				if (LINK->t->next != NULL && LINK->t->next->kind == tokvar &&
					LINK->t->next->UU.vp == lr.UU.U0.vp)
					j++;
				else
					i++;
			}
			if (LINK->t->kind == toknext)
			{
				if (LINK->t->next != NULL && LINK->t->next->kind == tokvar &&
					LINK->t->next->UU.vp == lr.UU.U0.vp)
					j--;
				else
					i--;
			}
			LINK->t = LINK->t->next;
		}
		while (i >= 0 && j >= 0);
		skiptoeos(LINK);
		return;
	}
	l = (looprec *) PhreeqcPtr->PHRQ_calloc(1, sizeof(looprec));
	if (l == NULL)
	{
		PhreeqcPtr->malloc_error();
	}
	else
	{
		*l = lr;
		loopbase = l;
	}
}

void PBasic::
cmdnext(struct LOC_exec *LINK)
{
	varrec *v;
	bool found;
	looprec *l, *WITH;

	if (!iseos(LINK))
		v = findvar(LINK);
	else
		v = NULL;
	do
	{
		if (loopbase == NULL || loopbase->kind == gosubloop)
		{
			if (phreeqci_gui)
			{
				_ASSERTE(nIDErrPrompt == 0);
				nIDErrPrompt = IDS_ERR_NEXT_WO_FOR;
			}
			errormsg("NEXT without FOR");
		}
		found = (bool) (loopbase->kind == forloop &&
						   (v == NULL || loopbase->UU.U0.vp == v));
		if (!found)
		{
			l = loopbase->next;
			PhreeqcPtr->PHRQ_free(loopbase);
			loopbase = l;
		}
	}
	while (!found);
	WITH = loopbase;
	*WITH->UU.U0.vp->UU.U0.val += WITH->UU.U0.step;
	if ((WITH->UU.U0.step < 0
		 || *WITH->UU.U0.vp->UU.U0.val <= WITH->UU.U0.max)
		&& (WITH->UU.U0.step > 0
			|| *WITH->UU.U0.vp->UU.U0.val >= WITH->UU.U0.max))
	{
		stmtline = WITH->homeline;
		LINK->t = WITH->hometok;
		return;
	}
	l = loopbase->next;
	PhreeqcPtr->PHRQ_free(loopbase);
	loopbase = l;
}

void PBasic::
cmdwhile(struct LOC_exec *LINK)
{
	looprec *l;

	l = (looprec *) PhreeqcPtr->PHRQ_calloc(1, sizeof(looprec));
	if (l == NULL)
	{
		PhreeqcPtr->malloc_error();
		return;
	}
	l->next = loopbase;
	loopbase = l;
	l->kind = whileloop;
	l->homeline = stmtline;
	l->hometok = LINK->t;
	if (iseos(LINK))
		return;
	if (realexpr(LINK) != 0)
		return;
	if (phreeqci_gui)
	{
		if (parse_whole_program == TRUE)
		{
			if (!skiploop(tokwhile, tokwend, LINK))
			{
				_ASSERTE(nIDErrPrompt == 0);
				nIDErrPrompt = IDS_ERR_WHILE_WO_WEND;
				errormsg("WHILE without WEND");
			}
			l = loopbase->next;
			PhreeqcPtr->PHRQ_free(loopbase);
			loopbase = l;
			skiptoeos(LINK);
		}
	}
	else
	{
		if (!skiploop(tokwhile, tokwend, LINK))
		{
			errormsg("WHILE without WEND");
		}
		l = loopbase->next;
		PhreeqcPtr->PHRQ_free(loopbase);
		loopbase = l;
		skiptoeos(LINK);
	}
}

void PBasic::
cmdwend(struct LOC_exec *LINK)
{
	tokenrec *tok;
	linerec *tokline;
	looprec *l;
	bool found;
	if (phreeqci_gui && !parse_whole_program)
	{
		return;
	}
	do
	{
		if (loopbase == NULL || loopbase->kind == gosubloop)
		{
			if (phreeqci_gui)
			{
				_ASSERTE(nIDErrPrompt == 0);
				nIDErrPrompt = IDS_ERR_WEND_WO_WHILE;
			}
			errormsg("WEND without WHILE");
		}
		found = (bool) (loopbase->kind == whileloop);
		if (!found)
		{
			l = loopbase->next;
			PhreeqcPtr->PHRQ_free(loopbase);
			loopbase = l;
		}
	}
	while (!found);
	if (!iseos(LINK))
	{
		if (realexpr(LINK) != 0)
			found = false;
	}
	tok = LINK->t;
	tokline = stmtline;
	if (found)
	{
		stmtline = loopbase->homeline;
		LINK->t = loopbase->hometok;
		if (!iseos(LINK))
		{
			if (realexpr(LINK) == 0)
				found = false;
		}
	}
	if (found)
		return;
	LINK->t = tok;
	stmtline = tokline;
	l = loopbase->next;
	PhreeqcPtr->PHRQ_free(loopbase);
	loopbase = l;
}

void PBasic::
cmdgosub(struct LOC_exec *LINK)
{
	looprec *l;

	l = (looprec *) PhreeqcPtr->PHRQ_calloc(1, sizeof(looprec));
	if (l == NULL)
	{
		PhreeqcPtr->malloc_error();
		return;
	}
	l->next = loopbase;
	loopbase = l;
	l->kind = gosubloop;
	l->homeline = stmtline;
	l->hometok = LINK->t;
	cmdgoto(LINK);
}

void PBasic::
cmdreturn(struct LOC_exec *LINK)
{
	looprec *l;
	bool found;

	if (phreeqci_gui && !parse_whole_program)
	{
		return;
	}

	do
	{
		if (loopbase == NULL)
		{
			if (phreeqci_gui)
			{
				_ASSERTE(nIDErrPrompt == 0);
				nIDErrPrompt = IDS_ERR_RETURN_WO_GOSUB;
			}
			errormsg("RETURN without GOSUB");
		}
		found = (bool) (loopbase->kind == gosubloop);
		if (!found)
		{
			l = loopbase->next;
			PhreeqcPtr->PHRQ_free(loopbase);
			loopbase = l;
		}
	}
	while (!found);
	stmtline = loopbase->homeline;
	LINK->t = loopbase->hometok;
	l = loopbase->next;
	PhreeqcPtr->PHRQ_free(loopbase);
	loopbase = l;
	skiptoeos(LINK);
}

void PBasic::
cmdread(struct LOC_exec *LINK)
{
	varrec *v;
	tokenrec *tok;
	bool found;

	do
	{
		v = findvar(LINK);
		tok = LINK->t;
		LINK->t = datatok;
		if (phreeqci_gui) 
		{
			if (parse_whole_program)
			{
				if (dataline == NULL)
				{
					dataline = linebase;
					LINK->t = dataline->txt;
				}
				if (LINK->t == NULL || LINK->t->kind != tokcomma)
				{
					do
					{
						while (LINK->t == NULL)
						{
							if (dataline == NULL || dataline->next == NULL)
							{
								_ASSERTE(nIDErrPrompt == 0);
								nIDErrPrompt = IDS_ERR_OUT_OF_DATA;
								errormsg("Out of Data");
							}
							dataline = dataline->next;
							LINK->t = dataline->txt;
						}
						found = (bool) (LINK->t->kind == tokdata);
						LINK->t = LINK->t->next;
					}
					while (!found || iseos(LINK));
				}
				else
					LINK->t = LINK->t->next;
				if (v->stringvar)
				{
					if (*v->UU.U1.sval != NULL)
						*v->UU.U1.sval = (char *) PhreeqcPtr->free_check_null(*v->UU.U1.sval);
					*v->UU.U1.sval = strexpr(LINK);
				}
				else
					*v->UU.U0.val = realexpr(LINK);
			}
		}
		else
		{
			if (dataline == NULL)
			{
				dataline = linebase;
				LINK->t = dataline->txt;
			}
			if (LINK->t == NULL || LINK->t->kind != tokcomma)
			{
				do
				{
					while (LINK->t == NULL)
					{
						if (dataline == NULL || dataline->next == NULL)
							errormsg("Out of Data");
						dataline = dataline->next;
						LINK->t = dataline->txt;
					}
					found = (bool) (LINK->t->kind == tokdata);
					LINK->t = LINK->t->next;
				}
				while (!found || iseos(LINK));
			}
			else
				LINK->t = LINK->t->next;
			if (v->stringvar)
			{
				if (*v->UU.U1.sval != NULL)
					*v->UU.U1.sval = (char *) PhreeqcPtr->free_check_null(*v->UU.U1.sval);
				*v->UU.U1.sval = strexpr(LINK);
			}
			else
				*v->UU.U0.val = realexpr(LINK);
		}
		datatok = LINK->t;
		LINK->t = tok;
		if (!iseos(LINK))
			require(tokcomma, LINK);

	}
	while (!iseos(LINK));
}

void PBasic::
cmddata(struct LOC_exec *LINK)
{
	skiptoeos(LINK);
}

void PBasic::
cmdrestore(struct LOC_exec *LINK)
{
	if (iseos(LINK))
		restoredata();
	else
	{
		dataline = mustfindline(intexpr(LINK));
		if (phreeqci_gui)
		{
			if (parse_whole_program)
			{
				datatok = dataline->txt;
			}
		}
		else
		{
			datatok = dataline->txt;
		}
	}
}

void PBasic::
cmdgotoxy(struct LOC_exec *LINK)
{
	intexpr(LINK);
	require(tokcomma, LINK);
}

void PBasic::
cmdon(struct LOC_exec *LINK)
{
	long i;
	looprec *l;

	i = intexpr(LINK);
	if (LINK->t != NULL && LINK->t->kind == tokgosub)
	{
		l = (looprec *) PhreeqcPtr->PHRQ_calloc(1, sizeof(looprec));
		if (l == NULL)
		{
			PhreeqcPtr->malloc_error();
		}
		else
		{
			l->next = loopbase;
			loopbase = l;
			l->kind = gosubloop;
			l->homeline = stmtline;
			l->hometok = LINK->t;
			LINK->t = LINK->t->next;
		}
	}
	else
		require(tokgoto, LINK);
	if (i < 1)
	{
		skiptoeos(LINK);
		return;
	}
	while (i > 1 && !iseos(LINK))
	{
		require(toknum, LINK);
		if (!iseos(LINK))
			require(tokcomma, LINK);
		i--;
	}
	if (!iseos(LINK))
		cmdgoto(LINK);
}

void PBasic::
cmddim(struct LOC_exec *LINK)
{
	long i, j, k;
	varrec *v;
	bool done;

	do
	{
		if (LINK->t == NULL || LINK->t->kind != tokvar)
			snerr(": error in DIM command");
		v = LINK->t->UU.vp;
		LINK->t = LINK->t->next;
		if (v->numdims != 0)
		{
			if (phreeqci_gui)
			{
				_ASSERTE(nIDErrPrompt == 0);
				nIDErrPrompt = IDS_ERR_ARRAY_ALREADY;
			}
			errormsg("Array already dimensioned before");
		}
		j = 1;
		i = 0;
		require(toklp, LINK);
		do
		{
			k = intexpr(LINK) + 1;
			if (k < 1)
				badsubscr();
			if (i >= maxdims)
				badsubscr();
			i++;
			v->dims[i - 1] = k;
			j *= k;
			done = (bool) (LINK->t != NULL && LINK->t->kind == tokrp);
			if (!done)
				require(tokcomma, LINK);
		}
		while (!done);
		LINK->t = LINK->t->next;
		v->numdims = (char) i;
		if (v->stringvar)
		{
			v->UU.U1.sarr = (char **) PhreeqcPtr->PHRQ_malloc(j * sizeof(char *));
			if (!v->UU.U1.sarr)
			{
				PhreeqcPtr->malloc_error();
#if !defined(R_SO)
				exit(4);
#endif
			}
			if (v->UU.U1.sarr == NULL)
				PhreeqcPtr->malloc_error();
			for (i = 0; i < j; i++)
				v->UU.U1.sarr[i] = NULL;
		}
		else
		{
			v->UU.U0.arr = (LDBLE *) PhreeqcPtr->PHRQ_malloc(j * sizeof(LDBLE));
			if (v->UU.U0.arr == NULL)
			{
				PhreeqcPtr->malloc_error();
			}
			else
			{
				for (i = 0; i < j; i++)
					v->UU.U0.arr[i] = 0.0;
			}
		}
		if (!iseos(LINK))
			require(tokcomma, LINK);
	}
	while (!iseos(LINK));
}

void PBasic::
cmderase(struct LOC_exec *LINK)
{
	varrec *v = NULL;
	do
	{
		if (LINK->t == NULL || LINK->t->kind != tokvar)
		{
			snerr(": error in DIM command");
		}
		else
		{
			v = LINK->t->UU.vp;
			LINK->t = LINK->t->next;
			clearvar(v);
			if (!iseos(LINK)) require(tokcomma, LINK);
		}
	}
	while (!iseos(LINK));
}


void PBasic::
cmdpoke(struct LOC_exec *LINK)
{
	union
	{
		long i;
		char *c;
	} trick;

/* p2c: basic.p, line 2073: Note: Range checking is OFF [216] */
	trick.i = intexpr(LINK);
	require(tokcomma, LINK);
	*trick.c = (char) intexpr(LINK);
/* p2c: basic.p, line 2077: Note: Range checking is ON [216] */
}

void PBasic::
exec(void)
{
	struct LOC_exec V;
	V.gotoflag = false;
	V.elseflag = false;
	V.t = NULL;
	char STR1[256] = {0};

	try
	{
		do
		{
			do
			{
				V.gotoflag = false;
				V.elseflag = false;
				while (stmttok != NULL && stmttok->kind == tokcolon)
					stmttok = stmttok->next;
				V.t = stmttok;
				if (V.t != NULL)
				{
					V.t = V.t->next;
#if defined(PHREEQCI_GUI)
					if (phreeqci_gui)
					{
						if (WaitForSingleObject(hInfiniteLoop, 0) == WAIT_OBJECT_0)
						{
							nIDErrPrompt = IDS_ERR_INFINITE_LOOP;
							errormsg("Possible infinite loop");
						}
					}
#endif
					switch (stmttok->kind)
					{

					case tokrem:
						/* blank case */
						break;

					case toklist:
						cmdlist(&V);
						break;

					case tokrun:
						cmdrun(&V);
						break;

					case toknew:
						cmdnew(&V);
						break;

					case tokload:
						cmdload(false, stringexpr(STR1, &V), &V);
						break;

					case tokmerge:
						cmdload(true, stringexpr(STR1, &V), &V);
						break;

					case toksave:
						cmdsave(&V);
						break;

					case tokbye:
						cmdbye();
						break;

					case tokdel:
						cmddel(&V);
						break;

					case tokrenum:
						cmdrenum(&V);
						break;

					case toklet:
						cmdlet(false, &V);
						break;

					case tokvar:
						cmdlet(true, &V);
						break;

					case tokprint:
						cmdprint(&V);
						break;

					case tokpunch:
						cmdpunch(&V);
						break;

					case tokput:
						cmdput(&V);
						break;

					case tokput_:
						cmdput_(&V);
						break;

					case tokchange_por:
						cmdchange_por(&V);
						break;

					case tokchange_surf:
						cmdchange_surf(&V);
						break;

#if defined PHREEQ98 || defined MULTICHART
					case tokgraph_x:
						cmdgraph_x(&V);
						break;

					case tokgraph_y:
						cmdgraph_y(&V);
						break;

					case tokgraph_sy:
						cmdgraph_sy(&V);
						break;
#endif
#if defined MULTICHART
					case tokplot_xy:
						cmdplot_xy(&V);
						break;
#endif

					case tokinput:
						if (phreeqci_gui)
						{
							_ASSERTE(nIDErrPrompt == 0);
							nIDErrPrompt = IDS_ERR_INPUT_NOTLEGAL;
							errormsg ("Basic command INPUT is not a legal command in PHREEQC.");
						}
						else
						{
							error_msg ("Basic command INPUT is not a legal command in PHREEQC.", STOP);
						}
						break;

					case tokgoto:
						cmdgoto(&V);
						break;

					case tokif:
						cmdif(&V);
						break;

					case tokelse:
						cmdelse(&V);
						break;

					case tokend:
						cmdend(&V);
						break;

					case tokstop:
						P_escapecode = -20;
						throw PBasicStop();
						//goto _Ltry1;
						break;

					case tokfor:
						cmdfor(&V);
						break;

					case toknext:
						cmdnext(&V);
						break;

					case tokwhile:
						cmdwhile(&V);
						break;

					case tokwend:
						cmdwend(&V);
						break;

					case tokgosub:
						cmdgosub(&V);
						break;

					case tokreturn:
						cmdreturn(&V);
						break;

					case tokread:
						cmdread(&V);
						break;

					case tokdata:
						cmddata(&V);
						break;

					case tokrestore:
						cmdrestore(&V);
						break;

					case tokgotoxy:
						cmdgotoxy(&V);
						break;

					case tokon:
						cmdon(&V);
						break;

					case tokdim:
						cmddim(&V);
						break;

					case tokerase:
						cmderase(&V);
						break;

					case tokpoke:
						cmdpoke(&V);
						break;

					default:
						if (phreeqci_gui)
						{
							_ASSERTE(nIDErrPrompt == 0);
							nIDErrPrompt = IDS_ERR_ILLEGAL;
						}
						Utilities::strcat_safe(STR1, MAX_LENGTH, "Illegal command in line: ");
						if (strcmp(inbuf, "run"))
							Utilities::strcat_safe(STR1, MAX_LENGTH, inbuf);
						errormsg(STR1);
						break;
					}
				}
				if (!V.elseflag && !iseos(&V))
					checkextra(&V);
				stmttok = V.t;
			}
			while (V.t != NULL);
			if (stmtline != NULL)
			{
				if (!V.gotoflag)
					stmtline = stmtline->next;
				if (stmtline != NULL)
					stmttok = stmtline->txt;
			}
		}
		while (stmtline != NULL);
	}
	catch (const PBasicStop&)
	{
		//_Ltry1:
		if (P_escapecode == -20)
			PhreeqcPtr->warning_msg("Break");
		/* printf("Break"); */
		else if (P_escapecode != 42)
		{
			switch (P_escapecode)
			{

			case -4:
				{
					char * error_string = PhreeqcPtr->sformatf( "Integer overflow in BASIC line\n %ld %s", stmtline->num, stmtline->inbuf);
					PhreeqcPtr->warning_msg(error_string);
				}
				break;

			case -5:
				{
					char * error_string = PhreeqcPtr->sformatf( "Divide by zero in BASIC line\n %ld %s", stmtline->num, stmtline->inbuf);
					PhreeqcPtr->warning_msg(error_string);
				}
				break;

			case -6:
				{
					char * error_string = PhreeqcPtr->sformatf( "Real math overflow in BASIC line\n %ld %s", stmtline->num, stmtline->inbuf);
					PhreeqcPtr->warning_msg(error_string);
				}
				break;

			case -7:
				{
					char * error_string = PhreeqcPtr->sformatf( "Real math underflow in BASIC line\n %ld %s", stmtline->num, stmtline->inbuf);
					PhreeqcPtr->warning_msg(error_string);
				}
				break;

			case -8:
			case -19:
			case -18:
			case -17:
			case -16:
			case -15:
				{
				char * error_string = PhreeqcPtr->sformatf( "Value range error in BASIC line\n %ld %s", stmtline->num, stmtline->inbuf);
				PhreeqcPtr->warning_msg(error_string);
				}
				break;

			case -10:
				{
					char * error_string = PhreeqcPtr->sformatf("I/O Error %d", (int) P_ioresult);
					PhreeqcPtr->warning_msg(error_string);
				}
				break;

			default:
				if (EXCP_LINE != -1)
				{
					char * error_string = PhreeqcPtr->sformatf( "%12ld\n", EXCP_LINE);
					PhreeqcPtr->warning_msg(error_string);
				}
				_Escape(P_escapecode);
				break;
			}
		}
		if (stmtline != NULL)
		{
			if (phreeqci_gui)
			{
				_ASSERTE(nErrLineNumber == 0);
				nErrLineNumber = stmtline->num;
			}
			else
			{
				char * error_string = PhreeqcPtr->sformatf( " in BASIC line\n %ld %s", stmtline->num, stmtline->inbuf);
				error_msg(error_string, CONTINUE);
			}
		}
	} // end catch
}								/*exec */

int PBasic::
free_dim_stringvar(varrec *l_varbase)
{
	int i, k;
	if (l_varbase->numdims > 0)
	{
		k = 1;
		for (i = 0; i < l_varbase->numdims; i++)
		{
			k = k * (l_varbase->dims[i]);
		}
		for (i = 0; i < k; i++)
		{
			PhreeqcPtr->free_check_null(l_varbase->UU.U1.sarr[i]);
		}
		l_varbase->UU.U1.sarr = (char **) PhreeqcPtr->free_check_null(l_varbase->UU.U1.sarr);
	}
	return (OK);
}

#if defined MULTICHART
void PBasic::
cmdplot_xy(struct LOC_exec *LINK)
{
	bool semiflag;
	valrec n[2];

	char STR[2][256];
	int i = 0;
	semiflag = false;

	while (!iseos(LINK) && i < 2)
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			i++;
			continue;
		}
		n[i] = expr(LINK);
		if (n[i].stringval)
		{
			Utilities::strcpy_safe(STR[i], MAX_LENGTH, n[i].UU.sval);
			PhreeqcPtr->PHRQ_free(n[i].UU.sval);
		}
		else
			numtostr(STR[i], n[i].UU.val);
	}

	ChartObject *chart = PhreeqcPtr->chart_handler.Get_current_chart();
	if (chart == NULL) return;

	std::string x_str(STR[0]), y_str(STR[1]);

// Code formerly in PlotXY, included here
	{

		bool new_sim = false, new_trans = false;

		if (chart->Get_FirstCallToUSER_GRAPH() && chart->Get_colnr() == 0)
			chart->Set_prev_sim_no(PhreeqcPtr->simulation);
		else
			chart->Set_prev_sim_no(PhreeqcPtr->simulation);

		// Add a curve if necessary
		if ((int) chart->Get_Curves().size() == chart->Get_colnr())
		{
			std::string head("");
			if (chart->Get_new_headings().size() > 0)
			{
				head = chart->Get_new_headings()[0];
				chart->Get_new_headings().erase(chart->Get_new_headings().begin());
			}
			if (chart->Get_new_plotxy_curves().size() > 0)
			{
				// find plotxy curve definitions
				chart->Add_curve(true, head,
					chart->Get_new_plotxy_curves()[0].Get_line_w(),
					chart->Get_new_plotxy_curves()[0].Get_symbol(),
					chart->Get_new_plotxy_curves()[0].Get_symbol_size(),
					chart->Get_new_plotxy_curves()[0].Get_y_axis(),
					chart->Get_new_plotxy_curves()[0].Get_color()					
					);
				// pop plotxy curve definition
				chart->Get_new_plotxy_curves().erase(chart->Get_new_plotxy_curves().begin());
			}
			else
			{
				chart->Add_curve(true, head);
			}
			chart->Set_curve_added(true);
		}

		if (x_str.size() > 0 && y_str.size() > 0)
		{
			chart->Get_Curves()[chart->Get_colnr()]->Get_x().push_back(atof(x_str.c_str()));
			chart->Get_Curves()[chart->Get_colnr()]->Get_y().push_back(atof(y_str.c_str()));
			chart->Set_point_added(true);

			// Mark added curve for first point, might have been invisible in DefineCurves
			if (chart->Get_Curves()[chart->Get_colnr()]->Get_x().size() == 1)
				chart->Set_curve_added(true);
		}
	}
	chart->Set_colnr(chart->Get_colnr() + 1);
}

void PBasic::
cmdgraph_x(struct LOC_exec *LINK)
{
	bool semiflag;
	valrec n;

	semiflag = false;

	ChartObject *chart = PhreeqcPtr->chart_handler.Get_current_chart();
	if (chart == NULL) return;

	while (!iseos(LINK))
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);

		if (n.stringval)
		{
			if (strlen(n.UU.sval) > 0)
			{
				chart->Set_graph_x(atof(n.UU.sval));
			}
			PhreeqcPtr->PHRQ_free(n.UU.sval);
		}
		else
			chart->Set_graph_x(n.UU.val);
	}
	if ((int) chart->Get_Curves().size() == chart->Get_colnr())
	{
		if (chart->Get_new_headings().size() > 0)
		{
			// remove x heading
			//if (chart->Get_colnr() == chart->Get_ColumnOffset())
			{
				chart->Get_new_headings().erase(chart->Get_new_headings().begin());
			}
		}
	}
}

void PBasic::
cmdgraph_y(struct LOC_exec *LINK)
{
	bool semiflag;
	valrec n;

	semiflag = false;

	ChartObject *chart = PhreeqcPtr->chart_handler.Get_current_chart();
	if (chart == NULL) return;

	while (!iseos(LINK))
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);

		// wait until x value is known before storing in curve
		if (n.stringval)
		{
			if (strlen(n.UU.sval) > 0)
			{
				chart->Get_graph_y()[chart->Get_colnr()] = atof(n.UU.sval);
			}	
			PhreeqcPtr->PHRQ_free(n.UU.sval);
		}
		else
		{
			chart->Get_graph_y()[chart->Get_colnr()] = n.UU.val;
		}
		chart->Set_point_added(true);

		// Add a new curve if necessary
		if ((int) chart->Get_Curves().size() == chart->Get_colnr())
		{
			if (chart->Get_new_headings().size() > 0)
			{
				chart->Add_curve(false, chart->Get_new_headings()[0]);
				chart->Get_new_headings().erase(chart->Get_new_headings().begin());
			}
			else
			{
				chart->Add_curve(false);
			}
			chart->Set_curve_added(true);
		}
		chart->Set_colnr(chart->Get_colnr() + 1);
	}
}

void PBasic::
cmdgraph_sy(struct LOC_exec *LINK)
{
	bool semiflag;
	valrec n;

	semiflag = false;

	ChartObject *chart = PhreeqcPtr->chart_handler.Get_current_chart();
	if (chart == NULL) return;

	while (!iseos(LINK))
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);

		chart->Get_secondary_y()[chart->Get_colnr()] = true;

		// wait until x value is known before storing in curve
		if (n.stringval)
		{
			if (strlen(n.UU.sval) > 0)
			{
				chart->Get_graph_y()[chart->Get_colnr()] = atof(n.UU.sval);
			}
			PhreeqcPtr->PHRQ_free(n.UU.sval);
		}
		else
			chart->Get_graph_y()[chart->Get_colnr()] = n.UU.val;
		chart->Set_point_added(true);
		// Add a new curve if necessary
		if ((int) chart->Get_Curves().size() == chart->Get_colnr())
		{
			if (chart->Get_new_headings().size() > 0)
			{
				chart->Add_curve(false, chart->Get_new_headings()[0]);
				chart->Get_new_headings().erase(chart->Get_new_headings().begin());
			}
			else
			{
				chart->Add_curve(false);
			}
			chart->Set_curve_added(true);
		}
		chart->Set_colnr(chart->Get_colnr() + 1);
	}
}
#endif // MULTICHART

/* In case your system lacks these... */

long PBasic::
my_labs(long l_x)
{
	return ((l_x > 0) ? l_x : -l_x);
}


/* #define __STDC__  */ /* PHREEQ98 */

void *  PBasic::
my_memmove(void * d, Const void * l_s, size_t n)
{
	char *dd = (char *) d, *ss = (char *) l_s;
	if (dd < ss || (unsigned int) (dd - ss) >= n)
	{
		memcpy(dd, ss, n);
	}
	else if (n > 0)
	{
		dd += n;
		ss += n;
		while (n-- > 0)
			*--dd = *--ss;
	}
	return d;
}

void *  PBasic::
my_memcpy(void * d, Const void * l_s, size_t n)
{
	char *ss = (char *) l_s, *dd = (char *) d;
	while (n-- > 0)
		*dd++ = *ss++;
	return d;
}

int PBasic::
my_memcmp(Const void * s1, Const void * s2, size_t n)
{
	char *a = (char *) s1, *b = (char *) s2;
	int i;
	while (n-- > 0)
		if ((i = (*a++) - (*b++)) != 0)
			return i;
	return 0;
}

void * PBasic::
my_memset(void * d, int c, size_t n)
{
	char *dd = (char *) d;
	while (n-- > 0)
		*dd++ = (char) c;
	return d;
}

int PBasic::
my_toupper(int c)
{
	if (islower(c))
		return _toupper(c);
	else
		return c;
}

int PBasic::
my_tolower(int c)
{
	if (isupper(c))
		return _tolower(c);
	else
		return c;
}

long PBasic::
ipow(long a, long b)
{
	long v;

	if (a == 0 || a == 1)
		return a;
	if (a == -1)
		return (b & 1) ? -1 : 1;
	if (b < 0)
		return 0;
	if (a == 2)
		return 1L << b;
	v = (b & 1) ? a : 1;
	while ((b >>= 1) > 0)
	{
		a *= a;
		if (b & 1)
			v *= a;
	}
	return v;
}

/* Common string functions: */

/* Store in "ret" the substring of length "len" starting from "pos" (1-based).
   Store a shorter or null string if out-of-range.  Return "ret". */
char * PBasic::
strsub(char *ret, char *l_s, int pos,
	   int len)
{
	char *s2;

	if (--pos < 0 || len <= 0)
	{
		*ret = 0;
		return ret;
	}
	while (pos > 0)
	{
		if (!*l_s++)
		{
			*ret = 0;
			return ret;
		}
		pos--;
	}
	s2 = ret;
	while (--len >= 0)
	{
		if (!(*s2++ = *l_s++))
			return ret;
	}
	*s2 = 0;
	return ret;
}

/* Return the index of the first occurrence of "pat" as a substring of "s",
   starting at index "pos" (1-based).  Result is 1-based, 0 if not found. */

int PBasic::
strpos2(char *l_s, char *pat, int pos)
{
	char *cp, ch;
	int slen;

	if (--pos < 0)
		return 0;
	slen = (int) strlen(l_s) - pos;
	cp = l_s + pos;
	if (!(ch = *pat++))
		return 0;
	pos = (int) strlen(pat);
	slen -= pos;
	while (--slen >= 0)
	{
		if (*cp++ == ch && !strncmp(cp, pat, pos))
			return (int) (cp - l_s);
	}
	return 0;
}

/* Case-insensitive version of strcmp. */
int PBasic::
strcicmp(char *s1, char *s2)
{
	unsigned char c1, c2;

	while (*s1)
	{
		if (*s1++ != *s2++)
		{
			if (!s2[-1])
				return 1;
			c1 = (unsigned char) toupper(s1[-1]);
			c2 = (unsigned char) toupper(s2[-1]);
			if (c1 != c2)
				return c1 - c2;
		}
	}
	if (*s2)
		return -1;
	return 0;
}

/* HP and Turbo Pascal string functions: */

/* Trim blanks at left end of string. */
char *  PBasic::
strltrim(char *l_s)
{
	while (Isspace((int) *l_s++));
	return l_s - 1;
}

/* Trim blanks at right end of string. */
char * PBasic::
strrtrim(char *l_s)
{
	char *s2 = l_s;

	if (!*l_s)
		return l_s;
	while (*++s2);
	while (s2 > l_s && Isspace((int) *--s2))
		*s2 = 0;
	return l_s;
}

/* Store in "ret" string "s" with enough pad chars added to reach "size". */

/* Copy the substring of length "len" from index "spos" of "s" (1-based)
   to index "dpos" of "d", lengthening "d" if necessary.  Length and
   indices must be in-range. */
void PBasic::
strmove(int len, char *l_s, int spos,
		char *d, int dpos)
{
	l_s += (size_t)spos - 1;
	d += (size_t)dpos - 1;
	while (*d && --len >= 0)
		*d++ = *l_s++;
	if (len > 0)
	{
		while (--len >= 0)
			*d++ = *l_s++;
		*d = 0;
	}
}

/* Insert string "src" at index "pos" of "dst". */
void PBasic::
strinsert(char *src, char *dst, int pos)
{
	int slen, dlen;

	if (--pos < 0)
		return;
	dlen = (int) strlen(dst);
	dst += dlen;
	dlen -= pos;
	if (dlen <= 0)
	{
		strcpy(dst, src);
		return;
	}
	slen = (int) strlen(src);
	do
	{
		dst[slen] = *dst;
		--dst;
	}
	while (--dlen >= 0);
	dst++;
	while (--slen >= 0)
		*dst++ = *src++;
}

/* File functions */

/* Peek at next character of input stream; return EOF at end-of-file. */
int PBasic::
P_peek(FILE * f)
{
	int ch;

	ch = getc(f);
	if (ch == EOF)
		return EOF;
	ungetc(ch, f);
	return (ch == '\n') ? ' ' : ch;
}

/* Check if at end of file, using Pascal "eof" semantics.  End-of-file for
   stdin is broken; remove the special case for it to be broken in a
   different way. */
/*int P_eof(FILE *f)*/
int PBasic::
P_eof(void)
{
	return 0;
}

/* Check if at end of line (or end of entire file). */
int PBasic::
P_eoln(FILE * f)
{
	int ch;

	ch = getc(f);
	if (ch == EOF)
		return 1;
	ungetc(ch, f);
	return (ch == '\n');
}

/* Read a packed array of characters from a file. */
void PBasic::
P_readpaoc(FILE * f, char *l_s, int len)
{
	int ch;

	for (;;)
	{
		if (len <= 0)
			return;
		ch = getc(f);
		if (ch == EOF || ch == '\n')
			break;
		*l_s++ = (char) ch;
		--len;
	}
	while (--len >= 0)
		*l_s++ = ' ';
	if (ch != EOF)
		ungetc(ch, f);
}

void PBasic::
P_readlnpaoc(FILE * f, char *l_s, int len)
{
	int ch;

	for (;;)
	{
		ch = getc(f);
		if (ch == EOF || ch == '\n')
			break;
		if (len > 0)
		{
			*l_s++ = (char) ch;
			--len;
		}
	}
	while (--len >= 0)
		*l_s++ = ' ';
}

/* Compute maximum legal "seek" index in file (0-based). */
long PBasic::
P_maxpos(FILE * f)
{
	long savepos = ftell(f);
	long val;

	if (fseek(f, 0L, SEEK_END))
		return -1;
	val = ftell(f);
	if (fseek(f, savepos, SEEK_SET))
		return -1;
	return val;
}

/* Use packed array of char for a file name. */
char * PBasic::
P_trimname(char * fn, int len)
{
	char *cp = fnbuf;

	while (--len >= 0 && *fn && !isspace((int) *fn))
		*cp++ = *fn++;
	*cp = 0;
	return fnbuf;
}

/* Pascal's "memavail" doesn`t make much sense in Unix with virtual memory.
   We fix memory size as 10Meg as a reasonable compromise. */

long PBasic::
memavail(void)
{
	return 10000000;			/* worry about this later! */
}

long  PBasic::
maxavail(void)
{
	return memavail();
}

/* Sets are stored as an array of longs.  S[0] is the size of the set;
   S[N] is the n`th 32-bit chunk of the set.  S[0] equals the maximum
   I such that S[I] is nonzero.  S[0] is zero for an empty set.  Within
   each long, bits are packed from lsb to msb.  The first bit of the
   set is the element with ordinal value 0.  (Thus, for a "set of 5..99",
   the lowest five bits of the first long are unused and always zero.) */

/* (Sets with 32 or fewer elements are normally stored as plain longs.) */
long * PBasic::
P_setunion(long *d, long *s1, long *s2)	/* d := s1 + s2 */
{
	long *dbase = d++;
	int sz1 = *s1++, sz2 = *s2++;
	while (sz1 > 0 && sz2 > 0)
	{
		*d++ = *s1++ | *s2++;
		sz1--, sz2--;
	}
	while (--sz1 >= 0)
		*d++ = *s1++;
	while (--sz2 >= 0)
		*d++ = *s2++;
	*dbase = (int) (d - dbase - 1);
	return dbase;
}

long * PBasic::
P_setint(long *d, long *s1, long *s2)	/* d := s1 * s2 */
{
	long *dbase = d++;
	int sz1 = *s1++, sz2 = *s2++;
	while (--sz1 >= 0 && --sz2 >= 0)
		*d++ = *s1++ & *s2++;
	while (--d > dbase && !*d);
	*dbase = (int) (d - dbase);
	return dbase;
}

long * PBasic::
P_setdiff(long *d, long *s1, long *s2)	/* d := s1 - s2 */
{
	long *dbase = d++;
	int sz1 = *s1++, sz2 = *s2++;
	while (--sz1 >= 0 && --sz2 >= 0)
		*d++ = *s1++ & ~*s2++;
	if (sz1 >= 0)
	{
		while (sz1-- >= 0)
			*d++ = *s1++;
	}
	while (--d > dbase && !*d);
	*dbase = (int) (d - dbase);
	return dbase;
}

long * PBasic::
P_setxor(long *d, long *s1, long *s2)	/* d := s1 / s2 */
{
	long *dbase = d++;
	int sz1 = *s1++, sz2 = *s2++;
	while (sz1 > 0 && sz2 > 0)
	{
		*d++ = *s1++ ^ *s2++;
		sz1--, sz2--;
	}
	while (--sz1 >= 0)
		*d++ = *s1++;
	while (--sz2 >= 0)
		*d++ = *s2++;
	while (--d > dbase && !*d);
	*dbase = (int) (d - dbase);
	return dbase;
}

long * PBasic::
P_addset(long *l_s, unsigned val)	/* s := s + [val] */
{
	long *sbase = l_s;
	int bit, size;
	bit = val % SETBITS;
	val /= SETBITS;
	size = *l_s;
	if ((long) ++val > size)
	{
		l_s += size;
		while ((long) val > size)
			*++l_s = 0, size++;
		*sbase = size;
	}
	else
		l_s += val;
	*l_s |= 1L << bit;
	return sbase;
}

//long * PBasic::
//P_addsetr(long *l_s, unsigned v1, unsigned v2)	/* s := s + [v1..v2] */
//{
//	long *sbase = l_s;
//	int b1, b2, size;
//	if ((int) v1 > (int) v2)
//		return sbase;
//	b1 = v1 % SETBITS;
//	v1 /= SETBITS;
//	b2 = v2 % SETBITS;
//	v2 /= SETBITS;
//	size = *l_s;
//	v1++;
//	if ((int) ++v2 > size)
//	{
//		while ((int) v2 > size)
//			l_s[++size] = 0;
//		l_s[v2] = 0;
//		*l_s = v2;
//	}
//	l_s += v1;
//	if (v1 == v2)
//	{
//		*l_s |= (~((-2L) << (b2 - b1))) << b1;
//	}
//	else
//	{
//		*l_s++ |= (-1L) << b1;
//		while (++v1 < v2)
//			*l_s++ = -1;
//		*l_s |= ~((-2L) << b2);
//	}
//	return sbase;
//}

long *  PBasic::
P_remset(long *l_s, unsigned val)	/* s := s - [val] */
{
	int bit;
	bit = val % SETBITS;
	val /= SETBITS;
	if ((long) ++val <= *l_s)
	{
		if (!(l_s[val] &= ~(1L << bit)))
			while (*l_s && !l_s[*l_s])
				(*l_s)--;
	}
	return l_s;
}

int PBasic::
P_setequal(long *s1, long *s2)	/* s1 = s2 */
{
	int size = *s1++;
	if (*s2++ != size)
		return 0;
	while (--size >= 0)
	{
		if (*s1++ != *s2++)
			return 0;
	}
	return 1;
}

int PBasic::
P_subset(long *s1, long *s2)	/* s1 <= s2 */
{
	int sz1 = *s1++, sz2 = *s2++;
	if (sz1 > sz2)
		return 0;
	while (--sz1 >= 0)
	{
		if (*s1++ & ~*s2++)
			return 0;
	}
	return 1;
}

long * PBasic::
P_setcpy(long *d, long *l_s)	/* d := s */
{
	long *save_d = d;

#ifdef SETCPY_MEMCPY
	memcpy(d, l_s, (*l_s + 1) * sizeof(long));
#else
	int i = *l_s + 1;
	while (--i >= 0)
		*d++ = *l_s++;
#endif
	return save_d;
}

/* s is a "smallset", i.e., a 32-bit or less set stored
   directly in a long. */
long * PBasic::
P_expset(long *d, long l_s)	/* d := s */
{
	if (l_s)
	{
		d[1] = l_s;
		*d = 1;
	}
	else
		*d = 0;
	return d;
}

long PBasic::
P_packset(long *l_s)		/* convert s to a small-set */
{
	if (*l_s++)
		return *l_s;
	else
		return 0;
}

int PBasic::
_OutMem(void)
{
	return _Escape(-2);
}

int  PBasic::
_CaseCheck(void)
{
	return _Escape(-9);
}

int  PBasic::
_NilCheck(void)
{
	return _Escape(-3);
}

#ifdef NPP
/* The following is suitable for the HP Pascal operating system.
   It might want to be revised when emulating another system. */

char * PBasic::
_ShowEscape(char *buf, int code, int ior, char *prefix)
{
	char *bufp;

	if (prefix && *prefix)
	{
		strcpy(buf, prefix);
		strcat(buf, ": ");
		bufp = buf + strlen(buf);
	}
	else
	{
		bufp = buf;
	}
	if (code == -10)
	{
		snprintf(bufp, sizeof(bufp), "Pascal system I/O error %d", ior); // FIXME -- replace sizeof
		switch (ior)
		{
		case 3:
			strcat(buf, " (illegal I/O request)");
			break;
		case 7:
			strcat(buf, " (bad file name)");
			break;
		case FileNotFound:		/*10 */
			strcat(buf, " (file not found)");
			break;
		case FileNotOpen:		/*13 */
			strcat(buf, " (file not open)");
			break;
		case BadInputFormat:	/*14 */
			strcat(buf, " (bad input format)");
			break;
		case 24:
			strcat(buf, " (not open for reading)");
			break;
		case 25:
			strcat(buf, " (not open for writing)");
			break;
		case 26:
			strcat(buf, " (not open for direct access)");
			break;
		case 28:
			strcat(buf, " (string subscript out of range)");
			break;
		case EndOfFile:		/*30 */
			strcat(buf, " (end-of-file)");
			break;
		case FileWriteError:	/*38 */
			strcat(buf, " (file write error)");
			break;
		}
	}
	else
	{
		snprintf(bufp, sizeof(bufp), "Pascal system error %d", code); // FIXME -- replace sizeof
		switch (code)
		{
		case -2:
			strcat(buf, " (out of memory)");
			break;
		case -3:
			strcat(buf, " (reference to NIL pointer)");
			break;
		case -4:
			strcat(buf, " (integer overflow)");
			break;
		case -5:
			strcat(buf, " (divide by zero)");
			break;
		case -6:
			strcat(buf, " (real math overflow)");
			break;
		case -8:
			strcat(buf, " (value range error)");
			break;
		case -9:
			strcat(buf, " (CASE value range error)");
			break;
		case -12:
			strcat(buf, " (bus error)");
			break;
		case -20:
			strcat(buf, " (stopped by user)");
			break;
		}
	}
	return buf;
}
#endif
int PBasic::
_Escape(int code)
{
	P_escapecode = code;
	throw PBasicStop();

	// following not used
}

int PBasic::
_EscIO(int code)
{
	P_ioresult = code;
	return _Escape(-10);
}

const std::map<const std::string, PBasic::BASIC_TOKEN>::value_type temp_tokens[] = {
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("+",                  PBasic::tokplus),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("-",                  PBasic::tokminus),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("*",                  PBasic::toktimes),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("/",                  PBasic::tokdiv),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("^",                  PBasic::tokup),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("( or [",             PBasic::toklp),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type(") or ]",             PBasic::tokrp),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("]",                  PBasic::tokcomma),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type(";",                  PBasic::toksemi),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type(":",                  PBasic::tokcolon),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("=",                  PBasic::tokeq),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("<",                  PBasic::toklt),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("<=",                 PBasic::tokle),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type(">",                  PBasic::tokgt),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type(">=",                 PBasic::tokge),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("and",                PBasic::tokand),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("or",                 PBasic::tokor),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("xor",                PBasic::tokxor),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("not",                PBasic::toknot),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("mod",                PBasic::tokmod),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("sqr",                PBasic::toksqr),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("sqrt",               PBasic::toksqrt),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("ceil",               PBasic::tokceil),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("floor",              PBasic::tokfloor),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("sin",                PBasic::toksin),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("cos",                PBasic::tokcos),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("tan",                PBasic::toktan),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("arctan",             PBasic::tokarctan),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("log",                PBasic::toklog),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("exp",                PBasic::tokexp),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("abs",                PBasic::tokabs),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("sgn",                PBasic::toksgn),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("str$",               PBasic::tokstr_),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("val",                PBasic::tokval),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("chr$",               PBasic::tokchr_),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("eol$",               PBasic::tokeol_),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("eol_notab$",         PBasic::tokeol_notab_),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("no_newline$",        PBasic::tokno_newline_),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("asc",                PBasic::tokasc),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("len",                PBasic::toklen),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("mid$",               PBasic::tokmid_),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("peek",               PBasic::tokpeek),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("let",                PBasic::toklet),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("input",              PBasic::tokinput),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("goto",               PBasic::tokgoto),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("go to",              PBasic::tokgoto),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("if",                 PBasic::tokif),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("end",                PBasic::tokend),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("stop",               PBasic::tokstop),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("for",                PBasic::tokfor),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("next",               PBasic::toknext),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("while",              PBasic::tokwhile),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("wend",               PBasic::tokwend),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("gosub",              PBasic::tokgosub),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("return",             PBasic::tokreturn),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("read",               PBasic::tokread),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("data",               PBasic::tokdata),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("restore",            PBasic::tokrestore),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("gotoxy",             PBasic::tokgotoxy),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("on",                 PBasic::tokon),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("dim",                PBasic::tokdim),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("erase",              PBasic::tokerase),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("poke",               PBasic::tokpoke),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("list",               PBasic::toklist),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("run",                PBasic::tokrun),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("new",                PBasic::toknew),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("load",               PBasic::tokload),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("merge",              PBasic::tokmerge),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("save",               PBasic::toksave),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("bye",                PBasic::tokbye),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("quit",               PBasic::tokbye),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("del",                PBasic::tokdel),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("renum",              PBasic::tokrenum),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("then",               PBasic::tokthen),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("else",               PBasic::tokelse),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("to",                 PBasic::tokto),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("step",               PBasic::tokstep),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("act",                PBasic::tokact),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("add_heading",        PBasic::tokadd_heading),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("alk",                PBasic::tokalk),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("aphi",               PBasic::tokaphi),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("calc_value",         PBasic::tokcalc_value),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("callback",           PBasic::tokcallback),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("cell_no",            PBasic::tokcell_no),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("change_por",         PBasic::tokchange_por),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("change_surf",        PBasic::tokchange_surf),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("charge_balance",     PBasic::tokcharge_balance),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("current_a",          PBasic::tokcurrent_a),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("debye_length",       PBasic::tokdebye_length),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("delta_h_phase",      PBasic::tokdelta_h_phase),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("delta_h_species",    PBasic::tokdelta_h_species),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("description",        PBasic::tokdescription),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("dh_a0",              PBasic::tokdh_a0),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("dh_a",               PBasic::tokdh_a),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("dh_av",              PBasic::tokdh_av),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("dh_b",               PBasic::tokdh_b),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("dh_bdot",            PBasic::tokdh_bdot),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("diff_c",             PBasic::tokdiff_c),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("dist",               PBasic::tokdist),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("edl",                PBasic::tokedl),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("edl_species",        PBasic::tokedl_species),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("eps_r",              PBasic::tokeps_r),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("eq_frac",            PBasic::tokeq_frac),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("equi",               PBasic::tokequi),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("equi_delta",         PBasic::tokequi_delta),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("equiv_frac",         PBasic::tokeq_frac),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("exists",             PBasic::tokexists),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("gamma",              PBasic::tokgamma),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("gas",                PBasic::tokgas),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("gas_p",              PBasic::tokgas_p),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("gas_vm",             PBasic::tokgas_vm),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("get",                PBasic::tokget),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("get$",               PBasic::tokget_),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("get_por",            PBasic::tokget_por),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("gfw",                PBasic::tokgfw),
#if defined (PHREEQ98) || defined (MULTICHART)
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("graph_x", PBasic::tokgraph_x),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("graph_y", PBasic::tokgraph_y),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("graph_sy", PBasic::tokgraph_sy),
#endif
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("instr",              PBasic::tokinstr),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("iso",                PBasic::tokiso),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("iso_unit",           PBasic::tokiso_unit),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("iterations",         PBasic::tokiterations),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("kappa",              PBasic::tokkappa),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("kin",                PBasic::tokkin),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("kin_delta",          PBasic::tokkin_delta),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("kin_time",           PBasic::tokkin_time),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("kinetics_formula",   PBasic::tokkinetics_formula),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("kinetics_formula$",  PBasic::tokkinetics_formula),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("la",                 PBasic::tokla),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("lg",                 PBasic::toklg),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("list_s_s",           PBasic::toklist_s_s),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("lk_named",           PBasic::toklk_named),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("lk_phase",           PBasic::toklk_phase),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("lk_species",         PBasic::toklk_species),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("lm",                 PBasic::toklm),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("log10",              PBasic::toklog10),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("ltrim",              PBasic::tokltrim),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("m0",                 PBasic::tokm0),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("m",                  PBasic::tokm),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("mcd_jtot",           PBasic::tokmcd_jtot),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("mcd_jconc",          PBasic::tokmcd_jconc),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("misc1",              PBasic::tokmisc1),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("misc2",              PBasic::tokmisc2),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("mol",                PBasic::tokmol),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("mu",                 PBasic::tokmu),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("osmotic",            PBasic::tokosmotic),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("pad",                PBasic::tokpad),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("pad$",               PBasic::tokpad_),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("parm",               PBasic::tokparm),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("rate_pk",            PBasic::tokrate_pk),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("rate_svd",           PBasic::tokrate_svd),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("rate_hermanska",     PBasic::tokrate_hermanska),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("meang",              PBasic::tokmeang),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("percent_error",      PBasic::tokpercent_error),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("phase_formula",      PBasic::tokphase_formula),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("phase_formula$",     PBasic::tokphase_formula_),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("phase_vm",           PBasic::tokphase_vm),
#if defined MULTICHART
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("plot_xy", PBasic::tokplot_xy),
#endif
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("porevolume",         PBasic::tokporevolume),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("pot_v",              PBasic::tokpot_v),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("pr_p",               PBasic::tokpr_p),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("pr_phi",             PBasic::tokpr_phi),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("pressure",           PBasic::tokpressure),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("print", PBasic::tokprint),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("punch", PBasic::tokpunch),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("put",                PBasic::tokput),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("put$",               PBasic::tokput_),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("qbrn",               PBasic::tokqbrn),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("rem",                PBasic::tokrem),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("rho",                PBasic::tokrho),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("rho_0",              PBasic::tokrho_0),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("rtrim",              PBasic::tokrtrim),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("rxn",                PBasic::tokrxn),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("s_s",                PBasic::toks_s),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("sc",                 PBasic::toksc),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("setdiff_c",          PBasic::toksetdiff_c),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("si",                 PBasic::toksi),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("sim_no",             PBasic::toksim_no),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("sim_time",           PBasic::toksim_time),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("soln_vol",           PBasic::toksoln_vol),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("species_formula",    PBasic::tokspecies_formula),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("species_formula$",   PBasic::tokspecies_formula_),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("phase_equation",     PBasic::tokphase_equation),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("phase_equation$",    PBasic::tokphase_equation_),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("species_equation",   PBasic::tokspecies_equation),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("species_equation$",  PBasic::tokspecies_equation_),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("sr",                 PBasic::toksr),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("step_no",            PBasic::tokstep_no),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("str_e$",             PBasic::tokstr_e_),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("str_f$",             PBasic::tokstr_f_),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("sum_gas",            PBasic::toksum_gas),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("sum_s_s",            PBasic::toksum_s_s),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("sum_species",        PBasic::toksum_species),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("surf",               PBasic::toksurf),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("sys",                PBasic::toksys),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("t_sc",               PBasic::tokt_sc),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("f_visc",             PBasic::tokf_visc),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("tc",                 PBasic::toktc),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("time",               PBasic::toktime),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("title",              PBasic::toktitle),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("tk",                 PBasic::toktk),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("tot",                PBasic::toktot),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("total_time",         PBasic::toktotal_time),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("totmol",             PBasic::toktotmol),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("totmole",            PBasic::toktotmole),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("totmoles",           PBasic::toktotmoles),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("trim",               PBasic::toktrim),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("viscos",             PBasic::tokviscos),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("viscos_0",           PBasic::tokviscos_0),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("vm",                 PBasic::tokvm),
	/* PHAST */
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("cell_pore_volume",   PBasic::tokcell_pore_volume),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("cell_porosity",      PBasic::tokcell_porosity),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("cell_saturation",    PBasic::tokcell_saturation),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("cell_volume",        PBasic::tokcell_volume),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("transport_cell_no",  PBasic::toktransport_cell_no),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("velocity_x",         PBasic::tokvelocity_x),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("velocity_y",         PBasic::tokvelocity_y),
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("velocity_z",         PBasic::tokvelocity_z),
	/* Undocumented */
	std::map<const std::string, PBasic::BASIC_TOKEN>::value_type("sa_declercq",        PBasic::toksa_declercq)
};
std::map<const std::string, PBasic::BASIC_TOKEN> PBasic::command_tokens(temp_tokens, temp_tokens + sizeof temp_tokens / sizeof temp_tokens[0]);

/* End. */

