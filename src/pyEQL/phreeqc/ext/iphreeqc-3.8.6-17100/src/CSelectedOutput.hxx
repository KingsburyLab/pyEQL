// SelectedOutput.h: interface for the CSelectedOutput class.
//
//////////////////////////////////////////////////////////////////////

#if !defined _INC_SELECTEDOUTPUT_H
#define _INC_SELECTEDOUTPUT_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string>
#include <map>
#include <list>
#include <vector>
#include "CVar.hxx"

#include "PHRQ_exports.h"

class IPQ_DLL_EXPORT CSelectedOutput
{
public:
	CSelectedOutput(void);
	virtual ~CSelectedOutput(void);

	int EndRow(void);
	void Clear(void);

	size_t GetRowCount(void)const;
	size_t GetColCount(void)const;

	CVar Get(int nRow, int nCol)const;
	VRESULT Get(int nRow, int nCol, VAR* pVAR)const;

	int PushBack(const char* key, const CVar& var);

	int PushBackDouble(const char* key, double dVal);
	int PushBackLong(const char* key, long lVal);
	int PushBackString(const char* key, const char* sVal);
	int PushBackEmpty(const char* key);

	// Serialize
	void Serialize(
		int row,
		std::vector<int> &types,
		std::vector<long> &longs,
		std::vector<double> &doubles,
		std::string &strings);
	void DeSerialize(
		std::vector<int> &types,        
		std::vector<long> &longs,       
		std::vector<double> &doubles,   
		std::string &strings);
	void Doublize(
		int &nrow,
		int &ncol,
		std::vector < double > &doubles);

#if defined(_DEBUG)
	void Dump(const char* heading);
	void AssertValid(void)const;
#endif

protected:
	friend std::ostream& operator<< (std::ostream &os, const CSelectedOutput &a);

	size_t m_nRowCount;

	std::vector< std::vector<CVar> > m_arrayVar;
	std::vector<CVar> m_vecVarHeadings;
	std::map< std::string, size_t > m_mapHeadingToCol;

private:
	static CSelectedOutput* s_instance;
};

#endif // !defined(_INC_SELECTEDOUTPUT_H)
