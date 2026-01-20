#if !defined(__ERROR_REPORTER_HXX_INC)
#define __ERROR_REPORTER_HXX_INC


#include <iosfwd>          // std::ostream
#include <cstdio>          // std::fprintf

class IErrorReporter
{
public:
	virtual size_t AddError(const char* error_msg) = 0;
	virtual void Clear(void) = 0;
	virtual ~IErrorReporter(void) = 0;
};

// Note: this is req'd in order for subclass dtors to be called
//
IErrorReporter::~IErrorReporter(void)
{
}

template <typename OS>
class CErrorReporter : public IErrorReporter
{
public:
	CErrorReporter(void);
	virtual ~CErrorReporter(void);

	virtual size_t AddError(const char *error_msg);
	virtual void Clear(void);
	OS* GetOS(void) { return m_pOS; }
protected:
	OS* m_pOS;
	size_t m_error_count;
};

template<typename OS>
CErrorReporter<OS>::CErrorReporter(void)
: m_pOS(0)
, m_error_count(0)
{
	this->m_pOS = new OS;
}

template<typename OS>
CErrorReporter<OS>::~CErrorReporter(void)
{
	delete this->m_pOS;
}

template<typename OS>
size_t CErrorReporter<OS>::AddError(const char* error_msg)
{
	++this->m_error_count;
	(*this->m_pOS) << error_msg;
	return this->m_error_count;
}

template<typename OS>
void CErrorReporter<OS>::Clear(void)
{
	this->m_error_count = 0;
	if (this->m_pOS->tellp() != std::ios::pos_type(-1))
	{
		delete this->m_pOS;
		this->m_pOS = new OS;
	}
}

#endif // __ERROR_REPORTER_HXX_INC
