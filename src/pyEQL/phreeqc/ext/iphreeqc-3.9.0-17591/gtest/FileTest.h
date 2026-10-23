#if !defined(FILETEST_H_INCLUDED)
#define FILETEST_H_INCLUDED

#include <string>

#if defined(_WIN32) || defined(__CYGWIN32__)
// DeleteFile defined in <windows.h>
#else
int DeleteFile(const char* szPathName);
#endif

bool FileExists(const char *szPathName);
size_t FileSize(const char *szPathName);

class FileTest
{
public:
	FileTest(std::string fn);
	~FileTest(void);

	std::string GetName(void);
	bool RemoveExisting(void);
	bool VerifyMissing(void);
	bool VerifyExists(void);
	bool Exists(void);
	int Delete(void);
	size_t Size(void);
	size_t LineCount(void);

protected:
	std::string _fn;
};

#endif // FILETEST_H_INCLUDED