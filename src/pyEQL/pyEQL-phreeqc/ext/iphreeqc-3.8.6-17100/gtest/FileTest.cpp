#include "FileTest.h"

#if defined(_WIN32) || defined(__CYGWIN32__)
#include <windows.h>
#else
#include <stdio.h>
#endif

#include "Phreeqc.h" /* snprintf */

#if defined(_WIN32) || defined(__CYGWIN32__)
bool FileExists(const char *szPathName)
{
	SECURITY_ATTRIBUTES sa;
	sa.nLength = sizeof(SECURITY_ATTRIBUTES);
	sa.lpSecurityDescriptor = NULL;
	sa.bInheritHandle = TRUE;
	HANDLE fileHandle = ::CreateFile(szPathName, GENERIC_READ, 0, &sa, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);

	bool retValue;
	if (fileHandle == INVALID_HANDLE_VALUE)
	{
		char buffer[100];
		snprintf(buffer, sizeof(buffer), "Could not open file (error %d)\n", GetLastError());
		retValue = false;
	}
	else
	{
		retValue = true;
		::CloseHandle(fileHandle);
	}
	return retValue;
}
#else
bool FileExists(const char *szPathName)
{
  FILE* fp;
  fp = fopen(szPathName, "r");
  if (fp == NULL) {
    return false;
  } else {
    fclose(fp);
    return true;
  }
}
#endif

#if defined(_WIN32) || defined(__CYGWIN32__)
// DeleteFile defined in <windows.h>
#else
int DeleteFile(const char* szPathName)
{
  if (remove(szPathName) == 0)
  {
    return 1;
  }
  return 0; // failure
}
#endif

#if defined(_WIN32) || defined(__CYGWIN32__)
size_t FileSize(const char *szPathName)
{
	HANDLE hFile = ::CreateFile(
		szPathName,            // file to open
		GENERIC_READ,          // open for reading
		FILE_SHARE_READ,       // share for reading
		NULL,                  // default security
		OPEN_EXISTING,         // existing file only
		FILE_ATTRIBUTE_NORMAL, // normal file
		NULL);                 // no attr. template

	if (hFile != INVALID_HANDLE_VALUE)
	{
		// read file size
		LARGE_INTEGER liFileSize;
		::GetFileSizeEx(hFile, &liFileSize);
		::CloseHandle(hFile);
		return (size_t) liFileSize.QuadPart;
	}
	return 0;
}
#else
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

size_t FileSize(const char *szPathName)
{
	struct stat s;
	stat(szPathName, &s);
	return (size_t) s.st_size;
}
#endif

FileTest::FileTest(std::string fn) : _fn(fn)
{
}

FileTest::~FileTest(void)
{
	int max_tries = 100;
	while(::FileExists(_fn.c_str()) && --max_tries)
	{
		::DeleteFile(_fn.c_str());
	}
}

std::string FileTest::GetName(void)
{
	return _fn;
}

bool FileTest::RemoveExisting(void)
{
	int max_tries = 100;
	while(::FileExists(_fn.c_str()) && --max_tries)
	{
		::DeleteFile(_fn.c_str());
	}
	return !::FileExists(_fn.c_str());
}

bool FileTest::VerifyMissing(void)
{
	int max_tries = 100;
	while(::FileExists(_fn.c_str()) && --max_tries);
	return !::FileExists(_fn.c_str());
}

bool FileTest::VerifyExists(void)
{
	int max_tries = 100;
	while(!::FileExists(_fn.c_str()) && --max_tries);
	return ::FileExists(_fn.c_str());
}

bool FileTest::Exists(void)
{
	return ::FileExists(_fn.c_str());
}

int FileTest::Delete(void)
{
	if (::FileExists(_fn.c_str()))
	{
		return ::DeleteFile(_fn.c_str());
	}
	return 1;
}

size_t FileTest::Size(void)
{
	return ::FileSize(_fn.c_str());
}

size_t FileTest::LineCount(void)
{
	size_t nlines = 0;
	if (::FileExists(_fn.c_str()))
	{
		std::ifstream ifs(_fn.c_str(), std::ifstream::in);
		std::string line;
		while (std::getline(ifs, line))
		{
			++nlines;
		}
		ifs.close();
	}
	return nlines;
}
