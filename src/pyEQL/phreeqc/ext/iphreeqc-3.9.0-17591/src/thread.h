#if defined(WIN32)

	#include <windows.h> 

/*	This is CriticalSection version */

/*
	typedef CRITICAL_SECTION mutex_t;

	#define mutex_init(m)		InitializeCriticalSection(m)
	#define mutex_delete(m)		DeleteCriticalSection(m)
	#define mutex_lock(m)		EnterCriticalSection(m)
	#define mutex_trylock(m)	TryEnterCriticalSection(m)
	#define mutex_unlock(m)		LeaveCriticalSection(m)
	#define atomic_xchg_1(m)	InterlockedExchange(m, 1)
	mutex_t map_lock;

	#define MUTEX_INITILIZER(m) \
		static unsigned int volatile m_init = 0; \
		static bool volatile m_wait = true; \
		if (!InterlockedExchange(&m_init, 1)) \
		{ mutex_init(&m); m_wait = false; } \
		while(m_wait);  
*/

/*	This is simplified spinlock version */

	typedef LONG volatile mutex_t;

	#define mutex_init(m)		InterlockedExchange(m, 0L)
	#define mutex_delete(m)
	#define mutex_lock(m)		while(InterlockedExchange(m, 1L)) Sleep(0)
	#define mutex_trylock(m)	InterlockedExchange(m, 1L)
	#define mutex_unlock(m)		InterlockedExchange(m, 0L)
	#define MUTEX_INITIALIZER	0L

#else

	#include <pthread.h>

	typedef pthread_mutex_t mutex_t;

	#define mutex_init(m)		pthread_mutex_init(m, NULL)
	#define mutex_delete(m)		pthread_mutex_destroy(m)
	#define mutex_lock(m)		pthread_mutex_lock(m) 
	#define mutex_trylock(m)	pthread_mutex_trylock(m) 
	#define mutex_unlock(m)		pthread_mutex_unlock(m)
	#define MUTEX_INITIALIZER	PTHREAD_MUTEX_INITIALIZER

#endif

#if !defined (_INC_PHREEQC_H)  || defined (PHREEQC) || defined (PHREEQC_PARALLEL)
	mutex_t map_lock = MUTEX_INITIALIZER;
	mutex_t qsort_lock = MUTEX_INITIALIZER;
#else
	extern mutex_t map_lock;
	extern mutex_t qsort_lock;
#endif

#ifndef WIN32
	#define qsort(a, b, c, d) mutex_lock(&qsort_lock); \
				qsort(a, b, c, d); \
				mutex_unlock(&qsort_lock)
#endif
