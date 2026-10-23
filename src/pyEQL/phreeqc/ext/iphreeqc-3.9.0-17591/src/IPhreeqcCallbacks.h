#ifndef _INC_IPHREEQC_CALLBACKS_H
#define _INC_IPHREEQC_CALLBACKS_H


#if defined(__cplusplus)
extern "C" {
#endif


typedef int (*PFN_PRERUN_CALLBACK)(void *cookie);
typedef int (*PFN_POSTRUN_CALLBACK)(void *cookie);
typedef int (*PFN_CATCH_CALLBACK)(void *cookie);


#if defined(__cplusplus)
}
#endif

#endif /* _INC_IPHREEQC_CALLBACKS_H */
