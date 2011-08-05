#ifndef __CONFIG_H__
#define __CONFIG_H__

#ifdef _MSC_VER
	#include "winsock2.h"
	#include "vigra/windows.h"
#endif

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

#cmakedefine HAVE_UNIX_ISNAN
#cmakedefine HAVE_UNIX_ISINF

#ifdef _WIN32
	#define PSF_EXPORT __declspec( dllexport )
	/* Disable a template related MSVC warning.
	   See: http://www.unknownroad.com/rtfm/VisualStudio/warningC4251.html */
	#pragma warning( disable: 4251 )
#else
	#define PSF_EXPORT
#endif

/*From file qlobal.h from Qt: Use this to avoid unsued variable warnings*/
#define PSF_UNUSED(x) (void)x;

#endif

