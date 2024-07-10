#ifndef POLYNOMIAL_KIT__VISIBILITY_H_
#define POLYNOMIAL_KIT__VISIBILITY_H_

// This logic was borrowed (then namespaced) from the examples on the gcc wiki:
//     https://gcc.gnu.org/wiki/Visibility

#if defined _WIN32 || defined __CYGWIN__
  #ifdef __GNUC__
    #define POLYNOMIAL_KIT_EXPORT __attribute__ ((dllexport))
    #define POLYNOMIAL_KIT_IMPORT __attribute__ ((dllimport))
  #else
    #define POLYNOMIAL_KIT_EXPORT __declspec(dllexport)
    #define POLYNOMIAL_KIT_IMPORT __declspec(dllimport)
  #endif
  #ifdef POLYNOMIAL_KIT_BUILDING_LIBRARY
    #define POLYNOMIAL_KIT_PUBLIC POLYNOMIAL_KIT_EXPORT
  #else
    #define POLYNOMIAL_KIT_PUBLIC POLYNOMIAL_KIT_IMPORT
  #endif
  #define POLYNOMIAL_KIT_PUBLIC_TYPE POLYNOMIAL_KIT_PUBLIC
  #define POLYNOMIAL_KIT_LOCAL
#else
  #define POLYNOMIAL_KIT_EXPORT __attribute__ ((visibility("default")))
  #define POLYNOMIAL_KIT_IMPORT
  #if __GNUC__ >= 4
    #define POLYNOMIAL_KIT_PUBLIC __attribute__ ((visibility("default")))
    #define POLYNOMIAL_KIT_LOCAL  __attribute__ ((visibility("hidden")))
  #else
    #define POLYNOMIAL_KIT_PUBLIC
    #define POLYNOMIAL_KIT_LOCAL
  #endif
  #define POLYNOMIAL_KIT_PUBLIC_TYPE
#endif

#endif  // POLYNOMIAL_KIT__VISIBILITY_H_
