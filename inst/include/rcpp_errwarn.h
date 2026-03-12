#ifndef __rcpp_errwarn_h__
#define __rcpp_errwarn_h__
#ifdef __cplusplus
extern "C"
{
#endif
  void rcpp_error(const char* format, ...);
  void rcpp_warning(const char* format, ...);
#ifdef __cplusplus
}
#endif
#endif // __rcpp_errwarn_h__
