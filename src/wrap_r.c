#include <R_ext/RS.h>
#include <R_ext/Random.h>

void F77_SUB(getrngseed)(void) {
  GetRNGstate();
}
void F77_SUB(putrngseed)(void) {
  PutRNGstate();
}
double F77_SUB(unifrand)(void) {
  return(unif_rand());
}
double unifCrand(void) {
    double res;
    GetRNGstate();
    res=unif_rand();
    PutRNGstate();
    return res;
}
