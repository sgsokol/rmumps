#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* MUMPS' double-precision C entry point (dmumps_c), compiled into this
 * package from the bundled MUMPS sources. Declared here only so we can take
 * its address for R_RegisterCCallable() below; the full prototype, together
 * with DMUMPS_STRUC_C, lives in dmumps_c.h (shipped in inst/include and made
 * available to client packages via LinkingTo: rmumps). */
extern void dmumps_c(void);

void R_init_rmumps(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
    /* Expose dmumps_c to other R packages on every platform via
     * R_GetCCallable("rmumps", "dmumps_c"). This is required on Windows,
     * where R_FindSymbol() cannot locate it: a Windows DLL exports only the
     * symbols it explicitly registers, so an unregistered dmumps_c is
     * invisible to a dependent package even though it is present in the DLL.
     * The LinkingTo: rmumps + R_GetCCallable() pattern is used by e.g. the
     * multbxxc and r2sundials packages. */
    R_RegisterCCallable("rmumps", "dmumps_c", (DL_FUNC) dmumps_c);
}
