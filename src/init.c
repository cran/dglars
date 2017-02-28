#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "dglars.h"

static const R_FortranMethodDef FortEntries[] = {
    {"pc_gaussian_c", (DL_FUNC) &F77_SUB(pc_gaussian_c), 21},
    {"pc_gaussian_g", (DL_FUNC) &F77_SUB(pc_gaussian_g), 27},
    {"pc_pois_c", (DL_FUNC) &F77_SUB(pc_pois_c), 26},
    {"pc_pois_g", (DL_FUNC) &F77_SUB(pc_pois_g), 27},
    {"pc_bin_c", (DL_FUNC) &F77_SUB(pc_bin_c), 27},
    {"pc_bin_g", (DL_FUNC) &F77_SUB(pc_bin_g), 28},
    {"pc_gamma_c", (DL_FUNC) &F77_SUB(pc_gamma_c), 26},
    {"pc_gamma_g", (DL_FUNC) &F77_SUB(pc_gamma_g), 27},
    {"pc_invgaus_c", (DL_FUNC) &F77_SUB(pc_invgaus_c), 26},
    {"pc_invgaus_g", (DL_FUNC) &F77_SUB(pc_invgaus_g), 27},
    {"pc_cvdglars", (DL_FUNC) &F77_SUB(pc_cvdglars), 30},
    {"ccd_pois_c", (DL_FUNC) &F77_SUB(ccd_pois_c), 23},
    {"ccd_bin_c", (DL_FUNC) &F77_SUB(ccd_bin_c), 24},
    {"ccd_cvdglars", (DL_FUNC) &F77_SUB(ccd_cvdglars), 28},
    {"predict", (DL_FUNC) &F77_SUB(predict), 14},
    {NULL, NULL, 0}
};

void attribute_visible R_init_dglars(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
