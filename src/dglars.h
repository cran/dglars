#include <R_ext/RS.h>

void
F77_SUB(pc_gaussian_c)(int *n, int *p, double *X, double *y,int *nup,
                       double *w, double *b, double *phi,double *ru,
                       double *dev, int *A, int *nv, int *nav,
                       int *nnonzero, double *g_seq, int *mthd,
                       double *g0, double *g_hat, double *eps, int *np,
                       int *conv);

void
F77_SUB(pc_gaussian_g)(int *linkid, int *n, int *p, double *X, double *y,
                       int *nup, double *w, double *b, double *phi,
                       double *ru, double *dev, int *A, int *nv, int *nav,
                       int *nnonzero, double *g_seq, int *mthd, double *g0,
                       double *g_hat, double *dg_max, double *eps, int *np,
                       int *ncrct, double *cf, double *NReps, int *nNR,
                       int *conv);

void
F77_SUB(pc_pois_c)(int *n, int *p, double *X, double *y, int *nup, double *w,
                   double *b, double *phi, double *ru, double *dev, int *A,
                   int *nv, int *nav, int *nnonzero, double *g_seq, int *mthd,
                   double *g0, double *g_hat, double *dg_max, double *eps,
                   int *np, int *ncrct, double *cf, double *NReps, int *nNR,
                   int *conv);

void
F77_SUB(pc_pois_g)(int *linkid, int *n, int *p, double *X, double *y, int *nup, double *w,
                   double *b, double *phi, double *ru, double *dev, int *A,
                   int *nv, int *nav, int *nnonzero, double *g_seq, int *mthd,
                   double *g0, double *g_hat, double *dg_max, double *eps,
                   int *np, int *ncrct, double *cf, double *NReps, int *nNR,
                   int *conv);

void
F77_SUB(pc_bin_c)(int *n, int *p, double *X, double *y, double *mi, int *nup, double *w,
                  double *b, double *phi, double *ru, double *dev, int *A,
                  int *nv, int *nav, int *nnonzero, double *g_seq, int *mthd,
                  double *g0, double *g_hat, double *dg_max, double *eps,
                  int *np, int *ncrct, double *cf, double *NReps, int *nNR,
                  int *conv);

void
F77_SUB(pc_bin_g)(int *linkid, int *n, int *p, double *X, double *y, double *mi,
                  int *nup, double *w, double *b, double *phi, double *ru, double *dev,
                  int *A, int *nv, int *nav, int *nnonzero, double *g_seq, int *mthd,
                  double *g0, double *g_hat, double *dg_max, double *eps,
                  int *np, int *ncrct, double *cf, double *NReps, int *nNR,
                  int *conv);

void
F77_SUB(pc_gamma_c)(int *n, int *p, double *X, double *y, int *nup, double *w,
                   double *b, double *phi, double *ru, double *dev, int *A,
                   int *nv, int *nav, int *nnonzero, double *g_seq, int *mthd,
                   double *g0, double *g_hat, double *dg_max, double *eps,
                   int *np, int *ncrct, double *cf, double *NReps, int *nNR,
                   int *conv);

void
F77_SUB(pc_gamma_g)(int *linkid, int *n, int *p, double *X, double *y, int *nup, double *w,
                    double *b, double *phi, double *ru, double *dev, int *A,
                    int *nv, int *nav, int *nnonzero, double *g_seq, int *mthd,
                    double *g0, double *g_hat, double *dg_max, double *eps,
                    int *np, int *ncrct, double *cf, double *NReps, int *nNR,
                    int *conv);

void
F77_SUB(pc_invgaus_c)(int *n, int *p, double *X, double *y, int *nup, double *w,
                    double *b, double *phi, double *ru, double *dev, int *A,
                    int *nv, int *nav, int *nnonzero, double *g_seq, int *mthd,
                    double *g0, double *g_hat, double *dg_max, double *eps,
                    int *np, int *ncrct, double *cf, double *NReps, int *nNR,
                    int *conv);

void
F77_SUB(pc_invgaus_g)(int *linkid, int *n, int *p, double *X, double *y, int *nup, double *w,
                      double *b, double *phi, double *ru, double *dev, int *A,
                      int *nv, int *nav, int *nnonzero, double *g_seq, int *mthd,
                      double *g0, double *g_hat, double *dg_max, double *eps,
                      int *np, int *ncrct, double *cf, double *NReps, int *nNR,
                      int *conv);

void
F77_SUB(pc_cvdglars)(int *familyid, int *linkid, int *n, int *p, double *X,
                     double *y, double *mi, int *nup, int *A, double *w,
                     int *foldid, int *nfold, int *ng, double *g, double *b,
                     double *phi, double *dev_m, double *dev_v, double *g_hat, int *nv,
                     int *mthd, double *g0, double *dg_max, double *eps, int *np,
                     int *ncrct, double *cf, double *NReps, int *nNR, int *conv);
                     
void
F77_SUB(ccd_pois_c)(int *n, int *p, double *X, double *y, int *nup, double *w, int *np,
                    double *g0, double *g_hat, int *nstp, double *eps, double *NReps,
                    int *nNR, int *mthd, double *b, double *phi, double *ru, double *dev,
                    double *g_seq, int *A, int *nnonzero, int *nav, int *conv);

void
F77_SUB(ccd_bin_c)(int *n, int *p, double *X, double *y, double *mi, int *nup, double *w, int *np,
                    double *g0, double *g_hat, int *nstp, double *eps, double *NReps,
                    int *nNR, int *mthd, double *b, double *phi, double *ru, double *dev,
                    double *g_seq, int *A, int *nnonzero, int *nav, int *conv);

void
F77_SUB(ccd_cvdglars)(int *familyid, int *linkid, int *n, int *p, double *X, double *y,
                      double *mi, int *nup, int *A, double *w, int *foldid, int *nfold,
                      int *np, double *g_seq, double *g0, int *ng, double *g_cv, int *nstp,
                      double *eps, double *NReps, int *nNR, double *g_hat, int *mthd,
                      double *b, double *phi, double *dev_m, double *dev_v, int *conv);

void
F77_SUB(predict)(int *familyid, int *linkid, int *n, int *p, double *X, double *y, double *mi,
                 int *np, double *b, double *g_seq, int *ng, double *g, double *coef, double *phi);
