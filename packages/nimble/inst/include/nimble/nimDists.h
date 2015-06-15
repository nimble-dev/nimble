#ifndef __NIM_DISTS
#define __NIM_DISTS

#include "NimArr.h"
#include "dists.h"

bool R_IsNA(NimArr<1, double> &P); // We use ISNA which is a macro for R_IsNA, so now to overload for the vector case, we muse overload R_IsNA
bool R_isnancpp(NimArr<1, double> &P);

double nimArr_ddirch(NimArr<1, double> &x, NimArr<1, double> &alpha, int give_log);
void nimArr_rdirch(NimArr<1, double> &ans, NimArr<1, double> &alpha);

double nimArr_dwish_chol(NimArr<2, double> &x, NimArr<2, double> &chol, double df, double scale_param, int give_log);
void nimArr_rwish_chol(NimArr<2, double> &ans, NimArr<2, double> &chol, double df, double prec_param);

double nimArr_dmulti(NimArr<1, double> &x, double size, NimArr<1, double> &prob, int give_log);
void nimArr_rmulti(NimArr<1, double> &ans, double size, NimArr<1, double> &prob);

double nimArr_dcat(double x, NimArr<1, double> &prob, int give_log);
double nimArr_rcat(NimArr<1, double> &prob);

double nimArr_dmnorm_chol(NimArr<1, double> &x, NimArr<1, double> &mean, NimArr<2, double> &chol, double prec_param, int give_log );
void nimArr_rmnorm_chol(NimArr<1, double> &ans, NimArr<1, double> &mean, NimArr<2, double> &chol, double prec_param);

double nimArr_dinterval(double x, double t, NimArr<1, double> &c, int give_log);
int nimArr_rinterval(double t, NimArr<1, double> &c);

double nimArr_dinterval(double x, double t, double c, int give_log);
int nimArr_rinterval(double t, double c);


#endif
