#ifndef __NIM_DISTS
#define __NIM_DISTS

#include "NimArr.h"
#include "dists.h"

double nimArr_ddirch(NimArr<1, double> &x, NimArr<1, double> &alpha, int give_log);
void nimArr_rdirch(NimArr<1, double> &ans, NimArr<1, double> &alpha);

double nimArr_dwish_chol(NimArr<2, double> &x, NimArr<2, double> &chol, double df, int scale_param, int give_log);
void nimArr_rwish_chol(NimArr<2, double> &ans, NimArr<2, double> &chol, double df, int prec_param);

double nimArr_dmulti(NimArr<1, double> &x, int size, NimArr<1, double> &prob, int give_log);
void nimArr_rmulti(NimArr<1, double> &ans, int size, NimArr<1, double> &prob);

double nimArr_dcat(int x, NimArr<1, double> &prob, int give_log);
int nimArr_rcat(NimArr<1, double> &prob);

double nimArr_dmnorm_chol(NimArr<1, double> &x, NimArr<1, double> &mean, NimArr<2, double> &chol, int prec_param, int give_log );
void nimArr_rmnorm_chol(NimArr<1, double> &ans, NimArr<1, double> &mean, NimArr<2, double> &chol, int prec_param);

#endif
