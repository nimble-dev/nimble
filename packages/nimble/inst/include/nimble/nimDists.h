/*
 * NIMBLE: an R package for programming with BUGS models.
 * Copyright (C) 2014-2017 Perry de Valpine, Christopher Paciorek,
 * Daniel Turek, Clifford Anderson-Bergman, Nick Michaud, Fritz Obermeyer,
 * Duncan Temple Lang.
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.R-project.org/Licenses/
 */

#ifndef __NIM_DISTS
#define __NIM_DISTS

#include "NimArr.h"
#include "dists.h"

double nim_dnorm(double x, double mu, double sigma, int give_log);

bool nimIsNA(double x);
bool nimIsNaN(double x);

bool nimAnyNA(NimArr<1, double> &P); 
bool nimAnyNaN(NimArr<1, double> &P);

double nimArr_ddirch(NimArr<1, double> &x, NimArr<1, double> &alpha, int give_log);
void nimArr_rdirch(NimArr<1, double> &ans, NimArr<1, double> &alpha);

double nimArr_dwish_chol(NimArr<2, double> &x, NimArr<2, double> &chol, double df, double scale_param, int give_log, int overwrite_inputs);
void nimArr_rwish_chol(NimArr<2, double> &ans, NimArr<2, double> &chol, double df, double prec_param, int overwrite_inputs);

double nimArr_dinvwish_chol(NimArr<2, double> &x, NimArr<2, double> &chol, double df, double scale_param, int give_log, int overwrite_inputs);
void nimArr_rinvwish_chol(NimArr<2, double> &ans, NimArr<2, double> &chol, double df, double prec_param, int overwrite_inputs);

double nimArr_dmulti(NimArr<1, double> &x, double size, NimArr<1, double> &prob, int give_log);
void nimArr_rmulti(NimArr<1, double> &ans, double size, NimArr<1, double> &prob);

double nimArr_dcat(double x, NimArr<1, double> &prob, int give_log);
double nimArr_rcat(NimArr<1, double> &prob);

double nimArr_dmnorm_chol(NimArr<1, double> &x, NimArr<1, double> &mean, NimArr<2, double> &chol, double prec_param, int give_log, int overwrite_inputs);
void nimArr_rmnorm_chol(NimArr<1, double> &ans, NimArr<1, double> &mean, NimArr<2, double> &chol, double prec_param);

double nimArr_dmvt_chol(NimArr<1, double> &x, NimArr<1, double> &mean, NimArr<2, double> &chol, double df, double prec_param, int give_log, int overwrite_inputs);
void nimArr_rmvt_chol(NimArr<1, double> &ans, NimArr<1, double> &mean, NimArr<2, double> &chol, double df, double prec_param);

double nimArr_dlkj_corr_cholesky(NimArr<2, double> &x, double eta, int p, int give_log);
void nimArr_rlkj_corr_cholesky(NimArr<2, double> &ans, double eta, int p);

double nimArr_dinterval(double x, double t, NimArr<1, double> &c, int give_log);
int nimArr_rinterval(double t, NimArr<1, double> &c);

double nimArr_dinterval(double x, double t, double c, int give_log);
int nimArr_rinterval(double t, double c);

double nimArr_dcar_normal(NimArr<1, double> &x, NimArr<1, double> &adj, NimArr<1, double> &wgts, NimArr<1, double> &num, double tau, int c, int zero_mean, int give_log);
void nimArr_rcar_normal(NimArr<1, double> &ans, NimArr<1, double> &adj, NimArr<1, double> &wgts, NimArr<1, double> &num, double tau, int c, int zero_mean);

double nimArr_dcar_proper(NimArr<1, double> &x, NimArr<1, double> &mu, NimArr<1, double> &C, NimArr<1, double> &adj, NimArr<1, double> &num, NimArr<1, double> &M, double tau, double gamma, NimArr<1, double> &evs, int give_log);
void nimArr_rcar_proper(NimArr<1, double> &ans, NimArr<1, double> &mu, NimArr<1, double> &C, NimArr<1, double> &adj, NimArr<1, double> &num, NimArr<1, double> &M, double tau, double gamma, NimArr<1, double> &evs);

#endif
