#pragma once

#include <vector>
#include <numeric>
#include <complex>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_result.h>
#include <complex_bessel.h>

using namespace sp_bessel;

std::complex<double> operator*(const std::complex<double>& a, const double& b);
std::complex<double> operator*(const double& a,const std::complex<double>& b);
std::complex<double> operator/(const std::complex<double>& a, const double& b);
std::complex<double> operator/(const double& a,const std::complex<double>& b);
std::complex<double> operator-(const std::complex<double>& a, const double& b);
std::complex<double> operator-(const double& a,const std::complex<double>& b);
std::complex<double> operator+(const std::complex<double>& a, const double& b);
std::complex<double> operator+(const double& a,const std::complex<double>& b);

bool operator==(const std::complex<double> &z,const double &a);
bool operator!=(const std::complex<double> &z,const double &a);
bool operator==(const double &a,const std::complex<double> &z);
bool operator!=(const double &a,const std::complex<double> &z);

std::complex<double> operator+(const std::complex<double> &z,const int n);
std::complex<double> operator-(const std::complex<double> &z,const int n);
std::complex<double> operator*(const std::complex<double> &z,const int n);
std::complex<double> operator/(const std::complex<double> &z,const int n);
std::complex<double> operator+(const int n,const std::complex<double> &z);
std::complex<double> operator-(const int n,const std::complex<double> &z);
std::complex<double> operator*(const int n,const std::complex<double> &z);
std::complex<double> operator/(const int n,const std::complex<double> &z);
std::complex<double> operator+(const std::complex<double> &z,const unsigned int n);
std::complex<double> operator-(const std::complex<double> &z,const unsigned int n);
std::complex<double> operator*(const std::complex<double> &z,const unsigned int n);
std::complex<double> operator/(const std::complex<double> &z,const unsigned int n);
std::complex<double> operator+(const unsigned int n,const std::complex<double> &z);
std::complex<double> operator-(const unsigned int n,const std::complex<double> &z);
std::complex<double> operator*(const unsigned int n,const std::complex<double> &z);
std::complex<double> operator/(const unsigned int n,const std::complex<double> &z);
std::complex<double> operator*(const std::vector<std::complex<double>> &c1,
                                    const std::vector<std::complex<double>> &c2);

bool operator==(const std::complex<double> &z,const int n);
bool operator!=(const std::complex<double> &z,const int n);
bool operator==(const int n,const std::complex<double> &z);
bool operator!=(const int n,const std::complex<double> &z);
bool operator==(const std::complex<double> &z,const unsigned int n);
bool operator!=(const std::complex<double> &z,const unsigned int n);
bool operator==(const unsigned int n,const std::complex<double> &z);
bool operator!=(const unsigned int n,const std::complex<double> &z);

std::complex<double> pow(const double& a, const std::complex<double> &b);
std::complex<double> pow(const std::complex<double>& a, const double &b);
std::complex<double> pow(const std::complex<double>& a, const std::complex<double> &b);

double inf_norm (const std::complex<double> &z);
bool isfinite (const std::complex<double> &z);

double log_r(double x);
double dilog_r(double x);

const std::complex<double> II(0.0,1.0);
std::complex<double> log_c(std::complex<double> z);
std::complex<double> dilog_c(std::complex<double> z);

std::complex<double> log_c_angle(std::complex<double> z, double angle);
std::complex<double> expm1 (const std::complex<double> &z);
std::complex<double> log1p (const std::complex<double> &z);

std::complex<double> LogGamma(std::complex<double> z);
std::complex<double> Gamma_inv (const std::complex<double> &z);


std::complex<double> LBesselJ(double k, std::complex<double> z);
std::complex<double> LBesselK(double k, std::complex<double> z);
std::complex<double> CBesselK(std::complex<double> nu, std::complex<double> z,
        bool warn=true);

std::complex<double> Hyp2F1 (std::complex<double> a,std::complex<double> b,
        std::complex<double> c,std::complex<double> z);
