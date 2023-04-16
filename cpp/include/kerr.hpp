// kerr.hpp

#ifndef KERR_HPP
#define KERR_HPP

#include "specialfunc.hpp"

double kerr_metric_blc(int i, int j, const double &a, const double &r, const double &th);
double kerr_metric_blc_z(int i, int j, const double &a, const double &r, const double &z);
double partial_kerr_metric_blc(int i, int j, int k, const double &a, const double &r, const double &th);
double kerr_connection_blc(int i, int j, int k, const double &a, const double &r, const double &th);

#endif
