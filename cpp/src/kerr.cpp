// kerr.cpp

#include "kerr.hpp"

double kerr_metric_blc(int i, int j, const double &a, const double &r, const double &th){
  if(i == 0 && j == 0){
    return - (1. - 2.*r/(r*r + pow(a*cos(th), 2)));
  }else if(i == 0 && j == 3){
    return - 2.*a*r*pow(sin(th), 2)/(r*r + pow(a*cos(th), 2));
  }else if(i == 1 && j == 1){
    return (r*r + pow(a*cos(th), 2))/(r*r - 2*r + a*a);
  }else if(i == 2 && j == 2){
    return (r*r + pow(a*cos(th), 2));
  }else if(i == 3 && j == 0){
    return - 2.*a*r*pow(sin(th), 2)/(r*r + pow(a*cos(th), 2));
  }else if(i == 3 && j == 3){
    return pow(sin(th), 2)*(r*r + a*a + 2.*r*a*a*pow(sin(th), 2)/(r*r + pow(a*cos(th), 2)));
  }else{
    return 0.;
  }
}

double kerr_metric_blc_z(int i, int j, const double &a, const double &r, const double &z){
  if(i == 0 && j == 0){
    return - (1. - 2.*r/(r*r + pow(a*z, 2)));
  }else if(i == 0 && j == 3){
    return - 2.*a*r*(1. - z*z)/(r*r + pow(a*z, 2));
  }else if(i == 1 && j == 1){
    return (r*r + pow(a*z, 2))/(r*r - 2*r + a*a);
  }else if(i == 2 && j == 2){
    return (r*r + pow(a*z, 2))/(1. - z*z);
  }else if(i == 3 && j == 0){
    return - 2.*a*r*(1. - z*z)/(r*r + pow(a*z, 2));
  }else if(i == 3 && j == 3){
    return (1. - z*z)*(r*r + a*a + 2.*r*a*a*(1. - z*z)/(r*r + pow(a*z, 2)));
  }else{
    return 0.;
  }
}

double kerr_metric_inverse_blc(int i, int j, const double &a, const double &r, const double &th){
  if(i == 0 && j == 0){
    return - (pow(r*r + a*a, 2) - (r*r - 2.*r + a*a)*pow(a*sin(th), 2))/(r*r + pow(a*cos(th), 2))/(r*r - 2*r + a*a);
  }else if(i == 0 && j == 3){
    return - 2.*a*r/(r*r + pow(a*cos(th), 2))/(r*r - 2*r + a*a);
  }else if(i == 1 && j == 1){
    return (r*r - 2*r + a*a)/(r*r + pow(a*cos(th), 2));
  }else if(i == 2 && j == 2){
    return 1/(r*r + pow(a*cos(th), 2));
  }else if(i == 3 && j == 0){
    return - 2.*a*r/(r*r + pow(a*cos(th), 2))/(r*r - 2*r + a*a);
  }else if(i == 3 && j == 3){
    return (r*r - 2*r + a*a - pow(a*sin(th), 2))/(r*r - 2*r + a*a)/(r*r + pow(a*cos(th), 2))/pow(sin(th), 2);
  }else{
    return 0.;
  }
}

double partial_kerr_metric_blc(int i, int j, int k, const double &a, const double &r, const double &th){
  if(k <= 0 || k >= 3){
    return 0.;
  }

  if(i == 0 && j == 0){
    if(k == 1){
      return - 2.*(r*r - pow(a*cos(th), 2))/pow(r*r + pow(a*cos(th), 2), 2);
    }else{
      return 4.*a*a*r*cos(th)*sin(th)/pow(r*r + pow(a*cos(th), 2), 2);
    }
  }else if(i == 0 && j == 3){
    if(k == 1){
      return 2.*(r*r - pow(a*cos(th), 2))*a*pow(sin(th), 2)/pow(r*r + pow(a*cos(th), 2), 2);
    }else{
      return -4.*a*r*(r*r + a*a)*sin(th)*cos(th)/pow(r*r + pow(a*cos(th), 2), 2);
    }
  }else if(i == 1 && j == 1){
    if(k == 1){
      return -2.*(r*(r - a*a) + (r - 1.)*pow(a*cos(th), 2))/pow(r*r - 2*r + a*a, 2);
    }else{
      return -2.*a*a*cos(th)*sin(th)/(r*r - 2*r + a*a);
    }
  }else if(i == 2 && j == 2){
    if(k == 1){
      return 2.*r;
    }else{
      return -2.*a*a*cos(th)*sin(th);
    }
  }else if(i == 3 && j == 0){
    if(k == 1){
      return 2.*(r*r - pow(a*cos(th), 2))*a*pow(sin(th), 2)/pow(r*r + pow(a*cos(th), 2), 2);
    }else{
      return -4.*a*r*(r*r + a*a)*sin(th)*cos(th)/pow(r*r + pow(a*cos(th), 2), 2);
    }
  }else if(i == 3 && j == 3){
    if(k == 1){
      return (2.*r*(r*r + a*a)*(cos(2.*th)*a*a + r*r) + 2.*a*a*(r*(a*a - r) - (r - 1.)*pow(a*cos(th), 2))*pow(sin(th), 2))*pow(sin(th), 2)/pow(r*r + pow(a*cos(th), 2), 2);
    }else{
      return 2.*cos(th)*sin(th)*(r*r - 2.*r + a*a + 2.*r*pow(r*r + a*a, 2)/pow(r*r + pow(a*cos(th), 2), 2));
    }
  }else{
    return 0.;
  }
}

double kerr_connection_blc(int i, int j, int k, const double &a, const double &r, const double &th){
  if(j > k){ return kerr_connection_blc(i, k, j, a, r, th); }
  double sigma = r*r + pow(a*cos(th), 2);

  if(i == 0){
    if(j == 0 && k == 1){
      return (r*r + a*a)*(r*r - pow(a*cos(th), 2))/pow(sigma, 2)/(r*r - 2.*r + a*a);
    }else if(j == 0 && k == 2){
      return - 2.*a*a*r*sin(th)*cos(th)/pow(sigma, 2);
    }else if(j == 1 && k == 3){
      return -a*pow(sin(th), 2)*(pow(a*cos(th), 2)*(r*r - a*a) + r*r*(3.*r*r + a*a))/pow(sigma, 2)/(r*r - 2.*r + a*a);
    }else if(j == 2 && k == 3){
      return 2.*r*pow(a*sin(th), 3)*cos(th)/pow(sigma, 2);
    }else{
      return 0.;
    }
  }else if(i == 1){
    if(j == 0 && k == 0){
      return (r*r - 2.*r + a*a)*(r*r - pow(a*cos(th), 2))/pow(sigma, 3);
    }else if(j == 0 && k == 3){
      return - (r*r - 2.*r + a*a)*a*pow(sin(th), 2)*(r*r - pow(a*cos(th), 2))/pow(sigma, 3);
    }else if(j == 1 && k == 1){
      return (r*pow(a*sin(th), 2) - (r*r - pow(a*cos(th), 2)))/sigma/(r*r - 2.*r + a*a);
    }else if(j == 1 && k == 2){
      return -a*a*cos(th)*sin(th)/sigma;
    }else if(j == 2 && k == 2){
      return -r*(r*r - 2.*r + a*a)/sigma;
    }else if(j == 3 && k == 3){
      return (r*r - 2.*r + a*a)*pow(sin(th), 2)*(-r*pow(sigma, 2) + pow(a*sin(th), 2)*(r*r - pow(a*cos(th), 2)))/pow(sigma, 3);
    }else{
      return 0.;
    }
  }else if(i == 2){
    if(j == 0 && k == 0){
      return - 2.*a*a*r*sin(th)*cos(th)/pow(sigma, 3);
    }else if(j == 0 && k == 3){
      return 2.*a*r*(r*r + a*a)*sin(th)*cos(th)/pow(sigma, 3);
    }else if(j == 1 && k == 1){
      return a*a*sin(th)*cos(th)/sigma/(r*r - 2.*r + a*a);
    }else if(j == 1 && k == 2){
      return r/sigma;
    }else if(j == 2 && k == 2){
      return -a*a*sin(th)*cos(th)/sigma;
    }else if(j == 3 && k == 3){
      return -sin(th)*cos(th)*((pow(r*r + a*a, 2) - (r*r - 2.*r + a*a)*pow(a*sin(th), 2))*sigma + 2.*(r*r + a*a)*r*pow(a*sin(th), 2))/pow(sigma, 3);
    }else{
      return 0.;
    }
  }else if(i == 3){
    if(j == 0 && k == 1){
      return a*(r*r - pow(a*cos(th), 2))/pow(sigma, 2)/(r*r - 2.*r + a*a);
    }else if(j == 0 && k == 2){
      return -2.*a*r*cot(th)/pow(sigma, 2);
    }else if(j == 1 && k == 3){
      return (r*pow(sigma, 2) + pow(a*sin(th), 2)*pow(a*cos(th), 2) - r*r*(sigma + r*r + a*a))/(r*r - 2.*r + a*a)/pow(sigma, 2);
    }else if(j == 2 && k == 3){
      return cot(th)*(pow(sigma, 2) + 2.*r*pow(a*sin(th), 2))/pow(sigma, 2);
    }else{
      return 0.;
    }
  }else{
    return 0.;
  }
}
