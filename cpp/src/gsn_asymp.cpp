// gsn_asymp.cpp

#include "gsn_asymp.hpp"

Complex gsn_asymptotic_derivative_initial_sum_term_1(const double &a, const int &, const int &, const double &omega, const double &, const double &){
	if(std::abs(a) < DBL_EPSILON){
		return 0.;
	}
	double kappa = sqrt(1. - pow(a, 2));
	return -6.*(-1. + pow(a,2) + pow(kappa,2))*pow(omega,3);
}

Complex gsn_asymptotic_derivative_initial_sum_term_2(const double &a, const int &, const int &mTemp, const double &omega, const double &lambda, const double &){
	double kappa = sqrt(1. - pow(a, 2));
	Complex m = Complex(mTemp);
	return 0.25*omega*(Complex(0.,-1.)*pow(lambda,2) + Complex(0.,2.)*lambda*(-1. + 2.*(-1. + pow(a,2) + pow(kappa,2))*pow(omega,2)) +
     12.*omega*(-1. - Complex(0.,1.)*a*m + omega*(Complex(0.,1.) + 4.*(-2. + kappa)*omega)*pow(a,2) + 8.*pow(omega,2) - 4.*kappa*pow(omega,2) - 8.*pow(kappa,2)*pow(omega,2) + 4.*pow(kappa,3)*pow(omega,2)));
}

Complex gsn_asymptotic_derivative_initial_sum_term_3(const double &a, const int &, const int &mTemp, const double &omega, const double &lambda, const double &){
	double kappa = sqrt(1. - pow(a, 2));
	Complex m = Complex(mTemp);
	return 0.125*omega*(pow(lambda,3) - 2.*pow(lambda,2)*(-1. - Complex(0.,3.)*(-1. + kappa)*omega + (-1. + pow(a,2) + pow(kappa,2))*pow(omega,2)) +
     4.*lambda*omega*(7.*a*m + Complex(0.,1.)*(-6. + 3.*kappa + Complex(0.,1.)*omega*(3. + 4.*pow(a,2)) + 3.*omega*(Complex(0.,-1.) + 4.*omega)*pow(kappa,2) + 12.*(-1. + pow(a,2))*pow(omega,2))) +
     24.*pow(omega,2)*(a*m*(Complex(0.,-3.) + Complex(0.,3.)*kappa + omega*(-1. + pow(a,2)) + omega*pow(kappa,2)) + 5.*pow(a,4)*pow(omega,2) +
        3.*(-1. + kappa)*(1. - Complex(0.,1.)*(1. + kappa)*omega + 8.*(-1. + pow(kappa,2))*pow(omega,2)) +
        pow(a,2)*(2. + 3.*kappa*omega*(Complex(0.,-1.) + 8.*omega) - 29.*pow(omega,2) + 5.*pow(kappa,2)*pow(omega,2))));
}

Complex gsn_asymptotic_derivative_initial_sum_term_4(const double &a, const int &, const int &mTemp, const double &omega, const double &lambda, const double &){
	double kappa = sqrt(1. - pow(a, 2));
	Complex m = Complex(mTemp);
	return -0.0026041666666666665*omega*(2.*(Complex(0.,-6.) + 5.*(-1. + kappa)*omega)*pow(lambda,4) + Complex(0.,1.)*pow(lambda,5) +
     4.*pow(lambda,3)*(Complex(0.,-15.) + 46.*(-1. + kappa)*omega + Complex(0.,10.)*a*m*omega - Complex(0.,2.)*(-6. + 11.*pow(a,2) + 6.*pow(kappa,2))*pow(omega,2)) +
     Complex(0.,48.)*lambda*omega*(3.*omega*pow(a,2)*pow(m,2) - 2.*a*m*(10. + Complex(0.,47.)*(-1. + kappa)*omega + (13.*pow(a,2) + 10.*(-1. + pow(kappa,2)))*pow(omega,2)) +
        omega*(75. - Complex(0.,12.)*omega + 16.*omega*(Complex(0.,-1.) + 2.*omega)*pow(kappa,3) + pow(kappa,2)*(12. + Complex(0.,12.)*omega - 192.*pow(omega,2)) + 176.*pow(omega,2) + 31.*pow(a,4)*pow(omega,2) +
           16.*pow(kappa,4)*pow(omega,2) - 4.*kappa*(21. - Complex(0.,4.)*omega + 8.*pow(omega,2)) +
           pow(a,2)*(-24. + Complex(0.,2.)*(-31. + 29.*kappa)*omega + 4.*(-51. + 8.*kappa + 11.*pow(kappa,2))*pow(omega,2)))) -
     8.*pow(lambda,2)*(Complex(0.,8.) + (41. - 41.*kappa + Complex(0.,86.)*a*m)*omega - 2.*(25.*a*(-1. + kappa)*m + Complex(0.,25.)*pow(a,2) + Complex(0.,12.)*(5. - 6.*kappa + pow(kappa,2)))*pow(omega,2) +
        2.*(-1. + kappa)*(49.*pow(a,2) + 24.*(-1. + pow(kappa,2)))*pow(omega,3)) + 288.*pow(omega,2)*
      (Complex(0.,-8.) + omega*(15. - 43.*kappa + 8.*(-2. + 3.*kappa)*pow(a,2) + 28.*pow(kappa,2)) + (Complex(0.,-12.) + 5.*(-1. + kappa)*omega)*pow(a,2)*pow(m,2) -
        Complex(0.,4.)*(4.*pow(a,4) + pow(a,2)*(-7. - 8.*kappa + 7.*pow(kappa,2)) + 4.*(3. - 1.*kappa - 3.*pow(kappa,2) + pow(kappa,3)))*pow(omega,2) -
        2.*a*m*(2. - Complex(0.,2.)*omega*(11. - 12.*kappa + pow(a,2) + pow(kappa,2)) + ((3. + 5.*kappa)*pow(a,2) + 8.*(-1. + pow(kappa,2)))*pow(omega,2)) +
        ((-149. + 37.*kappa)*pow(a,4) + 32.*pow(a,2)*(13. - 12.*kappa - 4.*pow(kappa,2) + pow(kappa,3)) + 16.*(-17. + 22.*kappa + 16.*pow(kappa,2) - 22.*pow(kappa,3) + pow(kappa,4)))*pow(omega,3)));
}

Complex gsn_asymptotic_derivative_initial_sum_term_5(const double &a, const int &, const int &mTemp, const double &omega, const double &lambda, const double &){
	double kappa = sqrt(1. - pow(a, 2));
	Complex m = Complex(mTemp);
	return -0.0026041666666666665*omega*((3. - Complex(0.,2.)*(-3. + 2.*kappa)*omega)*pow(lambda,5) +
     2.*pow(lambda,4)*(6. + (Complex(0.,-13.) + Complex(0.,16.)*kappa + a*m)*omega + 4.*(-9. + 13.*kappa + pow(a,2) - 4.*pow(kappa,2))*pow(omega,2)) -
     Complex(0.,4.)*pow(lambda,3)*(Complex(0.,3.) + 4.*(5. - 11.*kappa + Complex(0.,11.)*a*m)*omega +
        2.*(10.*a*(-3. + 2.*kappa)*m - Complex(0.,1.)*(78. - 124.*kappa + 11.*pow(a,2) + 46.*pow(kappa,2)))*pow(omega,2) -
        4.*((-33. + 16.*kappa)*pow(a,2) + 6.*(3. - 1.*kappa - 3.*pow(kappa,2) + pow(kappa,3)))*pow(omega,3)) -
     8.*omega*pow(lambda,2)*(Complex(0.,-23.) + 120.*omega + 2.*omega*pow(a,2) - 10.*omega*pow(a,2)*pow(m,2) - Complex(0.,300.)*pow(omega,2) - Complex(0.,142.)*pow(a,2)*pow(omega,2) +
        12.*(Complex(0.,-5.) + 22.*omega)*pow(kappa,3)*pow(omega,2) - 4.*omega*pow(kappa,2)*(-19. + Complex(0.,45.)*omega + (36. + 34.*pow(a,2))*pow(omega,2)) +
        a*m*(-43. - Complex(0.,2.)*(-125. + 126.*kappa)*omega - 2.*(-174. + 260.*kappa + 9.*pow(a,2) - 86.*pow(kappa,2))*pow(omega,2)) + 180.*pow(omega,3) - 600.*pow(a,2)*pow(omega,3) +
        100.*pow(a,4)*pow(omega,3) - 36.*pow(kappa,4)*pow(omega,3) + 4.*kappa*(Complex(0.,-6.) - 49.*omega + Complex(0.,15.)*(9. + 2.*pow(a,2))*pow(omega,2) + 2.*(-33. + 98.*pow(a,2))*pow(omega,3))) -
     Complex(0.,48.)*lambda*omega*(-8. - Complex(0.,1.)*omega*(-33. + 8.*pow(a,2)) + 3.*omega*(Complex(0.,13.) + (-6. + 4.*kappa)*omega)*pow(a,2)*pow(m,2) -
        2.*(97. - 156.*kappa + (-73. + 36.*kappa)*pow(a,2) + 64.*pow(kappa,2) - 2.*pow(kappa,3))*pow(omega,2) +
        Complex(0.,1.)*(5.*pow(a,4) + 4.*pow(a,2)*(43. - 38.*kappa + 27.*pow(kappa,2)) + 4.*(-9. - 20.*kappa + 6.*pow(kappa,2) + 20.*pow(kappa,3) + 3.*pow(kappa,4)))*pow(omega,3) -
        2.*a*m*(Complex(0.,-2.) + (-9. + 30.*kappa)*omega + Complex(0.,4.)*(35. - 55.*kappa + 3.*pow(a,2) + 20.*pow(kappa,2))*pow(omega,2) +
           ((-78. + 40.*kappa)*pow(a,2) + 4.*(15. - 7.*kappa - 15.*pow(kappa,2) + 7.*pow(kappa,3)))*pow(omega,3)) +
        2.*(5.*(-25. + 2.*kappa)*pow(a,4) + 4.*pow(a,2)*(99. - 21.*kappa - 39.*pow(kappa,2) + pow(kappa,3)) - 40.*(7. - 2.*kappa - 8.*pow(kappa,2) + 2.*pow(kappa,3) + pow(kappa,4)))*pow(omega,4)) +
     288.*pow(omega,2)*(8. + Complex(0.,1.)*omega*(-25. + 24.*kappa + 16.*pow(a,2)) + omega*pow(a,3)*pow(m,3) -
        4.*(pow(a,2)*(13. - 22.*kappa + 12.*pow(kappa,2)) + 3.*(-2. + 6.*kappa - 7.*pow(kappa,2) + 3.*pow(kappa,3)))*pow(omega,2) +
        pow(a,2)*pow(m,2)*(4. + Complex(0.,1.)*(-37. + 36.*kappa)*omega - 2.*(16. - 26.*kappa + pow(a,2) + 10.*pow(kappa,2))*pow(omega,2)) +
        Complex(0.,1.)*(3.*(-39. + 8.*kappa)*pow(a,4) + 4.*pow(a,2)*(43. + 7.*kappa - 53.*pow(kappa,2) + 7.*pow(kappa,3)) - 4.*(33. - 26.*kappa - 36.*pow(kappa,2) + 26.*pow(kappa,3) + 3.*pow(kappa,4)))*
         pow(omega,3) + a*m*(Complex(0.,4.) + omega*(-11. + 12.*kappa - 4.*pow(a,2)) + Complex(0.,2.)*(52. - 86.*kappa + (5. - 10.*kappa)*pow(a,2) + 28.*pow(kappa,2) + 6.*pow(kappa,3))*pow(omega,2) +
           (pow(a,4) + 4.*pow(a,2)*(7. - 24.*kappa + 13.*pow(kappa,2)) + 4.*(9. - 2.*kappa - 12.*pow(kappa,2) + 2.*pow(kappa,3) + 3.*pow(kappa,4)))*pow(omega,3)) +
        4.*(10.*pow(a,6) + pow(a,4)*(-175. + 87.*kappa + 4.*pow(kappa,2)) - 1.*pow(a,2)*(-341. + 346.*kappa + 132.*pow(kappa,2) - 74.*pow(kappa,3) + pow(kappa,4)) +
           8.*(-23. + 34.*kappa + 20.*pow(kappa,2) - 34.*pow(kappa,3) + 3.*pow(kappa,4)))*pow(omega,4)));
}

Complex gsn_asymptotic_derivative_initial_sum_term_6(const double &a, const int &, const int &mTemp, const double &omega, const double &lambda, const double &){
	double kappa = sqrt(1. - pow(a, 2));
	Complex m = Complex(mTemp);
	return 0.0026041666666666665*pow(omega,2)*(Complex(0.,1.)*(Complex(0.,11.) + a*m - 26.*omega + kappa*(Complex(0.,-13.) + 28.*omega) + 4.*omega*pow(a,2) - 6.*omega*pow(kappa,2))*pow(lambda,5) +
     2.*pow(lambda,4)*(-14. + Complex(0.,37.)*omega + 2.*a*m*(Complex(0.,-1.) + (-6. + 5.*kappa)*omega) - Complex(0.,13.)*omega*pow(a,2) + 6.*omega*(Complex(0.,2.) + 25.*omega)*pow(kappa,2) + 174.*pow(omega,2) -
        50.*pow(a,2)*pow(omega,2) - 18.*pow(kappa,3)*pow(omega,2) + kappa*(26. - Complex(0.,55.)*omega + (-306. + 32.*pow(a,2))*pow(omega,2))) +
     Complex(0.,4.)*pow(lambda,3)*(Complex(0.,-5.) + 80.*omega - 14.*omega*pow(a,2) + 10.*omega*pow(a,2)*pow(m,2) - Complex(0.,576.)*pow(omega,2) + Complex(0.,62.)*pow(a,2)*pow(omega,2) -
        24.*(Complex(0.,-4.) + 7.*omega)*pow(kappa,3)*pow(omega,2) + 6.*omega*pow(kappa,2)*(7. - Complex(0.,92.)*omega + 6.*(8. + pow(a,2))*pow(omega,2)) +
        a*m*(-7. - Complex(0.,2.)*(-85. + 93.*kappa)*omega + 2.*(-124. + 140.*kappa + 9.*pow(a,2) - 36.*pow(kappa,2))*pow(omega,2)) - 300.*pow(omega,3) + 596.*pow(a,2)*pow(omega,3) - 76.*pow(a,4)*pow(omega,3) +
        12.*pow(kappa,4)*pow(omega,3) + kappa*(Complex(0.,-13.) - 122.*omega + Complex(0.,6.)*(172. + 7.*pow(a,2))*pow(omega,2) - 56.*(-3. + 8.*pow(a,2))*pow(omega,3))) +
     8.*pow(lambda,2)*(8. + Complex(0.,1.)*omega*(-55. + 33.*kappa + 27.*pow(a,2) + 24.*pow(kappa,2)) + 2.*omega*(Complex(0.,-17.) + 10.*(-6. + 5.*kappa)*omega)*pow(a,2)*pow(m,2) +
        (474. - 726.*kappa + (-268. + 62.*kappa)*pow(a,2) + 330.*pow(kappa,2) - 78.*pow(kappa,3))*pow(omega,2) +
        Complex(0.,2.)*(51.*pow(a,4) + pow(a,2)*(-157. + 303.*kappa - 120.*pow(kappa,2)) + 6.*(-71. + 128.*kappa - 54.*pow(kappa,2) - 8.*pow(kappa,3) + 5.*pow(kappa,4)))*pow(omega,3) +
        2.*a*m*(Complex(0.,-2.) + 5.*(-7. + 16.*kappa)*omega - Complex(0.,1.)*(-409. + 543.*kappa + 76.*pow(a,2) - 132.*pow(kappa,2))*pow(omega,2) +
           2.*(423. - 753.*kappa + (-83. + 43.*kappa)*pow(a,2) + 387.*pow(kappa,2) - 57.*pow(kappa,3))*pow(omega,3)) -
        4.*((-227. + 122.*kappa)*pow(a,4) - 3.*pow(a,2)*(-231. + 357.*kappa - 111.*pow(kappa,2) + pow(kappa,3)) + 12.*(-13. + 22.*kappa + 8.*pow(kappa,2) - 22.*pow(kappa,3) + 5.*pow(kappa,4)))*pow(omega,4)) +
     288.*pow(omega,2)*(-16. + 24.*kappa + Complex(0.,65.)*omega - Complex(0.,111.)*kappa*omega + Complex(0.,48.)*omega*pow(kappa,2) + 2.*(Complex(0.,-2.) + (-6. + 5.*kappa)*omega)*pow(a,3)*pow(m,3) +
        6.*pow(omega,2) - 10.*kappa*pow(omega,2) - 162.*pow(kappa,2)*pow(omega,2) + 166.*pow(kappa,3)*pow(omega,2) -
        1.*pow(a,2)*pow(m,2)*(4. + Complex(0.,1.)*omega*(-113. + 13.*pow(a,2)) - 2.*omega*(Complex(0.,24.) + 83.*omega)*pow(kappa,2) + 2.*(-79. + 5.*pow(a,2))*pow(omega,2) + 34.*pow(kappa,3)*pow(omega,2) +
           kappa*(-12. + Complex(0.,159.)*omega + (290. + 4.*pow(a,2))*pow(omega,2))) + Complex(0.,336.)*pow(omega,3) - Complex(0.,352.)*kappa*pow(omega,3) +
        (Complex(0.,-53.) + (-410. + 88.*kappa)*omega)*pow(a,6)*pow(omega,3) - Complex(0.,384.)*pow(kappa,2)*pow(omega,3) + Complex(0.,352.)*pow(kappa,3)*pow(omega,3) +
        Complex(0.,48.)*pow(kappa,4)*pow(omega,3) + omega*pow(a,4)*(Complex(0.,8.) + 16.*(-1. + 2.*kappa)*omega - Complex(0.,1.)*(-557. + 191.*kappa + 44.*pow(kappa,2))*pow(omega,2) +
           2.*(1387. - 949.*kappa - 73.*pow(kappa,2) + 19.*pow(kappa,3))*pow(omega,3)) -
        2.*a*m*(kappa*(Complex(0.,-6.) + 19.*omega - Complex(0.,236.)*pow(omega,2)) + (Complex(0.,1.) + (-24. + 31.*kappa)*omega)*pow(a,4)*pow(omega,2) +
           4.*(Complex(0.,29.) - 16.*omega)*pow(kappa,2)*pow(omega,2) - Complex(0.,4.)*pow(kappa,3)*pow(omega,2) + 32.*pow(kappa,4)*pow(omega,3) +
           2.*(Complex(0.,4.) - 9.*omega + Complex(0.,62.)*pow(omega,2) + 16.*pow(omega,3)) +
           pow(a,2)*(Complex(0.,-4.) + 2.*(7. + 3.*kappa)*omega + Complex(0.,1.)*(-1. - 43.*kappa + 38.*pow(kappa,2))*pow(omega,2) - 2.*(-61. + 155.*kappa - 101.*pow(kappa,2) + 7.*pow(kappa,3))*pow(omega,3))) +
        1920.*pow(omega,4) - 3072.*kappa*pow(omega,4) - 1536.*pow(kappa,2)*pow(omega,4) + 3072.*pow(kappa,3)*pow(omega,4) - 384.*pow(kappa,4)*pow(omega,4) +
        pow(a,2)*(8. + Complex(0.,3.)*(-15. + 16.*kappa)*omega + (138. - 248.*kappa + 252.*pow(kappa,2) - 48.*pow(kappa,3))*pow(omega,2) -
           Complex(0.,16.)*(43. - 9.*kappa - 55.*pow(kappa,2) + 17.*pow(kappa,3))*pow(omega,3) + 64.*(-65. + 72.*kappa + 28.*pow(kappa,2) - 24.*pow(kappa,3) + pow(kappa,4))*pow(omega,4))) +
     Complex(0.,48.)*lambda*omega*(-24. + 24.*kappa + Complex(0.,81.)*omega - Complex(0.,111.)*kappa*omega + Complex(0.,24.)*omega*pow(kappa,2) + 3.*omega*pow(a,3)*pow(m,3) - 570.*pow(omega,2) +
        896.*kappa*pow(omega,2) - 486.*pow(kappa,2)*pow(omega,2) + 148.*pow(kappa,3)*pow(omega,2) +
        pow(a,4)*pow(omega,2)*(-98. + Complex(0.,1.)*(-63. + 17.*kappa)*omega + 2.*(-627. + 118.*kappa + 31.*pow(kappa,2))*pow(omega,2)) -
        1.*pow(a,2)*pow(m,2)*(8. + Complex(0.,1.)*(-141. + 143.*kappa)*omega + (50. - 84.*kappa + 22.*pow(a,2) + 46.*pow(kappa,2))*pow(omega,2)) - Complex(0.,180.)*pow(omega,3) -
        Complex(0.,136.)*kappa*pow(omega,3) + Complex(0.,192.)*pow(kappa,2)*pow(omega,3) + Complex(0.,136.)*pow(kappa,3)*pow(omega,3) - Complex(0.,12.)*pow(kappa,4)*pow(omega,3) +
        a*m*(Complex(0.,-4.) + omega*(73. + 6.*pow(a,2)) + Complex(0.,2.)*(-454. + 85.*pow(a,2))*pow(omega,2) - 8.*(Complex(0.,-16.) + 45.*omega)*pow(kappa,3)*pow(omega,2) +
           4.*omega*pow(kappa,2)*(18. - Complex(0.,215.)*omega + 8.*(15. + pow(a,2))*pow(omega,2)) + (-524. + 728.*pow(a,2) - 69.*pow(a,4))*pow(omega,3) + 44.*pow(kappa,4)*pow(omega,3) -
           2.*kappa*(Complex(0.,6.) + 69.*omega + Complex(0.,1.)*(-820. + 13.*pow(a,2))*pow(omega,2) + 12.*(-15. + 22.*pow(a,2))*pow(omega,3))) - 1632.*pow(omega,4) + 576.*kappa*pow(omega,4) +
        64.*pow(a,6)*pow(omega,4) + 1920.*pow(kappa,2)*pow(omega,4) - 576.*pow(kappa,3)*pow(omega,4) - 288.*pow(kappa,4)*pow(omega,4) +
        2.*pow(a,2)*(4. + Complex(0.,4.)*(-5. + 3.*kappa)*omega + (305. - 259.*kappa + 42.*pow(kappa,2))*pow(omega,2) - Complex(0.,4.)*(-49. + 37.*kappa - 31.*pow(kappa,2) + 5.*pow(kappa,3))*pow(omega,3) +
           2.*(689. - 182.*kappa - 360.*pow(kappa,2) + 38.*pow(kappa,3) + 7.*pow(kappa,4))*pow(omega,4))));
}

Complex gsn_asymptotic_derivative_initial_sum_term_7(const double &a, const int &, const int &mTemp, const double &omega, const double &lambda, const double &){
	double kappa = sqrt(1. - pow(a, 2));
	Complex m = Complex(mTemp);
	return 0.0026041666666666665*pow(omega,3)*(Complex(0.,-1.)*(Complex(0.,-36.) + Complex(0.,62.)*kappa + 2.*a*(-3. + 2.*kappa)*m + 96.*omega - 132.*kappa*omega + (Complex(0.,11.) - 34.*omega + 16.*kappa*omega)*pow(a,2) -
        Complex(0.,22.)*pow(kappa,2) + 48.*omega*pow(kappa,2) - 4.*omega*pow(kappa,3))*pow(lambda,5) +
     2.*pow(lambda,4)*(a*m*(Complex(0.,-8.) - 70.*omega + kappa*(Complex(0.,9.) + 92.*omega) + 9.*omega*pow(a,2) - 26.*omega*pow(kappa,2)) + pow(a,2)*pow(m,2) + 16.*pow(a,4)*pow(omega,2) +
        pow(a,2)*(13. + Complex(0.,1.)*(-65. + 38.*kappa)*omega - 4.*(86. - 91.*kappa + 19.*pow(kappa,2))*pow(omega,2)) +
        2.*(-20. + Complex(0.,57.)*omega + kappa*(45. - Complex(0.,96.)*omega - 696.*pow(omega,2)) + 348.*pow(omega,2) - 88.*pow(kappa,3)*pow(omega,2) + 4.*pow(kappa,4)*pow(omega,2) +
           pow(kappa,2)*(-22. + Complex(0.,33.)*omega + 432.*pow(omega,2)))) - Complex(0.,4.)*pow(lambda,3)*
      (2.*(Complex(0.,1.) + 10.*(-3. + 2.*kappa)*omega)*pow(a,2)*pow(m,2) - 2.*(Complex(0.,45.) + 2.*(-163. + 58.*kappa)*omega)*pow(a,4)*pow(omega,2) +
        2.*a*m*(17. + Complex(0.,2.)*omega*(-152. + 41.*pow(a,2)) + 2.*omega*(Complex(0.,-77.) + 138.*omega)*pow(kappa,2) + (444. - 104.*pow(a,2))*pow(omega,2) - 32.*pow(kappa,3)*pow(omega,2) +
           3.*kappa*(-5. + Complex(0.,158.)*omega + 8.*(-27. + 2.*pow(a,2))*pow(omega,2))) -
        1.*pow(a,2)*(Complex(0.,7.) + (-88. + 52.*kappa)*omega + Complex(0.,4.)*(172. - 105.*kappa + 7.*pow(kappa,2))*pow(omega,2) + 8.*(291. - 270.*kappa + 33.*pow(kappa,2) + 4.*pow(kappa,3))*pow(omega,3)) -
        2.*(-1. + kappa)*(4.*(Complex(0.,-5.) + 12.*omega)*pow(kappa,3)*pow(omega,2) - 6.*omega*pow(kappa,2)*(1. - Complex(0.,48.)*omega + 56.*pow(omega,2)) +
           kappa*(Complex(0.,11.) + 96.*omega - Complex(0.,972.)*pow(omega,2) + 144.*pow(omega,3)) + 2.*(Complex(0.,7.) - 75.*omega + Complex(0.,504.)*pow(omega,2) + 264.*pow(omega,3)))) -
     Complex(0.,48.)*lambda*omega*(2.*(Complex(0.,1.) + (-9. + 6.*kappa)*omega)*pow(a,3)*pow(m,3) +
        pow(a,2)*pow(m,2)*(40. + Complex(0.,1.)*omega*(-500. + 137.*pow(a,2)) + 2.*omega*(Complex(0.,-101.) + 156.*omega)*pow(kappa,2) + 6.*(20. + 17.*pow(a,2))*pow(omega,2) - 68.*pow(kappa,3)*pow(omega,2) -
           2.*kappa*(17. - Complex(0.,353.)*omega + 2.*(85. + 8.*pow(a,2))*pow(omega,2))) + (Complex(0.,13.) + (-686. + 72.*kappa)*omega)*pow(a,6)*pow(omega,3) +
        2.*omega*pow(a,4)*(Complex(0.,-9.) + (309. - 136.*kappa)*omega - Complex(0.,1.)*(-204. - 53.*kappa + 11.*pow(kappa,2))*pow(omega,2) +
           (2604. - 706.*kappa - 324.*pow(kappa,2) + 6.*pow(kappa,3))*pow(omega,3)) -
        2.*a*m*(Complex(0.,-14.) + 143.*omega - Complex(0.,1478.)*pow(omega,2) + (Complex(0.,32.) + (-317. + 122.*kappa)*omega)*pow(a,4)*pow(omega,2) +
           2.*(Complex(0.,-9.) + 76.*omega)*pow(kappa,4)*pow(omega,2) - 2.*omega*pow(kappa,3)*(11. - Complex(0.,198.)*omega + 408.*pow(omega,2)) - 952.*pow(omega,3) +
           2.*pow(kappa,2)*(Complex(0.,3.) + 75.*omega - Complex(0.,876.)*pow(omega,2) + 400.*pow(omega,3)) + 2.*kappa*(Complex(0.,5.) - 132.*omega + Complex(0.,1426.)*pow(omega,2) + 408.*pow(omega,3)) +
           pow(a,2)*(Complex(0.,3.) + 9.*(-1. + kappa)*omega + Complex(0.,2.)*(343. - 277.*kappa + 46.*pow(kappa,2))*pow(omega,2) + 4.*(366. - 323.*kappa + 18.*pow(kappa,2) + 17.*pow(kappa,3))*pow(omega,3))) +
        pow(a,2)*(-48. + Complex(0.,193.)*omega - 2098.*pow(omega,2) + 8.*pow(kappa,3)*(12. + Complex(0.,23.)*omega - 124.*pow(omega,2))*pow(omega,2) - Complex(0.,988.)*pow(omega,3) -
           4.*(Complex(0.,-3.) + 52.*omega)*pow(kappa,4)*pow(omega,3) + 4.*omega*pow(kappa,2)*(Complex(0.,12.) - 201.*omega + Complex(0.,60.)*pow(omega,2) + 1360.*pow(omega,3)) - 8816.*pow(omega,4) +
           4.*kappa*(6. - Complex(0.,47.)*omega + 602.*pow(omega,2) - Complex(0.,38.)*pow(omega,3) + 696.*pow(omega,4))) +
        2.*(36. - Complex(0.,114.)*omega + 764.*pow(omega,2) + 8.*pow(kappa,4)*pow(omega,2)*(5. + Complex(0.,8.)*omega + 56.*pow(omega,2)) + Complex(0.,320.)*pow(omega,3) +
           2.*omega*pow(kappa,3)*(Complex(0.,12.) - 145.*omega - Complex(0.,48.)*pow(omega,2) + 448.*pow(omega,3)) +
           kappa*(-48. + Complex(0.,213.)*omega - 1306.*pow(omega,2) + Complex(0.,96.)*pow(omega,3) - 896.*pow(omega,4)) + 2240.*pow(omega,4) -
           3.*pow(kappa,2)*(-4. + Complex(0.,39.)*omega - 268.*pow(omega,2) + Complex(0.,128.)*pow(omega,3) + 896.*pow(omega,4)))) +
     8.*pow(lambda,2)*(10.*omega*pow(a,3)*pow(m,3) + pow(a,2)*pow(m,2)*(1. + Complex(0.,2.)*(-82. + 73.*kappa)*omega + 20.*(-35. + 46.*kappa + 4.*pow(a,2) - 13.*pow(kappa,2))*pow(omega,2)) +
        2.*pow(a,4)*pow(omega,2)*(80. - Complex(0.,1.)*(-449. + 224.*kappa)*omega + (2734. - 2504.*kappa + 362.*pow(kappa,2))*pow(omega,2)) +
        a*m*(Complex(0.,-20.) + omega*(-182. + 57.*pow(a,2)) - Complex(0.,2.)*(-1356. + 499.*pow(a,2))*pow(omega,2) - 4.*(Complex(0.,29.) + 506.*omega)*pow(kappa,3)*pow(omega,2) +
           2.*omega*pow(kappa,2)*(-115. + Complex(0.,792.)*omega + (4392. - 268.*pow(a,2))*pow(omega,2)) - 2.*(-3390. + 1244.*pow(a,2) + pow(a,4))*pow(omega,3) + 116.*pow(kappa,4)*pow(omega,3) +
           kappa*(Complex(0.,17.) + 412.*omega + Complex(0.,2.)*(-2094. + 235.*pow(a,2))*pow(omega,2) + 8.*(-1707. + 307.*pow(a,2))*pow(omega,3))) - 304.*pow(a,6)*pow(omega,4) +
        pow(a,2)*(-9. - Complex(0.,1.)*(-155. + 74.*kappa)*omega - 4.*(346. - 311.*kappa + 40.*pow(kappa,2))*pow(omega,2) +
           Complex(0.,4.)*(-261. + 612.*kappa - 357.*pow(kappa,2) + 80.*pow(kappa,3))*pow(omega,3) + 4.*(-2715. + 4734.*kappa - 1980.*pow(kappa,2) + 98.*pow(kappa,3) + 7.*pow(kappa,4))*pow(omega,4)) -
        2.*(-16. + Complex(0.,63.)*omega - 804.*pow(omega,2) + 4.*pow(kappa,4)*pow(omega,2)*(-4. - Complex(0.,27.)*omega + 132.*pow(omega,2)) +
           4.*omega*pow(kappa,3)*(Complex(0.,1.) + 37.*omega + Complex(0.,18.)*pow(omega,2) - 456.*pow(omega,3)) + Complex(0.,1164.)*pow(omega,3) +
           3.*omega*pow(kappa,2)*(Complex(0.,-3.) - 228.*omega + Complex(0.,320.)*pow(omega,2) + 160.*pow(omega,3)) - 1008.*pow(omega,4) +
           kappa*(17. - Complex(0.,60.)*omega + 1356.*pow(omega,2) - Complex(0.,2088.)*pow(omega,3) + 1824.*pow(omega,4)))) +
     288.*pow(omega,2)*(-40. + Complex(0.,210.)*omega + pow(a,3)*(Complex(0.,-20.) + kappa*(Complex(0.,17.) + 92.*omega) + 7.*omega*(-10. + pow(a,2)) - 26.*omega*pow(kappa,2))*pow(m,3) + pow(a,4)*pow(m,4) +
        304.*pow(omega,2) - 16.*pow(kappa,4)*pow(omega,2)*(3. - Complex(0.,10.)*omega + 80.*pow(omega,2)) -
        1.*pow(a,6)*pow(omega,2)*(17. - Complex(0.,1.)*(-505. + 106.*kappa)*omega + 4.*(611. - 257.*kappa + 8.*pow(kappa,2))*pow(omega,2)) +
        pow(a,2)*pow(m,2)*(-3. + Complex(0.,358.)*omega - 12.*omega*(Complex(0.,3.) + 22.*omega)*pow(kappa,3) + 636.*pow(omega,2) - 25.*pow(a,4)*pow(omega,2) + 20.*pow(kappa,4)*pow(omega,2) +
           6.*pow(kappa,2)*(-2. + Complex(0.,45.)*omega + 152.*pow(omega,2)) - 2.*kappa*(-7. + Complex(0.,294.)*omega + 652.*pow(omega,2)) +
           pow(a,2)*(3. + Complex(0.,3.)*(-35. + 16.*kappa)*omega - 4.*(30. - 23.*kappa + 9.*pow(kappa,2))*pow(omega,2))) + Complex(0.,800.)*pow(omega,3) +
        8.*omega*pow(kappa,3)*(Complex(0.,-7.) + 52.*omega + Complex(0.,120.)*pow(omega,2) + 1024.*pow(omega,3)) -
        2.*omega*pow(a,4)*(Complex(0.,-24.) + 46.*omega - Complex(0.,1131.)*pow(omega,2) - 4.*(Complex(0.,9.) + 61.*omega)*pow(kappa,3)*pow(omega,2) +
           omega*pow(kappa,2)*(64. + Complex(0.,253.)*omega + 360.*pow(omega,2)) - 4954.*pow(omega,3) + 2.*pow(kappa,4)*pow(omega,3) +
           kappa*(Complex(0.,12.) - 85.*omega + Complex(0.,536.)*pow(omega,2) + 4100.*pow(omega,3))) +
        a*m*(Complex(0.,-56.) + 78.*omega - Complex(0.,572.)*pow(omega,2) - 4.*(Complex(0.,5.) + 64.*omega)*pow(kappa,4)*pow(omega,2) + 4.*omega*pow(kappa,3)*(5. + Complex(0.,50.)*omega + 32.*pow(omega,2)) +
           omega*pow(a,4)*(30. - Complex(0.,1.)*(-70. + 83.*kappa)*omega + 2.*(253. - 350.*kappa + 87.*pow(kappa,2))*pow(omega,2)) +
           kappa*(Complex(0.,65.) - 88.*omega + Complex(0.,1144.)*pow(omega,2) - 128.*pow(omega,3)) - 64.*pow(omega,3) - 7.*pow(a,6)*pow(omega,3) +
           2.*pow(kappa,2)*(Complex(0.,-6.) - 7.*omega - Complex(0.,376.)*pow(omega,2) + 160.*pow(omega,3)) +
           pow(a,2)*(Complex(0.,36.) - 151.*omega - Complex(0.,28.)*pow(omega,2) + 4.*(Complex(0.,13.) + 72.*omega)*pow(kappa,3)*pow(omega,2) -
              4.*omega*pow(kappa,2)*(-9. + Complex(0.,85.)*omega + 528.*pow(omega,2)) - 1232.*pow(omega,3) + 16.*pow(kappa,4)*pow(omega,3) +
              4.*kappa*(Complex(0.,-6.) + 20.*omega + Complex(0.,97.)*pow(omega,2) + 744.*pow(omega,3)))) + 4864.*pow(omega,4) + 72.*pow(a,8)*pow(omega,4) -
        2.*pow(kappa,2)*(12. - Complex(0.,129.)*omega + 16.*pow(omega,2) + Complex(0.,480.)*pow(omega,3) + 1792.*pow(omega,4)) -
        2.*kappa*(-31. + Complex(0.,204.)*omega + 320.*pow(omega,2) + Complex(0.,480.)*pow(omega,3) + 4096.*pow(omega,4)) +
        pow(a,2)*(39. - Complex(0.,157.)*omega + 192.*pow(omega,2) - 8.*pow(kappa,3)*pow(omega,2)*(40. + Complex(0.,173.)*omega + 768.*pow(omega,2)) - Complex(0.,2260.)*pow(omega,3) +
           4.*(Complex(0.,1.) + 112.*omega)*pow(kappa,4)*pow(omega,3) + 4.*omega*pow(kappa,2)*(Complex(0.,-18.) + 215.*omega + Complex(0.,740.)*pow(omega,2) + 1360.*pow(omega,3)) - 12032.*pow(omega,4) +
           2.*kappa*(-12. + Complex(0.,107.)*omega - 266.*pow(omega,2) + Complex(0.,500.)*pow(omega,3) + 7168.*pow(omega,4)))));
}

Complex gsn_asymptotic_derivative_initial_sum_term_8(const double &a, const int &, const int &mTemp, const double &omega, const double &lambda, const double &){
	double kappa = sqrt(1. - pow(a, 2));
	Complex m = Complex(mTemp);
	return -0.0026041666666666665*pow(omega,4)*(Complex(0.,1.)*(Complex(0.,-106.) + Complex(0.,226.)*kappa + 321.*omega - 516.*kappa*omega + 10.*omega*pow(a,4) + a*m*(-26. + 28.*kappa + 5.*pow(a,2) - 6.*pow(kappa,2)) -
        Complex(0.,130.)*pow(kappa,2) + 246.*omega*pow(kappa,2) + pow(a,2)*(Complex(0.,61.) - Complex(0.,49.)*kappa - 184.*omega + 152.*kappa*omega - 24.*omega*pow(kappa,2)) + Complex(0.,18.)*pow(kappa,3) -
        36.*omega*pow(kappa,3) + omega*pow(kappa,4))*pow(lambda,5) + 2.*pow(lambda,4)*
      (102. - 282.*kappa - Complex(0.,372.)*omega + Complex(0.,702.)*kappa*omega + 2.*omega*(Complex(0.,-12.) + (-106. + 49.*kappa)*omega)*pow(a,4) + 204.*pow(kappa,2) - Complex(0.,330.)*omega*pow(kappa,2) -
        36.*pow(kappa,3) + Complex(0.,22.)*omega*pow(kappa,3) + a*m*(Complex(0.,31.) + 312.*omega - 1.*kappa*(Complex(0.,49.) + 516.*omega) + 2.*(Complex(0.,-5.) + (-51. + 31.*kappa)*omega)*pow(a,2) +
           16.*(Complex(0.,1.) + 15.*omega)*pow(kappa,2) - 28.*omega*pow(kappa,3)) + Complex(0.,2.)*omega*pow(kappa,4) + (-7. + 5.*kappa)*pow(a,2)*pow(m,2) - 2481.*pow(omega,2) + 5445.*kappa*pow(omega,2) -
        3930.*pow(kappa,2)*pow(omega,2) + 1050.*pow(kappa,3)*pow(omega,2) - 85.*pow(kappa,4)*pow(omega,2) + pow(kappa,5)*pow(omega,2) +
        pow(a,2)*(-53. + Complex(0.,289.)*omega + 2.*omega*(Complex(0.,19.) + 444.*omega)*pow(kappa,2) + kappa*(59. - Complex(0.,281.)*omega - 2436.*pow(omega,2)) + 1776.*pow(omega,2) -
           76.*pow(kappa,3)*pow(omega,2))) + Complex(0.,4.)*pow(lambda,3)*(2.*pow(a,2)*(Complex(0.,7.) + 5.*kappa*(Complex(0.,-1.) + 28.*omega) + 5.*omega*(-26. + 5.*pow(a,2)) - 30.*omega*pow(kappa,2))*
         pow(m,2) + 2.*omega*pow(a,4)*(-3. + Complex(0.,3.)*(-9. + 47.*kappa)*omega + 2.*(895. - 554.*kappa + 51.*pow(kappa,2))*pow(omega,2)) +
        m*(2.*pow(a,5)*pow(omega,2) - 1.*pow(a,3)*(35. + Complex(0.,2.)*(-497. + 361.*kappa)*omega + 8.*(154. - 134.*kappa + 27.*pow(kappa,2))*pow(omega,2)) +
           2.*a*(70. - Complex(0.,992.)*omega - 4.*omega*(Complex(0.,-31.) + 66.*omega)*pow(kappa,3) + kappa*(-91. + Complex(0.,1880.)*omega - 2496.*pow(omega,2)) + 1455.*pow(omega,2) +
              11.*pow(kappa,4)*pow(omega,2) + pow(kappa,2)*(25. - Complex(0.,980.)*omega + 1374.*pow(omega,2)))) - 172.*pow(a,6)*pow(omega,3) +
        (-1. + kappa)*(Complex(0.,-2.)*pow(kappa,4)*pow(omega,2) - 3.*kappa*omega*(303. - Complex(0.,2480.)*omega + 160.*pow(omega,2)) + pow(kappa,3)*(omega + Complex(0.,288.)*pow(omega,2) - 480.*pow(omega,3)) +
           3.*pow(kappa,2)*(Complex(0.,6.) + 43.*omega - Complex(0.,876.)*pow(omega,2) + 800.*pow(omega,3)) - 3.*(Complex(0.,38.) - 369.*omega + Complex(0.,2222.)*pow(omega,2) + 1120.*pow(omega,3))) +
        pow(a,2)*(Complex(0.,-77.) + 554.*omega - Complex(0.,3956.)*pow(omega,2) + 4.*(Complex(0.,23.) - 84.*omega)*pow(kappa,3)*pow(omega,2) -
           4.*omega*pow(kappa,2)*(-19. + Complex(0.,281.)*omega + 327.*pow(omega,2)) - 8238.*pow(omega,3) + 26.*pow(kappa,4)*pow(omega,3) +
           kappa*(Complex(0.,29.) - 514.*omega + Complex(0.,4220.)*pow(omega,2) + 8736.*pow(omega,3)))) +
     Complex(0.,48.)*lambda*omega*(200. - 312.*kappa - Complex(0.,606.)*omega + Complex(0.,1374.)*kappa*omega + 120.*pow(kappa,2) - Complex(0.,1038.)*omega*pow(kappa,2) - 8.*pow(kappa,3) +
        Complex(0.,318.)*omega*pow(kappa,3) - Complex(0.,24.)*omega*pow(kappa,4) + pow(a,3)*(Complex(0.,14.) + 2.*kappa*(Complex(0.,-5.) + 42.*omega) + 3.*omega*(-26. + 5.*pow(a,2)) - 18.*omega*pow(kappa,2))*
         pow(m,3) + 4011.*pow(omega,2) - 7324.*kappa*pow(omega,2) + 4930.*pow(kappa,2)*pow(omega,2) - 1948.*pow(kappa,3)*pow(omega,2) + 379.*pow(kappa,4)*pow(omega,2) +
        pow(a,6)*pow(omega,2)*(-196. - Complex(0.,1.)*(37. + 147.*kappa)*omega + (-4244. + 880.*kappa + 44.*pow(kappa,2))*pow(omega,2)) +
        pow(a,2)*pow(m,2)*(166. - Complex(0.,1642.)*omega + 2.*(Complex(0.,69.) - 250.*omega)*omega*pow(kappa,3) + 263.*pow(omega,2) - 112.*pow(a,4)*pow(omega,2) + 31.*pow(kappa,4)*pow(omega,2) -
           2.*kappa*(105. - Complex(0.,1401.)*omega + 578.*pow(omega,2)) + 2.*pow(kappa,2)*(28. - Complex(0.,645.)*omega + 705.*pow(omega,2)) -
           1.*pow(a,2)*(40. + Complex(0.,1.)*(-821. + 533.*kappa)*omega + 8.*(-52. + 13.*kappa + 15.*pow(kappa,2))*pow(omega,2))) + Complex(0.,1952.)*pow(omega,3) + Complex(0.,192.)*kappa*pow(omega,3) -
        Complex(0.,2496.)*pow(kappa,2)*pow(omega,3) - Complex(0.,192.)*pow(kappa,3)*pow(omega,3) + Complex(0.,544.)*pow(kappa,4)*pow(omega,3) +
        m*(2.*omega*pow(a,5)*(26. + Complex(0.,1.)*(227. + 11.*kappa)*omega + (1849. - 1134.*kappa + 87.*pow(kappa,2))*pow(omega,2)) - 125.*pow(a,7)*pow(omega,3) +
           pow(a,3)*(Complex(0.,-78.) + 367.*omega - Complex(0.,6924.)*pow(omega,2) - 8.*(Complex(0.,-22.) + 153.*omega)*pow(kappa,3)*pow(omega,2) -
              12.*omega*pow(kappa,2)*(-10. + Complex(0.,203.)*omega + 35.*pow(omega,2)) - 10726.*pow(omega,3) + 122.*pow(kappa,4)*pow(omega,3) +
              6.*kappa*(Complex(0.,7.) - 68.*omega + Complex(0.,1328.)*pow(omega,2) + 1780.*pow(omega,3))) +
           2.*a*(Complex(0.,60.) - 557.*omega + omega*pow(kappa,4)*(-6. + Complex(0.,139.)*omega - 704.*pow(omega,2)) + Complex(0.,4651.)*pow(omega,2) - Complex(0.,1.)*pow(kappa,5)*pow(omega,2) +
              pow(kappa,2)*(Complex(0.,26.) - 681.*omega + Complex(0.,6410.)*pow(omega,2) - 2400.*pow(omega,3)) + 3104.*pow(omega,3) +
              2.*pow(kappa,3)*(Complex(0.,1.) + 72.*omega - Complex(0.,855.)*pow(omega,2) + 1536.*pow(omega,3)) - 1.*kappa*(Complex(0.,92.) - 1086.*omega + Complex(0.,9489.)*pow(omega,2) + 3072.*pow(omega,3))))\
         + 11776.*pow(omega,4) - 5120.*kappa*pow(omega,4) + 110.*pow(a,8)*pow(omega,4) - 14336.*pow(kappa,2)*pow(omega,4) + 5120.*pow(kappa,3)*pow(omega,4) + 2560.*pow(kappa,4)*pow(omega,4) +
        pow(a,4)*(24. + Complex(0.,6.)*(-27. + 17.*kappa)*omega + (2944. - 2246.*kappa + 392.*pow(kappa,2))*pow(omega,2) + Complex(0.,2.)*(611. + 795.*kappa - 465.*pow(kappa,2) + 47.*pow(kappa,3))*pow(omega,3) +
           (19295. - 6452.*kappa - 3870.*pow(kappa,2) + 316.*pow(kappa,3) + 39.*pow(kappa,4))*pow(omega,4)) +
        pow(a,2)*(-192. + Complex(0.,699.)*omega - 6800.*pow(omega,2) - 10.*pow(kappa,4)*pow(omega,2)*(6. + Complex(0.,5.)*omega + 96.*pow(omega,2)) - Complex(0.,2306.)*pow(omega,3) -
           Complex(0.,2.)*pow(kappa,5)*pow(omega,3) - 4.*omega*pow(kappa,3)*(Complex(0.,14.) - 227.*omega + Complex(0.,59.)*pow(omega,2) + 1120.*pow(omega,3)) - 26496.*pow(omega,4) +
           4.*pow(kappa,2)*(-6. + Complex(0.,112.)*omega - 1047.*pow(omega,2) + Complex(0.,1133.)*pow(omega,3) + 4560.*pow(omega,4)) +
           kappa*(168. - Complex(0.,1015.)*omega + 9252.*pow(omega,2) - Complex(0.,3538.)*pow(omega,3) + 9600.*pow(omega,4)))) +
     288.*pow(omega,2)*(90. - 142.*kappa - Complex(0.,702.)*omega + Complex(0.,1494.)*kappa*omega + 64.*pow(kappa,2) - Complex(0.,1110.)*omega*pow(kappa,2) - 8.*pow(kappa,3) +
        Complex(0.,334.)*omega*pow(kappa,3) - Complex(0.,24.)*omega*pow(kappa,4) + pow(a,3)*
         (Complex(0.,83.) + 312.*omega - 3.*kappa*(Complex(0.,35.) + 172.*omega) + 4.*(Complex(0.,-5.) + (-22. + 13.*kappa)*omega)*pow(a,2) + 4.*(Complex(0.,7.) + 60.*omega)*pow(kappa,2) -
           28.*omega*pow(kappa,3))*pow(m,3) + (-7. + 5.*kappa)*pow(a,4)*pow(m,4) - 1617.*pow(omega,2) + 3653.*kappa*pow(omega,2) - 1498.*pow(kappa,2)*pow(omega,2) - 742.*pow(kappa,3)*pow(omega,2) +
        203.*pow(kappa,4)*pow(omega,2) + pow(kappa,5)*pow(omega,2) + pow(a,2)*pow(m,2)*
         (-17. - Complex(0.,1122.)*omega + omega*(Complex(0.,8.) + (137. - 93.*kappa)*omega)*pow(a,4) - 3.*omega*(Complex(0.,4.) + 55.*omega)*pow(kappa,4) - 2273.*pow(omega,2) + pow(kappa,5)*pow(omega,2) +
           2.*pow(kappa,3)*(-2. + Complex(0.,119.)*omega + 701.*pow(omega,2)) - 2.*pow(kappa,2)*(-2. + Complex(0.,591.)*omega + 2029.*pow(omega,2)) + kappa*(19. + Complex(0.,2070.)*omega + 5093.*pow(omega,2)) +
           pow(a,2)*(9. + Complex(0.,527.)*omega + 30.*omega*(Complex(0.,3.) + 16.*omega)*pow(kappa,2) + 792.*pow(omega,2) - 92.*pow(kappa,3)*pow(omega,2) -
              1.*kappa*(3. + Complex(0.,471.)*omega + 980.*pow(omega,2)))) - Complex(0.,1792.)*pow(omega,3) + Complex(0.,2304.)*kappa*pow(omega,3) +
        2.*(Complex(0.,-64.) + (-462. + 93.*kappa)*omega)*pow(a,8)*pow(omega,3) + Complex(0.,2304.)*pow(kappa,2)*pow(omega,3) - Complex(0.,2304.)*pow(kappa,3)*pow(omega,3) -
        Complex(0.,512.)*pow(kappa,4)*pow(omega,3) + omega*pow(a,6)*(Complex(0.,24.) + (57. + 25.*kappa)*omega - Complex(0.,1.)*(-2893. + 1073.*kappa + 46.*pow(kappa,2))*pow(omega,2) +
           4.*(2858. - 1711.*kappa + 88.*pow(kappa,2) + 3.*pow(kappa,3))*pow(omega,3)) -
        1.*a*m*(Complex(0.,-183.) + 108.*omega - Complex(0.,1232.)*pow(omega,2) + 2.*(Complex(0.,-6.) + (-121. + 75.*kappa)*omega)*pow(a,6)*pow(omega,2) -
           4.*omega*pow(kappa,4)*(3. + Complex(0.,44.)*omega + 224.*pow(omega,2)) + 128.*pow(omega,3) + 8.*pow(kappa,2)*(Complex(0.,-11.) - 21.*omega - Complex(0.,272.)*pow(omega,2) + 96.*pow(omega,3)) +
           4.*pow(kappa,3)*(Complex(0.,1.) + 31.*omega + Complex(0.,248.)*pow(omega,2) + 192.*pow(omega,3)) - 3.*kappa*(Complex(0.,-87.) + 20.*omega - Complex(0.,864.)*pow(omega,2) + 256.*pow(omega,3)) +
           pow(a,4)*(Complex(0.,-24.) + (234. - 90.*kappa)*omega + Complex(0.,1.)*(607. - 849.*kappa + 216.*pow(kappa,2))*pow(omega,2) -
              4.*(-774. + 1205.*kappa - 468.*pow(kappa,2) + 27.*pow(kappa,3))*pow(omega,3)) +
           2.*pow(a,2)*(Complex(0.,76.) - 267.*omega - Complex(0.,230.)*pow(omega,2) + (Complex(0.,4.) + 43.*omega)*pow(kappa,4)*pow(omega,2) +
              2.*omega*pow(kappa,3)*(-13. + Complex(0.,53.)*omega + 413.*pow(omega,2)) - 2513.*pow(omega,3) + pow(kappa,5)*pow(omega,3) -
              2.*pow(kappa,2)*(Complex(0.,-6.) - 8.*omega + Complex(0.,333.)*pow(omega,2) + 2269.*pow(omega,3)) + kappa*(Complex(0.,-66.) + 253.*omega + Complex(0.,906.)*pow(omega,2) + 6053.*pow(omega,3)))) -
        12032.*pow(omega,4) + 20992.*kappa*pow(omega,4) + 8192.*pow(kappa,2)*pow(omega,4) - 20992.*pow(kappa,3)*pow(omega,4) + 3840.*pow(kappa,4)*pow(omega,4) +
        pow(a,4)*(24. - Complex(0.,4.)*omega*(47. - 48.*kappa + 6.*pow(kappa,2)) + (290. - 668.*kappa + 756.*pow(kappa,2) - 104.*pow(kappa,3))*pow(omega,2) +
           Complex(0.,2.)*(-4109. + 2383.*kappa + 1365.*pow(kappa,2) - 413.*pow(kappa,3) + 2.*pow(kappa,4))*pow(omega,3) +
           (-32801. + 31013.*kappa + 2662.*pow(kappa,2) - 3270.*pow(kappa,3) + 91.*pow(kappa,4) + pow(kappa,5))*pow(omega,4)) -
        1.*pow(a,2)*(123. - Complex(0.,577.)*omega - 316.*pow(omega,2) + 4.*pow(kappa,4)*pow(omega,2)*(17. + Complex(0.,20.)*omega + 544.*pow(omega,2)) - Complex(0.,6672.)*pow(omega,3) -
           4.*omega*pow(kappa,3)*(Complex(0.,-16.) + 333.*omega + Complex(0.,1304.)*pow(omega,2) + 5312.*pow(omega,3)) - 33408.*pow(omega,4) +
           2.*pow(kappa,2)*(12. - Complex(0.,235.)*omega + 1080.*pow(omega,2) + Complex(0.,4448.)*pow(omega,3) + 7680.*pow(omega,4)) +
           3.*kappa*(-43. + Complex(0.,307.)*omega - 52.*pow(omega,2) + Complex(0.,1312.)*pow(omega,3) + 14080.*pow(omega,4)))) +
     8.*pow(lambda,2)*(-110. + 170.*kappa + Complex(0.,186.)*omega - Complex(0.,138.)*kappa*omega - 56.*pow(kappa,2) - Complex(0.,102.)*omega*pow(kappa,2) + Complex(0.,46.)*omega*pow(kappa,3) +
        10.*(-7. + 5.*kappa)*omega*pow(a,3)*pow(m,3) - 5121.*pow(omega,2) + 9381.*kappa*pow(omega,2) - 5274.*pow(kappa,2)*pow(omega,2) + 1146.*pow(kappa,3)*pow(omega,2) - 133.*pow(kappa,4)*pow(omega,2) +
        pow(kappa,5)*pow(omega,2) + pow(a,2)*pow(m,2)*(-7. - Complex(0.,2.)*omega*(-337. + 85.*pow(a,2)) + 4.*omega*(Complex(0.,61.) + 600.*omega)*pow(kappa,2) + (3120. - 950.*pow(a,2))*pow(omega,2) -
           280.*pow(kappa,3)*pow(omega,2) + kappa*(5. - Complex(0.,882.)*omega + 30.*(-172. + 19.*pow(a,2))*pow(omega,2))) + Complex(0.,6144.)*pow(omega,3) - Complex(0.,10944.)*kappa*pow(omega,3) -
        4.*(Complex(0.,-68.) + (-872. + 335.*kappa)*omega)*pow(a,6)*pow(omega,3) + Complex(0.,5184.)*pow(kappa,2)*pow(omega,3) + Complex(0.,192.)*pow(kappa,3)*pow(omega,3) -
        Complex(0.,576.)*pow(kappa,4)*pow(omega,3) + 2.*omega*pow(a,4)*(Complex(0.,50.) + (-772. + 361.*kappa)*omega - Complex(0.,1.)*(2397. - 2361.*kappa + 466.*pow(kappa,2))*pow(omega,2) +
           4.*(-3309. + 3966.*kappa - 1065.*pow(kappa,2) + 50.*pow(kappa,3))*pow(omega,3)) +
        a*m*(Complex(0.,83.) + 448.*omega - Complex(0.,8838.)*pow(omega,2) + 2.*(Complex(0.,-81.) + (-322. + 72.*kappa)*omega)*pow(a,4)*pow(omega,2) - 2.*(Complex(0.,9.) + 545.*omega)*pow(kappa,4)*pow(omega,2) +
           4.*omega*pow(kappa,3)*(-40. + Complex(0.,247.)*omega + 2889.*pow(omega,2)) - 24186.*pow(omega,3) + 10.*pow(kappa,5)*pow(omega,3) -
           4.*pow(kappa,2)*(Complex(0.,-7.) - 205.*omega + Complex(0.,1884.)*pow(omega,2) + 9921.*pow(omega,3)) + kappa*(Complex(0.,-105.) - 1108.*omega + Complex(0.,15420.)*pow(omega,2) + 53394.*pow(omega,3)) -
           2.*pow(a,2)*(Complex(0.,10.) + (69. - 91.*kappa)*omega - Complex(0.,4.)*(617. - 522.*kappa + 68.*pow(kappa,2))*pow(omega,2) +
              12.*(-565. + 747.*kappa - 279.*pow(kappa,2) + 29.*pow(kappa,3))*pow(omega,3))) - 6144.*pow(omega,4) + 11520.*kappa*pow(omega,4) + 2304.*pow(kappa,2)*pow(omega,4) -
        11520.*pow(kappa,3)*pow(omega,4) + 3840.*pow(kappa,4)*pow(omega,4) - 1.*pow(a,2)*
         (-69. + Complex(0.,495.)*omega - 5816.*pow(omega,2) + 4.*pow(kappa,3)*pow(omega,2)*(34. + Complex(0.,475.)*omega + 873.*pow(omega,2)) - Complex(0.,3786.)*pow(omega,3) +
           2.*(Complex(0.,-75.) + 79.*omega)*pow(kappa,4)*pow(omega,3) - 2.*omega*pow(kappa,2)*(Complex(0.,-29.) + 1042.*omega + Complex(0.,3588.)*pow(omega,2) + 18402.*pow(omega,3)) - 38298.*pow(omega,4) +
           10.*pow(kappa,5)*pow(omega,4) + kappa*(39. - Complex(0.,423.)*omega + 7028.*pow(omega,2) + Complex(0.,10188.)*pow(omega,3) + 72978.*pow(omega,4)))));
}

Complex gsn_asymptotic_derivative_initial_sum_term_9(const double &a, const int &, const int &mTemp, const double &omega, const double &lambda, const double &){
	double kappa = sqrt(1. - pow(a, 2));
	Complex m = Complex(mTemp);
	return -0.0026041666666666665*(Complex(0.,-1.)*(Complex(0.,291.) - Complex(0.,720.)*kappa - 1002.*omega + 1800.*kappa*omega + 5.*(Complex(0.,5.) + (-22. + 8.*kappa)*omega)*pow(a,4) + Complex(0.,534.)*pow(kappa,2) -
        1020.*omega*pow(kappa,2) + 4.*a*m*(24. - 33.*kappa + 5.*(-2. + kappa)*pow(a,2) + 12.*pow(kappa,2) - 1.*pow(kappa,3)) - Complex(0.,128.)*pow(kappa,3) + 200.*omega*pow(kappa,3) -
        2.*pow(a,2)*(Complex(0.,123.) - 402.*omega + 4.*kappa*(Complex(0.,-41.) + 111.*omega) + (Complex(0.,43.) - 126.*omega)*pow(kappa,2) + 8.*omega*pow(kappa,3)) + Complex(0.,7.)*pow(kappa,4) -
        10.*omega*pow(kappa,4))*pow(lambda,5) + 2.*pow(lambda,4)*(246. - 804.*kappa - Complex(0.,1233.)*omega + Complex(0.,2580.)*kappa*omega + 732.*pow(kappa,2) - Complex(0.,1542.)*omega*pow(kappa,2) -
        212.*pow(kappa,3) + Complex(0.,244.)*omega*pow(kappa,3) + 14.*pow(kappa,4) - Complex(0.,1.)*omega*pow(kappa,4) +
        a*m*(7.*(Complex(0.,16.) + 171.*omega) - 2.*kappa*(Complex(0.,107.) + 1146.*omega) + 35.*omega*pow(a,4) + 2.*(Complex(0.,56.) + 687.*omega)*pow(kappa,2) -
           2.*pow(a,2)*(Complex(0.,33.) + 329.*omega - 2.*kappa*(Complex(0.,12.) + 169.*omega) + 71.*omega*pow(kappa,2)) - 2.*(Complex(0.,7.) + 138.*omega)*pow(kappa,3) + 13.*omega*pow(kappa,4)) +
        2.*pow(a,2)*(-17. + 20.*kappa + 3.*pow(a,2) - 5.*pow(kappa,2))*pow(m,2) - 8172.*pow(omega,2) + 19260.*kappa*pow(omega,2) + 40.*pow(a,6)*pow(omega,2) - 15480.*pow(kappa,2)*pow(omega,2) +
        4920.*pow(kappa,3)*pow(omega,2) - 540.*pow(kappa,4)*pow(omega,2) + 12.*pow(kappa,5)*pow(omega,2) +
        pow(a,4)*(13. + Complex(0.,1.)*(-221. + 90.*kappa)*omega - 4.*(403. - 331.*kappa + 52.*pow(kappa,2))*pow(omega,2)) +
        2.*pow(a,2)*(-87. + Complex(0.,613.)*omega - 8.*omega*(Complex(0.,1.) + 59.*omega)*pow(kappa,3) + 3864.*pow(omega,2) + 16.*pow(kappa,4)*pow(omega,2) +
           pow(kappa,2)*(-53. + Complex(0.,223.)*omega + 3048.*pow(omega,2)) - 2.*kappa*(-80. + Complex(0.,395.)*omega + 3132.*pow(omega,2)))) -
     Complex(0.,4.)*pow(lambda,3)*(4.*pow(a,2)*(Complex(0.,-17.) + kappa*(Complex(0.,20.) - 330.*omega) + 240.*omega + (Complex(0.,3.) + 50.*(-2. + kappa)*omega)*pow(a,2) +
           5.*(Complex(0.,-1.) + 24.*omega)*pow(kappa,2) - 10.*omega*pow(kappa,3))*pow(m,2) - 2.*(Complex(0.,85.) + (-958. + 272.*kappa)*omega)*pow(a,6)*pow(omega,2) +
        4.*a*m*(-128. + Complex(0.,1524.)*omega + omega*(Complex(0.,98.) + (-79. + 26.*kappa)*omega)*pow(a,4) + (Complex(0.,24.) - 49.*omega)*omega*pow(kappa,4) +
           pow(kappa,2)*(-92. + Complex(0.,2184.)*omega - 2790.*pow(omega,2)) - 2241.*pow(omega,2) + 2.*pow(kappa,3)*(5. - Complex(0.,238.)*omega + 346.*pow(omega,2)) +
           kappa*(206. - Complex(0.,3288.)*omega + 4308.*pow(omega,2)) + pow(a,2)*(63. - Complex(0.,1115.)*omega + omega*(Complex(0.,-313.) + 582.*omega)*pow(kappa,2) + 1410.*pow(omega,2) -
              54.*pow(kappa,3)*pow(omega,2) - 3.*kappa*(13. - Complex(0.,440.)*omega + 558.*pow(omega,2)))) +
        pow(a,4)*(Complex(0.,-49.) + (268. - 84.*kappa)*omega + Complex(0.,8.)*(-247. - 4.*kappa + 42.*pow(kappa,2))*pow(omega,2) + 8.*(-1995. + 1632.*kappa - 261.*pow(kappa,2) + 2.*pow(kappa,3))*pow(omega,3)) +
        (-1. + kappa)*(Complex(0.,381.) - 3888.*omega + Complex(0.,20952.)*pow(omega,2) + Complex(0.,24.)*pow(kappa,4)*pow(omega,2) + 9984.*pow(omega,3) +
           3.*kappa*(Complex(0.,-57.) + 1312.*omega - Complex(0.,8736.)*pow(omega,2) + 128.*pow(omega,3)) + pow(kappa,3)*(Complex(0.,7.) + 32.*omega - Complex(0.,1440.)*pow(omega,2) + 1920.*pow(omega,3)) -
           3.*pow(kappa,2)*(Complex(0.,11.) + 304.*omega - Complex(0.,3504.)*pow(omega,2) + 2560.*pow(omega,3))) +
        2.*pow(a,2)*(Complex(0.,195.) - 1448.*omega + Complex(0.,8925.)*pow(omega,2) + (Complex(0.,21.) - 142.*omega)*pow(kappa,4)*pow(omega,2) +
           4.*omega*pow(kappa,3)*(7. - Complex(0.,163.)*omega + 254.*pow(omega,2)) + 13554.*pow(omega,3) + pow(kappa,2)*(Complex(0.,23.) - 536.*omega + Complex(0.,4830.)*pow(omega,2) + 2700.*pow(omega,3)) -
           4.*kappa*(Complex(0.,43.) - 449.*omega + Complex(0.,3045.)*pow(omega,2) + 3954.*pow(omega,3)))) -
     Complex(0.,48.)*lambda*omega*(-528. + 912.*kappa + Complex(0.,1545.)*omega - Complex(0.,4032.)*kappa*omega - 432.*pow(kappa,2) + Complex(0.,3666.)*omega*pow(kappa,2) + 48.*pow(kappa,3) -
        Complex(0.,1392.)*omega*pow(kappa,3) + Complex(0.,165.)*omega*pow(kappa,4) +
        4.*pow(a,3)*(Complex(0.,-17.) + kappa*(Complex(0.,20.) - 99.*omega) + 72.*omega + 3.*(Complex(0.,1.) + 5.*(-2. + kappa)*omega)*pow(a,2) + (Complex(0.,-5.) + 36.*omega)*pow(kappa,2) -
           3.*omega*pow(kappa,3))*pow(m,3) - 10366.*pow(omega,2) + 19992.*kappa*pow(omega,2) - 14388.*pow(kappa,2)*pow(omega,2) + 5976.*pow(kappa,3)*pow(omega,2) - 1310.*pow(kappa,4)*pow(omega,2) +
        pow(a,2)*pow(m,2)*(-608. + Complex(0.,5115.)*omega + omega*(Complex(0.,325.) + (806. - 224.*kappa)*omega)*pow(a,4) + (Complex(0.,47.) - 254.*omega)*omega*pow(kappa,4) +
           pow(kappa,2)*(-416. + Complex(0.,5862.)*omega - 5300.*pow(omega,2)) - 542.*pow(omega,2) + 4.*pow(kappa,3)*(11. - Complex(0.,276.)*omega + 598.*pow(omega,2)) +
           4.*kappa*(239. - Complex(0.,2484.)*omega + 902.*pow(omega,2)) - 2.*pow(a,2)*
            (-146. + Complex(0.,1859.)*omega + (Complex(0.,407.) - 570.*omega)*omega*pow(kappa,2) + 762.*pow(omega,2) + 124.*pow(kappa,3)*pow(omega,2) -
              4.*kappa*(-22. + Complex(0.,491.)*omega + 18.*pow(omega,2)))) - Complex(0.,5440.)*pow(omega,3) + (Complex(0.,11.) + 2.*(-733. + 84.*kappa)*omega)*pow(a,8)*pow(omega,3) +
        Complex(0.,7296.)*pow(kappa,2)*pow(omega,3) - Complex(0.,1856.)*pow(kappa,4)*pow(omega,3) -
        2.*omega*pow(a,6)*(Complex(0.,25.) + (-881. + 322.*kappa)*omega + Complex(0.,1.)*(-205. - 876.*kappa + 247.*pow(kappa,2))*pow(omega,2) + 2.*(-5103. + 1510.*kappa + 183.*pow(kappa,2))*pow(omega,3)) -
        2.*a*m*(Complex(0.,204.) - 2037.*omega + Complex(0.,14108.)*pow(omega,2) + 2.*(Complex(0.,2.) + (-395. + 121.*kappa)*omega)*pow(a,6)*pow(omega,2) - Complex(0.,12.)*pow(kappa,5)*pow(omega,2) -
           3.*omega*pow(kappa,4)*(19. - Complex(0.,228.)*omega + 896.*pow(omega,2)) + pow(kappa,2)*(Complex(0.,228.) - 2994.*omega + Complex(0.,21944.)*pow(omega,2) - 6720.*pow(omega,3)) + 9408.*pow(omega,3) +
           2.*pow(kappa,3)*(Complex(0.,-16.) + 391.*omega - Complex(0.,3244.)*pow(omega,2) + 5152.*pow(omega,3)) -
           2.*kappa*(Complex(0.,204.) - 2139.*omega + Complex(0.,15118.)*pow(omega,2) + 5152.*pow(omega,3)) +
           pow(a,4)*(Complex(0.,25.) + (-11. + 6.*kappa)*omega + Complex(0.,4.)*(707. - 368.*kappa + 17.*pow(kappa,2))*pow(omega,2) +
              2.*(4314. - 3373.*kappa + 378.*pow(kappa,2) + 55.*pow(kappa,3))*pow(omega,3)) +
           2.*pow(a,2)*(Complex(0.,-102.) + 643.*omega - Complex(0.,7238.)*pow(omega,2) + (Complex(0.,-12.) + 281.*omega)*pow(kappa,4)*pow(omega,2) -
              2.*omega*pow(kappa,3)*(21. - Complex(0.,263.)*omega + 882.*pow(omega,2)) - 9127.*pow(omega,3) - 1.*pow(kappa,2)*(Complex(0.,24.) - 369.*omega + Complex(0.,4074.)*pow(omega,2) + 250.*pow(omega,3)) +
              kappa*(Complex(0.,114.) - 891.*omega + Complex(0.,10030.)*pow(omega,2) + 9916.*pow(omega,3)))) - 29952.*pow(omega,4) + 13824.*kappa*pow(omega,4) + 36864.*pow(kappa,2)*pow(omega,4) -
        13824.*pow(kappa,3)*pow(omega,4) - 6912.*pow(kappa,4)*pow(omega,4) + 2.*pow(a,2)*
         (328. - Complex(0.,1125.)*omega + 10515.*pow(omega,2) + Complex(0.,2388.)*pow(omega,3) + Complex(0.,12.)*pow(kappa,5)*pow(omega,3) +
           omega*pow(kappa,4)*(Complex(0.,-12.) + 265.*omega + Complex(0.,132.)*pow(omega,2) + 1792.*pow(omega,3)) +
           pow(kappa,2)*(96. - Complex(0.,1257.)*omega + 8796.*pow(omega,2) - Complex(0.,12376.)*pow(omega,3) - 28224.*pow(omega,4)) + 37952.*pow(omega,4) +
           4.*pow(kappa,3)*(-1. + Complex(0.,68.)*omega - 602.*pow(omega,2) + Complex(0.,574.)*pow(omega,3) + 2096.*pow(omega,4)) -
           4.*kappa*(93. - Complex(0.,519.)*omega + 4047.*pow(omega,2) - Complex(0.,2335.)*pow(omega,3) + 3824.*pow(omega,4))) +
        pow(a,4)*(-168. + Complex(0.,855.)*omega - 12126.*pow(omega,2) - 4.*pow(kappa,3)*pow(omega,2)*(-103. + Complex(0.,282.)*omega + 682.*pow(omega,2)) - Complex(0.,2375.)*pow(omega,3) +
           (Complex(0.,45.) - 382.*omega)*pow(kappa,4)*pow(omega,3) + 2.*omega*pow(kappa,2)*(Complex(0.,126.) - 1898.*omega + Complex(0.,4509.)*pow(omega,2) + 8870.*pow(omega,3)) - 65822.*pow(omega,4) +
           4.*kappa*(18. - Complex(0.,234.)*omega + 3095.*pow(omega,2) - Complex(0.,2970.)*pow(omega,3) + 6374.*pow(omega,4)))) +
     288.*pow(omega,2)*(192. - 276.*kappa - Complex(0.,2325.)*omega + Complex(0.,5316.)*kappa*omega + 96.*pow(kappa,2) - Complex(0.,4350.)*omega*pow(kappa,2) - 4.*pow(kappa,3) +
        Complex(0.,1508.)*omega*pow(kappa,3) - Complex(0.,165.)*omega*pow(kappa,4) +
        pow(a,3)*(19.*(Complex(0.,16.) + 63.*omega) - 2.*kappa*(Complex(0.,239.) + 1146.*omega) + 23.*omega*pow(a,4) + 2.*(Complex(0.,104.) + 687.*omega)*pow(kappa,2) -
           2.*pow(a,2)*(Complex(0.,73.) + 295.*omega - 2.*kappa*(Complex(0.,22.) + 149.*omega) + 61.*omega*pow(kappa,2)) - 2.*(Complex(0.,11.) + 138.*omega)*pow(kappa,3) + 13.*omega*pow(kappa,4))*pow(m,3) +
        2.*pow(a,4)*(-17. + 20.*kappa + 3.*pow(a,2) - 5.*pow(kappa,2))*pow(m,4) - 6380.*pow(omega,2) + 15228.*kappa*pow(omega,2) - 9144.*pow(kappa,2)*pow(omega,2) - 264.*pow(kappa,3)*pow(omega,2) +
        548.*pow(kappa,4)*pow(omega,2) + 12.*pow(kappa,5)*pow(omega,2) - 1.*pow(a,8)*pow(omega,2)*(37. - Complex(0.,3.)*(-491. + 94.*kappa)*omega + 4.*(1661. - 647.*kappa + 38.*pow(kappa,2))*pow(omega,2)) -
        1.*pow(a,2)*pow(m,2)*(106. + Complex(0.,3453.)*omega + omega*(Complex(0.,93.) + 892.*omega)*pow(kappa,4) + 7500.*pow(omega,2) + 96.*pow(a,6)*pow(omega,2) - 12.*pow(kappa,5)*pow(omega,2) -
           4.*pow(kappa,3)*(5. + Complex(0.,293.)*omega + 1534.*pow(omega,2)) - 4.*kappa*(55. + Complex(0.,1749.)*omega + 4511.*pow(omega,2)) +
           2.*pow(kappa,2)*(65. + Complex(0.,2319.)*omega + 7900.*pow(omega,2)) + pow(a,4)*(13. + Complex(0.,1.)*(101. + 6.*kappa)*omega + (-430. + 532.*kappa - 42.*pow(kappa,2))*pow(omega,2)) -
           2.*pow(a,2)*(48. + Complex(0.,1155.)*omega - 4.*omega*(Complex(0.,11.) + 117.*omega)*pow(kappa,3) + 1987.*pow(omega,2) + 27.*pow(kappa,4)*pow(omega,2) +
              3.*pow(kappa,2)*(5. + Complex(0.,155.)*omega + 606.*pow(omega,2)) - 2.*kappa*(30. + Complex(0.,699.)*omega + 1546.*pow(omega,2)))) - Complex(0.,3776.)*pow(omega,3) +
        Complex(0.,4992.)*kappa*pow(omega,3) + Complex(0.,5376.)*pow(kappa,2)*pow(omega,3) - Complex(0.,4992.)*pow(kappa,3)*pow(omega,3) - Complex(0.,1600.)*pow(kappa,4)*pow(omega,3) +
        2.*omega*pow(a,6)*(Complex(0.,72.) + 27.*omega + Complex(0.,6571.)*pow(omega,2) + 2.*(Complex(0.,31.) + 80.*omega)*pow(kappa,3)*pow(omega,2) +
           omega*pow(kappa,2)*(-75. - Complex(0.,287.)*omega + 1272.*pow(omega,2)) + 23188.*pow(omega,3) + 4.*pow(kappa,4)*pow(omega,3) -
           4.*kappa*(Complex(0.,9.) - 17.*omega + Complex(0.,843.)*pow(omega,2) + 4340.*pow(omega,3))) +
        a*m*(Complex(0.,568.) + 69.*omega + Complex(0.,2464.)*pow(omega,2) + omega*pow(kappa,4)*(85. + Complex(0.,864.)*omega + 2880.*pow(omega,2)) +
           2.*omega*pow(a,6)*(61. - Complex(0.,2.)*(-106. + 63.*kappa)*omega + (1199. - 1106.*kappa + 181.*pow(kappa,2))*pow(omega,2)) +
           pow(kappa,2)*(Complex(0.,424.) + 1086.*omega + Complex(0.,5888.)*pow(omega,2) - 1792.*pow(omega,3)) - 1088.*pow(omega,3) - 13.*pow(a,8)*pow(omega,3) -
           2.*pow(kappa,3)*(Complex(0.,23.) + 306.*omega + Complex(0.,1856.)*pow(omega,2) + 1600.*pow(omega,3)) + kappa*(Complex(0.,-934.) - 612.*omega - Complex(0.,5504.)*pow(omega,2) + 3200.*pow(omega,3)) -
           2.*pow(a,2)*(Complex(0.,287.) - 745.*omega - Complex(0.,1389.)*pow(omega,2) + 3.*omega*pow(kappa,4)*(4. + Complex(0.,17.)*omega + 44.*pow(omega,2)) - 9068.*pow(omega,3) +
              12.*pow(kappa,5)*pow(omega,3) + 2.*pow(kappa,3)*(Complex(0.,-2.) - 59.*omega + Complex(0.,130.)*pow(omega,2) + 1852.*pow(omega,3)) +
              4.*kappa*(Complex(0.,-77.) + 223.*omega + Complex(0.,1001.)*pow(omega,2) + 5519.*pow(omega,3)) - 1.*pow(kappa,2)*(Complex(0.,-78.) + 67.*omega + Complex(0.,2590.)*pow(omega,2) + 17240.*pow(omega,3))
              ) + pow(a,4)*(Complex(0.,144.) - 1141.*omega - Complex(0.,3508.)*pow(omega,2) + 98.*(Complex(0.,1.) + 14.*omega)*pow(kappa,3)*pow(omega,2) -
              6.*omega*pow(kappa,2)*(8. + Complex(0.,330.)*omega + 2059.*pow(omega,2)) - 14827.*pow(omega,3) + 5.*pow(kappa,4)*pow(omega,3) +
              kappa*(Complex(0.,-72.) + 840.*omega + Complex(0.,5538.)*pow(omega,2) + 25436.*pow(omega,3)))) - 29184.*pow(omega,4) + 52224.*kappa*pow(omega,4) + 120.*pow(a,10)*pow(omega,4) +
        18432.*pow(kappa,2)*pow(omega,4) - 52224.*pow(kappa,3)*pow(omega,4) + 10752.*pow(kappa,4)*pow(omega,4) -
        2.*pow(a,2)*(169. - Complex(0.,1117.)*omega - 2012.*pow(omega,2) - Complex(0.,9136.)*pow(omega,3) +
           4.*omega*pow(kappa,4)*(Complex(0.,-3.) + 61.*omega + Complex(0.,52.)*pow(omega,2) + 1080.*pow(omega,3)) - 44896.*pow(omega,4) -
           4.*pow(kappa,3)*(1. - Complex(0.,66.)*omega + 528.*pow(omega,2) + Complex(0.,2096.)*pow(omega,3) + 8368.*pow(omega,4)) +
           pow(kappa,2)*(63. - Complex(0.,1195.)*omega + 1560.*pow(omega,2) + Complex(0.,12448.)*pow(omega,3) + 20608.*pow(omega,4)) +
           2.*kappa*(-102. + Complex(0.,983.)*omega + 1384.*pow(omega,2) + Complex(0.,3136.)*pow(omega,3) + 29792.*pow(omega,4))) +
        pow(a,4)*(131. - Complex(0.,737.)*omega + 272.*pow(omega,2) + pow(kappa,4)*pow(omega,2)*(8. + Complex(0.,75.)*omega + 836.*pow(omega,2)) - Complex(0.,27237.)*pow(omega,3) -
           4.*omega*pow(kappa,3)*(Complex(0.,-2.) + 222.*omega + Complex(0.,1255.)*pow(omega,2) + 4002.*pow(omega,3)) +
           2.*omega*pow(kappa,2)*(Complex(0.,-156.) + 1510.*omega + Complex(0.,5489.)*pow(omega,2) + 4100.*pow(omega,3)) - 102476.*pow(omega,4) + 12.*pow(kappa,5)*pow(omega,4) +
           2.*kappa*(-36. + Complex(0.,481.)*omega - 778.*pow(omega,2) + Complex(0.,9218.)*pow(omega,3) + 53566.*pow(omega,4)))) +
     8.*pow(lambda,2)*(-336. + 636.*kappa - Complex(0.,69.)*omega + Complex(0.,708.)*kappa*omega - 336.*pow(kappa,2) - Complex(0.,990.)*omega*pow(kappa,2) + 44.*pow(kappa,3) +
        Complex(0.,356.)*omega*pow(kappa,3) - Complex(0.,21.)*omega*pow(kappa,4) + 20.*omega*pow(a,3)*(-17. + 20.*kappa + 3.*pow(a,2) - 5.*pow(kappa,2))*pow(m,3) - 15516.*pow(omega,2) +
        30492.*kappa*pow(omega,2) - 18936.*pow(kappa,2)*pow(omega,2) + 4440.*pow(kappa,3)*pow(omega,2) - 492.*pow(kappa,4)*pow(omega,2) + 12.*pow(kappa,5)*pow(omega,2) +
        2.*pow(a,6)*pow(omega,2)*(316. - Complex(0.,1.)*(-1827. + 718.*kappa)*omega + 4.*(3035. - 2141.*kappa + 254.*pow(kappa,2))*pow(omega,2)) +
        2.*pow(a,2)*pow(m,2)*(-17. + Complex(0.,1232.)*omega - 2.*omega*(Complex(0.,49.) + 690.*omega)*pow(kappa,3) + 5985.*pow(omega,2) + 145.*pow(a,4)*pow(omega,2) + 65.*pow(kappa,4)*pow(omega,2) -
           2.*kappa*(-10. + Complex(0.,997.)*omega + 5730.*pow(omega,2)) + pow(kappa,2)*(-5. + Complex(0.,896.)*omega + 6870.*pow(omega,2)) +
           pow(a,2)*(3. + Complex(0.,10.)*(-61. + 38.*kappa)*omega - 60.*(52. - 53.*kappa + 11.*pow(kappa,2))*pow(omega,2))) + Complex(0.,15744.)*pow(omega,3) - Complex(0.,27840.)*kappa*pow(omega,3) +
        Complex(0.,13248.)*pow(kappa,2)*pow(omega,3) + Complex(0.,192.)*pow(kappa,3)*pow(omega,3) - Complex(0.,1344.)*pow(kappa,4)*pow(omega,3) -
        1.*a*m*(Complex(0.,-304.) - 1167.*omega + Complex(0.,28062.)*pow(omega,2) + omega*pow(kappa,4)*(-55. + Complex(0.,270.)*omega + 6456.*pow(omega,2)) +
           omega*pow(a,4)*(41. - Complex(0.,2.)*(-1145. + 278.*kappa)*omega + 12.*(605. - 404.*kappa + 61.*pow(kappa,2))*pow(omega,2)) +
           kappa*(Complex(0.,478.) + 2988.*omega - Complex(0.,54216.)*pow(omega,2) - 188952.*pow(omega,3)) + pow(kappa,3)*(Complex(0.,22.) + 748.*omega - Complex(0.,6056.)*pow(omega,2) - 52848.*pow(omega,3)) +
           79704.*pow(omega,3) + 166.*pow(a,6)*pow(omega,3) - 120.*pow(kappa,5)*pow(omega,3) + 2.*pow(kappa,2)*(Complex(0.,-104.) - 1257.*omega + Complex(0.,15954.)*pow(omega,2) + 77880.*pow(omega,3)) +
           pow(a,2)*(Complex(0.,146.) + 322.*omega - Complex(0.,21592.)*pow(omega,2) + 8.*(Complex(0.,37.) + 1039.*omega)*pow(kappa,3)*pow(omega,2) -
              2.*omega*pow(kappa,2)*(-107. + Complex(0.,3260.)*omega + 24042.*pow(omega,2)) - 61230.*pow(omega,3) - 334.*pow(kappa,4)*pow(omega,3) +
              4.*kappa*(Complex(0.,-22.) - 145.*omega + Complex(0.,6178.)*pow(omega,2) + 24270.*pow(omega,3)))) - 17856.*pow(omega,4) + 34176.*kappa*pow(omega,4) - 664.*pow(a,8)*pow(omega,4) +
        5376.*pow(kappa,2)*pow(omega,4) - 34176.*pow(kappa,3)*pow(omega,4) + 12480.*pow(kappa,4)*pow(omega,4) -
        1.*pow(a,4)*(37. + Complex(0.,1.)*(-567. + 262.*kappa)*omega + 4.*(2140. - 1777.*kappa + 283.*pow(kappa,2))*pow(omega,2) -
           Complex(0.,4.)*(-5528. + 7385.*kappa - 2492.*pow(kappa,2) + 225.*pow(kappa,3))*pow(omega,3) + 8.*(13875. - 19770.*kappa + 7368.*pow(kappa,2) - 694.*pow(kappa,3) + 13.*pow(kappa,4))*pow(omega,4)) -
        2.*pow(a,2)*(-159. + Complex(0.,551.)*omega - 10761.*pow(omega,2) + pow(kappa,4)*pow(omega,2)*(-25. - Complex(0.,483.)*omega + 132.*pow(omega,2)) +
           4.*pow(kappa,3)*pow(omega,2)*(194. + Complex(0.,1165.)*omega + 2526.*pow(omega,2)) - Complex(0.,6843.)*pow(omega,3) - 62700.*pow(omega,4) + 60.*pow(kappa,5)*pow(omega,4) +
           2.*kappa*(84. - Complex(0.,277.)*omega + 7680.*pow(omega,2) + Complex(0.,10050.)*pow(omega,3) + 63942.*pow(omega,4)) -
           1.*pow(kappa,2)*(33. - Complex(0.,101.)*omega + 6222.*pow(omega,2) + Complex(0.,16074.)*pow(omega,3) + 73560.*pow(omega,4)))))*pow(omega,5);
}

Complex gsn_asymptotic_derivative_initial_sum_term_10(const double &a, const int &, const int &mTemp, const double &omega, const double &lambda, const double &){
	double kappa = sqrt(1. - pow(a, 2));
	Complex m = Complex(mTemp);
	return 0.0026041666666666665*(Complex(0.,1.)*(Complex(0.,759.) - Complex(0.,2109.)*kappa - 2972.*omega + 5808.*kappa*omega + 20.*omega*pow(a,6) + Complex(0.,1854.)*pow(kappa,2) - 3720.*omega*pow(kappa,2) -
        5.*pow(a,4)*(Complex(0.,-37.) + kappa*(Complex(0.,23.) - 96.*omega) + 144.*omega + 12.*omega*pow(kappa,2)) - Complex(0.,594.)*pow(kappa,3) + 880.*omega*pow(kappa,3) + Complex(0.,59.)*pow(kappa,4) -
        60.*omega*pow(kappa,4) + a*m*(321. - 516.*kappa + 15.*pow(a,4) + 246.*pow(kappa,2) - 30.*pow(a,2)*(7. - 6.*kappa + pow(kappa,2)) - 36.*pow(kappa,3) + pow(kappa,4)) +
        2.*pow(a,2)*(Complex(0.,-423.) + 1542.*omega - 57.*kappa*(Complex(0.,-13.) + 36.*omega) + (Complex(0.,-339.) + 792.*omega)*pow(kappa,2) + (Complex(0.,37.) - 92.*omega)*pow(kappa,3) +
           2.*omega*pow(kappa,4)) - Complex(0.,1.)*pow(kappa,5))*pow(lambda,5) + 2.*pow(lambda,4)*
      (-570. + 2166.*kappa + Complex(0.,4029.)*omega - Complex(0.,9207.)*kappa*omega + 10.*omega*(Complex(0.,-5.) + (-61. + 22.*kappa)*omega)*pow(a,6) - 2328.*pow(kappa,2) +
        Complex(0.,6570.)*omega*pow(kappa,2) + 880.*pow(kappa,3) - Complex(0.,1590.)*omega*pow(kappa,3) - 102.*pow(kappa,4) + Complex(0.,105.)*omega*pow(kappa,4) + 2.*pow(kappa,5) -
        Complex(0.,3.)*omega*pow(kappa,5) + m*(3.*(Complex(0.,-11.) + 2.*(-74. + 35.*kappa)*omega)*pow(a,5) +
           pow(a,3)*(Complex(0.,326.) + 3264.*omega - 2.*kappa*(Complex(0.,193.) + 2196.*omega) + 4.*(Complex(0.,23.) + 402.*omega)*pow(kappa,2) - 144.*omega*pow(kappa,3)) +
           2.*a*(-6.*(Complex(0.,32.) + 347.*omega) + kappa*(Complex(0.,417.) + 4425.*omega) - 3.*(Complex(0.,93.) + 1040.*omega)*pow(kappa,2) + (Complex(0.,61.) + 830.*omega)*pow(kappa,3) -
              1.*(Complex(0.,3.) + 70.*omega)*pow(kappa,4) + omega*pow(kappa,5))) + 2.*pow(a,2)*(69. - 105.*kappa + 3.*(-9. + 5.*kappa)*pow(a,2) + 45.*pow(kappa,2) - 5.*pow(kappa,3))*pow(m,2) +
        25396.*pow(omega,2) - 63364.*kappa*pow(omega,2) + 55272.*pow(kappa,2)*pow(omega,2) - 19880.*pow(kappa,3)*pow(omega,2) + 2660.*pow(kappa,4)*pow(omega,2) - 84.*pow(kappa,5)*pow(omega,2) +
        pow(a,4)*(-69. + Complex(0.,1411.)*omega + 8.*omega*(Complex(0.,17.) + 369.*omega)*pow(kappa,2) + kappa*(69. - Complex(0.,1093.)*omega - 10248.*pow(omega,2)) + 9228.*pow(omega,2) -
           196.*pow(kappa,3)*pow(omega,2)) + pow(a,2)*(498. - Complex(0.,4935.)*omega - 5.*omega*(Complex(0.,1.) + 86.*omega)*pow(kappa,4) - 29982.*pow(omega,2) + 4.*pow(kappa,5)*pow(omega,2) -
           6.*pow(kappa,2)*(-121. + Complex(0.,540.)*omega + 5350.*pow(omega,2)) + pow(kappa,3)*(-94. + Complex(0.,350.)*omega + 6800.*pow(omega,2)) +
           6.*kappa*(-211. + Complex(0.,1285.)*omega + 9130.*pow(omega,2)))) + Complex(0.,4.)*pow(lambda,3)*
      (2.*pow(a,2)*(-30.*kappa*(Complex(0.,-7.) + 86.*omega) + 3.*(Complex(0.,-46.) + 535.*omega) + 75.*omega*pow(a,4) + 30.*(Complex(0.,-3.) + 41.*omega)*pow(kappa,2) -
           6.*pow(a,2)*(Complex(0.,-9.) + kappa*(Complex(0.,5.) - 150.*omega) + 175.*omega + 25.*omega*pow(kappa,2)) + (Complex(0.,10.) - 180.*omega)*pow(kappa,3) + 5.*omega*pow(kappa,4))*pow(m,2) +
        2.*omega*pow(a,6)*(-8. + Complex(0.,1.)*(-99. + 331.*kappa)*omega + 12.*(528. - 274.*kappa + 23.*pow(kappa,2))*pow(omega,2)) -
        1.*a*m*(1731. - Complex(0.,17934.)*omega + Complex(0.,14.)*omega*pow(kappa,5) + 26360.*pow(omega,2) + 70.*pow(a,6)*pow(omega,2) + 15.*pow(kappa,4)*(1. - Complex(0.,58.)*omega + 72.*pow(omega,2)) -
           4.*pow(kappa,3)*(88. - Complex(0.,2387.)*omega + 2920.*pow(omega,2)) + 6.*pow(kappa,2)*(309. - Complex(0.,5506.)*omega + 6680.*pow(omega,2)) -
           6.*kappa*(536. - Complex(0.,7089.)*omega + 9200.*pow(omega,2)) + pow(a,4)*(111. + Complex(0.,2.)*(-1621. + 893.*kappa)*omega + 12.*(251. - 178.*kappa + 31.*pow(kappa,2))*pow(omega,2)) -
           2.*pow(a,2)*(641. - Complex(0.,8628.)*omega + 4.*(Complex(0.,133.) - 293.*omega)*omega*pow(kappa,3) + 11151.*pow(omega,2) + 39.*pow(kappa,4)*pow(omega,2) -
              4.*kappa*(164. - Complex(0.,3288.)*omega + 4017.*pow(omega,2)) + pow(kappa,2)*(137. - Complex(0.,5376.)*omega + 7410.*pow(omega,2)))) - 320.*pow(a,8)*pow(omega,3) +
        (-1. + kappa)*(Complex(0.,1137.) - 12942.*omega + Complex(0.,63080.)*pow(omega,2) + Complex(0.,1.)*pow(kappa,4)*(-1. + Complex(0.,6.)*omega + 168.*pow(omega,2)) + 28224.*pow(omega,3) -
           6.*kappa*(Complex(0.,143.) - 2568.*omega + Complex(0.,14416.)*pow(omega,2) + 224.*pow(omega,3)) -
           12.*pow(kappa,2)*(Complex(0.,-4.) + 405.*omega - Complex(0.,3204.)*pow(omega,2) + 1904.*pow(omega,3)) + pow(kappa,3)*(Complex(0.,26.) + 384.*omega - Complex(0.,6112.)*pow(omega,2) + 6720.*pow(omega,3))
           ) + 2.*pow(a,2)*(Complex(0.,771.) - 6501.*omega + omega*pow(kappa,4)*(-11. + Complex(0.,279.)*omega - 900.*pow(omega,2)) + Complex(0.,35067.)*pow(omega,2) + Complex(0.,1.)*pow(kappa,5)*pow(omega,2) +
           42172.*pow(omega,3) + 3.*pow(kappa,2)*(Complex(0.,97.) - 1404.*omega + Complex(0.,9274.)*pow(omega,2) + 3320.*pow(omega,3)) +
           pow(kappa,3)*(Complex(0.,-17.) + 518.*omega - Complex(0.,5078.)*pow(omega,2) + 4720.*pow(omega,3)) - 3.*kappa*(Complex(0.,319.) - 3266.*omega + Complex(0.,18617.)*pow(omega,2) + 17648.*pow(omega,3)))\
         + pow(a,4)*(Complex(0.,-417.) + 2606.*omega - Complex(0.,16584.)*pow(omega,2) + 8.*(Complex(0.,-27.) + 8.*omega)*pow(kappa,3)*pow(omega,2) -
           4.*omega*pow(kappa,2)*(-59. + Complex(0.,42.)*omega + 3240.*pow(omega,2)) - 62592.*pow(omega,3) + 32.*pow(kappa,4)*pow(omega,3) +
           kappa*(Complex(0.,207.) - 1898.*omega + Complex(0.,10392.)*pow(omega,2) + 61056.*pow(omega,3)))) +
     Complex(0.,48.)*lambda*omega*(-1344. + 2496.*kappa + Complex(0.,3813.)*omega - Complex(0.,11127.)*kappa*omega - 1344.*pow(kappa,2) + Complex(0.,11514.)*omega*pow(kappa,2) + 192.*pow(kappa,3) -
        Complex(0.,5046.)*omega*pow(kappa,3) + Complex(0.,753.)*omega*pow(kappa,4) - Complex(0.,3.)*omega*pow(kappa,5) +
        pow(a,3)*(Complex(0.,-276.) + kappa*(Complex(0.,420.) - 1548.*omega) + 963.*omega + 45.*omega*pow(a,4) + 18.*(Complex(0.,-10.) + 41.*omega)*pow(kappa,2) -
           6.*pow(a,2)*(Complex(0.,-18.) + Complex(0.,10.)*kappa + 105.*omega - 90.*kappa*omega + 15.*omega*pow(kappa,2)) - 4.*(Complex(0.,-5.) + 27.*omega)*pow(kappa,3) + 3.*omega*pow(kappa,4))*pow(m,3) -
        26452.*pow(omega,2) + 53456.*kappa*pow(omega,2) - 40536.*pow(kappa,2)*pow(omega,2) + 17296.*pow(kappa,3)*pow(omega,2) - 3956.*pow(kappa,4)*pow(omega,2) -
        1.*pow(a,8)*pow(omega,2)*(352. + Complex(0.,1.)*(1. + 473.*kappa)*omega + 4.*(2710. - 588.*kappa + 7.*pow(kappa,2))*pow(omega,2)) -
        1.*pow(a,2)*pow(m,2)*(2052. - Complex(0.,15309.)*omega + Complex(0.,7.)*omega*pow(kappa,5) + pow(kappa,3)*(-388. + Complex(0.,5774.)*omega - 9360.*pow(omega,2)) + 1076.*pow(omega,2) +
           310.*pow(a,6)*pow(omega,2) + pow(kappa,4)*(16. - Complex(0.,457.)*omega + 1300.*pow(omega,2)) + 6.*pow(kappa,2)*(350. - Complex(0.,3791.)*omega + 2980.*pow(omega,2)) -
           3.*kappa*(1244. - Complex(0.,10921.)*omega + 3568.*pow(omega,2)) + 3.*pow(a,4)*(42. + Complex(0.,1.)*(-909. + 443.*kappa)*omega + 4.*(-331. + 138.*kappa + 14.*pow(kappa,2))*pow(omega,2)) -
           2.*pow(a,2)*(746. - Complex(0.,7317.)*omega + (Complex(0.,307.) - 1176.*omega)*omega*pow(kappa,3) + kappa*(-746. + Complex(0.,9975.)*omega - 240.*pow(omega,2)) - 2553.*pow(omega,2) +
              59.*pow(kappa,4)*pow(omega,2) + pow(kappa,2)*(152. - Complex(0.,3597.)*omega + 3318.*pow(omega,2)))) - Complex(0.,14272.)*pow(omega,3) + Complex(0.,640.)*kappa*pow(omega,3) +
        Complex(0.,19968.)*pow(kappa,2)*pow(omega,3) - Complex(0.,640.)*pow(kappa,3)*pow(omega,3) - Complex(0.,5696.)*pow(kappa,4)*pow(omega,3) +
        a*m*(Complex(0.,-1224.) + 13933.*omega + 6.*(-1. + Complex(0.,28.)*omega)*omega*pow(kappa,5) - Complex(0.,82696.)*pow(omega,2) +
           2.*omega*pow(a,6)*(31. - Complex(0.,1.)*(-967. + 85.*kappa)*omega + (5689. - 2938.*kappa + 203.*pow(kappa,2))*pow(omega,2)) - 54080.*pow(omega,3) - 187.*pow(a,8)*pow(omega,3) -
           8.*pow(kappa,3)*(Complex(0.,-65.) + 927.*omega - Complex(0.,5722.)*pow(omega,2) + 7984.*pow(omega,3)) +
           2.*pow(kappa,2)*(Complex(0.,-1044.) + 12123.*omega - Complex(0.,71400.)*pow(omega,2) + 17920.*pow(omega,3)) +
           pow(kappa,4)*(Complex(0.,-32.) + 741.*omega - Complex(0.,5672.)*pow(omega,2) + 18240.*pow(omega,3)) +
           2.*kappa*(Complex(0.,1428.) - 15693.*omega + Complex(0.,92612.)*pow(omega,2) + 31936.*pow(omega,3)) +
           pow(a,4)*(Complex(0.,-434.) + 2019.*omega - Complex(0.,35316.)*pow(omega,2) + 4.*(Complex(0.,17.) - 635.*omega)*pow(kappa,3)*pow(omega,2) -
              6.*omega*pow(kappa,2)*(-70. + Complex(0.,898.)*omega + 1381.*pow(omega,2)) - 70365.*pow(omega,3) + 227.*pow(kappa,4)*pow(omega,3) +
              2.*kappa*(Complex(0.,125.) - 927.*omega + Complex(0.,14886.)*pow(omega,2) + 32058.*pow(omega,3))) +
           2.*pow(a,2)*(Complex(0.,816.) - 6516.*omega + omega*pow(kappa,4)*(-39. + Complex(0.,345.)*omega - 3244.*pow(omega,2)) + Complex(0.,53965.)*pow(omega,2) - Complex(0.,3.)*pow(kappa,5)*pow(omega,2) +
              58580.*pow(omega,3) + pow(kappa,2)*(Complex(0.,504.) - 5469.*omega + Complex(0.,41290.)*pow(omega,2) + 1080.*pow(omega,3)) +
              2.*pow(kappa,3)*(Complex(0.,-24.) + 489.*omega - Complex(0.,3555.)*pow(omega,2) + 8056.*pow(omega,3)) -
              1.*kappa*(Complex(0.,1248.) - 10572.*omega + Complex(0.,84775.)*pow(omega,2) + 68112.*pow(omega,3)))) - 74240.*pow(omega,4) + 35840.*kappa*pow(omega,4) + 172.*pow(a,10)*pow(omega,4) +
        92160.*pow(kappa,2)*pow(omega,4) - 35840.*pow(kappa,3)*pow(omega,4) - 17920.*pow(kappa,4)*pow(omega,4) +
        2.*pow(a,6)*(24. + Complex(0.,23.)*(-11. + 7.*kappa)*omega + (5219. - 3481.*kappa + 538.*pow(kappa,2))*pow(omega,2) +
           Complex(0.,1.)*(233. + 6717.*kappa - 3147.*pow(kappa,2) + 253.*pow(kappa,3))*pow(omega,3) + 2.*(21167. - 7802.*kappa - 1452.*pow(kappa,2) + 90.*pow(kappa,3) + 13.*pow(kappa,4))*pow(omega,4)) +
        2.*pow(a,2)*(1020. - Complex(0.,3345.)*omega + 31383.*pow(omega,2) + 3.*(-1. + Complex(0.,28.)*omega)*pow(kappa,5)*pow(omega,2) + Complex(0.,3996.)*pow(omega,3) +
           omega*pow(kappa,4)*(Complex(0.,-124.) + 1399.*omega + Complex(0.,172.)*pow(omega,2) + 5920.*pow(omega,3)) +
           pow(kappa,2)*(468. - Complex(0.,5433.)*omega + 32802.*pow(omega,2) - Complex(0.,50952.)*pow(omega,3) - 82432.*pow(omega,4)) +
           kappa*(-1356. + Complex(0.,7311.)*omega - 53187.*pow(omega,2) + Complex(0.,37444.)*pow(omega,3) - 45888.*pow(omega,4)) + 104672.*pow(omega,4) +
           pow(kappa,3)*(-36. + Complex(0.,1559.)*omega - 10250.*pow(omega,2) + Complex(0.,13224.)*pow(omega,3) + 27968.*pow(omega,4))) +
        pow(a,4)*(-792. + Complex(0.,3543.)*omega - 45542.*pow(omega,2) + pow(kappa,4)*(-198. + Complex(0.,341.)*omega - 2132.*pow(omega,2))*pow(omega,2) - Complex(0.,359.)*pow(omega,3) -
           Complex(0.,11.)*pow(kappa,5)*pow(omega,3) - 2.*omega*pow(kappa,3)*(Complex(0.,138.) - 2032.*omega + Complex(0.,4931.)*pow(omega,2) + 7672.*pow(omega,3)) - 210868.*pow(omega,4) +
           pow(kappa,2)*(-72. + Complex(0.,2340.)*omega - 23184.*pow(omega,2) + Complex(0.,55234.)*pow(omega,3) + 69480.*pow(omega,4)) +
           kappa*(576. - Complex(0.,5253.)*omega + 56392.*pow(omega,2) - Complex(0.,63487.)*pow(omega,3) + 91344.*pow(omega,4)))) +
     288.*pow(omega,2)*(-396. + 444.*kappa + Complex(0.,7477.)*omega - Complex(0.,18135.)*kappa*omega + 36.*pow(kappa,2) + Complex(0.,15930.)*omega*pow(kappa,2) - 116.*pow(kappa,3) -
        Complex(0.,6038.)*omega*pow(kappa,3) + 16.*pow(kappa,4) + Complex(0.,801.)*omega*pow(kappa,4) - Complex(0.,3.)*omega*pow(kappa,5) +
        (3.*(Complex(0.,-21.) + 2.*(-56. + 25.*kappa)*omega)*pow(a,7) - 2.*pow(a,5)*(Complex(0.,-373.) - 1494.*omega + kappa*(Complex(0.,373.) + 1986.*omega) - 2.*(Complex(0.,38.) + 357.*omega)*pow(kappa,2) +
              62.*omega*pow(kappa,3)) + 2.*pow(a,3)*(-3.*(Complex(0.,171.) + 694.*omega) + kappa*(Complex(0.,933.) + 4425.*omega) - 15.*(Complex(0.,35.) + 208.*omega)*pow(kappa,2) +
              (Complex(0.,97.) + 830.*omega)*pow(kappa,3) - 2.*(Complex(0.,2.) + 35.*omega)*pow(kappa,4) + omega*pow(kappa,5)))*pow(m,3) +
        2.*pow(a,4)*(69. - 105.*kappa + 3.*(-9. + 5.*kappa)*pow(a,2) + 45.*pow(kappa,2) - 5.*pow(kappa,3))*pow(m,4) + 21940.*pow(omega,2) - 54788.*kappa*pow(omega,2) + 39272.*pow(kappa,2)*pow(omega,2) -
        5416.*pow(kappa,3)*pow(omega,2) - 924.*pow(kappa,4)*pow(omega,2) - 84.*pow(kappa,5)*pow(omega,2) +
        pow(a,2)*pow(m,2)*(414. + Complex(0.,10389.)*omega + 2.*omega*(Complex(0.,34.) + (416. - 181.*kappa)*omega)*pow(a,6) - 3.*omega*(Complex(0.,1.) + 28.*omega)*pow(kappa,5) + 23348.*pow(omega,2) +
           pow(kappa,4)*(16. + Complex(0.,513.)*omega + 3940.*pow(omega,2)) + 6.*pow(kappa,2)*(133. + Complex(0.,2815.)*omega + 9340.*pow(omega,2)) -
           2.*pow(kappa,3)*(111. + Complex(0.,2507.)*omega + 11860.*pow(omega,2)) - 1.*kappa*(1014. + Complex(0.,22743.)*omega + 59524.*pow(omega,2)) +
           pow(a,4)*(133. + Complex(0.,1179.)*omega + 6.*omega*(Complex(0.,10.) + 3.*omega)*pow(kappa,2) - 426.*pow(omega,2) - 110.*pow(kappa,3)*pow(omega,2) +
              kappa*(-89. - Complex(0.,633.)*omega + 1206.*pow(omega,2))) - 1.*pow(a,2)*
            (492. + Complex(0.,9207.)*omega + omega*(Complex(0.,33.) + 566.*omega)*pow(kappa,4) + 16950.*pow(omega,2) - 6.*pow(kappa,3)*(7. + Complex(0.,155.)*omega + 1028.*pow(omega,2)) +
              6.*pow(kappa,2)*(63. + Complex(0.,1006.)*omega + 3430.*pow(omega,2)) - 6.*kappa*(138. + Complex(0.,2251.)*omega + 5092.*pow(omega,2)))) + Complex(0.,7424.)*pow(omega,3) -
        Complex(0.,9728.)*kappa*pow(omega,3) + 2.*(Complex(0.,-131.) + 3.*(-307. + 58.*kappa)*omega)*pow(a,10)*pow(omega,3) - Complex(0.,12288.)*pow(kappa,2)*pow(omega,3) +
        Complex(0.,9728.)*pow(kappa,3)*pow(omega,3) + Complex(0.,4864.)*pow(kappa,4)*pow(omega,3) +
        omega*pow(a,8)*(Complex(0.,48.) + (217. - 41.*kappa)*omega + Complex(0.,3.)*(3309. - 1183.*kappa + 12.*pow(kappa,2))*pow(omega,2) -
           4.*(-9031. + 5042.*kappa - 526.*pow(kappa,2) + 17.*pow(kappa,3))*pow(omega,3)) +
        m*(3.*(Complex(0.,27.) + (200. - 82.*kappa)*omega)*pow(a,9)*pow(omega,2) + 2.*pow(a,7)*
            (Complex(0.,24.) + (-489. + 197.*kappa)*omega - Complex(0.,2.)*(979. - 876.*kappa + 122.*pow(kappa,2))*pow(omega,2) + 4.*(-1953. + 2263.*kappa - 635.*pow(kappa,2) + 29.*pow(kappa,3))*pow(omega,3)) +
           2.*a*(Complex(0.,-849.) - 626.*omega + omega*pow(kappa,5) - Complex(0.,2240.)*pow(omega,2) + kappa*(Complex(0.,1557.) + 2121.*omega + Complex(0.,5440.)*pow(omega,2) - 5632.*pow(omega,3)) +
              2304.*pow(omega,3) + pow(kappa,2)*(Complex(0.,-861.) - 2640.*omega - Complex(0.,7616.)*pow(omega,2) + 2048.*pow(omega,3)) -
              2.*pow(kappa,4)*(Complex(0.,2.) + 107.*omega + Complex(0.,832.)*pow(omega,2) + 2176.*pow(omega,3)) +
              pow(kappa,3)*(Complex(0.,145.) + 1342.*omega + Complex(0.,6080.)*pow(omega,2) + 5632.*pow(omega,3))) +
           pow(a,5)*(Complex(0.,-687.) + 4296.*omega + Complex(0.,16700.)*pow(omega,2) + 2.*(Complex(0.,-15.) + 56.*omega)*pow(kappa,4)*pow(omega,2) -
              2.*omega*pow(kappa,3)*(-24. + Complex(0.,501.)*omega + 4978.*pow(omega,2)) + 61272.*pow(omega,3) - 6.*pow(kappa,5)*pow(omega,3) +
              6.*pow(kappa,2)*(Complex(0.,-12.) + 126.*omega + Complex(0.,2033.)*pow(omega,2) + 10620.*pow(omega,3)) -
              2.*kappa*(Complex(0.,-252.) + 2187.*omega + Complex(0.,14273.)*pow(omega,2) + 56787.*pow(omega,3))) +
           2.*pow(a,3)*(Complex(0.,1015.) - 1662.*omega - Complex(0.,6197.)*pow(omega,2) + 3.*(Complex(0.,1.) + 28.*omega)*pow(kappa,5)*pow(omega,2) +
              omega*pow(kappa,4)*(74. + Complex(0.,351.)*omega + 220.*pow(omega,2)) - 30132.*pow(omega,3) + 2.*pow(kappa,3)*(Complex(0.,-15.) - 224.*omega + Complex(0.,283.)*pow(omega,2) + 7220.*pow(omega,3)) -
              2.*pow(kappa,2)*(Complex(0.,-209.) + 72.*omega + Complex(0.,4925.)*pow(omega,2) + 30036.*pow(omega,3)) +
              kappa*(Complex(0.,-1279.) + 2172.*omega + Complex(0.,15991.)*pow(omega,2) + 74436.*pow(omega,3)))) + 69632.*pow(omega,4) - 126976.*kappa*pow(omega,4) - 40960.*pow(kappa,2)*pow(omega,4) +
        126976.*pow(kappa,3)*pow(omega,4) - 28672.*pow(kappa,4)*pow(omega,4) + pow(a,6)*
         (48. - Complex(0.,6.)*omega*(107. - 96.*kappa + 12.*pow(kappa,2)) - 2.*(-16. + 313.*kappa - 563.*pow(kappa,2) + 57.*pow(kappa,3))*pow(omega,2) +
           Complex(0.,1.)*(-52051. + 33030.*kappa + 3288.*pow(kappa,2) - 1546.*pow(kappa,3) + 7.*pow(kappa,4))*pow(omega,3) +
           2.*(-85359. + 74846.*kappa - 7250.*pow(kappa,2) - 1496.*pow(kappa,3) - 23.*pow(kappa,4) + 2.*pow(kappa,5))*pow(omega,4)) +
        pow(a,2)*(4.*pow(kappa,5)*pow(omega,2) + omega*pow(kappa,4)*(Complex(0.,-229.) + 2098.*omega + Complex(0.,1408.)*pow(omega,2) + 30208.*pow(omega,3)) +
           2.*pow(kappa,2)*(153. - Complex(0.,5508.)*omega - 2802.*pow(omega,2) + Complex(0.,33152.)*pow(omega,3) + 53248.*pow(omega,4)) -
           3.*(-282. + Complex(0.,2861.)*omega + 7082.*pow(omega,2) + Complex(0.,15744.)*pow(omega,3) + 78336.*pow(omega,4)) -
           2.*pow(kappa,3)*(9. - Complex(0.,1487.)*omega + 5240.*pow(omega,2) + Complex(0.,24320.)*pow(omega,3) + 98816.*pow(omega,4)) +
           2.*kappa*(-507. + Complex(0.,8175.)*omega + 18558.*pow(omega,2) + Complex(0.,17664.)*pow(omega,3) + 162304.*pow(omega,4))) -
        1.*pow(a,4)*(491. - Complex(0.,3115.)*omega - 3140.*pow(omega,2) + pow(kappa,4)*pow(omega,2)*(228. + Complex(0.,639.)*omega + 5020.*pow(omega,2)) - Complex(0.,84117.)*pow(omega,3) +
           3.*(Complex(0.,1.) + 28.*omega)*pow(kappa,5)*pow(omega,3) - 2.*omega*pow(kappa,3)*(Complex(0.,-120.) + 2232.*omega + Complex(0.,11157.)*pow(omega,2) + 32620.*pow(omega,3)) - 305844.*pow(omega,4) +
           pow(kappa,2)*(72. - Complex(0.,2056.)*omega + 8932.*pow(omega,2) + Complex(0.,37510.)*pow(omega,3) + 21912.*pow(omega,4)) +
           kappa*(-415. + Complex(0.,4549.)*omega + 692.*pow(omega,2) + Complex(0.,64151.)*pow(omega,3) + 346116.*pow(omega,4)))) +
     8.*pow(lambda,2)*(948. - 2052.*kappa + Complex(0.,1909.)*omega - Complex(0.,5847.)*kappa*omega + 1380.*pow(kappa,2) + Complex(0.,5946.)*omega*pow(kappa,2) - 308.*pow(kappa,3) -
        Complex(0.,2198.)*omega*pow(kappa,3) + 16.*pow(kappa,4) + Complex(0.,225.)*omega*pow(kappa,4) - Complex(0.,3.)*omega*pow(kappa,5) +
        20.*omega*pow(a,3)*(69. - 105.*kappa + 3.*(-9. + 5.*kappa)*pow(a,2) + 45.*pow(kappa,2) - 5.*pow(kappa,3))*pow(m,3) + 45172.*pow(omega,2) - 94276.*kappa*pow(omega,2) + 63912.*pow(kappa,2)*pow(omega,2) -
        16616.*pow(kappa,3)*pow(omega,2) + 1892.*pow(kappa,4)*pow(omega,2) - 84.*pow(kappa,5)*pow(omega,2) +
        2.*pow(a,2)*pow(m,2)*(69. - Complex(0.,4167.)*omega + 30.*omega*(Complex(0.,-9.) + 5.*(-13. + 6.*kappa)*omega)*pow(a,4) - 1.*omega*(Complex(0.,37.) + 700.*omega)*pow(kappa,4) +
           pow(kappa,2)*(45. - Complex(0.,4512.)*omega - 31200.*pow(omega,2)) - 20820.*pow(omega,2) + 10.*pow(kappa,5)*pow(omega,2) + pow(kappa,3)*(-5. + Complex(0.,862.)*omega + 8300.*pow(omega,2)) +
           3.*kappa*(-35. + Complex(0.,2594.)*omega + 14750.*pow(omega,2)) + pow(a,2)*
            (-27. + Complex(0.,3100.)*omega + 10.*omega*(Complex(0.,67.) + 759.*omega)*pow(kappa,2) + 15630.*pow(omega,2) - 670.*pow(kappa,3)*pow(omega,2) -
              5.*kappa*(-3. + Complex(0.,638.)*omega + 4182.*pow(omega,2)))) - Complex(0.,39360.)*pow(omega,3) + Complex(0.,69120.)*kappa*pow(omega,3) +
        4.*(Complex(0.,192.) + (2329. - 706.*kappa)*omega)*pow(a,8)*pow(omega,3) - Complex(0.,32640.)*pow(kappa,2)*pow(omega,3) + Complex(0.,2880.)*pow(kappa,4)*pow(omega,3) +
        2.*omega*pow(a,6)*(Complex(0.,99.) + (-2941. + 1284.*kappa)*omega - Complex(0.,1.)*(13013. - 9893.*kappa + 1394.*pow(kappa,2))*pow(omega,2) +
           4.*(-16461. + 15960.*kappa - 3612.*pow(kappa,2) + 167.*pow(kappa,3))*pow(omega,3)) -
        1.*m*(4.*(Complex(0.,21.) + (-41. + 119.*kappa)*omega)*pow(a,7)*pow(omega,2) +
           pow(a,5)*(Complex(0.,63.) + (-448. + 278.*kappa)*omega - Complex(0.,2.)*(8509. - 5095.*kappa + 376.*pow(kappa,2))*pow(omega,2) +
              8.*(-6231. + 6240.*kappa - 1752.*pow(kappa,2) + 137.*pow(kappa,3))*pow(omega,3)) -
           2.*pow(a,3)*(Complex(0.,373.) + 468.*omega - Complex(0.,43068.)*pow(omega,2) - 2.*(Complex(0.,20.) + 1037.*omega)*pow(kappa,4)*pow(omega,2) +
              4.*omega*pow(kappa,3)*(-15. + Complex(0.,607.)*omega + 7433.*pow(omega,2)) - 122034.*pow(omega,3) + 10.*pow(kappa,5)*pow(omega,3) -
              4.*pow(kappa,2)*(Complex(0.,-19.) - 111.*omega + Complex(0.,5793.)*pow(omega,2) + 32685.*pow(omega,3)) +
              kappa*(Complex(0.,-373.) - 852.*omega + Complex(0.,60252.)*pow(omega,2) + 219858.*pow(omega,3))) +
           2.*a*(Complex(0.,513.) + 1635.*omega - Complex(0.,43253.)*pow(omega,2) + omega*pow(kappa,5)*(-4. + Complex(0.,15.)*omega + 420.*pow(omega,2)) +
              pow(kappa,2)*(Complex(0.,525.) + 3798.*omega - Complex(0.,61770.)*pow(omega,2) - 277512.*pow(omega,3)) +
              pow(kappa,4)*(Complex(0.,4.) + 167.*omega - Complex(0.,1137.)*pow(omega,2) - 15220.*pow(omega,3)) - 123908.*pow(omega,3) +
              pow(kappa,3)*(Complex(0.,-97.) - 1360.*omega + Complex(0.,15166.)*pow(omega,2) + 105160.*pow(omega,3)) +
              kappa*(Complex(0.,-933.) - 4236.*omega + Complex(0.,91011.)*pow(omega,2) + 311060.*pow(omega,3)))) + 49920.*pow(omega,4) - 96768.*kappa*pow(omega,4) - 12288.*pow(kappa,2)*pow(omega,4) +
        96768.*pow(kappa,3)*pow(omega,4) - 37632.*pow(kappa,4)*pow(omega,4) + pow(a,2)*
         (-1194. + Complex(0.,1137.)*omega - 73368.*pow(omega,2) + 10.*pow(kappa,5)*pow(omega,2)*(1. + Complex(0.,3.)*omega + 84.*pow(omega,2)) - Complex(0.,47482.)*pow(omega,3) -
           1.*omega*pow(kappa,4)*(Complex(0.,13.) + 584.*omega + Complex(0.,4722.)*pow(omega,2) + 1640.*pow(omega,3)) - 388360.*pow(omega,4) -
           6.*pow(kappa,2)*(105. + Complex(0.,132.)*omega + 9536.*pow(omega,2) + Complex(0.,21630.)*pow(omega,3) + 88472.*pow(omega,4)) +
           pow(kappa,3)*(54. + Complex(0.,286.)*omega + 9748.*pow(omega,2) + Complex(0.,41276.)*pow(omega,3) + 93200.*pow(omega,4)) +
           2.*kappa*(849. - Complex(0.,177.)*omega + 58689.*pow(omega,2) + Complex(0.,73827.)*pow(omega,3) + 418004.*pow(omega,4))) +
        pow(a,4)*(301. - Complex(0.,1901.)*omega + 38160.*pow(omega,2) - 4.*pow(kappa,3)*pow(omega,2)*(202. + Complex(0.,2401.)*omega + 10856.*pow(omega,2)) + Complex(0.,91386.)*pow(omega,3) +
           2.*(Complex(0.,167.) + 854.*omega)*pow(kappa,4)*pow(omega,3) + 4.*omega*pow(kappa,2)*(Complex(0.,-44.) + 2943.*omega + Complex(0.,16704.)*pow(omega,2) + 78090.*pow(omega,3)) + 420492.*pow(omega,4) -
           40.*pow(kappa,5)*pow(omega,4) - 1.*kappa*(161. - Complex(0.,1403.)*omega + 41484.*pow(omega,2) + Complex(0.,147540.)*pow(omega,3) + 677976.*pow(omega,4)))))*pow(omega,6);
}

Complex gsn_asymptotic_derivative_initial_sum_term_11(const double &a, const int &, const int &mTemp, const double &omega, const double &lambda, const double &){
	double kappa = sqrt(1. - pow(a, 2));
	Complex m = Complex(mTemp);
	return 0.0026041666666666665*(Complex(0.,-1.)*(5.*(Complex(0.,9.) + 2.*(-27. + 8.*kappa)*omega)*pow(a,6) -
        10.*pow(a,4)*(Complex(0.,90.) - 366.*omega + kappa*(Complex(0.,-99.) + 336.*omega) + (Complex(0.,21.) - 78.*omega)*pow(kappa,2) + 4.*omega*pow(kappa,3)) +
        2.*a*m*(-501. + 900.*kappa + 15.*(-5. + 2.*kappa)*pow(a,4) - 510.*pow(kappa,2) + 100.*pow(kappa,3) - 10.*pow(a,2)*(-45. + 51.*kappa - 15.*pow(kappa,2) + pow(kappa,3)) - 5.*pow(kappa,4)) +
        pow(a,2)*(Complex(0.,2631.) - 10818.*omega + 12.*kappa*(Complex(0.,-463.) + 1370.*omega) + (Complex(0.,3426.) - 7740.*omega)*pow(kappa,2) + 4.*(Complex(0.,-169.) + 310.*omega)*pow(kappa,3) +
           (Complex(0.,31.) - 50.*omega)*pow(kappa,4)) + Complex(0.,2.)*(-951. + kappa*(2913. + Complex(0.,8848.)*omega) - Complex(0.,4236.)*omega + (-2910. - Complex(0.,6216.)*omega)*pow(kappa,2) +
           10.*(113. + Complex(0.,168.)*omega)*pow(kappa,3) + (-155. - Complex(0.,140.)*omega)*pow(kappa,4) + 5.*pow(kappa,5)))*pow(lambda,5) +
     2.*pow(lambda,4)*(a*m*(-4.*(Complex(0.,315.) + 3377.*omega) + kappa*(Complex(0.,3009.) + 31080.*omega) + 95.*omega*pow(a,6) - 12.*(Complex(0.,197.) + 2050.*omega)*pow(kappa,2) -
           3.*pow(a,4)*(Complex(0.,100.) - Complex(0.,55.)*kappa + 1078.*omega - 900.*kappa*omega + 150.*omega*pow(kappa,2)) + (Complex(0.,694.) + 7840.*omega)*pow(kappa,3) -
           4.*(Complex(0.,16.) + 225.*omega)*pow(kappa,4) + pow(a,2)*(87.*(Complex(0.,16.) + 159.*omega) - 36.*kappa*(Complex(0.,59.) + 615.*omega) + 6.*(Complex(0.,146.) + 1805.*omega)*pow(kappa,2) -
              4.*(Complex(0.,22.) + 435.*omega)*pow(kappa,3) + 65.*omega*pow(kappa,4)) + (Complex(0.,1.) + 24.*omega)*pow(kappa,5)) +
        pow(a,2)*(501. - 900.*kappa + 21.*pow(a,4) + 510.*pow(kappa,2) - 12.*pow(a,2)*(26. - 25.*kappa + 5.*pow(kappa,2)) - 100.*pow(kappa,3) + 5.*pow(kappa,4))*pow(m,2) + 80.*pow(a,8)*pow(omega,2) -
        1.*pow(a,6)*(3. - Complex(0.,1.)*(-671. + 230.*kappa)*omega + 4.*(1311. - 875.*kappa + 110.*pow(kappa,2))*pow(omega,2)) +
        2.*pow(a,4)*(-121. + Complex(0.,3735.)*omega - 2.*omega*(Complex(0.,27.) + 740.*omega)*pow(kappa,3) + kappa*(259. - Complex(0.,4062.)*omega - 29880.*pow(omega,2)) + 22128.*pow(omega,2) +
           40.*pow(kappa,4)*pow(omega,2) + pow(kappa,2)*(-73. + Complex(0.,1053.)*omega + 11880.*pow(omega,2))) +
        pow(a,2)*(1305. - Complex(0.,18753.)*omega + 2.*omega*(Complex(0.,1.) + 30.*omega)*pow(kappa,5) + pow(kappa,4)*(41. - Complex(0.,161.)*omega - 3300.*pow(omega,2)) - 107060.*pow(omega,2) +
           4.*pow(kappa,3)*(-197. + Complex(0.,814.)*omega + 9310.*pow(omega,2)) - 6.*pow(kappa,2)*(-565. + Complex(0.,3041.)*omega + 23900.*pow(omega,2)) +
           6.*kappa*(-722. + Complex(0.,5621.)*omega + 35714.*pow(omega,2))) - 2.*(642. - Complex(0.,6400.)*omega + pow(kappa,4)*(246. - Complex(0.,440.)*omega - 5600.*pow(omega,2)) - 37728.*pow(omega,2) +
           pow(kappa,5)*(-9. + Complex(0.,18.)*omega + 224.*pow(omega,2)) - 8.*pow(kappa,2)*(-432. + Complex(0.,1599.)*omega + 11480.*pow(omega,2)) +
           2.*pow(kappa,3)*(-783. + Complex(0.,1990.)*omega + 18144.*pow(omega,2)) + kappa*(-2817. + Complex(0.,15730.)*omega + 98656.*pow(omega,2)))) -
     Complex(0.,4.)*pow(lambda,3)*(2.*(3.*(Complex(0.,7.) + 50.*(-5. + 2.*kappa)*omega)*pow(a,6) -
           4.*pow(a,4)*(Complex(0.,78.) - 1125.*omega + 75.*kappa*(Complex(0.,-1.) + 17.*omega) + (Complex(0.,15.) - 375.*omega)*pow(kappa,2) + 25.*omega*pow(kappa,3)) -
           1.*(Complex(0.,-1.) + 10.*omega)*pow(a,2)*(501. - 900.*kappa + 510.*pow(kappa,2) - 100.*pow(kappa,3) + 5.*pow(kappa,4)))*pow(m,2) -
        10.*(Complex(0.,23.) + 2.*(-219. + 52.*kappa)*omega)*pow(a,8)*pow(omega,2) -
        2.*a*m*(-2763. + Complex(0.,25514.)*omega + 20.*omega*(Complex(0.,-19.) + kappa*omega)*pow(a,6) + (1. - Complex(0.,74.)*omega)*pow(kappa,5) +
           pow(kappa,4)*(-79. + Complex(0.,2450.)*omega - 2360.*pow(omega,2)) - 37368.*pow(omega,2) + 2.*pow(kappa,3)*(497. - Complex(0.,9810.)*omega + 10800.*pow(omega,2)) -
           6.*pow(kappa,2)*(649. - Complex(0.,9430.)*omega + 11032.*pow(omega,2)) + kappa*(5709. - Complex(0.,65106.)*omega + 83680.*pow(omega,2)) +
           pow(a,4)*(-525. + Complex(0.,8922.)*omega - 24.*omega*(Complex(0.,-67.) + 115.*omega)*pow(kappa,2) - 9120.*pow(omega,2) + 220.*pow(kappa,3)*pow(omega,2) +
              3.*kappa*(85. - Complex(0.,2838.)*omega + 3060.*pow(omega,2))) + 2.*pow(a,2)*
            (110.*omega*(Complex(0.,-1.) + 2.*omega)*pow(kappa,4) + pow(kappa,3)*(-59. + Complex(0.,2632.)*omega - 3800.*pow(omega,2)) + 3.*pow(kappa,2)*(221. - Complex(0.,4948.)*omega + 6080.*pow(omega,2)) +
              3.*(457. - Complex(0.,5074.)*omega + 6668.*pow(omega,2)) - 3.*kappa*(609. - Complex(0.,9112.)*omega + 11000.*pow(omega,2)))) +
        pow(a,6)*(Complex(0.,-141.) + (904. - 268.*kappa)*omega + Complex(0.,4.)*(-1416. - 353.*kappa + 281.*pow(kappa,2))*pow(omega,2) +
           40.*(-1629. + 1158.*kappa - 177.*pow(kappa,2) + 4.*pow(kappa,3))*pow(omega,3)) +
        2.*(-1. + kappa)*(Complex(0.,-1569.) + 20516.*omega + pow(kappa,4)*(Complex(0.,3.) + 36.*omega - Complex(0.,448.)*pow(omega,2)) - Complex(0.,91584.)*pow(omega,2) - 38400.*pow(omega,3) +
           16.*kappa*(Complex(0.,96.) - 1731.*omega + Complex(0.,8468.)*pow(omega,2) + 288.*pow(omega,3)) - 8.*pow(kappa,3)*(Complex(0.,3.) + 158.*omega - Complex(0.,1464.)*pow(omega,2) + 1344.*pow(omega,3)) +
           6.*pow(kappa,2)*(Complex(0.,-47.) + 1796.*omega - Complex(0.,10944.)*pow(omega,2) + 5376.*pow(omega,3))) +
        pow(a,2)*(Complex(0.,-5283.) + 52104.*omega + 4.*(-1. + Complex(0.,5.)*omega)*omega*pow(kappa,5) - Complex(0.,251164.)*pow(omega,2) +
           pow(kappa,3)*(Complex(0.,452.) - 9080.*omega + Complex(0.,58600.)*pow(omega,2) - 37440.*pow(omega,3)) - 250992.*pow(omega,3) +
           pow(kappa,4)*(Complex(0.,-11.) + 472.*omega - Complex(0.,4460.)*pow(omega,2) + 8720.*pow(omega,3)) -
           6.*pow(kappa,2)*(Complex(0.,583.) - 8128.*omega + Complex(0.,43380.)*pow(omega,2) + 11312.*pow(omega,3)) +
           4.*kappa*(Complex(0.,2001.) - 22569.*omega + Complex(0.,111729.)*pow(omega,2) + 83504.*pow(omega,3))) -
        2.*pow(a,4)*(Complex(0.,-1108.) + 8316.*omega - Complex(0.,45501.)*pow(omega,2) + (Complex(0.,51.) + 230.*omega)*pow(kappa,4)*pow(omega,2) -
           4.*omega*pow(kappa,3)*(33. + Complex(0.,3.)*omega + 10.*pow(omega,2)) - 112266.*pow(omega,3) - 1.*pow(kappa,2)*(Complex(0.,169.) - 2304.*omega + Complex(0.,9582.)*pow(omega,2) + 31500.*pow(omega,3)) +
           kappa*(Complex(0.,967.) - 8736.*omega + Complex(0.,44436.)*pow(omega,2) + 123960.*pow(omega,3)))) -
     Complex(0.,48.)*lambda*omega*(2.*(3.*(Complex(0.,7.) + 15.*(-5. + 2.*kappa)*omega)*pow(a,7) -
           6.*pow(a,5)*(Complex(0.,52.) - 225.*omega + 5.*kappa*(Complex(0.,-10.) + 51.*omega) + (Complex(0.,10.) - 75.*omega)*pow(kappa,2) + 5.*omega*pow(kappa,3)) -
           1.*(Complex(0.,-1.) + 3.*omega)*pow(a,3)*(501. - 900.*kappa + 510.*pow(kappa,2) - 100.*pow(kappa,3) + 5.*pow(kappa,4)))*pow(m,3) +
        pow(m,2)*(omega*(Complex(0.,647.) - 170.*(-17. + 4.*kappa)*omega)*pow(a,8) -
           2.*pow(a,6)*(-600. + Complex(0.,7668.)*omega - 3.*omega*(Complex(0.,-359.) + 370.*omega)*pow(kappa,2) + 8130.*pow(omega,2) + 280.*pow(kappa,3)*pow(omega,2) -
              3.*kappa*(-95. + Complex(0.,2167.)*omega + 1300.*pow(omega,2))) + pow(a,4)*
            (-6384. + Complex(0.,52689.)*omega + (Complex(0.,233.) - 1210.*omega)*omega*pow(kappa,4) + 15894.*pow(omega,2) + 12.*kappa*(694. - Complex(0.,7077.)*omega + 410.*pow(omega,2)) +
              4.*pow(kappa,3)*(64. - Complex(0.,1595.)*omega + 3430.*pow(omega,2)) - 6.*pow(kappa,2)*(492. - Complex(0.,6837.)*omega + 5090.*pow(omega,2))) +
           2.*pow(a,2)*(3264. - Complex(0.,22185.)*omega + (-1. + Complex(0.,39.)*omega)*pow(kappa,5) + kappa*(-6609. + Complex(0.,51243.)*omega - 15344.*pow(omega,2)) + 1060.*pow(omega,2) +
              7.*pow(kappa,4)*(12. - Complex(0.,195.)*omega + 380.*pow(omega,2)) + 6.*pow(kappa,2)*(734. - Complex(0.,6675.)*omega + 4676.*pow(omega,2)) -
              2.*pow(kappa,3)*(547. - Complex(0.,6175.)*omega + 8120.*pow(omega,2)))) + (Complex(0.,-1.) + 10.*(-273. + 32.*kappa)*omega)*pow(a,10)*pow(omega,3) -
        2.*omega*pow(a,8)*(Complex(0.,61.) + (-2103. + 668.*kappa)*omega + Complex(0.,1.)*(152. - 3073.*kappa + 709.*pow(kappa,2))*pow(omega,2) +
           10.*(-3005. + 924.*kappa + 3.*pow(kappa,2) + 2.*pow(kappa,3))*pow(omega,3)) -
        2.*a*m*(Complex(0.,-1688.) + 22479.*omega - Complex(0.,117536.)*pow(omega,2) + (Complex(0.,-122.) + 15.*(-103. + 26.*kappa)*omega)*pow(a,8)*pow(omega,2) +
           Complex(0.,2.)*pow(kappa,5)*(1. + Complex(0.,18.)*omega + 224.*pow(omega,2)) - 74624.*pow(omega,3) +
           6.*pow(kappa,2)*(Complex(0.,-628.) + 7527.*omega - Complex(0.,37120.)*pow(omega,2) + 7680.*pow(omega,3)) -
           4.*pow(kappa,3)*(Complex(0.,-307.) + 3905.*omega - Complex(0.,19152.)*pow(omega,2) + 23360.*pow(omega,3)) +
           pow(kappa,4)*(Complex(0.,-128.) + 1935.*omega - Complex(0.,10720.)*pow(omega,2) + 28544.*pow(omega,3)) +
           kappa*(Complex(0.,4386.) - 53808.*omega + Complex(0.,273920.)*pow(omega,2) + 93440.*pow(omega,3)) +
           pow(a,6)*(Complex(0.,73.) + (-307. + 159.*kappa)*omega + Complex(0.,2.)*(5341. - 2251.*kappa + 42.*pow(kappa,2))*pow(omega,2) +
              10.*(3107. - 2111.*kappa + 237.*pow(kappa,2) + 15.*pow(kappa,3))*pow(omega,3)) +
           pow(a,4)*(Complex(0.,-1162.) + 8235.*omega - Complex(0.,85248.)*pow(omega,2) + (Complex(0.,8.) + 1305.*omega)*pow(kappa,4)*pow(omega,2) -
              4.*omega*pow(kappa,3)*(79. - Complex(0.,434.)*omega + 2225.*pow(omega,2)) - 130599.*pow(omega,3) -
              2.*pow(kappa,2)*(Complex(0.,119.) - 1725.*omega + Complex(0.,13188.)*pow(omega,2) + 9165.*pow(omega,3)) +
              2.*kappa*(Complex(0.,587.) - 4941.*omega + Complex(0.,46236.)*pow(omega,2) + 66270.*pow(omega,3))) +
           pow(a,2)*(Complex(0.,2805.) - 27723.*omega + 3.*(1. - Complex(0.,18.)*omega)*omega*pow(kappa,5) + Complex(0.,185758.)*pow(omega,2) +
              pow(kappa,4)*(Complex(0.,21.) - 489.*omega + Complex(0.,2646.)*pow(omega,2) - 14840.*pow(omega,3)) + 179208.*pow(omega,3) +
              6.*pow(kappa,2)*(Complex(0.,481.) - 5110.*omega + Complex(0.,29730.)*pow(omega,2) + 392.*pow(omega,3)) +
              4.*pow(kappa,3)*(Complex(0.,-129.) + 1763.*omega - Complex(0.,9327.)*pow(omega,2) + 15880.*pow(omega,3)) -
              3.*kappa*(Complex(0.,1732.) - 16851.*omega + Complex(0.,106906.)*pow(omega,2) + 73376.*pow(omega,3)))) +
        2.*pow(a,4)*(1548. - Complex(0.,6363.)*omega + 79999.*pow(omega,2) + (-2. + Complex(0.,75.)*omega)*pow(kappa,5)*pow(omega,2) - Complex(0.,12009.)*pow(omega,3) +
           omega*pow(kappa,4)*(Complex(0.,-57.) + 1013.*omega - Complex(0.,1413.)*pow(omega,2) + 4580.*pow(omega,3)) +
           kappa*(-1476. + Complex(0.,11565.)*omega - 114070.*pow(omega,2) + Complex(0.,140559.)*pow(omega,3) - 152176.*pow(omega,4)) +
           pow(kappa,2)*(324. - Complex(0.,6765.)*omega + 56728.*pow(omega,2) - Complex(0.,133410.)*pow(omega,3) - 122472.*pow(omega,4)) + 321444.*pow(omega,4) +
           2.*pow(kappa,3)*(-6. + Complex(0.,672.)*omega - 6342.*pow(omega,2) + Complex(0.,15307.)*pow(omega,3) + 17160.*pow(omega,4))) +
        pow(a,6)*(-416. + Complex(0.,2987.)*omega - 50862.*pow(omega,2) - 4.*pow(kappa,3)*pow(omega,2)*(-268. + Complex(0.,1599.)*omega + 1150.*pow(omega,2)) + Complex(0.,7885.)*pow(omega,3) +
           15.*(Complex(0.,11.) - 42.*omega)*pow(kappa,4)*pow(omega,3) + 2.*omega*pow(kappa,2)*(Complex(0.,382.) - 6334.*omega + Complex(0.,24723.)*pow(omega,2) + 16390.*pow(omega,3)) - 317478.*pow(omega,4) +
           4.*kappa*(36. - Complex(0.,803.)*omega + 11656.*pow(omega,2) - Complex(0.,20515.)*pow(omega,3) + 34090.*pow(omega,4))) +
        2.*(1664. - Complex(0.,4581.)*omega + Complex(0.,15.)*omega*pow(kappa,5) + 33380.*pow(omega,2) + Complex(0.,17920.)*pow(omega,3) +
           omega*pow(kappa,4)*(Complex(0.,-1425.) + 5540.*omega + Complex(0.,8192.)*pow(omega,2) + 22528.*pow(omega,3)) + 90112.*pow(omega,4) +
           2.*pow(kappa,3)*(-160. + Complex(0.,4095.)*omega - 12024.*pow(omega,2) + Complex(0.,512.)*pow(omega,3) + 22528.*pow(omega,4)) -
           1.*kappa*(3264. - Complex(0.,14691.)*omega + 70256.*pow(omega,2) + Complex(0.,1024.)*pow(omega,3) + 45056.*pow(omega,4)) -
           2.*pow(kappa,2)*(-960. + Complex(0.,8397.)*omega - 27788.*pow(omega,2) + Complex(0.,13056.)*pow(omega,3) + 56320.*pow(omega,4))) +
        pow(a,2)*(-5952. + Complex(0.,18757.)*omega - 181990.*pow(omega,2) - Complex(0.,4.)*omega*pow(kappa,5)*(1. + Complex(0.,18.)*omega + 224.*pow(omega,2)) +
           omega*pow(kappa,4)*(Complex(0.,1501.) - 11862.*omega + Complex(0.,2240.)*pow(omega,2) - 36096.*pow(omega,3)) - Complex(0.,6848.)*pow(omega,3) - 559872.*pow(omega,4) +
           16.*kappa*(552. - Complex(0.,2922.)*omega + 20835.*pow(omega,2) - Complex(0.,16448.)*pow(omega,3) + 16416.*pow(omega,4)) -
           4.*pow(kappa,3)*(-96. + Complex(0.,3473.)*omega - 19294.*pow(omega,2) + Complex(0.,28320.)*pow(omega,3) + 43136.*pow(omega,4)) +
           pow(kappa,2)*(-3648. + Complex(0.,40278.)*omega - 226068.*pow(omega,2) + Complex(0.,364544.)*pow(omega,3) + 460800.*pow(omega,4)))) +
     288.*pow(omega,2)*(pow(a,3)*(-4.*(Complex(0.,816.) + 3377.*omega) + kappa*(Complex(0.,6609.) + 31080.*omega) + 53.*omega*pow(a,6) - 12.*(Complex(0.,367.) + 2050.*omega)*pow(kappa,2) -
           15.*pow(a,4)*(Complex(0.,40.) + 174.*omega - 1.*kappa*(Complex(0.,19.) + 140.*omega) + 22.*omega*pow(kappa,2)) + 2.*(Complex(0.,547.) + 3920.*omega)*pow(kappa,3) -
           12.*(Complex(0.,7.) + 75.*omega)*pow(kappa,4) + pow(a,2)*(21.*(Complex(0.,152.) + 611.*omega) - 12.*kappa*(Complex(0.,347.) + 1695.*omega) + 18.*(Complex(0.,82.) + 545.*omega)*pow(kappa,2) -
              4.*(Complex(0.,32.) + 385.*omega)*pow(kappa,3) + 55.*omega*pow(kappa,4)) + (Complex(0.,1.) + 24.*omega)*pow(kappa,5))*pow(m,3) +
        pow(a,4)*(501. - 900.*kappa + 21.*pow(a,4) + 510.*pow(kappa,2) - 12.*pow(a,2)*(26. - 25.*kappa + 5.*pow(kappa,2)) - 100.*pow(kappa,3) + 5.*pow(kappa,4))*pow(m,4) -
        1.*pow(a,10)*pow(omega,2)*(69. - Complex(0.,1.)*(-3539. + 634.*kappa)*omega + 4.*(3868. - 1415.*kappa + 93.*pow(kappa,2))*pow(omega,2)) +
        pow(a,2)*pow(m,2)*(1357. + Complex(0.,30512.)*omega + 69504.*pow(omega,2) - 249.*pow(a,8)*pow(omega,2) - 2.*pow(kappa,5)*(1. + Complex(0.,18.)*omega + 224.*pow(omega,2)) -
           8.*pow(kappa,3)*(146. + Complex(0.,2435.)*omega + 10496.*pow(omega,2)) + pow(kappa,4)*(133. + Complex(0.,2400.)*omega + 15360.*pow(omega,2)) +
           6.*pow(kappa,2)*(553. + Complex(0.,9656.)*omega + 30912.*pow(omega,2)) - 2.*kappa*(1827. + Complex(0.,35634.)*omega + 92960.*pow(omega,2)) +
           pow(a,6)*(-53. - Complex(0.,1.)*(-205. + 204.*kappa)*omega + 8.*(549. - 430.*kappa + 45.*pow(kappa,2))*pow(omega,2)) +
           pow(a,4)*(787. + Complex(0.,7758.)*omega - 20.*omega*(Complex(0.,5.) + 71.*omega)*pow(kappa,3) + 4491.*pow(omega,2) + 75.*pow(kappa,4)*pow(omega,2) +
              2.*pow(kappa,2)*(101. + Complex(0.,825.)*omega + 1785.*pow(omega,2)) - 2.*kappa*(443. + Complex(0.,3510.)*omega + 1830.*pow(omega,2))) +
           pow(a,2)*(12.*pow(kappa,5)*pow(omega,2) - 1.*pow(kappa,4)*(21. + Complex(0.,381.)*omega + 3740.*pow(omega,2)) + 4.*pow(kappa,3)*(117. + Complex(0.,1583.)*omega + 8110.*pow(omega,2)) +
              12.*kappa*(326. + Complex(0.,4777.)*omega + 10821.*pow(omega,2)) - 6.*pow(kappa,2)*(395. + Complex(0.,5193.)*omega + 16180.*pow(omega,2)) - 3.*(655. + Complex(0.,11351.)*omega + 21604.*pow(omega,2))
              )) + 2.*omega*pow(a,8)*(Complex(0.,176.) + 337.*omega + Complex(0.,25729.)*pow(omega,2) + 2.*(Complex(0.,33.) - 200.*omega)*pow(kappa,3)*pow(omega,2) +
           omega*pow(kappa,2)*(-19. + Complex(0.,371.)*omega + 8520.*pow(omega,2)) + 82948.*pow(omega,3) + 20.*pow(kappa,4)*pow(omega,3) -
           1.*kappa*(Complex(0.,72.) + 205.*omega + Complex(0.,13030.)*pow(omega,2) + 58640.*pow(omega,3))) +
        a*m*(Complex(0.,-4928.) - 6212.*omega + (Complex(0.,1.) + 24.*omega)*pow(kappa,5) - Complex(0.,6976.)*pow(omega,2) +
           omega*pow(a,8)*(290. - Complex(0.,1.)*(-1674. + 635.*kappa)*omega + (6954. - 4780.*kappa + 570.*pow(kappa,2))*pow(omega,2)) +
           kappa*(Complex(0.,9873.) + 18984.*omega + Complex(0.,19584.)*pow(omega,2) - 35840.*pow(omega,3)) + 15872.*pow(omega,3) - 5.*pow(a,10)*pow(omega,3) +
           4.*pow(kappa,2)*(Complex(0.,-1581.) - 5430.*omega - Complex(0.,9536.)*pow(omega,2) + 2304.*pow(omega,3)) -
           4.*pow(kappa,4)*(Complex(0.,21.) + 465.*omega + Complex(0.,2800.)*pow(omega,2) + 6272.*pow(omega,3)) +
           2.*pow(kappa,3)*(Complex(0.,707.) + 5360.*omega + Complex(0.,18368.)*pow(omega,2) + 17920.*pow(omega,3)) -
           1.*pow(a,6)*(Complex(0.,-376.) + 5109.*omega + Complex(0.,24228.)*pow(omega,2) - 4.*(Complex(0.,44.) + 925.*omega)*pow(kappa,3)*pow(omega,2) +
              2.*omega*pow(kappa,2)*(202. + Complex(0.,3344.)*omega + 20505.*pow(omega,2)) + 80999.*pow(omega,3) + 15.*pow(kappa,4)*pow(omega,3) -
              4.*kappa*(Complex(0.,-36.) + 905.*omega + Complex(0.,6851.)*pow(omega,2) + 27425.*pow(omega,3))) +
           pow(a,2)*(Complex(0.,6840.) - 5095.*omega + omega*pow(kappa,4)*(673. + Complex(0.,3328.)*omega - 640.*pow(omega,2)) - Complex(0.,47136.)*pow(omega,2) +
              4.*omega*pow(kappa,5)*(1. + Complex(0.,18.)*omega + 224.*pow(omega,2)) - 188544.*pow(omega,3) +
              4.*pow(kappa,3)*(Complex(0.,-104.) - 949.*omega + Complex(0.,764.)*pow(omega,2) + 25664.*pow(omega,3)) -
              6.*pow(kappa,2)*(Complex(0.,-662.) - 469.*omega + Complex(0.,11952.)*pow(omega,2) + 65408.*pow(omega,3)) +
              4.*kappa*(Complex(0.,-2457.) + 1386.*omega + Complex(0.,29154.)*pow(omega,2) + 118240.*pow(omega,3))) +
           pow(a,4)*(Complex(0.,-2904.) + 13512.*omega + Complex(0.,69618.)*pow(omega,2) - 3.*(Complex(0.,1.) + 32.*omega)*pow(kappa,5)*pow(omega,2) +
              10.*omega*pow(kappa,4)*(-3. - Complex(0.,21.)*omega + 170.*pow(omega,2)) + pow(kappa,3)*(Complex(0.,24.) + 96.*omega - Complex(0.,6682.)*pow(omega,2) - 54800.*pow(omega,3)) + 228932.*pow(omega,3) +
              6.*pow(kappa,2)*(Complex(0.,-96.) + 755.*omega + Complex(0.,10268.)*pow(omega,2) + 46900.*pow(omega,3)) -
              3.*kappa*(Complex(0.,-887.) + 5596.*omega + Complex(0.,42281.)*pow(omega,2) + 150608.*pow(omega,3)))) + 188.*pow(a,12)*pow(omega,4) -
        2.*(pow(kappa,5)*(1. + Complex(0.,18.)*omega + 224.*pow(omega,2)) + 16.*pow(kappa,4)*(-4. - Complex(0.,105.)*omega - 10.*pow(omega,2) - Complex(0.,448.)*pow(omega,3) + 2304.*pow(omega,4)) +
           4.*pow(kappa,2)*(-111. - Complex(0.,6882.)*omega - 18032.*pow(omega,2) + Complex(0.,3456.)*pow(omega,3) + 11264.*pow(omega,4)) -
           4.*(-101. + Complex(0.,2902.)*omega + 8680.*pow(omega,2) + Complex(0.,1664.)*pow(omega,3) + 20480.*pow(omega,4)) -
           2.*pow(kappa,3)*(-187. - Complex(0.,5590.)*omega - 8352.*pow(omega,2) + Complex(0.,4096.)*pow(omega,3) + 75776.*pow(omega,4)) +
           kappa*(-255. + Complex(0.,29586.)*omega + 90080.*pow(omega,2) + Complex(0.,8192.)*pow(omega,3) + 151552.*pow(omega,4))) +
        pow(a,6)*(323. - Complex(0.,2795.)*omega - 1241.*pow(omega,2) + pow(kappa,4)*pow(omega,2)*(3. + Complex(0.,115.)*omega + 60.*pow(omega,2)) +
           2.*omega*pow(kappa,2)*(Complex(0.,-432.) + 2665.*omega + Complex(0.,6417.)*pow(omega,2) - 35060.*pow(omega,3)) - Complex(0.,187165.)*pow(omega,3) +
           2.*(Complex(0.,1.) + 30.*omega)*pow(kappa,5)*pow(omega,3) - 4.*omega*pow(kappa,3)*(Complex(0.,-6.) + 335.*omega + Complex(0.,2626.)*pow(omega,2) + 4530.*pow(omega,3)) - 584980.*pow(omega,4) +
           2.*kappa*(-72. + Complex(0.,1599.)*omega - 340.*pow(omega,2) + Complex(0.,69455.)*pow(omega,3) + 288998.*pow(omega,4))) +
        pow(a,2)*(1995. - Complex(0.,31685.)*omega + 2.*omega*(Complex(0.,1.) + 30.*omega)*pow(kappa,5) - 88052.*pow(omega,2) - Complex(0.,116416.)*pow(omega,3) - 601600.*pow(omega,4) +
           pow(kappa,4)*(-21. - Complex(0.,1413.)*omega + 6748.*pow(omega,2) + Complex(0.,3520.)*pow(omega,3) + 96768.*pow(omega,4)) +
           2.*pow(kappa,2)*(93. - Complex(0.,23583.)*omega - 35028.*pow(omega,2) + Complex(0.,85120.)*pow(omega,3) + 133632.*pow(omega,4)) -
           4.*pow(kappa,3)*(-45. - Complex(0.,3578.)*omega + 4290.*pow(omega,2) + Complex(0.,32736.)*pow(omega,3) + 139008.*pow(omega,4)) +
           2.*kappa*(-1026. + Complex(0.,32343.)*omega + 86214.*pow(omega,2) + Complex(0.,45504.)*pow(omega,3) + 429568.*pow(omega,4))) -
        2.*pow(a,4)*(769. - Complex(0.,6783.)*omega - 13380.*pow(omega,2) + pow(kappa,5)*pow(omega,2)*(1. + Complex(0.,18.)*omega + 224.*pow(omega,2)) - Complex(0.,122824.)*pow(omega,3) +
           4.*omega*pow(kappa,4)*(Complex(0.,-9.) + 216.*omega + Complex(0.,416.)*pow(omega,2) + 2960.*pow(omega,3)) - 439808.*pow(omega,4) +
           pow(kappa,2)*(187. - Complex(0.,5673.)*omega + 8156.*pow(omega,2) + Complex(0.,57624.)*pow(omega,3) + 25536.*pow(omega,4)) -
           2.*pow(kappa,3)*(6. - Complex(0.,515.)*omega + 4233.*pow(omega,2) + Complex(0.,20666.)*pow(omega,3) + 58944.*pow(omega,4)) +
           kappa*(-745. + Complex(0.,10710.)*omega + 15529.*pow(omega,2) + Complex(0.,102898.)*pow(omega,3) + 530848.*pow(omega,4)))) +
     8.*pow(lambda,2)*(10.*omega*pow(a,3)*(501. - 900.*kappa + 21.*pow(a,4) + 510.*pow(kappa,2) - 12.*pow(a,2)*(26. - 25.*kappa + 5.*pow(kappa,2)) - 100.*pow(kappa,3) + 5.*pow(kappa,4))*pow(m,3) +
        2.*pow(a,8)*pow(omega,2)*(766. - Complex(0.,1.)*(-5959. + 1956.*kappa)*omega + 10.*(3765. - 2170.*kappa + 217.*pow(kappa,2))*pow(omega,2)) +
        pow(a,2)*pow(m,2)*(501. - Complex(0.,26628.)*omega + 10.*omega*(Complex(0.,1.) + 24.*omega)*pow(kappa,5) + pow(kappa,4)*(5. - Complex(0.,780.)*omega - 9000.*pow(omega,2)) - 135080.*pow(omega,2) +
           740.*pow(a,6)*pow(omega,2) + 20.*pow(kappa,3)*(-5. + Complex(0.,487.)*omega + 3920.*pow(omega,2)) - 30.*pow(kappa,2)*(-17. + Complex(0.,1264.)*omega + 8200.*pow(omega,2)) +
           30.*kappa*(-30. + Complex(0.,1843.)*omega + 10360.*pow(omega,2)) - 3.*pow(a,4)*(-7. - Complex(0.,10.)*(-170. + 83.*kappa)*omega + 20.*(487. - 400.*kappa + 65.*pow(kappa,2))*pow(omega,2)) +
           4.*pow(a,2)*(-78. + Complex(0.,6630.)*omega - 10.*omega*(Complex(0.,29.) + 410.*omega)*pow(kappa,3) + 33330.*pow(omega,2) + 150.*pow(kappa,4)*pow(omega,2) +
              15.*pow(kappa,2)*(-1. + Complex(0.,216.)*omega + 1720.*pow(omega,2)) - 15.*kappa*(-5. + Complex(0.,592.)*omega + 3540.*pow(omega,2)))) +
        a*m*(-1.*omega*pow(a,6)*(325. - Complex(0.,2.)*(-2071. + 189.*kappa)*omega + 20.*(483. - 124.*kappa + pow(kappa,2))*pow(omega,2)) +
           kappa*(Complex(0.,6609.) + 25440.*omega - Complex(0.,585192.)*pow(omega,2) - 1.938944e6*pow(omega,3)) + pow(kappa,5)*(Complex(0.,1.) + 64.*omega - Complex(0.,360.)*pow(omega,2) - 4480.*pow(omega,3)) -
           630.*pow(a,8)*pow(omega,3) + 4.*pow(kappa,4)*(Complex(0.,-21.) - 370.*omega + Complex(0.,3420.)*pow(omega,2) + 31120.*pow(omega,3)) +
           8.*(Complex(0.,-408.) - 1213.*omega + Complex(0.,32362.)*pow(omega,2) + 92088.*pow(omega,3)) -
           2.*pow(kappa,3)*(Complex(0.,-547.) - 4720.*omega + Complex(0.,65720.)*pow(omega,2) + 379968.*pow(omega,3)) +
           4.*pow(kappa,2)*(Complex(0.,-1101.) - 5940.*omega + Complex(0.,111072.)*pow(omega,2) + 460544.*pow(omega,3)) +
           pow(a,2)*(Complex(0.,3192.) + 3873.*omega + omega*pow(kappa,4)*(41. - Complex(0.,1794.)*omega - 30720.*pow(omega,2)) - Complex(0.,319050.)*pow(omega,2) +
              10.*(Complex(0.,1.) + 36.*omega)*pow(kappa,5)*pow(omega,2) - 889824.*pow(omega,3) + 4.*pow(kappa,3)*(Complex(0.,-32.) - 183.*omega + Complex(0.,10333.)*pow(omega,2) + 81660.*pow(omega,3)) -
              6.*pow(kappa,2)*(Complex(0.,-246.) - 613.*omega + Complex(0.,42270.)*pow(omega,2) + 199440.*pow(omega,3)) +
              6.*kappa*(Complex(0.,-694.) - 1114.*omega + Complex(0.,86123.)*pow(omega,2) + 294204.*pow(omega,3))) +
           pow(a,4)*(Complex(0.,-600.) + 2194.*omega + Complex(0.,96036.)*pow(omega,2) - 8.*(Complex(0.,63.) + 2215.*omega)*pow(kappa,3)*pow(omega,2) +
              2.*omega*pow(kappa,2)*(299. + Complex(0.,8754.)*omega + 66090.*pow(omega,2)) + 266430.*pow(omega,3) + 510.*pow(kappa,4)*pow(omega,3) -
              1.*kappa*(Complex(0.,-285.) + 2524.*omega + Complex(0.,87072.)*pow(omega,2) + 338040.*pow(omega,3)))) - 1220.*pow(a,10)*pow(omega,4) -
        1.*pow(a,6)*(93. + Complex(0.,1.)*(-1285. + 474.*kappa)*omega + 4.*(8734. - 6495.*kappa + 972.*pow(kappa,2))*pow(omega,2) -
           Complex(0.,4.)*(-35481. + 37578.*kappa - 9831.*pow(kappa,2) + 590.*pow(kappa,3))*pow(omega,3) + 40.*(15237. - 18150.*kappa + 5832.*pow(kappa,2) - 530.*pow(kappa,3) + 11.*pow(kappa,4))*pow(omega,4)) -
        2.*(pow(kappa,5)*(1. + Complex(0.,18.)*omega + 224.*pow(omega,2)) + pow(kappa,3)*(694. + Complex(0.,5420.)*omega + 29568.*pow(omega,2) - Complex(0.,384.)*pow(omega,3) - 132096.*pow(omega,4)) +
           16.*pow(kappa,4)*(-4. - Complex(0.,45.)*omega - 230.*pow(omega,2) - Complex(0.,180.)*pow(omega,3) + 3360.*pow(omega,4)) +
           4.*pow(kappa,2)*(-591. - Complex(0.,3426.)*omega - 25600.*pow(omega,2) + Complex(0.,9792.)*pow(omega,3) + 3456.*pow(omega,4)) -
           4.*(315. + Complex(0.,1222.)*omega + 15912.*pow(omega,2) - Complex(0.,12048.)*pow(omega,3) + 16896.*pow(omega,4)) +
           kappa*(3009. + Complex(0.,13842.)*omega + 139936.*pow(omega,2) - Complex(0.,84096.)*pow(omega,3) + 132096.*pow(omega,4))) +
        pow(a,2)*(-3957. - Complex(0.,4677.)*omega - 235832.*pow(omega,2) - Complex(0.,156848.)*pow(omega,3) +
           2.*omega*pow(kappa,5)*(Complex(0.,1.) + 70.*omega + Complex(0.,180.)*pow(omega,2) + 2240.*pow(omega,3)) - 1.151424e6*pow(omega,4) -
           1.*pow(kappa,4)*(21. + Complex(0.,261.)*omega + 3944.*pow(omega,2) + Complex(0.,20880.)*pow(omega,3) + 16960.*pow(omega,4)) +
           4.*pow(kappa,3)*(141. + Complex(0.,858.)*omega + 11898.*pow(omega,2) + Complex(0.,41740.)*pow(omega,3) + 93216.*pow(omega,4)) -
           2.*pow(kappa,2)*(1731. + Complex(0.,5919.)*omega + 114240.*pow(omega,2) + Complex(0.,240000.)*pow(omega,3) + 888832.*pow(omega,4)) +
           2.*kappa*(3390. + Complex(0.,6807.)*omega + 205686.*pow(omega,2) + Complex(0.,253716.)*pow(omega,3) + 1.295104e6*pow(omega,4))) -
        2.*pow(a,4)*(-779. + Complex(0.,1893.)*omega - 74970.*pow(omega,2) - 1.*pow(kappa,4)*pow(omega,2)*(154. + Complex(0.,1911.)*omega + 7620.*pow(omega,2)) - Complex(0.,173271.)*pow(omega,3) +
           10.*(Complex(0.,1.) + 30.*omega)*pow(kappa,5)*pow(omega,3) + 2.*omega*pow(kappa,3)*(Complex(0.,11.) + 2252.*omega + Complex(0.,16704.)*pow(omega,2) + 63100.*pow(omega,3)) - 738388.*pow(omega,4) -
           1.*pow(kappa,2)*(137. - Complex(0.,111.)*omega + 36372.*pow(omega,2) + Complex(0.,176010.)*pow(omega,3) + 701160.*pow(omega,4)) +
           kappa*(731. - Complex(0.,1434.)*omega + 96288.*pow(omega,2) + Complex(0.,319422.)*pow(omega,3) + 1.307868e6*pow(omega,4)))))*pow(omega,7);
}

Complex gsn_asymptotic_derivative_initial_sum_term_12(const double &a, const int &, const int &mTemp, const double &omega, const double &lambda, const double &){
	double kappa = sqrt(1. - pow(a, 2));
	Complex m = Complex(mTemp);
	return -0.0026041666666666665*(Complex(0.,1.)*(35.*omega*pow(a,8) - 5.*pow(a,6)*(Complex(0.,-83.) + kappa*(Complex(0.,43.) - 232.*omega) + 416.*omega + 24.*omega*pow(kappa,2)) +
        10.*pow(a,4)*(Complex(0.,-357.) + kappa*(Complex(0.,537.) - 1800.*omega) + 1593.*omega + 3.*(Complex(0.,-69.) + 194.*omega)*pow(kappa,2) + (Complex(0.,19.) - 56.*omega)*pow(kappa,3) +
           omega*pow(kappa,4)) + m*(35.*pow(a,7) - 30.*pow(a,5)*(31. - 22.*kappa + 3.*pow(kappa,2)) + 5.*pow(a,3)*(681. - 924.*kappa + 366.*pow(kappa,2) - 44.*pow(kappa,3) + pow(kappa,4)) -
           4.*a*(743. - 1452.*kappa + 930.*pow(kappa,2) - 220.*pow(kappa,3) + 15.*pow(kappa,4))) +
        pow(a,2)*(Complex(0.,7617.) - 35528.*omega + 3.*kappa*(Complex(0.,-6195.) + 19904.*omega) - 30.*(Complex(0.,-467.) + 1080.*omega)*pow(kappa,2) + 10.*(Complex(0.,-381.) + 640.*omega)*pow(kappa,3) -
           5.*(Complex(0.,-65.) + 72.*omega)*pow(kappa,4) - Complex(0.,5.)*pow(kappa,5)) +
        Complex(0.,4.)*(-1153. + kappa*(3851. + Complex(0.,12896.)*omega) - Complex(0.,5848.)*omega + (-4266. - Complex(0.,9744.)*omega)*pow(kappa,2) + 2.*(955. + Complex(0.,1456.)*omega)*pow(kappa,3) +
           (-325. - Complex(0.,280.)*omega)*pow(kappa,4) + 15.*pow(kappa,5)))*pow(lambda,5) +
     2.*pow(lambda,4)*(5.*omega*(Complex(0.,-22.) + (-281. + 83.*kappa)*omega)*pow(a,8) +
        a*m*(Complex(0.,3979.) + 41568.*omega - 1.*kappa*(Complex(0.,10233.) + 101864.*omega) + 2.*(Complex(0.,-44.) + (-687. + 265.*kappa)*omega)*pow(a,6) + 42.*(Complex(0.,215.) + 2096.*omega)*pow(kappa,2) -
           30.*(Complex(0.,107.) + 1064.*omega)*pow(kappa,3) + pow(a,4)*(-15.*kappa*(Complex(0.,115.) + 1344.*omega) + 3.*(Complex(0.,603.) + 5972.*omega) + 30.*(Complex(0.,11.) + 202.*omega)*pow(kappa,2) -
              440.*omega*pow(kappa,3)) + 5.*(Complex(0.,83.) + 896.*omega)*pow(kappa,4) - 1.*(Complex(0.,13.) + 168.*omega)*pow(kappa,5) +
           2.*pow(a,2)*(-3.*(Complex(0.,901.) + 8783.*omega) + kappa*(Complex(0.,4866.) + 47745.*omega) - 6.*(Complex(0.,446.) + 4695.*omega)*pow(kappa,2) + (Complex(0.,486.) + 6170.*omega)*pow(kappa,3) -
              1.*(Complex(0.,21.) + 425.*omega)*pow(kappa,4) + 5.*omega*pow(kappa,5))) +
        pow(a,2)*(-1683. + 3405.*kappa + 21.*(-11. + 5.*kappa)*pow(a,4) - 2310.*pow(kappa,2) + 610.*pow(kappa,3) - 12.*pow(a,2)*(-121. + 155.*kappa - 55.*pow(kappa,2) + 5.*pow(kappa,3)) - 55.*pow(kappa,4) +
           pow(kappa,5))*pow(m,2) + pow(a,6)*(53. + Complex(0.,5357.)*omega + 60.*omega*(Complex(0.,7.) + 123.*omega)*pow(kappa,2) + 33492.*pow(omega,2) - 400.*pow(kappa,3)*pow(omega,2) -
           5.*kappa*(-5. + Complex(0.,715.)*omega + 6216.*pow(omega,2))) + 2.*pow(a,4)*
         (336. - Complex(0.,17367.)*omega - 2.*omega*(Complex(0.,12.) + 325.*omega)*pow(kappa,4) + pow(kappa,2)*(699. - Complex(0.,8949.)*omega - 71580.*pow(omega,2)) - 93906.*pow(omega,2) +
           5.*pow(kappa,5)*pow(omega,2) + pow(kappa,3)*(-77. + Complex(0.,1009.)*omega + 12530.*pow(omega,2)) + 3.*kappa*(-432. + Complex(0.,7833.)*omega + 48635.*pow(omega,2))) +
        pow(a,2)*(-3225. + Complex(0.,67391.)*omega + kappa*(13683. - Complex(0.,134901.)*omega - 772408.*pow(omega,2)) + pow(kappa,5)*(7. - Complex(0.,41.)*omega - 504.*pow(omega,2)) + 358912.*pow(omega,2) +
           7.*pow(kappa,4)*(-59. + Complex(0.,245.)*omega + 2720.*pow(omega,2)) - 2.*pow(kappa,3)*(-2139. + Complex(0.,10425.)*omega + 85960.*pow(omega,2)) +
           6.*pow(kappa,2)*(-2223. + Complex(0.,14485.)*omega + 95312.*pow(omega,2))) +
        2.*(1421. - Complex(0.,19690.)*omega + pow(kappa,4)*(985. - Complex(0.,2450.)*omega - 21168.*pow(omega,2)) - 108144.*pow(omega,2) + 3.*pow(kappa,5)*(-17. + Complex(0.,42.)*omega + 336.*pow(omega,2)) +
           14.*pow(kappa,3)*(-365. + Complex(0.,1210.)*omega + 8784.*pow(omega,2)) - 6.*pow(kappa,2)*(-1639. + Complex(0.,7686.)*omega + 48144.*pow(omega,2)) +
           kappa*(-7175. + Complex(0.,51382.)*omega + 294192.*pow(omega,2)))) + Complex(0.,4.)*pow(lambda,3)*
      (2.*pow(a,2)*(Complex(0.,1683.) - 14860.*omega + 15.*kappa*(Complex(0.,-227.) + 1936.*omega) + 175.*omega*pow(a,6) - 30.*(Complex(0.,-77.) + 620.*omega)*pow(kappa,2) -
           3.*pow(a,4)*(Complex(0.,-77.) + 1550.*omega - 5.*kappa*(Complex(0.,-7.) + 220.*omega) + 150.*omega*pow(kappa,2)) + 10.*(Complex(0.,-61.) + 440.*omega)*pow(kappa,3) +
           (Complex(0.,55.) - 300.*omega)*pow(kappa,4) + pow(a,2)*(-60.*kappa*(Complex(0.,-31.) + 385.*omega) + 3.*(Complex(0.,-484.) + 5675.*omega) + 30.*(Complex(0.,-22.) + 305.*omega)*pow(kappa,2) -
              20.*(Complex(0.,-3.) + 55.*omega)*pow(kappa,3) + 25.*omega*pow(kappa,4)) - Complex(0.,1.)*pow(kappa,5))*pow(m,2) +
        5.*omega*pow(a,8)*(-17. + Complex(0.,4.)*(2. + 57.*kappa)*omega + 4.*(1703. - 760.*kappa + 57.*pow(kappa,2))*pow(omega,2)) -
        1.*m*(240.*pow(a,9)*pow(omega,2) + pow(a,7)*(281. + Complex(0.,2.)*(-3989. + 1785.*kappa)*omega + 20.*(197. - 118.*kappa + 21.*pow(kappa,2))*pow(omega,2)) -
           2.*pow(a,5)*(3204. - Complex(0.,40410.)*omega + 2.*(Complex(0.,713.) - 1510.*omega)*omega*pow(kappa,3) + 44085.*pow(omega,2) + 85.*pow(kappa,4)*pow(omega,2) +
              3.*pow(kappa,2)*(155. - Complex(0.,5794.)*omega + 7370.*pow(omega,2)) - 3.*kappa*(905. - Complex(0.,17146.)*omega + 18340.*pow(omega,2))) +
           pow(a,3)*(Complex(0.,70.)*omega*pow(kappa,5) + pow(kappa,4)*(99. - Complex(0.,4950.)*omega + 5880.*pow(omega,2)) - 4.*pow(kappa,3)*(651. - Complex(0.,16055.)*omega + 19080.*pow(omega,2)) +
              6.*pow(kappa,2)*(2699. - Complex(0.,44290.)*omega + 51160.*pow(omega,2)) - 6.*kappa*(5554. - Complex(0.,67845.)*omega + 81392.*pow(omega,2)) +
              3.*(7009. - Complex(0.,67146.)*omega + 89192.*pow(omega,2))) + 2.*a*(-8437. + Complex(0.,70676.)*omega + 13.*(1. - Complex(0.,36.)*omega)*pow(kappa,5) - 102848.*pow(omega,2) -
              5.*pow(kappa,4)*(101. - Complex(0.,2196.)*omega + 1792.*pow(omega,2)) - 6.*pow(kappa,2)*(2435. - Complex(0.,30188.)*omega + 34272.*pow(omega,2)) +
              pow(kappa,3)*(4530. - Complex(0.,71560.)*omega + 73024.*pow(omega,2)) + kappa*(18945. - Complex(0.,191268.)*omega + 243136.*pow(omega,2)))) - 530.*pow(a,10)*pow(omega,3) -
        Complex(0.,8.)*(-1. + kappa)*(1019. + Complex(0.,15593.)*omega + 64488.*pow(omega,2) + pow(kappa,4)*(-3. + Complex(0.,63.)*omega + 504.*pow(omega,2)) +
           24.*pow(kappa,2)*(13. + Complex(0.,439.)*omega + 2212.*pow(omega,2) + Complex(0.,912.)*pow(omega,3)) +
           2.*pow(kappa,3)*(1. - Complex(0.,791.)*omega - 5208.*pow(omega,2) - Complex(0.,4032.)*pow(omega,3)) - Complex(0.,25344.)*pow(omega,3) +
           Complex(0.,6.)*kappa*(Complex(0.,195.) - 3875.*omega + Complex(0.,16952.)*pow(omega,2) + 768.*pow(omega,3))) +
        Complex(0.,1.)*pow(a,2)*(-16401. - Complex(0.,191638.)*omega - 840664.*pow(omega,2) + pow(kappa,5)*(1. + Complex(0.,82.)*omega + 408.*pow(omega,2)) +
           pow(kappa,4)*(-149. - Complex(0.,4510.)*omega - 27000.*pow(omega,2) - Complex(0.,35840.)*pow(omega,3)) +
           2.*pow(kappa,3)*(1437. + Complex(0.,28050.)*omega + 140440.*pow(omega,2) + Complex(0.,66752.)*pow(omega,3)) +
           kappa*(28389. + Complex(0.,368874.)*omega + 1.62924e6*pow(omega,2) - Complex(0.,1.00544e6)*pow(omega,3)) + Complex(0.,720256.)*pow(omega,3) +
           Complex(0.,6.)*pow(kappa,2)*(Complex(0.,2559.) - 38930.*omega + Complex(0.,177736.)*pow(omega,2) + 36288.*pow(omega,3))) +
        2.*pow(a,4)*(Complex(0.,4683.) - 42429.*omega + Complex(0.,205095.)*pow(omega,2) + Complex(0.,15.)*pow(kappa,5)*pow(omega,2) - 1.*omega*pow(kappa,4)*(63. + Complex(0.,205.)*omega + 1800.*pow(omega,2)) +
           376792.*pow(omega,3) + pow(kappa,3)*(Complex(0.,-131.) + 2498.*omega - Complex(0.,7210.)*pow(omega,2) + 1120.*pow(omega,3)) +
           3.*pow(kappa,2)*(Complex(0.,569.) - 7076.*omega + Complex(0.,28690.)*pow(omega,2) + 43920.*pow(omega,3)) -
           3.*kappa*(Complex(0.,1821.) - 18546.*omega + Complex(0.,84295.)*pow(omega,2) + 152224.*pow(omega,3))) +
        pow(a,6)*(Complex(0.,-1351.) + 10042.*omega + 12.*omega*pow(kappa,2)*(64. + Complex(0.,387.)*omega - 4370.*pow(omega,2)) - Complex(0.,56676.)*pow(omega,2) +
           4.*(Complex(0.,-263.) + 500.*omega)*pow(kappa,3)*pow(omega,2) - 287100.*pow(omega,3) + 20.*pow(kappa,4)*pow(omega,3) +
           kappa*(Complex(0.,595.) - 6454.*omega + Complex(0.,22044.)*pow(omega,2) + 249840.*pow(omega,3)))) +
     Complex(0.,48.)*lambda*omega*((105.*omega*pow(a,9) - 6.*pow(a,7)*(Complex(0.,-77.) + 465.*omega - 5.*kappa*(Complex(0.,-7.) + 66.*omega) + 45.*omega*pow(kappa,2)) +
           3.*pow(a,5)*(Complex(0.,-968.) + 3405.*omega - 20.*kappa*(Complex(0.,-62.) + 231.*omega) + 10.*(Complex(0.,-44.) + 183.*omega)*pow(kappa,2) + (Complex(0.,40.) - 220.*omega)*pow(kappa,3) +
              5.*omega*pow(kappa,4)) - Complex(0.,2.)*pow(a,3)*(-1683. + kappa*(3405. + Complex(0.,8712.)*omega) - Complex(0.,4458.)*omega + (-2310. - Complex(0.,5580.)*omega)*pow(kappa,2) +
              10.*(61. + Complex(0.,132.)*omega)*pow(kappa,3) + (-55. - Complex(0.,90.)*omega)*pow(kappa,4) + pow(kappa,5)))*pow(m,3) -
        1.*pow(a,10)*pow(omega,2)*(616. + Complex(0.,1.)*(-39. + 991.*kappa)*omega + 4.*(5867. - 1290.*kappa + 43.*pow(kappa,2))*pow(omega,2)) +
        pow(a,2)*pow(m,2)*(-665.*pow(a,8)*pow(omega,2) - 1.*pow(a,6)*(316. + Complex(0.,5.)*(-1395. + 551.*kappa)*omega + 20.*(-849. + 332.*kappa + 5.*pow(kappa,2))*pow(omega,2)) +
           2.*pow(a,4)*(3669. - Complex(0.,35559.)*omega + (Complex(0.,869.) - 3260.*omega)*omega*pow(kappa,3) - 29580.*pow(omega,2) + 140.*pow(kappa,4)*pow(omega,2) +
              3.*pow(kappa,2)*(170. - Complex(0.,4031.)*omega + 2760.*pow(omega,2)) + 3.*kappa*(-1015. + Complex(0.,13469.)*omega + 4860.*pow(omega,2))) +
           2.*(9923. - Complex(0.,62598.)*omega + (-13. + Complex(0.,258.)*omega)*pow(kappa,5) + pow(kappa,3)*(-4970. + Complex(0.,46740.)*omega - 51968.*pow(omega,2)) +
              kappa*(-21849. + Complex(0.,153866.)*omega - 42880.*pow(omega,2)) + 2160.*pow(omega,2) + pow(kappa,4)*(535. - Complex(0.,6430.)*omega + 9520.*pow(omega,2)) +
              2.*pow(kappa,2)*(8235. - Complex(0.,65886.)*omega + 41776.*pow(omega,2))) +
           pow(a,2)*(-24432. + Complex(0.,177981.)*omega - Complex(0.,37.)*omega*pow(kappa,5) + pow(kappa,4)*(-104. + Complex(0.,2785.)*omega - 7440.*pow(omega,2)) + 46544.*pow(omega,2) +
              3.*kappa*(12648. - Complex(0.,108035.)*omega + 8672.*pow(omega,2)) - 6.*pow(kappa,2)*(3004. - Complex(0.,31595.)*omega + 20320.*pow(omega,2)) +
              pow(kappa,3)*(2824. - Complex(0.,40690.)*omega + 63200.*pow(omega,2)))) -
        1.*m*(2.*omega*pow(a,9)*(28. + Complex(0.,2.)*(-1626. + 221.*kappa)*omega - 5.*(2687. - 1186.*kappa + 71.*pow(kappa,2))*pow(omega,2)) + 245.*pow(a,11)*pow(omega,3) +
           pow(a,7)*(Complex(0.,1406.) - 10397.*omega + Complex(0.,143140.)*pow(omega,2) + 220.*(Complex(0.,1.) + 19.*omega)*pow(kappa,3)*pow(omega,2) +
              2.*omega*pow(kappa,2)*(-750. + Complex(0.,6142.)*omega + 16395.*pow(omega,2)) + 286885.*pow(omega,3) - 355.*pow(kappa,4)*pow(omega,3) -
              2.*kappa*(Complex(0.,345.) - 4200.*omega + Complex(0.,49006.)*pow(omega,2) + 116270.*pow(omega,3))) +
           4.*a*(Complex(0.,-2183.) + 34559.*omega - Complex(0.,162656.)*pow(omega,2) + Complex(0.,9.)*pow(kappa,5)*(1. + Complex(0.,14.)*omega + 112.*pow(omega,2)) - 99712.*pow(omega,3) +
              pow(kappa,4)*(Complex(0.,-315.) + 4335.*omega - Complex(0.,19008.)*pow(omega,2) + 42112.*pow(omega,3)) +
              pow(kappa,2)*(Complex(0.,-5934.) + 78330.*omega - Complex(0.,335072.)*pow(omega,2) + 57600.*pow(omega,3)) -
              2.*pow(kappa,3)*(Complex(0.,-1145.) + 14936.*omega - Complex(0.,61408.)*pow(omega,2) + 65408.*pow(omega,3)) +
              kappa*(Complex(0.,6165.) - 87114.*omega + Complex(0.,392912.)*pow(omega,2) + 130816.*pow(omega,3))) +
           2.*pow(a,5)*(Complex(0.,-4920.) + 45933.*omega - Complex(0.,355407.)*pow(omega,2) + Complex(0.,5.)*pow(kappa,5)*pow(omega,2) + omega*pow(kappa,4)*(156. - Complex(0.,419.)*omega + 9010.*pow(omega,2)) +
              pow(kappa,3)*(Complex(0.,214.) - 4514.*omega + Complex(0.,19082.)*pow(omega,2) - 47880.*pow(omega,3)) - 451774.*pow(omega,3) -
              3.*pow(kappa,2)*(Complex(0.,806.) - 9859.*omega + Complex(0.,56490.)*pow(omega,2) + 23820.*pow(omega,3)) +
              kappa*(Complex(0.,6672.) - 65700.*omega + Complex(0.,455553.)*pow(omega,2) + 497304.*pow(omega,3))) +
           pow(a,3)*(Complex(0.,17394.) - 210699.*omega + 6.*pow(kappa,5)*(Complex(0.,-1.) + 18.*omega - Complex(0.,164.)*pow(omega,2)) + Complex(0.,1.204968e6)*pow(omega,2) +
              pow(kappa,4)*(Complex(0.,474.) - 7395.*omega + Complex(0.,30600.)*pow(omega,2) - 117440.*pow(omega,3)) + 1.054016e6*pow(omega,3) +
              2.*pow(kappa,2)*(Complex(0.,12522.) - 143805.*omega + Complex(0.,692424.)*pow(omega,2) + 5824.*pow(omega,3)) +
              4.*pow(kappa,3)*(Complex(0.,-1551.) + 19675.*omega - Complex(0.,83980.)*pow(omega,2) + 113088.*pow(omega,3)) -
              2.*kappa*(Complex(0.,18447.) - 210288.*omega + Complex(0.,1.121788e6)*pow(omega,2) + 677504.*pow(omega,3)))) + 253.*pow(a,12)*pow(omega,4) +
        pow(a,8)*(80. + Complex(0.,2.)*(-667. + 381.*kappa)*omega + (30395. - 17986.*kappa + 2388.*pow(kappa,2))*pow(omega,2) +
           Complex(0.,2.)*(-3909. + 25247.*kappa - 9815.*pow(kappa,2) + 661.*pow(kappa,3))*pow(omega,3) + 10.*(28007. - 10784.*kappa - 270.*pow(kappa,2) + 7.*pow(kappa,4))*pow(omega,4)) +
        Complex(0.,4.)*(Complex(0.,-2016.) - 5379.*omega + 45.*omega*pow(kappa,5) - Complex(0.,41704.)*pow(omega,2) +
           omega*pow(kappa,4)*(-2415. - Complex(0.,7400.)*omega + 11264.*pow(omega,2) - Complex(0.,27648.)*pow(omega,3)) + 21760.*pow(omega,3) +
           2.*pow(kappa,3)*(Complex(0.,240.) + 6177.*omega + Complex(0.,16240.)*pow(omega,2) + 512.*pow(omega,3) - Complex(0.,27648.)*pow(omega,4)) +
           kappa*(Complex(0.,4128.) + 18753.*omega + Complex(0.,90976.)*pow(omega,2) - 1024.*pow(omega,3) + Complex(0.,55296.)*pow(omega,4)) - Complex(0.,107520.)*pow(omega,4) +
           Complex(0.,6.)*pow(kappa,2)*(-432. + Complex(0.,3877.)*omega - 12424.*pow(omega,2) + Complex(0.,5504.)*pow(omega,3) + 22528.*pow(omega,4))) +
        2.*pow(a,4)*(5400. - Complex(0.,20724.)*omega + 266814.*pow(omega,2) + Complex(0.,1.)*omega*pow(kappa,5)*(3. + Complex(0.,41.)*omega + 594.*pow(omega,2)) - Complex(0.,75126.)*pow(omega,3) +
           omega*pow(kappa,4)*(Complex(0.,-597.) + 6490.*omega - Complex(0.,9710.)*pow(omega,2) + 16880.*pow(omega,3)) +
           kappa*(-6120. + Complex(0.,43902.)*omega - 424317.*pow(omega,2) + Complex(0.,551066.)*pow(omega,3) - 479232.*pow(omega,4)) +
           pow(kappa,2)*(1800. - Complex(0.,30867.)*omega + 242040.*pow(omega,2) - Complex(0.,556092.)*pow(omega,3) - 400288.*pow(omega,4)) + 941488.*pow(omega,4) +
           pow(kappa,3)*(-120. + Complex(0.,8067.)*omega - 63370.*pow(omega,2) + Complex(0.,152180.)*pow(omega,3) + 132992.*pow(omega,4))) +
        pow(a,6)*(-2304. + Complex(0.,13653.)*omega - 218956.*pow(omega,2) + pow(kappa,4)*(-464. + Complex(0.,2127.)*omega - 4120.*pow(omega,2))*pow(omega,2) + Complex(0.,72179.)*pow(omega,3) -
           Complex(0.,35.)*pow(kappa,5)*pow(omega,3) - 2.*omega*pow(kappa,3)*(Complex(0.,394.) - 6466.*omega + Complex(0.,25983.)*pow(omega,2) + 16480.*pow(omega,3)) - 1.104952e6*pow(omega,4) +
           6.*pow(kappa,2)*(-24. + Complex(0.,1314.)*omega - 15178.*pow(omega,2) + Complex(0.,49445.)*pow(omega,3) + 25240.*pow(omega,4)) +
           kappa*(1392. - Complex(0.,19533.)*omega + 247492.*pow(omega,2) - Complex(0.,422399.)*pow(omega,3) + 532224.*pow(omega,4))) +
        pow(a,2)*(-16576. + Complex(0.,50207.)*omega + 3.*omega*pow(kappa,5)*(Complex(0.,-17.) + 168.*omega - Complex(0.,1344.)*pow(omega,2)) - 514880.*pow(omega,2) +
           omega*pow(kappa,4)*(Complex(0.,7035.) - 44480.*omega + Complex(0.,19584.)*pow(omega,2) - 103936.*pow(omega,3)) + Complex(0.,22016.)*pow(omega,3) - 1.459712e6*pow(omega,4) -
           2.*pow(kappa,3)*(-800. + Complex(0.,26615.)*omega - 133976.*pow(omega,2) + Complex(0.,208896.)*pow(omega,3) + 251392.*pow(omega,4)) +
           2.*pow(kappa,2)*(-6240. + Complex(0.,67443.)*omega - 367584.*pow(omega,2) + Complex(0.,597440.)*pow(omega,3) + 622080.*pow(omega,4)) +
           kappa*(26688. - Complex(0.,139551.)*omega + 1.005976e6*pow(omega,2) - Complex(0.,852544.)*pow(omega,3) + 723968.*pow(omega,4)))) +
     288.*pow(omega,2)*(pow(a,3)*(Complex(0.,9923.) + 41568.*omega - 1.*kappa*(Complex(0.,21849.) + 101864.*omega) + 2.*(Complex(0.,-79.) + 8.*(-57. + 20.*kappa)*omega)*pow(a,6) +
           6.*(Complex(0.,2745.) + 14672.*omega)*pow(kappa,2) - 70.*(Complex(0.,71.) + 456.*omega)*pow(kappa,3) +
           pow(a,4)*(Complex(0.,3669.) + 15012.*omega - 15.*kappa*(Complex(0.,203.) + 1096.*omega) + 30.*(Complex(0.,17.) + 158.*omega)*pow(kappa,2) - 320.*omega*pow(kappa,3)) +
           5.*(Complex(0.,107.) + 896.*omega)*pow(kappa,4) - 1.*(Complex(0.,13.) + 168.*omega)*pow(kappa,5) +
           4.*pow(a,2)*(-3.*(Complex(0.,1018.) + 4111.*omega) + 3.*kappa*(Complex(0.,1581.) + 7390.*omega) - 3.*(Complex(0.,751.) + 4310.*omega)*pow(kappa,2) + (Complex(0.,353.) + 2780.*omega)*pow(kappa,3) -
              1.*(Complex(0.,13.) + 185.*omega)*pow(kappa,4) + 2.*omega*pow(kappa,5)))*pow(m,3) +
        pow(a,4)*(-1683. + 3405.*kappa + 21.*(-11. + 5.*kappa)*pow(a,4) - 2310.*pow(kappa,2) + 610.*pow(kappa,3) - 12.*pow(a,2)*(-121. + 155.*kappa - 55.*pow(kappa,2) + 5.*pow(kappa,3)) - 55.*pow(kappa,4) +
           pow(kappa,5))*pow(m,4) + pow(a,2)*pow(m,2)*(-4033. - Complex(0.,87508.)*omega + 4.*omega*(Complex(0.,49.) + (698. - 235.*kappa)*omega)*pow(a,8) - 199648.*pow(omega,2) +
           pow(kappa,5)*(19. + Complex(0.,252.)*omega + 2016.*pow(omega,2)) - 5.*pow(kappa,4)*(137. + Complex(0.,2004.)*omega + 10976.*pow(omega,2)) +
           3.*kappa*(3869. + Complex(0.,71972.)*omega + 185376.*pow(omega,2)) + pow(kappa,3)*(4710. + Complex(0.,70232.)*omega + 278208.*pow(omega,2)) -
           2.*pow(kappa,2)*(5793. + Complex(0.,94500.)*omega + 290912.*pow(omega,2)) -
           1.*pow(a,4)*(3603. + Complex(0.,39930.)*omega + omega*(Complex(0.,36.) + 935.*omega)*pow(kappa,4) + 39171.*pow(omega,2) + 9.*pow(kappa,5)*pow(omega,2) -
              2.*pow(kappa,3)*(101. + Complex(0.,925.)*omega + 5935.*pow(omega,2)) - 9.*kappa*(589. + Complex(0.,5262.)*omega + 6035.*pow(omega,2)) +
              6.*pow(kappa,2)*(349. + Complex(0.,2759.)*omega + 6105.*pow(omega,2))) +
           pow(a,6)*(523. + Complex(0.,1147.)*omega - 32.*omega*(Complex(0.,6.) + 115.*omega)*pow(kappa,2) - 17288.*pow(omega,2) + 20.*pow(kappa,3)*pow(omega,2) +
              3.*kappa*(-95. + Complex(0.,97.)*omega + 6220.*pow(omega,2))) + pow(a,2)*
            (6837. + Complex(0.,118505.)*omega + 228928.*pow(omega,2) - 3.*pow(kappa,5)*(1. + Complex(0.,5.)*omega + 56.*pow(omega,2)) + 3.*pow(kappa,4)*(79. + Complex(0.,935.)*omega + 6560.*pow(omega,2)) +
              6.*pow(kappa,2)*(1837. + Complex(0.,23155.)*omega + 67312.*pow(omega,2)) - 2.*pow(kappa,3)*(1461. + Complex(0.,17135.)*omega + 73560.*pow(omega,2)) -
              1.*kappa*(15267. + Complex(0.,221571.)*omega + 497384.*pow(omega,2)))) + (Complex(0.,-484.) + (-3357. + 591.*kappa)*omega)*pow(a,12)*pow(omega,3) +
        omega*pow(a,10)*(Complex(0.,80.) + (539. - 205.*kappa)*omega + Complex(0.,1.)*(27433. - 9431.*kappa + 328.*pow(kappa,2))*pow(omega,2) +
           (96004. - 50616.*kappa + 6148.*pow(kappa,2) - 224.*pow(kappa,3))*pow(omega,3)) -
        1.*a*m*(Complex(0.,-13955.) - 23712.*omega + (Complex(0.,13.) + 168.*omega)*pow(kappa,5) - Complex(0.,7424.)*pow(omega,2) + 2.*(Complex(0.,-123.) + 2.*(-259. + 75.*kappa)*omega)*pow(a,10)*pow(omega,2) +
           kappa*(Complex(0.,30105.) + 71144.*omega + Complex(0.,30208.)*pow(omega,2) - 106496.*pow(omega,3)) + 49152.*pow(omega,3) +
           pow(kappa,2)*(Complex(0.,-21654.) - 79968.*omega - Complex(0.,93184.)*pow(omega,2) + 20480.*pow(omega,3)) +
           2.*pow(kappa,3)*(Complex(0.,2965.) + 19800.*omega + Complex(0.,52480.)*pow(omega,2) + 53248.*pow(omega,3)) -
           1.*pow(kappa,4)*(Complex(0.,535.) + 7360.*omega + Complex(0.,34560.)*pow(omega,2) + 69632.*pow(omega,3)) +
           pow(a,8)*(Complex(0.,-80.) + (2738. - 978.*kappa)*omega + Complex(0.,1.)*(15809. - 10505.*kappa + 1006.*pow(kappa,2))*pow(omega,2) -
              4.*(-12847. + 11780.*kappa - 2555.*pow(kappa,2) + 90.*pow(kappa,3))*pow(omega,3)) +
           2.*pow(a,6)*(Complex(0.,1051.) - 10533.*omega - Complex(0.,60710.)*pow(omega,2) + 15.*(Complex(0.,2.) - 21.*omega)*pow(kappa,4)*pow(omega,2) +
              2.*omega*pow(kappa,3)*(35. + Complex(0.,836.)*omega + 8165.*pow(omega,2)) - 181127.*pow(omega,3) + 5.*pow(kappa,5)*pow(omega,3) -
              2.*pow(kappa,2)*(Complex(0.,-36.) + 1083.*omega + Complex(0.,13256.)*pow(omega,2) + 62215.*pow(omega,3)) +
              kappa*(Complex(0.,-636.) + 10103.*omega + Complex(0.,80320.)*pow(omega,2) + 275425.*pow(omega,3))) -
           1.*pow(a,4)*(Complex(0.,11373.) - 36174.*omega + omega*pow(kappa,4)*(174. + Complex(0.,657.)*omega - 13120.*pow(omega,2)) - Complex(0.,263451.)*pow(omega,2) +
              3.*omega*pow(kappa,5)*(2. + Complex(0.,23.)*omega + 280.*pow(omega,2)) - 794272.*pow(omega,3) +
              2.*pow(kappa,3)*(Complex(0.,-108.) + 50.*omega + Complex(0.,18397.)*pow(omega,2) + 127160.*pow(omega,3)) -
              6.*pow(kappa,2)*(Complex(0.,-577.) + 2884.*omega + Complex(0.,45425.)*pow(omega,2) + 185808.*pow(omega,3)) +
              3.*kappa*(Complex(0.,-4111.) + 17066.*omega + Complex(0.,168859.)*pow(omega,2) + 550360.*pow(omega,3))) +
           2.*pow(a,2)*(Complex(0.,11084.) + 431.*omega - Complex(0.,80980.)*pow(omega,2) + omega*pow(kappa,5)*(13. + Complex(0.,252.)*omega + 2016.*pow(omega,2)) +
              pow(kappa,4)*(Complex(0.,26.) + 1475.*omega + Complex(0.,6300.)*pow(omega,2) - 4704.*pow(omega,3)) - 281568.*pow(omega,3) +
              2.*pow(kappa,3)*(Complex(0.,-633.) - 4395.*omega + Complex(0.,2956.)*pow(omega,2) + 85344.*pow(omega,3)) -
              2.*pow(kappa,2)*(Complex(0.,-4293.) - 7143.*omega + Complex(0.,62020.)*pow(omega,2) + 304736.*pow(omega,3)) +
              kappa*(Complex(0.,-17790.) - 7095.*omega + Complex(0.,197548.)*pow(omega,2) + 716896.*pow(omega,3)))) +
        pow(a,8)*(80. - Complex(0.,12.)*omega*(151. - 112.*kappa + 12.*pow(kappa,2)) + (-1901. + 2227.*kappa + 390.*pow(kappa,2) - 14.*pow(kappa,3))*pow(omega,2) -
           Complex(0.,2.)*(113717. - 72635.*kappa + 4009.*pow(kappa,2) + 799.*pow(kappa,3) + 6.*pow(kappa,4))*pow(omega,3) +
           2.*(-339482. + 284545.*kappa - 52220.*pow(kappa,2) + 2850.*pow(kappa,3) - 290.*pow(kappa,4) + 5.*pow(kappa,5))*pow(omega,4)) +
        2.*(841. - Complex(0.,34826.)*omega - 103792.*pow(omega,2) + 9.*pow(kappa,5)*(1. + Complex(0.,14.)*omega + 112.*pow(omega,2)) - Complex(0.,10240.)*pow(omega,3) +
           pow(kappa,3)*(1570. + Complex(0.,38956.)*omega + 71264.*pow(omega,2) - Complex(0.,10240.)*pow(omega,3) - 356352.*pow(omega,4)) - 190464.*pow(omega,4) +
           6.*pow(kappa,2)*(-341. - Complex(0.,15078.)*omega - 40208.*pow(omega,2) + Complex(0.,5120.)*pow(omega,3) + 16384.*pow(omega,4)) +
           pow(kappa,4)*(-315. - Complex(0.,6450.)*omega - 5552.*pow(omega,2) - Complex(0.,20480.)*pow(omega,3) + 92160.*pow(omega,4)) +
           kappa*(-27. + Complex(0.,92598.)*omega + 278320.*pow(omega,2) + Complex(0.,10240.)*pow(omega,3) + 356352.*pow(omega,4))) +
        2.*pow(a,4)*(2166. - Complex(0.,28745.)*omega - 71257.*pow(omega,2) + 14.*pow(kappa,5)*pow(omega,2)*(1. + Complex(0.,9.)*omega + 72.*pow(omega,2)) - Complex(0.,342698.)*pow(omega,3) +
           omega*pow(kappa,4)*(Complex(0.,-394.) + 4235.*omega + Complex(0.,6510.)*pow(omega,2) + 47824.*pow(omega,3)) - 1.226736e6*pow(omega,4) +
           pow(kappa,2)*(429. - Complex(0.,28689.)*omega - 7962.*pow(omega,2) + Complex(0.,163996.)*pow(omega,3) + 50080.*pow(omega,4)) -
           1.*pow(kappa,3)*(7. - Complex(0.,6449.)*omega + 25260.*pow(omega,2) + Complex(0.,135572.)*pow(omega,3) + 391328.*pow(omega,4)) +
           kappa*(-2046. + Complex(0.,48843.)*omega + 106390.*pow(omega,2) + Complex(0.,309622.)*pow(omega,3) + 1.56216e6*pow(omega,4))) -
        1.*pow(a,6)*(1421. - Complex(0.,12877.)*omega - 16805.*pow(omega,2) + pow(kappa,5)*pow(omega,2)*(3. + Complex(0.,41.)*omega + 504.*pow(omega,2)) +
           pow(kappa,4)*pow(omega,2)*(331. + Complex(0.,957.)*omega + 2080.*pow(omega,2)) - Complex(0.,625743.)*pow(omega,3) -
           2.*omega*pow(kappa,3)*(Complex(0.,-288.) + 4241.*omega + Complex(0.,25863.)*pow(omega,2) + 43000.*pow(omega,3)) - 1.896192e6*pow(omega,4) -
           6.*pow(kappa,2)*(-24. + Complex(0.,1038.)*omega - 2923.*pow(omega,2) - Complex(0.,6235.)*pow(omega,3) + 50064.*pow(omega,4)) +
           kappa*(-987. + Complex(0.,16551.)*omega + 13991.*pow(omega,2) + Complex(0.,525173.)*pow(omega,3) + 2.058936e6*pow(omega,4))) -
        1.*pow(a,2)*(4567. - Complex(0.,111375.)*omega - 321600.*pow(omega,2) + pow(kappa,5)*(3. + Complex(0.,41.)*omega + 504.*pow(omega,2)) - Complex(0.,275200.)*pow(omega,3) +
           pow(kappa,3)*(1742. + Complex(0.,62450.)*omega + 8080.*pow(omega,2) - Complex(0.,332288.)*pow(omega,3) - 1.507328e6*pow(omega,4)) - 1.511424e6*pow(omega,4) +
           pow(kappa,4)*(-237. - Complex(0.,7235.)*omega + 16480.*pow(omega,2) + Complex(0.,5888.)*pow(omega,3) + 290816.*pow(omega,4)) +
           2.*pow(kappa,2)*(-1101. - Complex(0.,94287.)*omega - 187440.*pow(omega,2) + Complex(0.,212480.)*pow(omega,3) + 327680.*pow(omega,4)) +
           kappa*(-3201. + Complex(0.,241653.)*omega + 679736.*pow(omega,2) + Complex(0.,217600.)*pow(omega,3) + 2.220032e6*pow(omega,4)))) +
     8.*pow(lambda,2)*(10.*omega*pow(a,3)*(-1683. + 3405.*kappa + 21.*(-11. + 5.*kappa)*pow(a,4) - 2310.*pow(kappa,2) + 610.*pow(kappa,3) -
           12.*pow(a,2)*(-121. + 155.*kappa - 55.*pow(kappa,2) + 5.*pow(kappa,3)) - 55.*pow(kappa,4) + pow(kappa,5))*pow(m,3) +
        pow(a,2)*pow(m,2)*(-1683. + Complex(0.,81398.)*omega + 10.*omega*(Complex(0.,-137.) + (-1143. + 425.*kappa)*omega)*pow(a,6) + kappa*(3405. - Complex(0.,183642.)*omega - 1.01864e6*pow(omega,2)) +
           pow(kappa,5)*(1. - Complex(0.,130.)*omega - 1680.*pow(omega,2)) + 415680.*pow(omega,2) + 210.*pow(kappa,2)*(-11. + Complex(0.,678.)*omega + 4192.*pow(omega,2)) +
           5.*pow(kappa,4)*(-11. + Complex(0.,998.)*omega + 8960.*pow(omega,2)) - 10.*pow(kappa,3)*(-61. + Complex(0.,4442.)*omega + 31920.*pow(omega,2)) +
           2.*pow(a,2)*(726. - Complex(0.,50865.)*omega - 5.*omega*(Complex(0.,49.) + 795.*omega)*pow(kappa,4) - 255075.*pow(omega,2) + 45.*pow(kappa,5)*pow(omega,2) +
              10.*pow(kappa,3)*(-3. + Complex(0.,640.)*omega + 5865.*pow(omega,2)) - 30.*pow(kappa,2)*(-11. + Complex(0.,1319.)*omega + 9005.*pow(omega,2)) +
              15.*kappa*(-62. + Complex(0.,5400.)*omega + 30695.*pow(omega,2))) + pow(a,4)*
            (240.*omega*(Complex(0.,19.) + 225.*omega)*pow(kappa,2) - 3800.*pow(kappa,3)*pow(omega,2) - 15.*kappa*(-7. + Complex(0.,1766.)*omega + 12200.*pow(omega,2)) +
              3.*(-77. + Complex(0.,10370.)*omega + 54880.*pow(omega,2)))) + 2.*(Complex(0.,939.) - 5.*(-2029. + 511.*kappa)*omega)*pow(a,10)*pow(omega,3) +
        omega*pow(a,8)*(Complex(0.,300.) + (-16121. + 6079.*kappa)*omega - Complex(0.,2.)*(49277. - 31131.*kappa + 3552.*pow(kappa,2))*pow(omega,2) +
           40.*(-11517. + 9378.*kappa - 1821.*pow(kappa,2) + 76.*pow(kappa,3))*pow(omega,3)) +
        a*m*(Complex(0.,9923.) + 29612.*omega - Complex(0.,753544.)*pow(omega,2) + 10.*(Complex(0.,14.) + (473. - 235.*kappa)*omega)*pow(a,8)*pow(omega,2) +
           pow(kappa,4)*(Complex(0.,535.) + 5980.*omega - Complex(0.,66600.)*pow(omega,2) - 460992.*pow(omega,3)) - 2.11296e6*pow(omega,3) +
           pow(kappa,5)*(Complex(0.,-13.) - 348.*omega + Complex(0.,2520.)*pow(omega,2) + 20160.*pow(omega,3)) -
           6.*pow(kappa,2)*(Complex(0.,-2745.) - 12820.*omega + Complex(0.,250488.)*pow(omega,2) + 964928.*pow(omega,3)) +
           pow(kappa,3)*(Complex(0.,-4970.) - 32920.*omega + Complex(0.,511472.)*pow(omega,2) + 2.556288e6*pow(omega,3)) +
           kappa*(Complex(0.,-21849.) - 79244.*omega + Complex(0.,1.809336e6)*pow(omega,2) + 5.787072e6*pow(omega,3)) -
           2.*pow(a,6)*(Complex(0.,79.) + (-1403. + 755.*kappa)*omega - Complex(0.,2.)*(11221. - 4983.*kappa + 156.*pow(kappa,2))*pow(omega,2) +
              80.*(-666. + 471.*kappa - 90.*pow(kappa,2) + 5.*pow(kappa,3))*pow(omega,3)) +
           pow(a,4)*(-2.*(Complex(0.,83.) + 4170.*omega)*pow(kappa,4)*pow(omega,2) + 4.*omega*pow(kappa,3)*(134. + Complex(0.,3637.)*omega + 39960.*pow(omega,2)) -
              6.*pow(kappa,2)*(Complex(0.,-85.) + 842.*omega + Complex(0.,28360.)*pow(omega,2) + 146620.*pow(omega,3)) -
              3.*(Complex(0.,-1223.) + 2256.*omega + Complex(0.,152718.)*pow(omega,2) + 405452.*pow(omega,3)) +
              3.*kappa*(Complex(0.,-1015.) + 3796.*omega + Complex(0.,178404.)*pow(omega,2) + 604480.*pow(omega,3))) -
           2.*pow(a,2)*(Complex(0.,6108.) + 9027.*omega - Complex(0.,555968.)*pow(omega,2) + 5.*omega*pow(kappa,5)*(-1. + Complex(0.,28.)*omega + 336.*pow(omega,2)) +
              pow(kappa,4)*(Complex(0.,26.) + 215.*omega - Complex(0.,8960.)*pow(omega,2) - 87200.*pow(omega,3)) - 1.516448e6*pow(omega,3) +
              2.*pow(kappa,3)*(Complex(0.,-353.) - 1165.*omega + Complex(0.,63420.)*pow(omega,2) + 379280.*pow(omega,3)) -
              6.*pow(kappa,2)*(Complex(0.,-751.) - 1625.*omega + Complex(0.,98000.)*pow(omega,2) + 405216.*pow(omega,3)) +
              kappa*(Complex(0.,-9486.) - 16305.*omega + Complex(0.,1.00758e6)*pow(omega,2) + 3.245776e6*pow(omega,3)))) +
        2.*(-3191. - Complex(0.,18890.)*omega - 174576.*pow(omega,2) + 9.*pow(kappa,5)*(1. + Complex(0.,14.)*omega + 112.*pow(omega,2)) + Complex(0.,115968.)*pow(omega,3) +
           pow(kappa,3)*(2530. + Complex(0.,22828.)*omega + 99936.*pow(omega,2) - Complex(0.,1536.)*pow(omega,3) - 350208.*pow(omega,4)) - 178176.*pow(omega,4) +
           6.*pow(kappa,2)*(-1205. - Complex(0.,8998.)*omega - 52368.*pow(omega,2) + Complex(0.,15360.)*pow(omega,3) + 5120.*pow(omega,4)) +
           3.*pow(kappa,4)*(-105. - Complex(0.,1190.)*omega - 4624.*pow(omega,2) - Complex(0.,1792.)*pow(omega,3) + 49152.*pow(omega,4)) +
           3.*kappa*(2743. + Complex(0.,17810.)*omega + 133904.*pow(omega,2) - Complex(0.,67072.)*pow(omega,3) + 116736.*pow(omega,4))) +
        pow(a,6)*(883. - Complex(0.,4307.)*omega + 167304.*pow(omega,2) - 4.*pow(kappa,3)*pow(omega,2)*(710. + Complex(0.,8697.)*omega + 46810.*pow(omega,2)) + Complex(0.,658464.)*pow(omega,3) +
           4.*(Complex(0.,209.) + 1930.*omega)*pow(kappa,4)*pow(omega,3) + 12.*omega*pow(kappa,2)*(Complex(0.,-17.) + 3551.*omega + Complex(0.,26265.)*pow(omega,2) + 117460.*pow(omega,3)) +
           2.522568e6*pow(omega,4) - 100.*pow(kappa,5)*pow(omega,4) - 3.*kappa*(135. - Complex(0.,819.)*omega + 54372.*pow(omega,2) + Complex(0.,287044.)*pow(omega,3) + 1.16086e6*pow(omega,4))) +
        2.*pow(a,4)*(-3234. - Complex(0.,897.)*omega - 270690.*pow(omega,2) + 5.*pow(kappa,5)*pow(omega,2)*(7. + Complex(0.,41.)*omega + 504.*pow(omega,2)) - Complex(0.,613443.)*pow(omega,3) -
           1.*omega*pow(kappa,4)*(Complex(0.,34.) + 1850.*omega + Complex(0.,14215.)*pow(omega,2) + 49120.*pow(omega,3)) - 2.445248e6*pow(omega,4) +
           pow(kappa,3)*(113. + Complex(0.,849.)*omega + 29630.*pow(omega,2) + Complex(0.,184090.)*pow(omega,3) + 609040.*pow(omega,4)) -
           3.*pow(kappa,2)*(457. + Complex(0.,1323.)*omega + 58980.*pow(omega,2) + Complex(0.,264610.)*pow(omega,3) + 935648.*pow(omega,4)) +
           kappa*(4074. + Complex(0.,4779.)*omega + 391095.*pow(omega,2) + Complex(0.,1.248873e6)*pow(omega,3) + 4.666328e6*pow(omega,4))) +
        pow(a,2)*(12009. + Complex(0.,38991.)*omega + 724772.*pow(omega,2) + Complex(0.,494728.)*pow(omega,3) + 3.29568e6*pow(omega,4) -
           1.*pow(kappa,5)*(3. + Complex(0.,41.)*omega + 1044.*pow(omega,2) + Complex(0.,2520.)*pow(omega,3) + 20160.*pow(omega,4)) +
           pow(kappa,4)*(237. + Complex(0.,2435.)*omega + 20660.*pow(omega,2) + Complex(0.,85800.)*pow(omega,3) + 95424.*pow(omega,4)) -
           2.*pow(kappa,3)*(1671. + Complex(0.,12025.)*omega + 100580.*pow(omega,2) + Complex(0.,311800.)*pow(omega,3) + 676032.*pow(omega,4)) +
           6.*pow(kappa,2)*(2447. + Complex(0.,12805.)*omega + 138780.*pow(omega,2) + Complex(0.,275576.)*pow(omega,3) + 937280.*pow(omega,4)) -
           1.*kappa*(23487. + Complex(0.,94197.)*omega + 1.355428e6*pow(omega,2) + Complex(0.,1.648056e6)*pow(omega,3) + 7.691712e6*pow(omega,4)))))*pow(omega,8);
}
