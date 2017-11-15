#ifndef DEF_ORTH
#define DEF_ORTH
/**
 * @file Orthonoamality.h
 *
 * Interface de la classe Energy
 */
#include <armadillo>

/**
 * @class Orthonormality
 *
 * # Orthogonalité de la fonction d'onde
 * \f[\forall (m,n), \int \psi^*_m(z)\psi_n(z) dz = \delta_{mn}.\f]
 *
 * # Quadrature de Guass-Hermite
 * \f[\int_{a}^b \omega(x)g(x) dx \simeq \sum_{i=0}^{n-1}w^\omega_ig(x^\omega_i)\f]
 * | \f$a\f$       | \f$b\f$       | \f$\omega(x)\f$                                    | Quadrature      | Associated polynomial |
 * | :-----------: | :-----------: | :------------------------------------------------: | :-------------: | :-------------------: |
 * | \f$−\infty\f$ | \f$+\infty\f$ | \f$e^{-x^2}\f$                                     | Gauss–Hermite   | Hermite               |
 *
 */
class Orthonormality
{
public:
    static arma::mat gaussHermiteG(double n, double m, arma::mat X);
    static long double quadrature(double n, double m);

};

#endif