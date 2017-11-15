#ifndef DEF_HERMITE
#define DEF_HERMITE
/**
 * @file HermitePolynomials.h
 *
 * Interface de la classe HermitePolynomials
 */
#include <armadillo>

/**
 * @class HermitePolynomials
 *
 * Permet de calculer le polynôme d'Hermite
 *
 * # Definition
 * \f[\forall n\ge 0, H_n(z)\equiv (-1)^n e^{z^2} \frac{d^n}{dz^n}\left( e^{-z^2} \right).\f]
 * # Relation de récurrence
 * \f[H_0(z) = 1\f]
 * \f[H_1(z) = 2z\f]
 * \f[\forall n\ge 1, H_{n+1}(z) = 2zH_n(z)-2nH_{n-1}(z).\f]
 */
class HermitePolynomials
{
public:
    static arma::mat calculateHermite(int n, arma::mat Z);
    static double calculateHermite(int n, double z);
};

#endif