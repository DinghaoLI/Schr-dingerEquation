#ifndef DEF_UTILS
#define DEF_UTILS
/**
 * @file Utils.h
 *
 * Interface de la classe Utils
 *
 */
#include <armadillo>

/**
 * @class Utils
 *
 * #les fonctions de calculation mathématique
 */
class Utils
{
public:
    static arma::mat derivative(arma::mat F1, arma::mat F2, arma::mat Z1, arma::mat Z2);
    static int fac(int n);

};

#endif