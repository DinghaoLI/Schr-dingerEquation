/**
 * @file HermitePolynomials.cpp
 *
 * Fichier d'implémentation pour le calcul du polynôme d'Hermite
 */

#include <iostream>
#include <armadillo>
#include "HermitePolynomials.h"

/**
 * Calcul du polynôme Hermite par recurrence
 *
 * Cette fonction a pour but d'Calculation du polynôme Hermite par recurrence
 *
 * @param n
 * @param z
 *
 * @return la valeur de la fonction Hn(z)
 */
double HermitePolynomials::calculateHermite(int n, double z)
{
    if (n == 0)
        return 1;
    else if (n == 1)
        return (2*z);
    else
        return (2 * z * calculateHermite(n-1, z) - 2 * (n-1) * calculateHermite(n-2, z));
}



/**
 * Calcul du polynôme d'Hermite par récurrence
 *
 * Cette fonction a pour but de calculer le polynôme d'Hermite par récurrence avec des vecteurs
 *
 * @param n
 * @param Z un vecteur de valeurs
 *
 * @return un vecteur contenant les valeurs de la fonction Hermite Hn(z)
 */
arma::mat HermitePolynomials::calculateHermite(int n, arma::mat Z)
{
    if (n == 0)
    {
        arma::mat A = Z.ones(size(Z));
        return  A;
    }
    else
    {
        if (n == 1)
        {
            return Z.for_each([](arma::mat::elem_type& val)
            {
                val = 2 * val;
            });
        }
        else
        {
            return (2 * (Z % calculateHermite(n - 1, Z))) - (2 * (n-1) * calculateHermite(n - 2, Z));
        }
    }
}



