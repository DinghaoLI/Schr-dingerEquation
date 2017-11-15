/**
 * @file Utils.cpp
 *
 * Fichier d'implémentation pour la calculation mathématique
 */

#include <iostream>
#include <armadillo>
#include "Solution1DHO.h"
#include "Utils.h"


/** Calcul de la factorielle de n
*
* Cette fonction a pour but de calculer la factorielle de n par récurrence.
*
* @param n entier
*
* @return le résultat de la factorielle de n
*/
int Utils::fac(int n)
{
    int f;

    if (n==0 || n==1)
        f=1;
    else
        f=fac(n-1)*n;

    return f;
}

/**
 * Calcul de dérivée
 *
 * Cette fonction a pour but d'approximer la dérivée du premier ordre
 *
 * @param F1 vecteur contenant les n valeurs d'une fonction
 * @param F2 vecteur contenant les n valeurs d'une fonction décalée d'une valeur
 * @param Z1 vecteur contenant les n points où sont évaluée la fonction
 * @param Z2 vecteur contenant les n points où sont évaluée la fonction décalée d'une valeur
 *
 * @return res vecteur contenant les valeurs de la fonction dérivée
 */
arma::mat Utils::derivative(arma::mat F1, arma::mat F2, arma::mat Z1, arma::mat Z2)
{
    arma::mat res = (F2-F1)/(Z2-Z1);
    return res;
}

