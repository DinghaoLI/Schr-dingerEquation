/**
 * @file Orthonormality.cpp
 *
 * Fichier d'implémentation pour la vérification de la propriété orthonormal de la fonction d'onde
 */

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include "Orthonormality.h"
#include "Solution1DHO.h"
#include "HermitePolynomials.h"
#include "Utils.h"

/**
 * Calcul du résultat de g(x) de la méthode Gauss-Hermite
 *
 * Cette fonction a pour but de calculer la partie g(x) de la quadrature pour orthogonalité.
 *
 * @param n
 * @param m
 * @param X un vecteur des points Gauss-Hermite
 *
 * @return un vecteur contenant les valeurs de la fonction g(x)
 */
arma::mat Orthonormality::gaussHermiteG(double n, double m, arma::mat X)
{
    double f1;
    double f2;

    arma::mat Hnm;
    arma::mat res;
    f1 = 1 / std::sqrt(std::exp(n * std::log(2)) * Utils::fac(n)) / std::sqrt(std::exp(m * std::log(2)) * Utils::fac(m));
    f2 = std::sqrt(1 / M_PI);
    Hnm = HermitePolynomials::calculateHermite(n, X) % HermitePolynomials::calculateHermite(m, X);
    res = f1 * f2 * Hnm;

    return res;

}

/** Calcul d'orthogonalité de la solution 1D-HO
*
* Cette fonction a pour but de calculer l'orthogonalité de la solution 1D-HO par deux paramètres n et m.
*
* @param n
* @param m
*
*
* @return le résultat d'orthogonalité de la solution 1D-HO avec les paramètres n et m
*/
long double Orthonormality::quadrature(double n, double m)
{
    arma::mat gx;
    arma::mat gauss_point = {{
            -2.453407083009012499038365306336166239661e-1,
            2.453407083009012499038365306336166239661e-1,
            -7.374737285453943587056051442521042290772e-1,
            7.374737285453943587056051442521042290772e-1,
            1.234076215395323007885818346959410229585,
            -1.234076215395323007885818346959410229584,
            -1.738537712116586206780865662136406442958,
            1.738537712116586206780865662136406442953,
            2.254974002089275523082333344734565128082,
            -2.254974002089275523082333344734565128065,
            -2.788806058428130480525033756403185410695,
            2.788806058428130480525033756403185410655,
            3.347854567383216326914924522996463698566,
            -3.347854567383216326914924522996463698495,
            -3.94476404011562521037562880052441180715,
            3.944764040115625210375628800524411807067,
            4.603682449550744273077675248978347585171,
            -4.603682449550744273077675248978347585109,
            5.387480890011232862016900410681120753981,
            -5.387480890011232862016900410681120754003,
        }
    };

    arma::mat gauss_point_weight = {{
            4.622436696006100896503286398612081142142e-1,
            4.622436696006100896503286398612081142142e-1,
            2.866755053628341297196597062280879168236e-1,
            2.866755053628341297196597062280879168236e-1,
            1.090172060200233200137550335354255770852e-1,
            1.090172060200233200137550335354255770846e-1,
            2.481052088746361088216495255894039439922e-2,
            2.481052088746361088216495255894039440028e-2,
            3.24377334223786183218324713235370544232e-3,
            3.243773342237861832183247132353705443042e-3,
            2.283386360163539672571459179634955394906e-4,
            2.283386360163539672571459179634955393512e-4,
            7.802556478532063694145991999647569104495e-6,
            7.802556478532063694145991999647569095955e-6,
            1.086069370769281693999524563447163430255e-7,
            1.086069370769281693999524563447163432688e-7,
            4.399340992273180553628851455467928211995e-10,
            4.399340992273180553628851455467928212879e-10,
            2.229393645534151292522500616029095785758e-13,
            2.22939364553415129252250061602909578525e-13,
        }
    };

    gx = Orthonormality::gaussHermiteG(n, m, gauss_point);

    arma::mat res;
    res = gx * gauss_point_weight.t();
    long double resDouble = res(0, 0);
    return resDouble;
}