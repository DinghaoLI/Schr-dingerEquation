/**
 * @file Solution1DHO.cpp
 *
 * Calcul de l'équation schrödinger d'oscillateur harmonique à une dimension
 */
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include "HermitePolynomials.h"
#include "Solution1DHO.h"
#include "Utils.h"

/** Constructeur
 *
 * valeur m,w vaut 1 par défaut
 */

Solution1DHO::Solution1DHO()
{
    m=1;
    w=1;
    h=1;
}

/** Constructeur avec des valeurs données
 *
 * @param masse
 * @param pulsation
 */

Solution1DHO::Solution1DHO(double masse, double pulsation)
{
    m=masse;
    w=pulsation;
}

/** Initialisation des constantes
 *
 * Cette fonction permet d'initialiser les constantes
 *
 * @param m masse
 * @param w pulsation de l'onde
 *
 * @return rien
 */
void Solution1DHO::setConstData(double m, double w)
{
    this->w = w;
    this->m = m;
}

/** Calcul de la fonction d'onde Psi
*
* Cette fonction a pour but de calculer la fonction Psi solution de l'équation de Schrödinger 1D-HO.
*
* @param n
* @param Z
*
* @return le résultat de la solution de 1D-HO Pn(z)
*/
arma::mat Solution1DHO::solution(int n, arma::mat Z)
{
    double f1;
    double f2;
    arma::mat f3;
    arma::mat res;

    f1 = (1/std::sqrt(std::exp(n * std::log(2)) * Utils::fac(n)));
    f2 = std::sqrt(std::sqrt(m * w / M_PI / h));
    f3 =  arma::exp(-(m * w * (Z % Z) / 2 / h));

    res = f1 * f2 * (f3 % HermitePolynomials::calculateHermite(n, std::sqrt(m * w / h) * Z));

    return res;
}

/** Stockage des valeurs de la fonction d'onde
*
* Cette fonction a pour but de stocker les valeurs de la fonction d'onde dans un fichier plot.txt
*
* @param n
* @param Z
*
* @return rien
*/
void Solution1DHO::solutionToFile(int n, arma::mat Z)
{
    arma::mat res = Solution1DHO::solution(n,Z);

    std::ofstream out("plot.txt");
    if (out.is_open())
    {
        out << Z;
        out << res << std::endl;
        out.close();
    }
    std::cout << "La solution de l'équation est stockée dans plot.txt" << std::endl;

}

/** Tracer le graphe de la fonction d'onde Psi
*
* Cette fonction a pour but de tracer le graph de la fonction d'onde Psi
* avec python par appel système.
*
* @return rien
*/
void Solution1DHO::drawSolution()
{
    std::cout << "Tracer du graphe avec les données de plot.txt! " << std::endl;
    system("python plot.py");
}

/**
 * Calcul de l'énergie
 *
 * Cette fonction a pout but de calculer l'énergie en connaissant la fonction d'onde pour vérifier qu'elle est constante.
 * # Equation de Schrödinger
 * \f[\left(\frac{\hat{p}_{(z)}^2}{2m} + \frac{1}{2}m\omega^2\hat{z}^2\right)\psi_n = E_n\psi_n.\f]
 *
 * @param n mode de l'énergie
 * @param a limite supérieur de Z
 * @param b limite inférieur de Z
 * @param N nombre de points z à caculer
 *
 *
 *
 * @retval E un vecteur contenant les valeurs des énergies pour N points.
 */
arma::mat Solution1DHO::energy(int n, double a, double b, int N)
{
    double pas = (b-a)/N;
    arma::mat Z = arma::linspace(a,b+2*pas,N+2);
    arma::mat Psi = Solution1DHO::solution(n,Z);
    arma::mat dPsi = Utils::derivative(Psi.rows(0,N), Psi.rows(1,N+1),Z.rows(0,N), Z.rows(1,N+1));
    arma::mat d2Psi = Utils::derivative(dPsi.rows(0,N-1), dPsi.rows(1,N), Z.rows(0,N-1), Z.rows(1,N));

    arma::mat ZN = Z.rows(0,N-1);
    arma::mat PsiN = Psi.rows(0,N-1);
    arma::mat E;
    E = -1 * h * h * d2Psi / (2 * m * PsiN) + m * w * w * ZN % ZN / 2;
    return E;
}
