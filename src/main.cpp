/**
 *@file main.cpp
 */

#include <fstream>
#include <armadillo>
#include "HermitePolynomials.h"
#include "Solution1DHO.h"
#include "Orthonormality.h"


/**
 * @mainpage Calcul de l'équation schrödinger d'oscillateur harmonique à une dimension
 *
 * # Equation de Schrödinger
 * \f[\left(\frac{\hat{p}_{(z)}^2}{2m} + \frac{1}{2}m\omega^2\hat{z}^2\right)\psi_n = E_n\psi_n.\f]
 *
 * # Solution forme de la fonction d'onde
 * \f[\psi_n(z) = \frac{1}{\sqrt{2^n n!}}\left(\frac{m\omega}{\pi\hbar}\right)^{1/4}e^{-\frac{m\omega z^2}{2\hbar}}H_n\left(\sqrt{\frac{m\omega}{\hbar}} . z\right).\f]
 */

/**
 *La fonction principale
 *
 *Cette fonction a pour but de tester.
 *
 *@return la fonction toujours retourne 0.
 */

int main()
{
    int i,j;
    int a=-4, b=4, N=100;
    int n=0;

    Solution1DHO ex = Solution1DHO();

    arma::mat A = arma::linspace(-4,4,100);
    std::cout << "Tracer la fonction d\'onde pour n = " << "[" << a << " ," << b << "] pour N =" << N << std::endl<< std::endl;
    ex.solutionToFile(0, A.t());
    ex.drawSolution();


    n = 4;
    std::cout.precision(50);
    for (i = 1; i<n; i++)
    {
        for (j = 1; j<n; j++)
        {
            std::cout << "Test Orthonormality n = " << i << " m = " << j << " quadrature: " << Orthonormality::quadrature(i, j) << std::endl;
        }
    }
    std::cout << std::endl << std::endl;

    N = 1000000;
    n = 2;
    for (i = 0; i<=n; i++)
    {

        std::cout<< "La plage de valeur z: (-2,2) n: "<< i << std::endl;
        std::cout<< "Nombre d'échantillons: N = 1000000" << std::endl;
        std::cout<< "E" << i << " (1-5): ";
        arma::mat res;
        res = ex.energy(i, -2, 2, N);
        std::cout<< res.rows(0,4) << std::endl;

        std::cout<<  "E" << i << " (500100-500104): " ;
        std::cout<< res.rows(500099,500103) << std::endl;

        std::cout<< "E" << i << " (999995-1000000): " ;
        std::cout<< res.rows(999994,999999) << std::endl;

    }



    return 0;
}


