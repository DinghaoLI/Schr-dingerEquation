#ifndef DEF_SOLUTION1DHO
#define DEF_SOLUTION1DHO
/**
 * @file Solution1DHO.h
 *
 * Interface de la classe Solution1DHO
 */
#include <armadillo>

/**
 * @class Solution1DHO
 *
 * # Equation de Schrödinger
 * \f[\left(\frac{\hat{p}_{(z)}^2}{2m} + \frac{1}{2}m\omega^2\hat{z}^2\right)\psi_n = E_n\psi_n.\f]
 *
 * # Solution forme de la fonction d'onde
 * \f[\psi_n(z) = \frac{1}{\sqrt{2^n n!}}\left(\frac{m\omega}{\pi\hbar}\right)^{1/4}e^{-\frac{m\omega z^2}{2\hbar}}H_n\left(\sqrt{\frac{m\omega}{\hbar}} . z\right).\f]
 *
 */
class Solution1DHO
{
private:
    double h;
    double m,w;
public:
    Solution1DHO();
    Solution1DHO(double m, double w);
    void setConstData(double m, double w);
    arma::mat solution(int n, arma::mat Z);
    double solution(int n, double z);
    void solutionToFile(int n, arma::mat z);
    arma::mat energy(int n, double a, double b, int N);
    void drawSolution();
};

#endif