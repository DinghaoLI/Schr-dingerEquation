// MyTest.h
#include <cxxtest/TestSuite.h>
#include <armadillo>
#include "HermitePolynomials.h"
#include "Utils.h"
#include "Orthonormality.h"
#include "Solution1DHO.h"

class MyTestSuite1 : public CxxTest::TestSuite
{
public:
    void testCalculateHermite_0(void)
    {
        TS_TRACE("testCalculateHermite_0");
        arma::mat A = {{1,2,3}};
        arma::mat B = HermitePolynomials::calculateHermite(0,A);
        arma::mat C = {{1,1,1}} ;
        double diff = arma::norm(B-C,2);
        TS_ASSERT_EQUALS(diff,0);
    }

    void testCalculateHermite_1(void)
    {
        TS_TRACE("testCalculateHermite_1");
        arma::mat A = {{1,2,3}};
        arma::mat B = HermitePolynomials::calculateHermite(1,A);
        arma::mat C = {{2,4,6}};
        double diff = arma::norm(B-C,2);
        TS_ASSERT_EQUALS(diff,0);
    }

    void testCalculateFac(void)
    {
        TS_TRACE("CalculateFac");
        double res = Utils::fac(4);
        TS_ASSERT_EQUALS(res, 24);
    }

    void testOrthonormality(void)
    {
        TS_TRACE("CheckOrthonormality");
        double res1 = Orthonormality::quadrature(2, 2);
        double res2 = Orthonormality::quadrature(1, 1);
        double res3 = Orthonormality::quadrature(3, 3);
        double res4 = Orthonormality::quadrature(1, 3);
        double res5 = Orthonormality::quadrature(2, 1);
        double res6 = Orthonormality::quadrature(1, 4);

        TS_ASSERT_DELTA(res1, 1, 0.00000001);
        TS_ASSERT_DELTA(res2, 1, 0.00000001);
        TS_ASSERT_DELTA(res3, 1, 0.00000001);
        TS_ASSERT_DELTA(res4, 0, 0.00000001);
        TS_ASSERT_DELTA(res5, 0, 0.00000001);
        TS_ASSERT_DELTA(res6, 0, 0.00000001);
    }

    void testEnergy(void)
    {
        TS_TRACE("CheckEnergy");
        arma::mat res;
        Solution1DHO ex;
        ex.setConstData(1,1);
        res = ex.energy(0, -1, 1, 1000000);
        arma::mat n1 = {{0.5, 0.5, 0.5}};
        double diff = arma::norm(res.rows(4, 6).t()-n1, 2);
        TS_ASSERT_DELTA(diff, 0, 0.01);
        diff = arma::norm(res.rows(1004, 1006).t()-n1, 2);
        TS_ASSERT_DELTA(diff, 0, 0.01);
        diff = arma::norm(res.rows(10004, 10006).t()-n1, 2);
        TS_ASSERT_DELTA(diff, 0, 0.01);
        diff = arma::norm(res.rows(100004, 100006).t()-n1, 2);
        TS_ASSERT_DELTA(diff, 0, 0.01);

        res = ex.energy(1, -1, 1, 1000000);
        n1 = {{1.5, 1.5, 1.5}};
        diff = arma::norm(res.rows(4, 6).t()-n1, 2);
        TS_ASSERT_DELTA(diff, 0, 0.01);
        diff = arma::norm(res.rows(1004, 1006).t()-n1, 2);
        TS_ASSERT_DELTA(diff, 0, 0.01);
        diff = arma::norm(res.rows(10004, 10006).t()-n1, 2);
        TS_ASSERT_DELTA(diff, 0, 0.01);
        diff = arma::norm(res.rows(100004, 100006).t()-n1, 2);
        TS_ASSERT_DELTA(diff, 0, 0.01);

        res = ex.energy(2, -1, 1, 1000000);
        n1 = {{2.5, 2.5, 2.5}};
        diff = arma::norm(res.rows(4, 6).t()-n1, 2);
        TS_ASSERT_DELTA(diff, 0, 0.01);
        diff = arma::norm(res.rows(1004, 1006).t()-n1, 2);
        TS_ASSERT_DELTA(diff, 0, 0.01);
        diff = arma::norm(res.rows(10004, 10006).t()-n1, 2);
        TS_ASSERT_DELTA(diff, 0, 0.01);
        diff = arma::norm(res.rows(100004, 100006).t()-n1, 2);
        TS_ASSERT_DELTA(diff, 0, 0.01);

    }

};