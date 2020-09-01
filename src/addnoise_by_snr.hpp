#pragma once

#include "../inc/TJU_SS_HEAD.h"

/*============================================================================*/
template <typename T>
float Norm2Sq(const arma::Mat<T>& A)
{
    return powf(arma::norm(A, "fro"), 2);
}
/*============================================================================*/
template <typename T>
float Norm2Sq(const arma::Cube<T>& A)
{
    int nsli = A.n_slices;
    float ans = 0;
    for (int i = 0; i < nsli; i++)
    {
        ans += Norm2Sq(A.slice(i));
    }
    return ans;
}

/*============================================================================*/
/**
 * Assume SNR_db = 10*lg(power_s / power_n)
 */
template<typename T>
void addGaussianNoiseBySNR(arma::Mat<T> & data, float SNR_db)
{
    arma::Mat<T>  noise(size(data),arma::fill::randn); //Gauss Noise
    //arma::Mat<T>  noise(size(data),fill::randu); //Non-gauss Noise
    float   power_n      = Norm2Sq(noise);
    float   power_s      = Norm2Sq(data );
    float   SNR_no_db    = pow( 10, (SNR_db / 10.0));
    float   power_n_need = power_s / SNR_no_db;
    float   factor       = sqrt(power_n_need / power_n);

    noise *= factor;
    data  += noise;

    return;
}
