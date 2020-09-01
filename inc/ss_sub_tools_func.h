#ifndef SS_SUB_TOOLS_FUNC_H
#define SS_SUB_TOOLS_FUNC_H

#include "TJU_SS_HEAD.h"
#include "mycomplex.h"



//------------------------------傅里叶变换------------------------------------//
extern "C"
{
void FFTCC1(int isign, fcomplex *z, int n);

void FFTCC2(int isign1, int isign2, fcomplex *zz, int n1, int n2,fcomplex *work);

void FFTCC3(int isign1   , int isign2 , int isign3 ,
            int n1       , int n2     , int n3     ,
            fcomplex *zz , int LevelVerbose);

void fft1d_(fcomplex *z, int *nn, int *ii);
}


//--------------------针对二维傅里叶变换，将0点移到中心-----------------------//
void fft_shift(arma::Mat<fcomplex> & data_in );


//--------------------------按照信噪比添加噪声--------------------------------//
template<typename T>

void addGaussianNoiseBySNR(arma::Mat<T> & data, float SNR_db);


//------------------矩阵求解(稀疏群LASSO、LASSO、共轭梯度)--------------------//
void group_lasso_func(fcomplex *In_A_raw, fcomplex *Out_x_raw, fcomplex *In_d_raw,
        int *Nhx,  int *Np, float *errr);

//----------------------------------------------------------------------------//













#endif
