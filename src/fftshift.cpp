//----------------------------------------------------------------------------//
//文件名    ：  fft_shitf.cpp
//CopyRight(c) 2020-
//创建人    :   盛燊(Sheng Shen)
//单位      :   同济大学 海洋与地球科学学院 波现象与智能反演成像研究组(WPI)
//日期      :   2020.06
//描述      :   常用的一些子函数集合系列:   二维傅里叶变换的时移
//----------------------------------------------------------------------------//



#include "../inc/ss_sub_tools_func.h"



void fft_shift(arma::Mat<fcomplex> & data_in )
{
    //---------------------------针对二维FFT，将0点移到中心-----------------------//
    // Data_in ：复数矩阵,引用传参
    // nx  : 慢维
    // nt  : 快维
    int nx = data_in.n_cols;
    int nt = data_in.n_rows;


    arma::Mat<fcomplex> sei_fk_temp_1(size(data_in));
    arma::Mat<fcomplex> sei_fk_temp_2(size(data_in));
    //arma::Mat<fcomplex> & sei_fk_1 = data_in;

    /*Lfet & RIght*/
    for(int ix=0; ix<nx/2; ix++)
    {
        for(int iz=0; iz<nt; iz++)
        {
            sei_fk_temp_1(iz,ix) = data_in(iz,ix+nx/2);
        }
    }
    for(int ix=nx/2; ix<nx; ix++)
    {
        for(int iz=0; iz<nt; iz++)
        {
            sei_fk_temp_1(iz,ix) = data_in(iz,ix-nx/2);
        }
    }
    /*Up & Down*/
    for(int iz=0; iz<nt/2; iz++)
    {
        for(int ix=0; ix<nx; ix++)
        {
            sei_fk_temp_2(iz,ix) = sei_fk_temp_1(iz+nt/2,ix);
        }
    }
    for(int iz=nt/2; iz<nt; iz++)
    {
        for(int ix=0; ix<nx; ix++)
        {
            sei_fk_temp_2(iz,ix) = sei_fk_temp_1(iz-nt/2,ix);
        }
    }
    /*UP TO DOWM*/

    for(int iz=0; iz<nt; iz++)
    {
        for(int ix=0; ix<nx; ix++)
        {
            //data_in(iz,ix) = sei_fk_temp_2(nt-iz-1,ix);
            data_in(iz,ix) = sei_fk_temp_2(iz,ix);
        }
    }


}
