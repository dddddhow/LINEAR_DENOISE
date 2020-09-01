//----------------------------------------------------------------------------//
//程序目的： 1、线性去噪
//
//程序原理： 对原数据进行FK变换，由于FK域和t-x域斜率代表视速度，在FK域根据视速度
//           的概念压制线性噪声，然后反变换回t-x域。
//
//程序输入输出参数说明：
//          nxori：横向样点数（个）
//          ntori：时间样点数（个）
//          sei_in ：读入去噪前地震数据
//          sei_out：输出去噪后地震数据
//          cur：视速度，将小于视速度的FK成分按照比例压缩(m/s)
//          cur_1：视速度，将此视速度的数据在FK域进行压缩(m/s)
//          注意：其他参数为中间参数，会在程序中进行说明。
//
//Copyright：2020-
//          WPI TONGJI University
//Author  ：ShengShen
//Time    ：2020 07 22
//----------------------------------------------------------------------------//

#include "../inc/TJU_SS_HEAD.h"
#include "../inc/alloc.h"
#include "../inc/mycomplex.h"
#include "../inc/suhdr240.h"
#include "../inc/ss_sub_tools_func.h"
#include "./addnoise_by_snr.hpp"

using namespace std;
using namespace arma;

int main()
{
    //========================================================================//
    //                        生成地震剖面sei   FKDENOISE ANALYSIS            //
    //========================================================================//
    cout<<"--------------------------------------------------"<<endl;
    cout<<"START"<<endl;
    cout<<"--------------------------------------------------"<<endl;


    cout<<"--------------------------------------------------"<<endl;
    cout<<"Define Matrix & Request Memory"<<endl;
    //为了使频谱泄露的能量不影响去噪后结果，将nx_ori, nt_ori扩大一倍,即nx,nt
    int nx_ori=1441;
    int nt_ori=2000;

    int nx=nx_ori*2;
    int nt=nt_ori*2;

    arma::Mat<float> sei_in(nt,nx,fill::zeros);
    arma::Mat<float> sei_ou(nt,nx,fill::zeros);
    arma::Mat<cx_float> sei_fk(nt,nx,fill::zeros);

    arma::Mat<float> sei_in_ori(nt_ori,nx_ori,fill::zeros);
    arma::Mat<float> sei_ou_ori(nt_ori,nx_ori,fill::zeros);
    arma::Mat<cx_float> sei_fk_ori(nt_ori,nx_ori,fill::zeros);
    arma::Mat<float> sei_dif_ori(nt_ori,nx_ori,fill::zeros);
    arma::Mat<float> sei_fk_abs_ori(nt_ori,nx_ori,fill::zeros);

    cout<<"--------------------------------------------------"<<endl;


    //-----------------------------0:Load Data--------------------------------//
    cout<<"--------------------------------------------------"<<endl;
    cout<<"First Step : Load seismic data"<<endl;

    //sei_in_ori.load("../file/shot_vz_newsample_nt2000_nx2880.dat",raw_binary);
    //sei_in_ori.load("../xrw/seis_newsample_nt2000_nx1441.dat",raw_binary);
    sei_in_ori.load("../xrw_1/seis_newsample_nt2000_nx1441.dat",raw_binary);
    sei_in_ori.reshape(nt_ori,nx_ori);

    for(int ix=0; ix<nx_ori;ix++)
    {
        for(int it=0; it<nt_ori;it++)
        {
            sei_in(it,ix)=sei_in_ori(it,ix);
        }
    }

    cout<<"Load seismic data finish"<<endl;
    cout<<"--------------------------------------------------"<<endl;


    //-------------------------------1:FK-------------------------------------//
    cout<<"--------------------------------------------------"<<endl;
    cout<<"Second Step : FK"<<endl;

    {
        sei_fk.set_real(sei_in);
        arma::Col<cx_float> working(nx);
        FFTCC2(-1, 1, &sei_fk(0,0), nt, nx, &working(0));
        fft_shift(sei_fk);
    }

    for(int ix=0; ix<nx_ori;ix++)
    {
        for(int it=0; it<nt_ori;it++)
        {
            sei_fk_ori(it,ix)=sei_fk(it*2,ix*2);
        }
    }
    sei_fk_abs_ori=abs(sei_fk_ori);

    cout<<"FK finish"<<endl;
    cout<<"--------------------------------------------------"<<endl;


    //----------------------------------2:DENOISE-----------------------------//
    cout<<"--------------------------------------------------"<<endl;
    cout<<"Third Step : Denoise"<<endl;
    //参数说明：
    //cur：视速度，将小于此视速度的FK域成分按比例压缩
    //cur_1：视速度：将此视速度的FK成分按比例压缩
    //cur_v1：FK域三角窗左侧曲率
    //cur_v2：FK域三角窗右侧曲率

    float cur=1200;
    float cur_1=1000;


    float real_median_seis_fk=max(max(abs(real(sei_fk))));
    float imag_median_seis_fk=max(max(abs(imag(sei_fk))));
    cout<<"     "<<real_median_seis_fk<<endl;
    cout<<"     "<<imag_median_seis_fk<<endl;

    arma::Col<cx_float> cur_sum(nt,fill::zeros);
    float curv_1_1=(nt*1.0/2-cur_1)*1.0/(nx*1.0/2);

    for(int it=0; it<nt; it++)
    {
        for(int ix=0; ix<nx; ix++)
        {
            int i_it=int(ix*curv_1_1)+it;
            if(i_it<nt)
            {
                cur_sum(it)=cur_sum(it)+sei_fk(i_it,ix);
            }
        }
        float real_sum=real(cur_sum(it));
        float imag_sum=imag(cur_sum(it));

        //cout<<real_sum<<endl;
        //cout<<imag_sum<<endl;

        if(abs(real_sum) > real_median_seis_fk*10 &&
           abs(imag_sum) > imag_median_seis_fk*10)
        {
            for(int ix=0; ix<nx; ix++)
            {
                int i_it=int(ix*curv_1_1);
                if(i_it<nt)
                {
                    float rea =real( sei_fk(i_it,ix));
                    float ima =imag( sei_fk(i_it,ix));
                    sei_fk(i_it,ix).real(rea*0.0001);
                    sei_fk(i_it,ix).imag(ima*0.0001);
                }
            }

        }
    }

/*
    float curv_1=(nt*1.0/2-cur)*1.0/(nx*1.0/2);
    int  curv_d1=int((nt*1.0/2)-curv_1*(nx*1.0/2));

    float curv_2=(cur-nt*1.0/2)*1.0/(nx*1.0/2);
    int  curv_d2=int((nt*1.0/2)-curv_2*(nx*1.0/2));

    cout<<"     Here we Suppress data with apparent-velocity less than "<<cur<<endl;
    cout<<"     And here we Suppress data with apparent-velocity near "<<cur_1<<endl;
    cout<<"     cur_v1 = "<<curv_1 <<"   d1 = "<<curv_d1<<endl;
    cout<<"     cur_v2 = "<<curv_2 <<"   d2 = "<<curv_d2<<endl;


    for(int ix=0; ix<nx; ix++)
    {
        for(int it=0; it<nt; it++)
        {
            if((nt/2-it)*1.0/(nx/2-ix)<curv_1 &&
                (nt/2-it)*1.0/(nx/2-ix)>curv_2)
            {
                float rea =real( sei_fk(it,ix));
                float ima =imag( sei_fk(it,ix));
                float distance_1=abs(it-(ix*curv_1+curv_d1));
                float distance_2=abs(it-(ix*curv_2+curv_d2));
                //cout<<"dus 1 =" <<distance_1<<endl;
                //cout<<"dus 2 =" <<distance_2<<endl;
                float distance = min(distance_1, distance_2);
                //cout<<"dus  =" <<distance<<endl;
                //cout<<endl;
                //float ratio = (1-distance*1.0/cur) *0.1;
                float ratio = 1.0/(distance+1);
                //cout<<"Ratio is :"<<ratio<<endl;

                sei_fk(it,ix).real(rea*ratio);
                sei_fk(it,ix).imag(ima*ratio);
            }
        }
    }

*/
    for(int ix=0; ix<nx_ori;ix++)
    {
        for(int it=0; it<nt_ori;it++)
        {
            sei_fk_ori(it,ix)=sei_fk(it*2,ix*2);
        }
    }
    arma::Mat<float> sei_fk_abs_ori_denoise=abs(sei_fk_ori);

    cout<<"Denoise finish"<<endl;
    cout<<"--------------------------------------------------"<<endl;


    //-----------------------------3:IFK--------------------------------------//
    cout<<"--------------------------------------------------"<<endl;
    cout<<"Fouth Step : IFK"<<endl;

    {
        fft_shift(sei_fk);
        arma::Col<cx_float> working(nx);
        FFTCC2(1, -1, &sei_fk(0,0), nt, nx, &working(0));
        sei_ou=real(sei_fk);
    }

    for(int ix=0; ix<nx_ori;ix++)
    {
        for(int it=0; it<nt_ori;it++)
        {
            sei_ou_ori(it,ix)=sei_ou(it,ix);
        }
    }

    cout<<"IFK finish"<<endl;
    cout<<"--------------------------------------------------"<<endl;


    //--------------------------4:Cal Dif-------------------------------------//
    sei_dif_ori = sei_in_ori - sei_ou_ori;
    //sei_dif_ori.save("../file/shot_vz_newsample_dif.dat",raw_binary);


    //--------------------------5:Save Files----------------------------------//
    cout<<"--------------------------------------------------"<<endl;
    cout<<"Fifth Save Files"<<endl;

    string fn= "../xrw_1/";

    string fn_sei_fk_abs_ori=fn+"shot_vz_newsample_fk_abs.dat";
    sei_fk_abs_ori.save(fn_sei_fk_abs_ori,raw_binary);

    string fn_sei_fk_abs_ori_denoise=fn+"shot_vz_newsample_fk_denoise_abs.dat";
    sei_fk_abs_ori_denoise.save(fn_sei_fk_abs_ori_denoise,raw_binary);

    string fn_sei_ou_ori=fn+"shot_vz_newsample_ou.dat";
    sei_ou_ori.save(fn_sei_ou_ori,raw_binary);

    string fn_dif_ori=fn+"shot_vz_newsample_dif.dat";
    sei_dif_ori.save(fn_dif_ori,raw_binary);


    cout<<"Save FIles finish"<<endl;
    cout<<"--------------------------------------------------"<<endl;

    cout<<"--------------------------------------------------"<<endl;
    cout<<"FINISH"<<endl;
    cout<<"--------------------------------------------------"<<endl;



    return 0.0;


    /*



       int nt=1000;
       int nx=100;
       int nw=100; //length of wavelet
       double dt=0.001;
       double dx=5;
       double fre=25; // frequence of wavelet


       arma::Mat<float> ref(nt,nx); // reflection coefficient
       arma::Mat<float> sei(nt,nx); // seismic matrix
       Col<float> wavelet(nw);// wavelet vecter

    //========================================================================//
    //                        生成地震剖面sei 或者 加载地震数据               //
    //========================================================================//
    //reflection coefficition//
    ref.zeros();
    //ref.row(100).fill(1);
    float cur = 3;
    for(int i=0; i<nx; i++)
    {
    ref(cur*i,i)=1;
    if(i%6==0 && i>nx/5 && i<nx/2)
    //if(i%6==0 )
    {
    ref(cur*i+10,i)=1;
    }
    }


    //Non Linear TX relationship
    for(int i=0; i<nx; i++)
    {
    //ref(floor(sqrt(250000*dt*dt+i*i*dx*dx/1500/500)/dt),i)=1;
    }

    //Non Linear Amp
    //for(int i=0; i<nx; i++)
    //{
    //ref(2*i+500,i)=1;
    //if(i%6==0 && i>nx/3 && i<nx/2)
    //{
    //ref(2*i+500,i)=3;
    //}
    //}

    //Non Linear Amp
    //for(int i=0; i<nx; i++)
    //{
    //ref(2*i+500,i)=1;
    //if(i%6==0)
    //{
    //ref(2*i+500+10,i)=1;
    //}
    //}


    //rickey wavelet//
    wavelet.zeros();
    for (int t=0; t<nw; t++)
    {
    float tt=t*dt;
    wavelet(t)=(1-2*pow((PI*fre*(tt-1.0/(1.5*fre))),2))
     *exp(-pow((PI*fre*(tt-1.0/(1.5*fre))),2));
     }


    //convolution//
    arma::Mat<float> sei_temp(nt+nw-1,nx);
    sei_temp.zeros();
    for(int i=0; i<nx; i++)
    {
        sei_temp.col(i)=conv(ref.col(i),wavelet);
    }

    sei=sei_temp(span(0,nt-1),span(0,nx-1));

    //========================================================================//
    //                              FK Filter Analysis                        //
    //========================================================================//

    //--------------------------------1 : 2D FFT-----------------------------//

    int flag_ft=0;
    //  flag = 0  FK变换（频率-波数域）
    //  flag = 1  FFT2变换（Kx-Ky域）


    arma::Mat<cx_float> sei_fk(nt,nx,fill::zeros);
    if(flag_ft == 0)
    {
        sei_fk.set_real(sei);
        arma::Col<cx_float> working(nx);
        FFTCC2(-1, 1, &sei_fk(0,0), nt, nx, &working(0));
        fft_shift(sei_fk);
    }


    if(flag_ft == 1)
    {
        arma::Mat<cx_float> sei_fk_temp_1d_wx(nt,nx,fill::zeros);

        for(int ix=0; ix<nx; ix++)
        {
            sei_fk_temp_1d_wx.col(ix) = fft(sei.col(ix));
        }
        for(int iz=0; iz<nt; iz++)
        {
            sei_fk.row(iz) = fft(sei_fk_temp_1d_wx.row(iz));
        }
        fft_shift(sei_fk);
    }


    arma::Mat<float> sei_fk_imag  = imag(sei_fk);
    arma::Mat<float> sei_fk_real  = real(sei_fk);
    arma::Mat<float> sei_fk_amp  = abs(sei_fk);

    sei.save("../file/sei.dat",raw_binary);
    sei_fk.save("../file/sei_fk.dat",raw_binary);
    sei_fk_amp.save("../file/sei_fk_amp.dat",raw_binary);



    //--------------------------------2 : 2D IFFT-----------------------------//
    arma::Mat<float> sei_ift(nt,nx,fill::zeros);

    if(flag_ft ==0)
    {
        fft_shift(sei_fk);
        arma::Col<cx_float> working(nx);
        FFTCC2(1, -1, &sei_fk(0,0), nt, nx, &working(0));
        sei_ift=real(sei_fk);
    }



    if(flag_ft == 1)
    {
    }

    arma::Mat<float> sei_ift_amp = abs(sei_ift);
    sei_ift.save("../file/sei_ift.dat",raw_binary);
    sei_ift_amp.save("../file/sei_ift_amp.dat",raw_binary);

    return 0;

    //========================================================================//
    //                              Denoise分析                               //
    //========================================================================//
    //1 : Gauss noise analysis
    float snr = -1;
    addGaussianNoiseBySNR(sei,snr);
    //peppersalt noise

    arma::uvec peppersalt_x = randperm(nx*nt);
    arma::Col<float> peppersalt_1d(nx*nt,fill::zeros);
    for(int ix=0; ix<nx*nt/70;ix++)
    {
        peppersalt_1d(peppersalt_x(ix)) = 1;
    }
    for(int ix=nx*nt/70; ix<nx*nt/35;ix++)
    {
        peppersalt_1d(peppersalt_x(ix)) = -1;
    }
    arma::Mat<float> peppersalt = reshape (peppersalt_1d,nt,nx);

    peppersalt.save("../file/peppersalt.dat",raw_binary);

    sei=sei+peppersalt;

    //2 : NonLinear  t-x relationship


    //3 : change the Matrix ref


    //4 : NonLinearAmp


    //========================================================================//
    //               Aniostropic  Gauss Filter Denoise Analysis               //
    //========================================================================//

    //测方向


    //Anios Gauss kernel
    int nx_gauss_kernel = 5;
    float cita_x = 2;
    float cita_z = 1;
    float cita   = 0*1.0/180*PI;
    //float cita = atan(cur);

    arma::Mat<float> gauss_kernel(nx_gauss_kernel,nx_gauss_kernel,fill::zeros);
    for(int ix=0; ix<nx_gauss_kernel; ix++)
    {
        for(int iz=0; iz<nx_gauss_kernel; iz++)
        {
            int z_t = ix-(nx_gauss_kernel-1)/2;
            int x_t = iz-(nx_gauss_kernel-1)/2;
            gauss_kernel(ix,iz)=1.0/(2.0*PI*cita_x*cita_z)
                *exp(-1.0/2*
                        (
                         (pow((x_t*cos(cita)+z_t*sin(cita)),2)/cita_x/cita_x)
                         +
                         (pow((-x_t*sin(cita)+z_t*cos(cita)),2)/cita_z/cita_z)
                        )
                    );
        }
    }
    float max_g = max(max(gauss_kernel));
    gauss_kernel = gauss_kernel / max_g *1;
    //gauss_kernel.print();


    // Gauss Filter Denoise
    arma::Mat<float> sei_gauss_fiter(nt,nx,fill::zeros);
    sei_gauss_fiter = sei;
    for (int ix=0; ix<nx-nx_gauss_kernel; ix++)
    {

        for (int iz=0; iz<nt-nx_gauss_kernel; iz++)
        {
            sei_gauss_fiter(span(iz,iz+nx_gauss_kernel-1),span(ix,ix+nx_gauss_kernel-1))
                =
                sei(span(iz,iz+nx_gauss_kernel-1),span(ix,ix+nx_gauss_kernel-1))
                *
                gauss_kernel;
        }
    }
    sei_gauss_fiter.save("../file/sei_gauss.dat",raw_binary);
    arma::Mat <float> sei_gauss_dif =sei-sei_gauss_fiter;
    sei_gauss_dif.save("../file/sei_gauss_dif.dat",raw_binary);



    //  Medium Filter Denoise
    int nx_me_filter = 5;
    arma::Mat<float> sei_me_fiter(nt,nx,fill::zeros);
    //sei_me_fiter = sei;
    arma::Col<float> temp_me(nx_me_filter * nx_me_filter,fill::zeros);

    for(int ix=0; ix<nx-nx_me_filter; ix++)
    {
        for(int iz=0 ; iz<nt-nx_me_filter; iz++)
        {
            int ime=0;
            for(int ixx=ix; ixx<ix+nx_me_filter; ixx++)
            {
                for(int izz=iz; izz<iz+nx_me_filter;izz++)
                {
                    temp_me(ime)=sei(izz,ixx);
                    ime++;
                }
            }
            sei_me_fiter(iz,ix) = median(temp_me);

        }
    }


    sei_me_fiter.save("../file/sei_me.dat",raw_binary);
    arma::Mat <float> sei_me_dif =sei-sei_me_fiter;
    sei_me_dif.save("../file/sei_me_dif.dat",raw_binary);


    //========================================================================//
    //                              线性RADON变换                             //
    //========================================================================//

    //-----------------------define source s------------------------//
    const int nf = nt;
    const int np = 200;

    arma::Mat<cx_float> s(nf,nx);

    for(int i=0; i<nx; i++)
    {
        Col<cx_float> s_temp_1 = fft(sei.col(i));
        for(int j=0;j <nf; j++)
        {
            s(j,i)= s_temp_1(j);
        }
    }

    //----------------- define x as omega-P data -----------------//
    arma::Mat<cx_float> x(nf,np,fill::zeros);

    //-------------------define transform matrix--------------------//
    arma::Mat<cx_float> A(nx,np,fill::zeros);

    const float df = 1.0 / (nf*dt);
    const float x0=0;
    //const float z0=0;
    const float k_pmax = 1.0 / 1500.0;
    auto p_list = arma::linspace<fvec>(-k_pmax,+k_pmax,np);
    const fcomplex unit_imag = fcomplex(0.0f,1.0f);

    arma::Mat<cx_float> inv_radon_sei_fre(nf,nx,fill::zeros);

    for (int w=0; w<nf; w++)
    {
        const float freqHz = (w >= nf/2)? (w-nf)*df : w*df;
        //gen A mat
        for(int ix=0; ix<nx ;ix++)
        {
            for (int ip=0; ip<np; ip++)
            {
                float omega   = 2.0*PI* freqHz;
                float delta_x = (ix*dx)-x0;
                float p       = p_list[ip];
                float phase   = omega*p*delta_x;
                A(ix,ip)      = exp( unit_imag*phase);
            }
        }


        //---------------------solve x----------------------------//
        Col<fcomplex> x_temp(np,fill::zeros);
        Col<fcomplex> s_temp = s.row(w).st();
        float err=0.001;
        int np_temp = np;
        //-----------LASSO----------//
        group_lasso_func(&A(0,0), &x_temp(0), &s_temp(0), &nx, &np_temp, &err);

        x.row(w) = x_temp.st();
        Col<cx_float> inv_s_temp = A*x_temp;
        inv_radon_sei_fre.row(w)=inv_s_temp.st();
    }

    fmat x_amp = arma::abs(x);


    arma::Mat<float> m_radon(nt,np);
    for(int i=0; i<np; i++)
    {
        arma::Col<cx_float> m_radon_temp = ifft(x.col(i));
        m_radon.col(i)                   = arma::abs(m_radon_temp);
    }





    //========================================================================//
    //                              逆RADON变换                               //
    //========================================================================//

    arma::Mat<float> invradon_sei(nt,nx,fill::zeros);
    arma::Mat<cx_float> inv_temp(nt,nx,fill::zeros);
    for(int i=0; i<nx;i++)
    {
        inv_temp.col(i)=ifft(inv_radon_sei_fre.col(i));

        for(int j=0; j<nt;j++)
        {
            invradon_sei(j,i)=inv_temp(j,i).real();
        }
    }
    //difference of ori_seis and inv_seis//
    arma::Mat<float> sei_diff(nt,nx,fill::zeros);
    sei_diff = sei - invradon_sei;







    //========================================================================//
    //                              文件保存                                  //
    //========================================================================//

    //ref.save("../file/ref.dat",raw_binary);             // Reflection coefficition
    //wavelet.save("../file/wavelet.dat",raw_binary);     // Ricker Wavelet
    //sei.save("../file/sei.dat",raw_binary);               // Synthetic seismogram
    //x_amp.save("../file/x.dat",raw_binary);             // Model in omega-p domain
    //m_radon.save("../file/m_radon.dat",raw_binary);       // Model in t-p domain
    //invradon_sei.save("../file/inv_sei.dat",raw_binary);  // The result of I_Radon_T
    //sei_diff.save("../file/sei_diff.dat",raw_binary);     // The difference




    return 0.0;

    */


}
