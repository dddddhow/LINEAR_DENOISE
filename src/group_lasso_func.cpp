#include "../inc/group_lasso_func.h"
#include "../inc/mycomplex.h"

using namespace arma;
using namespace std;

/*=============================================================================*/
// In_A  : mat
// x     : solution
// f     : f = Ax
// In_d  : right vector = f
// Nhx   :
// Np    :
// R     : residual
// err   : allowed length of residual
/*============================================================================*/


void lasso_func(const void *In_A_raw, void *Out_x_raw,const void *In_d_raw,
                const int *Np, const float *errr);
void cg_func( const void *In_A_raw, void *Out_x_raw, const void *In_d_raw,
                const int *Np, const float *errr);


//=============================================================================//

void group_lasso_func(fcomplex *In_A_raw, fcomplex *Out_x_raw, fcomplex *In_d_raw,
        int *Nhx,  int *Np, float *errr)
{
    complex<float> * In_A        = (complex<float> * )In_A_raw;
    complex<float> * const Out_x = Out_x_raw;
    complex<float> * In_d        = (complex<float> * )In_d_raw;
    const int np                 = *Np;
    const int nhx                = *Nhx;
    float err                    = *errr;


    //=============================Definition=================================//

    arma::Mat<cx_float> A(nhx, np);
    cx_fcolvec    X(np);
    cx_fcolvec    f(nhx);
    cx_fcolvec    R(nhx);

    float lmada    = 0.00;
    float lmada_max=50;
    float lamada_m = 0.0;
    float lamada_m_max = 50;
    int iter_num   = 5;
    int group_size = 5;
    int group_num  = np/group_size;
    arma::Mat<cx_float>  A_temp(nhx, group_size);
    arma::Mat<cx_float>  A_temp_in(group_size, group_size);
    cx_fcolvec     X_temp(group_size);
    cx_fcolvec     X_temp_out(group_size);
    cx_fcolvec     A_temp_f(nhx);
    cx_fcolvec     temp_f(nhx);
    cx_fcolvec     temp_vec_1(nhx);
    cx_fcolvec     temp_vec_2(nhx);

    memcpy( &A(0,0), In_A, sizeof(fcomplex)*nhx*np );

    for(int i=0;i<nhx;i++)
    {
        f(i)=In_d[i];
    }

    fvec vc_tmp = {1.0E-2, 5*err};

    //init
    X.fill(fcomplex(0.0f));
    arma::Mat<cx_float> ATA(nhx, np,fill::zeros);
    cx_fcolvec    ATD(nhx,fill::zeros);
    cx_fcolvec    X_L2(np,fill::zeros);

    ATA          = A.t() * A;
    ATA.diag()  += lmada ;
    ATD          = A.t() * f;

    // computation for the initial X (here define as X_L2 , L2+L2)
    cg_func(&ATA(0,0), &X_L2(0),&ATD(0),&np,&err);
    //lasso_func(&ATA(0,0), &X_L2(0),&ATD(0),&np,&err);


    //lamada_m_max = abs(X_L2.max());

    for(int iter = 0; iter < iter_num;  iter++)
    {
        //自适应lmada_m
        lamada_m = lamada_m_max * 0.01;
        lmada    = lmada_max *0.01;

        for(int group = 0; group < group_num; group++)
        {
            int p = group * group_size;
            int q = (group+1) * group_size-1;
            A_temp = A.cols(p,q);
            X_temp = X_L2(span(p,q));


            //x_temp_abs
            float x_temp_abs = 0;
            for (int i = 0; i < group_size; i++)
            {
                x_temp_abs = x_temp_abs + abs(X_temp(i));
            }

            //lamada_m_max
            if(x_temp_abs >= lamada_m_max )
            {
                lamada_m_max = x_temp_abs;
            }

            if( x_temp_abs<=lamada_m)
            {
                X_temp_out.fill(fcomplex(0.0f));
            }
            else if( x_temp_abs >lamada_m)
            {
                //computate AkXk
                //temp_vec_1 = temp_f =AKXK
                //temp_vec_2 = Ak(f-AKXK)

                temp_vec_1.fill(fcomplex(0.0f));
                temp_vec_2.fill(fcomplex(0.0f));
                X_temp_out.fill(fcomplex(0.0f));

                temp_vec_1 = f - (A * X_L2 - A.cols(p,q) * X_L2(span(p,q)));
                temp_vec_2 = A_temp.t() * temp_vec_1;

                float lmada_max_temp=abs(temp_vec_2.max());
                if(lmada_max <= lmada_max_temp)
                {
                    lmada_max = lmada_max_temp;
                }

                double temp = 0;
                for (int i = 0; i < group_size; i++)
                {
                    temp = temp + abs(temp_vec_2(i));
                }

                if (temp < lmada)
                {
                    X_temp_out.fill(fcomplex(0.0f));
                }
                else if (temp >= lmada)
                {
                    A_temp_in          = A_temp.t() * A_temp;
                    A_temp_in.diag()  += lmada * 1.0 / x_temp_abs;
                    temp_f             = f - (A * X_L2 - A.cols(p,q) * X_L2(span(p,q)));
                    A_temp_f           = A_temp.t() * temp_f;

                    lasso_func( &A_temp_in(0,0), &X_temp_out(0), &A_temp_f(0), &group_size, &err);
                }
            }

            if(iter == iter_num - 1)
            {
                X_L2(span(p,q)) = X_temp_out;
            }
        }
    }

    X = X_L2 ;

    //=============================== Output =================================
    for(int i=0;i<np;i++)
    {
        Out_x[i] = X(i);
    }
    return;
}




/*============================================================================*/
void cg_func( const void *In_A_raw, void *Out_x_raw, const void *In_d_raw, const int *Np, const float *errr)
{
    typedef complex<float> fcomplex;

    complex<float> * In_A  = (complex<float> * )In_A_raw;
    complex<float> * Out_x = (complex<float> * )Out_x_raw;
    complex<float> * In_d  = (complex<float> * )In_d_raw;
    int np=*Np;
    float err=*errr;

    //=============================Definition=================================

    arma::Mat<cx_float> A(np, np);
    cx_fcolvec X(np);
    cx_fcolvec f(np);

    cx_fcolvec R(np);
    cx_fcolvec AR(np);


    int iter_num = 100;
    double norm2min, errs;
    double alpha   , beta;
    double rs_old  ,rs_new;
    double RAR;

    memcpy( &A(0,0), In_A, sizeof(fcomplex)*np*np );

    for(int i=0;i<np;i++)
    {
        f(i)=In_d[i];
    }

    fvec vc_tmp = {1.0E-2, 5*err};
    norm2min = arma::min(vc_tmp);
    errs     = double(err) * double(err);

    //Zeros
    R.fill(fcomplex(0.0f));
    AR.fill(fcomplex(0.0f));
    X.fill(fcomplex(0.0f));

    //Initial
    //float lmada = 0.005;
    //A.diag() += lmada *1.0;
    R = A*X - f;
    rs_new = real(cdot(R,R));
    rs_old = rs_new;

    //==========================Begin GD Process==============================
    for (int iteration = 0; iteration<iter_num; iteration++)
    {

        if(abs(rs_old) <= norm2min)
        {
            break;
        }

        //compute alpha dot(r,r)/dot(AP,P)
        AR  = A*R;
        RAR = std::real(cdot(R,AR));

        if(abs(RAR) <= norm2min)
        {
            break;
        }
        //compute alpha
        alpha  = rs_old*1.0/RAR;

        //compute X   X=X-alpha*R
        X      -= alpha*R;

        //=========================== L1  Process== ==============================
        float lamda_r ;
        float lamda_i ;
        lamda_r = alpha * 0.000;
        lamda_i = alpha * 0.000;


        for(int i=0; i<np; i++)
        {
            if(X(i).real() == abs(X(i).real()))
            {
                X(i).real(abs(X(i).real()) - lamda_r);
            }
            else if(X(i).real() == -abs(X(i).real()))
            {
                X(i).real(-(abs(X(i).real()) -lamda_r));
            }


            if(X(i).imag() == abs(X(i).imag()))
            {
                X(i).imag(abs(X(i).imag()) - lamda_i);
            }
            else if(X(i).imag() == -abs(X(i).imag()))
            {
                X(i).imag(-(abs(X(i).imag()) -lamda_i));
            }

            if(abs(X(i).real()) < lamda_r)
            {
                X(i).real(0);
            }
            if(abs(X(i).imag()) < lamda_i)
            {
                X(i).imag(0);
            }

        }

        //compute R   R = AX - f
        R      = A*X - f;
        rs_new = real(cdot(R,R));
        rs_old = rs_new;
    }


    //=============================== Output =================================
    for(int i=0;i<np;i++)
    {
        Out_x[i] =X(i);
    }


}

void lasso_func(const void *In_A_raw, void *Out_x_raw,const void *In_d_raw,const int *Np, const float *errr)
{
    typedef complex<float> fcomplex;

    complex<float> * In_A  = (complex<float> * )In_A_raw;
    complex<float> * Out_x = (complex<float> * )Out_x_raw;
    complex<float> * In_d  = (complex<float> * )In_d_raw;
    int np=*Np;
    float err=*errr;

    //=============================Definition=================================

    arma::Mat<cx_float> A(np, np);
    cx_fcolvec X(np);
    cx_fcolvec f(np);

    cx_fcolvec R(np);
    cx_fcolvec AR(np);


    int iter_num = np;
    double norm2min, errs;
    double alpha   , beta;
    double rs_old  ,rs_new;
    double RAR;

    memcpy( &A(0,0), In_A, sizeof(fcomplex)*np*np );

    for(int i=0;i<np;i++)
    {
        f(i)=In_d[i];
    }

    fvec vc_tmp = {1.0E-2, 5*err};
    norm2min = arma::min(vc_tmp);
    errs     = double(err) * double(err);

    //Zeros
    R.fill(fcomplex(0.0f));
    AR.fill(fcomplex(0.0f));
    X.fill(fcomplex(0.0f));

    //Initial
    R = A*X - f;
    rs_new = real(cdot(R,R));
    rs_old = rs_new;

    //==========================Begin GD Process==============================
    for (int iteration = 0; iteration<iter_num; iteration++)
    {

        if(abs(rs_old) <= norm2min)
        {
            break;
        }

        //compute alpha dot(r,r)/dot(AP,P)
        AR  = A*R;
        RAR = std::real(cdot(R,AR));

        if(abs(RAR) <= norm2min)
        {
            break;
        }
        //compute alpha
        alpha  = rs_old*1.0/RAR;

        //compute X   X=X-alpha*R
        X      -= alpha*R;

        //=========================== L1  Process== ==============================
        float lamda_r ;
        float lamda_i ;
        lamda_r = alpha * 50;
        lamda_i = alpha * 50;


        for(int i=0; i<np; i++)
        {
            if(X(i).real() == abs(X(i).real()))
            {
                X(i).real(abs(X(i).real()) - lamda_r);
            }
            else if(X(i).real() == -abs(X(i).real()))
            {
                X(i).real(-(abs(X(i).real()) -lamda_r));
            }


            if(X(i).imag() == abs(X(i).imag()))
            {
                X(i).imag(abs(X(i).imag()) - lamda_i);
            }
            else if(X(i).imag() == -abs(X(i).imag()))
            {
                X(i).imag(-(abs(X(i).imag()) -lamda_i));
            }

            if(abs(X(i).real()) < lamda_r)
            {
                X(i).real(0);
            }
            if(abs(X(i).imag()) < lamda_i)
            {
                X(i).imag(0);
            }

        }

        //compute R   R = AX - f
        R      = A*X - f;
        rs_new = real(cdot(R,R));
        rs_old = rs_new;
    }


    //=============================== Output =================================
    for(int i=0;i<np;i++)
    {
        Out_x[i] =X(i);
    }


}

