#include "../inc/ss_sub_tools_func.h"

extern "C"
{

    void FFTCC1(int isign, fcomplex *z, int n)
        // isign : sign of i when doing fft
        // z     : fcomplex data array
        // n     : length of z
    {
        int nn = n;         // tmp for fortran subroutine
        int ii = isign;     // tmp for fortran subroutine
        // if isign==0 ,then DO NOTHING!!!
        if ( isign != 0)
        { fft1d_( z, &nn, &ii); }
    }

    void FFTCC2(int isign1, int isign2, fcomplex *zz, int n1, int n2,fcomplex *work)
        // isign1 : isign  of fast dim
        // isign2 : isign  of slow dim
        // zz     : fcomplex data array (saved as 1d array) length of zz = n1*n2
        // n1     : length of fast dim
        // n2     : length of slow dim
        // work   : working array ,length of work = n2
    {
        int i,j;
        // STEP #1
        // DO 1d-fft on fast dim
        // fft each sub-array
        for (i=0;i<n2;i++)
        {FFTCC1( isign1, &zz[0+i*n1], n1);}

        // STEP #2
        // DO 2d-fft on slow dim

        for (j=0;j<n1;j++){
            // get j-th slow sub-array into a successive array
            for (i=0;i<n2;i++)
            {work[i] = zz[j+i*n1];}
            // do fft
            FFTCC1( isign2, work, n2);
            // put back
            for (i=0;i<n2;i++)
            {zz[j+i*n1] = work[i];}
        }


    }

    /*============================================================================*/
    void FFTCC3(int isign1   , int isign2 , int isign3 ,
            int n1       , int n2     , int n3     ,
            fcomplex *zz , int LevelVerbose)
    {
        if (LevelVerbose > 0)
        {
            printf("FFTCC3 n1,n2,n3=%d %d %d BEGIN\n",n1,n2,n3);
        }
        //auto t0 = GetTimer();
        int num2d = n2*n1;
        int work_len = MAX(n2,n3);
        //#pragma omp parallel
        {

            fcomplex * work = (fcomplex*) calloc(work_len,sizeof(fcomplex));
            //#pragma omp for
            for (int I3=0;I3<n3;I3++)
            {
                fcomplex * p = &zz[I3*num2d];
                FFTCC2(isign1, isign2, p, n1, n2,work);
            }

            //#pragma omp barrier

            // STEP #2
            // DO 1d-fft on slow dim

            //#pragma omp for schedule(dynamic,2)
            for (int I=0;I<num2d;I++)
            {
                // get j-th slow sub-array into a successive array
                for (int i=0;i<n3;i++)
                { work[i] = zz[I+i*num2d]; }
                // do fft
                FFTCC1( isign3, work, n3);
                // put back
                for (int i=0;i<n3;i++)
                {zz[I+i*num2d] = work[i];}
            }
            free(work);
        }

        //auto t1 = GetTimer();
        //double t_used= GetTimeDiff(t0,t1);

        //if (LevelVerbose > 0)
        //{
        //printf("FFTCC3 n1,n2,n3=%d %d %d END, time_used=[%g] sec\n",n1,n2,n3,t_used);
        //}


    }



}



