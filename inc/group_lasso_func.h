#ifndef GROUP_LASSO_FUNC_H
#define GROUP_LASSO_FUNC_H

#include "TJU_SS_HEAD.h"
#include "mycomplex.h"

void group_lasso_func(fcomplex *In_A_raw, fcomplex *Out_x_raw, fcomplex *In_d_raw,
        int *Nhx,  int *Np, float *errr);

#endif
