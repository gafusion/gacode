/*
 *
 * =====================================================
 * Direct solve routines using the CUDA CUSOLVER library
 * =====================================================
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <cusolverSp.h>
#include <cusparse.h>

#define USE_DOUBLE 1

#if USE_DOUBLE
int cusolve_sparse_(int *n, int *nnz, double *csrValA_H, int *csrRowPtrA_H, int *csrColIndA_H, double *b_H, double *x_H)
{
    int m = *n;
    int nnzA = *nnz;
    int i,j,k;
    int singularity = -1;
    cusolverSpHandle_t cusolverH = NULL;
    cudaStream_t stream = NULL;

    cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
    cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;

    cudaError_t cudaStat1 = cudaSuccess;
    cudaError_t cudaStat2 = cudaSuccess;
    cudaError_t cudaStat3 = cudaSuccess;
    cudaError_t cudaStat4 = cudaSuccess;
    cudaError_t cudaStat5 = cudaSuccess;

    if (x_H == NULL) x_H = b_H;
    assert(csrValA_H != NULL);
    assert(csrRowPtrA_H != NULL);
    assert(csrColIndA_H != NULL);
    assert(b_H != NULL);

/* step 1: create cusolver handle, bind a stream */
    status = cusolverSpCreate(&cusolverH);
    assert(CUSOLVER_STATUS_SUCCESS == status);

    cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    assert(cudaSuccess == cudaStat1);

    status = cusolverSpSetStream(cusolverH, stream);
    assert(CUSOLVER_STATUS_SUCCESS == status);

/* step 2: copy A to device */
    double *csrValA_D    = NULL;
    int    *csrRowPtrA_D = NULL;
    int    *csrColIndA_D = NULL;
    double *b_D          = NULL;
    double *x_D          = NULL;

    cudaStat1 = cudaMalloc ((void**)&csrValA_D, sizeof(double) * nnzA);
    cudaStat2 = cudaMalloc ((void**)&csrRowPtrA_D, sizeof(int) * (m+1));
    cudaStat3 = cudaMalloc ((void**)&csrColIndA_D, sizeof(int) * nnzA);
    cudaStat4 = cudaMalloc ((void**)&b_D, sizeof(double) * m);
    cudaStat5 = cudaMalloc ((void**)&x_D, sizeof(double) * m);

    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);
    assert(cudaSuccess == cudaStat3);
    assert(cudaSuccess == cudaStat4);
    assert(cudaSuccess == cudaStat5);

    cudaStat1 = cudaMemcpy(csrValA_D, csrValA_H, sizeof(double)*nnzA, cudaMemcpyHostToDevice);
    cudaStat2 = cudaMemcpy(csrRowPtrA_D, csrRowPtrA_H, sizeof(int)*(m+1), cudaMemcpyHostToDevice);
    cudaStat3 = cudaMemcpy(csrColIndA_D, csrColIndA_H, sizeof(int)*nnzA, cudaMemcpyHostToDevice);
    cudaStat4 = cudaMemcpy(b_D, b_H, sizeof(double)*m, cudaMemcpyHostToDevice);

    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);
    assert(cudaSuccess == cudaStat3);
    assert(cudaSuccess == cudaStat4);

/* step 5: solve A*X = B */
    cusparseMatDescr_t descrA = NULL;
    cusparse_status = cusparseCreateMatDescr(&descrA);
    assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);
    cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descrA,CUSPARSE_INDEX_BASE_ONE);

    status = cusolverSpDcsrlsvqr( cusolverH, m, nnzA, descrA, csrValA_D, csrRowPtrA_D, csrColIndA_D, b_D, 1.0E-12, 1, x_D, &singularity);

    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == status);
    //printf("Singularity = %d\n",singularity);
    assert(singularity == -1);

    cudaStat1 = cudaMemcpy(x_H , x_D, sizeof(double)*m, cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cudaStat1);

/* free resources */
    if (csrValA_D   ) cudaFree(csrValA_D   );
    if (csrRowPtrA_D) cudaFree(csrRowPtrA_D);
    if (csrColIndA_D) cudaFree(csrColIndA_D);
    if (b_D         ) cudaFree(b_D         );
    if (x_D         ) cudaFree(x_D         );

    if (cusolverH   ) cusolverSpDestroy(cusolverH);
    if (stream      ) cudaStreamDestroy(stream);

    cudaDeviceReset();

    return 0;
}

#else

int cusolve_sparse_(int *n, int *nnz, double *csrValA_H_D, int *csrRowPtrA_H, int *csrColIndA_H, double *b_H_D, double *x_H_D)
{
    int m = *n;
    int nnzA = *nnz;
    int i,j,k;
    int singularity = -1;
    cusolverSpHandle_t cusolverH = NULL;
    cudaStream_t stream = NULL;

    cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
    cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;

    cudaError_t cudaStat1 = cudaSuccess;
    cudaError_t cudaStat2 = cudaSuccess;
    cudaError_t cudaStat3 = cudaSuccess;
    cudaError_t cudaStat4 = cudaSuccess;
    cudaError_t cudaStat5 = cudaSuccess;

    if (x_H_D == NULL) x_H_D = b_H_D;
    assert(csrValA_H_D != NULL);
    assert(csrRowPtrA_H != NULL);
    assert(csrColIndA_H != NULL);
    assert(b_H_D != NULL);

/* step 0: create float arrays */
    float *csrValA_H = (float*)malloc(sizeof(float) * nnzA);
    float *b_H       = (float*)malloc(sizeof(float) * m);
    float *x_H       = (float*)malloc(sizeof(float) * m);

    for (i = 0; i < nnzA; i++) csrValA_H[i] = (float)csrValA_H_D[i];
    for (i = 0; i < m; i++)    b_H[i]       = (float)b_H_D[i];

/* step 1: create cusolver handle, bind a stream */
    status = cusolverSpCreate(&cusolverH);
    assert(CUSOLVER_STATUS_SUCCESS == status);

    cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    assert(cudaSuccess == cudaStat1);

    status = cusolverSpSetStream(cusolverH, stream);
    assert(CUSOLVER_STATUS_SUCCESS == status);

/* step 2: copy A to device */
    float *csrValA_D     = NULL;
    int    *csrRowPtrA_D = NULL;
    int    *csrColIndA_D = NULL;
    float *b_D           = NULL;
    float *x_D           = NULL;

    cudaStat1 = cudaMalloc ((void**)&csrValA_D, sizeof(float) * nnzA);
    cudaStat2 = cudaMalloc ((void**)&csrRowPtrA_D, sizeof(int) * (m+1));
    cudaStat3 = cudaMalloc ((void**)&csrColIndA_D, sizeof(int) * nnzA);
    cudaStat4 = cudaMalloc ((void**)&b_D, sizeof(float) * m);
    cudaStat5 = cudaMalloc ((void**)&x_D, sizeof(float) * m);

    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);
    assert(cudaSuccess == cudaStat3);
    assert(cudaSuccess == cudaStat4);
    assert(cudaSuccess == cudaStat5);

    cudaStat1 = cudaMemcpy(csrValA_D, csrValA_H, sizeof(float)*nnzA, cudaMemcpyHostToDevice);
    cudaStat2 = cudaMemcpy(csrRowPtrA_D, csrRowPtrA_H, sizeof(int)*(m+1), cudaMemcpyHostToDevice);
    cudaStat3 = cudaMemcpy(csrColIndA_D, csrColIndA_H, sizeof(int)*nnzA, cudaMemcpyHostToDevice);
    cudaStat4 = cudaMemcpy(b_D, b_H, sizeof(float)*m, cudaMemcpyHostToDevice);

    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);
    assert(cudaSuccess == cudaStat3);
    assert(cudaSuccess == cudaStat4);

/* step 5: solve A*X = B */
    cusparseMatDescr_t descrA = NULL;
    cusparse_status = cusparseCreateMatDescr(&descrA);
    assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);
    cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descrA,CUSPARSE_INDEX_BASE_ONE);

    status = cusolverSpScsrlsvqr( cusolverH, m, nnzA, descrA, csrValA_D, csrRowPtrA_D, csrColIndA_D, b_D, 1.0E-12, 1, x_D, &singularity);

    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == status);
    //printf("Singularity = %d\n",singularity);
    assert(singularity == -1);

    cudaStat1 = cudaMemcpy(x_H , x_D, sizeof(float)*m, cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cudaStat1);
    for (i = 0; i < m; i++) x_H_D[i] = (double)x_H[i];

/* free resources */
    if (csrValA_D   ) cudaFree(csrValA_D   );
    if (csrRowPtrA_D) cudaFree(csrRowPtrA_D);
    if (csrColIndA_D) cudaFree(csrColIndA_D);
    if (b_D         ) cudaFree(b_D         );
    if (x_D         ) cudaFree(x_D         );

    if (cusolverH   ) cusolverSpDestroy(cusolverH);
    if (stream      ) cudaStreamDestroy(stream);

    cudaDeviceReset();

    free(csrValA_H);
    free(b_H);
    free(x_H);

    return 0;
}

#endif
