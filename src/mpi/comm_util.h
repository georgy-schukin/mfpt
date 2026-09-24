#pragma once

#include <mpi.h>

#include "defs.h"
#include "block_decomp.h"

MPI_Datatype makeColType(const DArray2 &array);
MPI_Datatype makeRowType(const DArray2 &array);
MPI_Datatype makeDataVectorType(int num_of_blocks, int block_size, int stride);
MPI_Datatype makeDataVectorType(const DArray2 &array);

void syncShadowsKPrev(DArray2 &arr, MPI_Datatype col_type, int rank, int size);
void syncShadowsKNext(DArray2 &arr, MPI_Datatype col_type, int rank, int size);
void syncShadowsIPrev(DArray2 &arr, MPI_Datatype row_type, int rank, int size);
void syncShadowsINext(DArray2 &arr, MPI_Datatype row_type, int rank, int size);
void syncShadowsK(DArray2 &arr, MPI_Datatype col_type, int rank, int size);
void syncShadowsI(DArray2 &arr, MPI_Datatype row_type, int rank, int size);

DArray2 gatherArrayK(const DArray2 &local_data, const BlockDecomposition &k_decomp, int rank, int size, int root = 0);
DArray2 scatterArrayK(const DArray2 &data, int size_x, const BlockDecomposition &k_decomp, int shadow_x, int shadow_y, int rank, int size, int root = 0);

void combineIFromK(const DArray2 &src, DArray2 &dst, const BlockDecomposition &i_decomp, const BlockDecomposition &k_decomp, int rank, int size);
void combineKFromI(const DArray2 &src, DArray2 &dst, const BlockDecomposition &k_decomp, const BlockDecomposition &i_decomp, int rank, int size);

MPI_Op makeVectorSumOp();
