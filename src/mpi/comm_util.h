#pragma once

#include <mpi.h>

#include "defs.h"
#include "block_decomp.h"

MPI_Datatype makeColType(const DArray2 &array);
MPI_Datatype makeRowType(const DArray2 &array);
MPI_Datatype makeRowType(int block_size, int extent);
MPI_Datatype makeDataVectorType(int num_of_blocks, int block_size, int stride, int extent = -1);
MPI_Datatype makeDataVectorType(const DArray2 &array);

void syncShadowsColsPrev(DArray2 &arr, MPI_Datatype col_type, int rank, int size);
void syncShadowsColsNext(DArray2 &arr, MPI_Datatype col_type, int rank, int size);
void syncShadowsRowsPrev(DArray2 &arr, MPI_Datatype row_type, int rank, int size);
void syncShadowsRowsNext(DArray2 &arr, MPI_Datatype row_type, int rank, int size);
void syncShadowsCols(DArray2 &arr, MPI_Datatype col_type, int rank, int size);
void syncShadowsRows(DArray2 &arr, MPI_Datatype row_type, int rank, int size);

DArray2 gatherArrayCols(const DArray2 &local_data, const BlockDecomposition &cols_decomp, int rank, int size, int root = 0);
DArray2 scatterArrayCols(const DArray2 &data, int size_x, const BlockDecomposition &cols_decomp, int shadow_x, int shadow_y, int rank, int size, int root = 0);

void combineRowsFromCols(const DArray2 &src, DArray2 &dst, const BlockDecomposition &rows_decomp, const BlockDecomposition &cols_decomp, int rank, int size);
void combineColsFromRows(const DArray2 &src, DArray2 &dst, const BlockDecomposition &cols_decomp, const BlockDecomposition &rows_decomp, int rank, int size);

MPI_Op makeVectorSumOp();
