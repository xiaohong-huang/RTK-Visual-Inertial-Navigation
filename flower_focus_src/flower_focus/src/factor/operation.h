
#pragma once
#include <eigen3/Eigen/Dense>
//this file is from Ceres-Solver https://github.com/ceres-solver/ceres-solver


#define GEMM_STORE_SINGLE(p, index, value) \
  if (kOperation > 0) {                          \
    p[index] += value;                           \
  } else if (kOperation < 0) {                   \
    p[index] -= value;                           \
  } else {                                       \
    p[index] = value;                            \
  }

#define GEMM_STORE_PAIR(p, index, v1, v2) \
  if (kOperation > 0) {                         \
    p[index] += v1;                             \
    p[index + 1] += v2;                         \
  } else if (kOperation < 0) {                  \
    p[index] -= v1;                             \
    p[index + 1] -= v2;                         \
  } else {                                      \
    p[index] = v1;                              \
    p[index + 1] = v2;                          \
  }

#define GEMM_OPT_NAIVE_HEADER \
  double c0 = 0.0;                  \
  double c1 = 0.0;                  \
  double c2 = 0.0;                  \
  double c3 = 0.0;                  \
  const double* pa = a;             \
  const double* pb = b;             \
  const int span = 4;               \
  int col_r = col_a & (span - 1);   \
  int col_m = col_a - col_r;


#define GEMM_OPT_STORE_MAT1X4 \
  if (kOperation > 0) {             \
    *c++ += c0;                     \
    *c++ += c1;                     \
    *c++ += c2;                     \
    *c++ += c3;                     \
  } else if (kOperation < 0) {      \
    *c++ -= c0;                     \
    *c++ -= c1;                     \
    *c++ -= c2;                     \
    *c++ -= c3;                     \
  } else {                          \
    *c++ = c0;                      \
    *c++ = c1;                      \
    *c++ = c2;                      \
    *c++ = c3;                      \
  }

static inline void MMM_mat1x4(const int col_a,
                              const double* a,
                              const double* b,
                              const int col_stride_b,
                              double* c,
                              const int kOperation) {
    GEMM_OPT_NAIVE_HEADER
    double av = 0.0;
    int bi = 0;

#define GEMM_OPT_MMM_MAT1X4_MUL \
  av = pa[k];                         \
  pb = b + bi;                        \
  c0 += av * *pb++;                   \
  c1 += av * *pb++;                   \
  c2 += av * *pb++;                   \
  c3 += av * *pb++;                   \
  bi += col_stride_b;                 \
  k++;

    for (int k = 0; k < col_m;) {
        GEMM_OPT_MMM_MAT1X4_MUL
        GEMM_OPT_MMM_MAT1X4_MUL
        GEMM_OPT_MMM_MAT1X4_MUL
        GEMM_OPT_MMM_MAT1X4_MUL
    }

    for (int k = col_m; k < col_a;) {
        GEMM_OPT_MMM_MAT1X4_MUL
    }

    GEMM_OPT_STORE_MAT1X4

#undef GEMM_OPT_MMM_MAT1X4_MUL
}


static inline void MTV_mat4x1(const int col_a,
                              const double* a,
                              const int col_stride_a,
                              const double* b,
                              double* c,
                              const int kOperation) {
    GEMM_OPT_NAIVE_HEADER
    double bv = 0.0;

    // clang-format off
#define GEMM_OPT_MTV_MAT4X1_MUL \
  bv = *pb;                           \
  c0 += *(pa    ) * bv;               \
  c1 += *(pa + 1) * bv;               \
  c2 += *(pa + 2) * bv;               \
  c3 += *(pa + 3) * bv;               \
  pa += col_stride_a;                 \
  pb++;
    // clang-format on

    for (int k = 0; k < col_m; k += span) {
        GEMM_OPT_MTV_MAT4X1_MUL
        GEMM_OPT_MTV_MAT4X1_MUL
        GEMM_OPT_MTV_MAT4X1_MUL
        GEMM_OPT_MTV_MAT4X1_MUL
    }

    for (int k = col_m; k < col_a; k++) {
        GEMM_OPT_MTV_MAT4X1_MUL
    }

    GEMM_OPT_STORE_MAT1X4

#undef GEMM_OPT_MTV_MAT4X1_MUL
}


static inline void MVM_mat4x1(const int col_a,
                              const double* a,
                              const int col_stride_a,
                              const double* b,
                              double* c,
                              const int kOperation) {
    GEMM_OPT_NAIVE_HEADER
    double bv = 0.0;

    // clang-format off
#define GEMM_OPT_MVM_MAT4X1_MUL  \
  bv = *pb;                            \
  c0 += *(pa                   ) * bv; \
  c1 += *(pa + col_stride_a    ) * bv; \
  c2 += *(pa + col_stride_a * 2) * bv; \
  c3 += *(pa + col_stride_a * 3) * bv; \
  pa++;                                \
  pb++;
    // clang-format on

    for (int k = 0; k < col_m; k += span) {
        GEMM_OPT_MVM_MAT4X1_MUL
        GEMM_OPT_MVM_MAT4X1_MUL
        GEMM_OPT_MVM_MAT4X1_MUL
        GEMM_OPT_MVM_MAT4X1_MUL
    }

    for (int k = col_m; k < col_a; k++) {
        GEMM_OPT_MVM_MAT4X1_MUL
    }

    GEMM_OPT_STORE_MAT1X4

#undef GEMM_OPT_MVM_MAT4X1_MUL
}


static inline void MTM_mat1x4(const int col_a,
                              const double* a,
                              const int col_stride_a,
                              const double* b,
                              const int col_stride_b,
                              double* c,
                              const int kOperation) {
    GEMM_OPT_NAIVE_HEADER
    double av = 0.0;
    int ai = 0;
    int bi = 0;

#define GEMM_OPT_MTM_MAT1X4_MUL \
  av = pa[ai];                        \
  pb = b + bi;                        \
  c0 += av * *pb++;                   \
  c1 += av * *pb++;                   \
  c2 += av * *pb++;                   \
  c3 += av * *pb++;                   \
  ai += col_stride_a;                 \
  bi += col_stride_b;

    for (int k = 0; k < col_m; k += span) {
        GEMM_OPT_MTM_MAT1X4_MUL
        GEMM_OPT_MTM_MAT1X4_MUL
        GEMM_OPT_MTM_MAT1X4_MUL
        GEMM_OPT_MTM_MAT1X4_MUL
    }

    for (int k = col_m; k < col_a; k++) {
        GEMM_OPT_MTM_MAT1X4_MUL
    }

    GEMM_OPT_STORE_MAT1X4

#undef GEMM_OPT_MTM_MAT1X4_MUL
}


template <int kRowA, int kColA, int kRowB, int kColB, int kOperation>
inline void MatrixMatrixMultiply(const double* A,
                                 const int num_row_a,
                                 const int num_col_a,
                                 const double* B,
                                 const int num_row_b,
                                 const int num_col_b,
                                 double* C,
                                 const int start_row_c,
                                 const int start_col_c,
                                 const int row_stride_c,
                                 const int col_stride_c,
                                 bool is_sym) {
    const int NUM_ROW_A = (kRowA != Eigen::Dynamic ? kRowA : num_row_a);
    \
    const int NUM_COL_A = (kColA != Eigen::Dynamic ? kColA : num_col_a);
    \
    const int NUM_COL_B = (kColB != Eigen::Dynamic ? kColB : num_col_b);
    const int NUM_ROW_C = NUM_ROW_A;
    const int NUM_COL_C = NUM_COL_B;
    const int span = 4;

    // Calculate the remainder part first.

    // Process the last odd column if present.
    if (NUM_COL_C & 1) {
        int col = NUM_COL_C - 1;
        const double* pa = &A[0];
        for (int row = 0; row < NUM_ROW_C; ++row, pa += NUM_COL_A) {
            const double* pb = &B[col];
            double tmp = 0.0;
            for (int k = 0; k < NUM_COL_A; ++k, pb += NUM_COL_B)
                tmp += pa[k] * pb[0];

            const int index = (row + start_row_c) * col_stride_c + start_col_c + col;
            GEMM_STORE_SINGLE(C, index, tmp);
        }

        // Return directly for efficiency of extremely small matrix multiply.
        if (NUM_COL_C == 1)
            return;
    }

    // Process the couple columns in remainder if present.
    if (NUM_COL_C & 2) {
        int col = NUM_COL_C & (int)(~(span - 1));
        const double* pa = &A[0];
        for (int row = 0; row < NUM_ROW_C; ++row, pa += NUM_COL_A) {
            const double* pb = &B[col];
            double tmp1 = 0.0, tmp2 = 0.0;
            for (int k = 0; k < NUM_COL_A; ++k, pb += NUM_COL_B) {
                double av = pa[k];
                tmp1 += av * pb[0];
                tmp2 += av * pb[1];
            }

            const int index = (row + start_row_c) * col_stride_c + start_col_c + col;
            GEMM_STORE_PAIR(C, index, tmp1, tmp2);
        }

        // Return directly for efficiency of extremely small matrix multiply.
        if (NUM_COL_C < span)
            return;
    }

    // Calculate the main part with multiples of 4.
    int col_m = NUM_COL_C & (int)(~(span - 1));
    if (is_sym) {
        for (int row = 0; row < NUM_ROW_C; ++row) {
            for (int col = row / span * span; col < col_m; col += span) {
                const int index = (row + start_row_c) * col_stride_c + start_col_c + col;
                // clang-format off
                MMM_mat1x4(NUM_COL_A, &A[row * NUM_COL_A],
                           &B[col], NUM_COL_B, &C[index], kOperation);
                // clang-format on
            }
        }
    } else {
        for (int row = 0; row < NUM_ROW_C; ++row) {
            for (int col = 0; col < col_m; col += span) {
                const int index = (row + start_row_c) * col_stride_c + start_col_c + col;
                // clang-format off
                MMM_mat1x4(NUM_COL_A, &A[row * NUM_COL_A],
                           &B[col], NUM_COL_B, &C[index], kOperation);
                // clang-format on
            }
        }
    }

}





inline flower_focus::MatrixXd InvertPSDMatrix(const flower_focus::MatrixXd& m) {
    return m.template selfadjointView<Eigen::Upper>().llt().solve(flower_focus::MatrixXd::Identity(m.rows(), m.rows()));
}






















template <int num_rows = Eigen::Dynamic, int num_cols = Eigen::Dynamic>
struct EigenTypes {
    using Matrix =
        Eigen::Matrix < double,
        num_rows,
        num_cols,
        num_cols == 1 ? Eigen::ColMajor : Eigen::RowMajor >;

    using MatrixRef = Eigen::Map<Matrix>;
    using ConstMatrixRef = Eigen::Map<const Matrix>;
    using Vector = Eigen::Matrix<double, num_rows, 1>;
    using VectorRef = Eigen::Map<Eigen::Matrix<double, num_rows, 1>>;
    using ConstVectorRef = Eigen::Map<const Eigen::Matrix<double, num_rows, 1>>;
};

// The following three macros are used to share code and reduce
// template junk across the various GEMM variants.
#define GEMM_BEGIN(name)                                          \
  template <int kRowA, int kColA, int kRowB, int kColB, int kOperation> \
  inline void name(const double* A,                                     \
                   const int num_row_a,                                 \
                   const int num_col_a,                                 \
                   const double* B,                                     \
                   const int num_row_b,                                 \
                   const int num_col_b,                                 \
                   double* C,                                           \
                   const int start_row_c,                               \
                   const int start_col_c,                               \
                   const int row_stride_c,                              \
                   const int col_stride_c)

#define GEMM_NAIVE_HEADER                                        \
  const int NUM_ROW_A = (kRowA != Eigen::Dynamic ? kRowA : num_row_a); \
  const int NUM_COL_A = (kColA != Eigen::Dynamic ? kColA : num_col_a); \
  const int NUM_COL_B = (kColB != Eigen::Dynamic ? kColB : num_col_b);

#define GEMM_EIGEN_HEADER                                 \
  const typename EigenTypes<kRowA, kColA>::ConstMatrixRef Aref( \
      A, num_row_a, num_col_a);                                 \
  const typename EigenTypes<kRowB, kColB>::ConstMatrixRef Bref( \
      B, num_row_b, num_col_b);                                 \
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> Cref(C, row_stride_c, col_stride_c);

// clang-format off
#define CALL_GEMM(name)                                           \
  name<kRowA, kColA, kRowB, kColB, kOperation>(                         \
      A, num_row_a, num_col_a,                                          \
      B, num_row_b, num_col_b,                                          \
      C, start_row_c, start_col_c, row_stride_c, col_stride_c);
// clang-format on


// For the matrix-matrix functions below, there are three variants for
// each functionality. Foo, FooNaive and FooEigen. Foo is the one to
// be called by the user. FooNaive is a basic loop based
// implementation and FooEigen uses Eigen's implementation. Foo
// chooses between FooNaive and FooEigen depending on how many of the
// template arguments are fixed at compile time. Currently, FooEigen
// is called if all matrix dimensions are compile time
// constants. FooNaive is called otherwise. This leads to the best
// performance currently.
//
// The MatrixMatrixMultiply variants compute:
//
//   C op A * B;
//
// The MatrixTransposeMatrixMultiply variants compute:
//
//   C op A' * B
//
// where op can be +=, -=, or =.
//
// The template parameters (kRowA, kColA, kRowB, kColB) allow
// specialization of the loop at compile time. If this information is
// not available, then Eigen::Dynamic should be used as the template
// argument.
//
//   kOperation =  1  -> C += A * B
//   kOperation = -1  -> C -= A * B
//   kOperation =  0  -> C  = A * B
//
// The functions can write into matrices C which are larger than the
// matrix A * B. This is done by specifying the true size of C via
// row_stride_c and col_stride_c, and then indicating where A * B
// should be written into by start_row_c and start_col_c.
//
// Graphically if row_stride_c = 10, col_stride_c = 12, start_row_c =
// 4 and start_col_c = 5, then if A = 3x2 and B = 2x4, we get
//
//   ------------
//   ------------
//   ------------
//   ------------
//   -----xxxx---
//   -----xxxx---
//   -----xxxx---
//   ------------
//   ------------
//   ------------
//
GEMM_BEGIN(MatrixMatrixMultiplyEigen) {
    GEMM_EIGEN_HEADER
    Eigen::Block<Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>, kRowA, kColB> block(
                                                                                                     Cref, start_row_c, start_col_c, num_row_a, num_col_b);

    if (kOperation > 0)
        block.noalias() += Aref * Bref;
    else if (kOperation < 0)
        block.noalias() -= Aref * Bref;
    else
        block.noalias() = Aref * Bref;
}

GEMM_BEGIN(MatrixMatrixMultiplyNaive) {
    GEMM_NAIVE_HEADER

    const int NUM_ROW_C = NUM_ROW_A;
    const int NUM_COL_C = NUM_COL_B;
    const int span = 4;

    // Calculate the remainder part first.

    // Process the last odd column if present.
    if (NUM_COL_C & 1) {
        int col = NUM_COL_C - 1;
        const double* pa = &A[0];
        for (int row = 0; row < NUM_ROW_C; ++row, pa += NUM_COL_A) {
            const double* pb = &B[col];
            double tmp = 0.0;
            for (int k = 0; k < NUM_COL_A; ++k, pb += NUM_COL_B)
                tmp += pa[k] * pb[0];

            const int index = (row + start_row_c) * col_stride_c + start_col_c + col;
            GEMM_STORE_SINGLE(C, index, tmp);
        }

        // Return directly for efficiency of extremely small matrix multiply.
        if (NUM_COL_C == 1)
            return;
    }

    // Process the couple columns in remainder if present.
    if (NUM_COL_C & 2) {
        int col = NUM_COL_C & (int)(~(span - 1));
        const double* pa = &A[0];
        for (int row = 0; row < NUM_ROW_C; ++row, pa += NUM_COL_A) {
            const double* pb = &B[col];
            double tmp1 = 0.0, tmp2 = 0.0;
            for (int k = 0; k < NUM_COL_A; ++k, pb += NUM_COL_B) {
                double av = pa[k];
                tmp1 += av * pb[0];
                tmp2 += av * pb[1];
            }

            const int index = (row + start_row_c) * col_stride_c + start_col_c + col;
            GEMM_STORE_PAIR(C, index, tmp1, tmp2);
        }

        // Return directly for efficiency of extremely small matrix multiply.
        if (NUM_COL_C < span)
            return;
    }

    // Calculate the main part with multiples of 4.
    int col_m = NUM_COL_C & (int)(~(span - 1));
    for (int col = 0; col < col_m; col += span) {
        for (int row = 0; row < NUM_ROW_C; ++row) {
            const int index = (row + start_row_c) * col_stride_c + start_col_c + col;
            // clang-format off
            MMM_mat1x4(NUM_COL_A, &A[row * NUM_COL_A],
                       &B[col], NUM_COL_B, &C[index], kOperation);
            // clang-format on
        }
    }
}

GEMM_BEGIN(MatrixMatrixMultiply) {


    // if (kRowA != Eigen::Dynamic && kColA != Eigen::Dynamic &&
    //         kRowB != Eigen::Dynamic && kColB != Eigen::Dynamic) {
    //     CALL_GEMM(MatrixMatrixMultiplyEigen)
    // } else
    {
        CALL_GEMM(MatrixMatrixMultiplyNaive)
    }

}

GEMM_BEGIN(MatrixTransposeMatrixMultiplyEigen) {
    GEMM_EIGEN_HEADER
    // clang-format off
    Eigen::Block<Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>, kColA, kColB> block(Cref,
                                                                                                 start_row_c, start_col_c,
                                                                                                 num_col_a, num_col_b);
    // clang-format on
    if (kOperation > 0)
        block.noalias() += Aref.transpose() * Bref;
    else if (kOperation < 0)
        block.noalias() -= Aref.transpose() * Bref;
    else
        block.noalias() = Aref.transpose() * Bref;
}

GEMM_BEGIN(MatrixTransposeMatrixMultiplyNaive) {
    GEMM_NAIVE_HEADER
    const int NUM_ROW_C = NUM_COL_A;
    const int NUM_COL_C = NUM_COL_B;
    const int span = 4;

    // Process the remainder part first.

    // Process the last odd column if present.
    if (NUM_COL_C & 1) {
        int col = NUM_COL_C - 1;
        for (int row = 0; row < NUM_ROW_C; ++row) {
            const double* pa = &A[row];
            const double* pb = &B[col];
            double tmp = 0.0;
            for (int k = 0; k < NUM_ROW_A; ++k) {
                tmp += pa[0] * pb[0];
                pa += NUM_COL_A;
                pb += NUM_COL_B;
            }

            const int index = (row + start_row_c) * col_stride_c + start_col_c + col;
            GEMM_STORE_SINGLE(C, index, tmp);
        }

        // Return directly for efficiency of extremely small matrix multiply.
        if (NUM_COL_C == 1)
            return;
    }

    // Process the couple columns in remainder if present.
    if (NUM_COL_C & 2) {
        int col = NUM_COL_C & (int)(~(span - 1));
        for (int row = 0; row < NUM_ROW_C; ++row) {
            const double* pa = &A[row];
            const double* pb = &B[col];
            double tmp1 = 0.0, tmp2 = 0.0;
            for (int k = 0; k < NUM_ROW_A; ++k) {
                double av = *pa;
                tmp1 += av * pb[0];
                tmp2 += av * pb[1];
                pa += NUM_COL_A;
                pb += NUM_COL_B;
            }

            const int index = (row + start_row_c) * col_stride_c + start_col_c + col;
            GEMM_STORE_PAIR(C, index, tmp1, tmp2);
        }

        // Return directly for efficiency of extremely small matrix multiply.
        if (NUM_COL_C < span)
            return;
    }

    // Process the main part with multiples of 4.
    int col_m = NUM_COL_C & (int)(~(span - 1));
    for (int col = 0; col < col_m; col += span) {
        for (int row = 0; row < NUM_ROW_C; ++row) {
            const int index = (row + start_row_c) * col_stride_c + start_col_c + col;
            // clang-format off
            MTM_mat1x4(NUM_ROW_A, &A[row], NUM_COL_A,
                       &B[col], NUM_COL_B, &C[index], kOperation);
            // clang-format on
        }
    }
}

GEMM_BEGIN(MatrixTransposeMatrixMultiply) {


    // if (kRowA != Eigen::Dynamic && kColA != Eigen::Dynamic &&
    //         kRowB != Eigen::Dynamic && kColB != Eigen::Dynamic) {
    //     CALL_GEMM(MatrixTransposeMatrixMultiplyEigen)
    // } else
    {
        CALL_GEMM(MatrixTransposeMatrixMultiplyNaive)
    }

}

// Matrix-Vector multiplication
//
// c op A * b;
//
// where op can be +=, -=, or =.
//
// The template parameters (kRowA, kColA) allow specialization of the
// loop at compile time. If this information is not available, then
// Eigen::Dynamic should be used as the template argument.
//
// kOperation =  1  -> c += A' * b
// kOperation = -1  -> c -= A' * b
// kOperation =  0  -> c  = A' * b
template <int kRowA, int kColA, int kOperation>
inline void MatrixVectorMultiply(const double* A,
                                 const int num_row_a,
                                 const int num_col_a,
                                 const double* b,
                                 double* c) {


    const int NUM_ROW_A = (kRowA != Eigen::Dynamic ? kRowA : num_row_a);
    const int NUM_COL_A = (kColA != Eigen::Dynamic ? kColA : num_col_a);
    const int span = 4;

    // Calculate the remainder part first.

    // Process the last odd row if present.
    if (NUM_ROW_A & 1) {
        int row = NUM_ROW_A - 1;
        const double* pa = &A[row * NUM_COL_A];
        const double* pb = &b[0];
        double tmp = 0.0;
        for (int col = 0; col < NUM_COL_A; ++col)
            tmp += (*pa++) * (*pb++);
        GEMM_STORE_SINGLE(c, row, tmp);

        // Return directly for efficiency of extremely small matrix multiply.
        if (NUM_ROW_A == 1)
            return;
    }

    // Process the couple rows in remainder if present.
    if (NUM_ROW_A & 2) {
        int row = NUM_ROW_A & (int)(~(span - 1));
        const double* pa1 = &A[row * NUM_COL_A];
        const double* pa2 = pa1 + NUM_COL_A;
        const double* pb = &b[0];
        double tmp1 = 0.0, tmp2 = 0.0;
        for (int col = 0; col < NUM_COL_A; ++col) {
            double bv = *pb++;
            tmp1 += *(pa1++) * bv;
            tmp2 += *(pa2++) * bv;
        }
        GEMM_STORE_PAIR(c, row, tmp1, tmp2);

        // Return directly for efficiency of extremely small matrix multiply.
        if (NUM_ROW_A < span)
            return;
    }

    // Calculate the main part with multiples of 4.
    int row_m = NUM_ROW_A & (int)(~(span - 1));
    for (int row = 0; row < row_m; row += span) {
        // clang-format off
        MVM_mat4x1(NUM_COL_A, &A[row * NUM_COL_A], NUM_COL_A,
                   &b[0], &c[row], kOperation);
        // clang-format on
    }
}

// Similar to MatrixVectorMultiply, except that A is transposed, i.e.,
//
// c op A' * b;
template <int kRowA, int kColA, int kOperation>
inline void MatrixTransposeVectorMultiply(const double* A,
                                          const int num_row_a,
                                          const int num_col_a,
                                          const double* b,
                                          double* c) {



    const int NUM_ROW_A = (kRowA != Eigen::Dynamic ? kRowA : num_row_a);
    const int NUM_COL_A = (kColA != Eigen::Dynamic ? kColA : num_col_a);
    const int span = 4;

    // Calculate the remainder part first.

    // Process the last odd column if present.
    if (NUM_COL_A & 1) {
        int row = NUM_COL_A - 1;
        const double* pa = &A[row];
        const double* pb = &b[0];
        double tmp = 0.0;
        for (int col = 0; col < NUM_ROW_A; ++col) {
            tmp += *pa * (*pb++);
            pa += NUM_COL_A;
        }
        GEMM_STORE_SINGLE(c, row, tmp);

        // Return directly for efficiency of extremely small matrix multiply.
        if (NUM_COL_A == 1)
            return;
    }

    // Process the couple columns in remainder if present.
    if (NUM_COL_A & 2) {
        int row = NUM_COL_A & (int)(~(span - 1));
        const double* pa = &A[row];
        const double* pb = &b[0];
        double tmp1 = 0.0, tmp2 = 0.0;
        for (int col = 0; col < NUM_ROW_A; ++col) {
            // clang-format off
            double bv = *pb++;
            tmp1 += *(pa    ) * bv;
            tmp2 += *(pa + 1) * bv;
            pa += NUM_COL_A;
            // clang-format on
        }
        GEMM_STORE_PAIR(c, row, tmp1, tmp2);

        // Return directly for efficiency of extremely small matrix multiply.
        if (NUM_COL_A < span)
            return;
    }

    // Calculate the main part with multiples of 4.
    int row_m = NUM_COL_A & (int)(~(span - 1));
    for (int row = 0; row < row_m; row += span) {
        // clang-format off
        MTV_mat4x1(NUM_ROW_A, &A[row], NUM_COL_A,
                   &b[0], &c[row], kOperation);
        // clang-format on
    }

}

#undef GEMM_BEGIN
#undef GEMM_EIGEN_HEADER
#undef GEMM_NAIVE_HEADER
#undef CALL_GEMM
#undef GEMM_STORE_SINGLE
#undef GEMM_STORE_PAIR
