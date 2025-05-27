#pragma once
#include <eigen3/Eigen/Dense>
#include "../parameter/parameters.h"

class MatirxInfo {
  public:
    MatirxInfo(int idx1_, int idx2_, flower_focus::MatrixXd  m_, bool is_identity_, int b_row_, int b_col_): idx1(idx1_), idx2(idx2_), m(m_), is_identity(is_identity_),
        b_row(b_row_), b_col(b_col_) {}
    MatirxInfo(int idx1_, int idx2_, bool is_identity_, int b_row_, int b_col_): idx1(idx1_), idx2(idx2_),  is_identity(is_identity_),
        b_row(b_row_), b_col(b_col_) {}
    int idx1;
    int idx2;
    flower_focus::MatrixXd  m;
    bool is_identity = false;
    int b_row;
    int b_col;
};

//B is a sparse matrix
inline void matrix_update(flower_focus::MatrixXd& A, flower_focus::MatrixXd B, flower_focus::MatrixXd& C, int a, int b,  bool is_left, bool is_identity, int B_ROW, int B_COL) {

    if (is_left) {
        //C=A*B
        if (is_identity) {
            if (B_ROW == 1 && B_COL == 1) C.col(b) += A.col( a);
            else C.block(0, b, A.cols(), B_COL) += A.block(0, a, A.cols(), B_ROW);
        } else {
            if (B_ROW == 1 && B_COL == 1) C.col(b) += A.col( a) * B(0);
            else C.block(0, b, A.cols(), B_COL) += A.block(0, a, A.cols(), B_ROW) * B;
        }

    } else {
        if (is_identity)
            C.block(b, b, B_COL, A.rows() - b) += A.block(a, b, B_ROW, A.rows() - b);
        else
            C.block(b, b, B_COL, A.rows() - b) += B.transpose() * A.block(a, b, B_ROW, A.rows() - b);
    }
}


//B is a sparse matrix
inline void matrix_update(Eigen::MatrixXd& A, Eigen::MatrixXd B, Eigen::MatrixXd& C, int a, int b,  bool is_left, bool is_identity, int B_ROW, int B_COL) {

    if (is_left) {
        //C=A*B
        if (is_identity) {
            if (B_ROW == 1 && B_COL == 1) C.col(b) += A.col( a);
            else C.block(0, b, A.cols(), B_COL) += A.block(0, a, A.cols(), B_ROW);
        } else {
            if (B_ROW == 1 && B_COL == 1) C.col(b) += A.col( a) * B(0);
            else C.block(0, b, A.cols(), B_COL) += A.block(0, a, A.cols(), B_ROW) * B;
        }

    } else {
        if (is_identity)
            C.block(b, b, B_COL, A.rows() - b) += A.block(a, b, B_ROW, A.rows() - b);
        else
            C.block(b, b, B_COL, A.rows() - b) += B.transpose() * A.block(a, b, B_ROW, A.rows() - b);
    }
}

inline void vector_update(Eigen::VectorXd& A, flower_focus::MatrixXd B, Eigen::VectorXd& C, int a, int b,  bool is_left, bool is_identity, int B_ROW, int B_COL) {

    if (is_left) {
        // C=B*A
        if (is_identity)
            C.segment(a, B_ROW) += A.segment(b, B_COL);
        else
            C.segment(a,  B_ROW) += B * A.segment(b, B_COL);
    } else {
        //C=B.transpose()*A;
        if (is_identity)
            C.segment(b, B_COL) += A.segment(a, B_ROW);
        else
            C.segment(b,  B_COL) += B.transpose() * A.segment(a, B_ROW);
    }
}