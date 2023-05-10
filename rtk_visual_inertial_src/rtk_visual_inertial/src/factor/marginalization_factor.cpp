

#include "marginalization_factor.h"
#include<assert.h>
#include "../parameter/parameters.h"
#include "../utility/utility.h"
#include "../utility/tic_toc.h"

void ResidualBlockInfo::Evaluate() {
    residuals.resize(cost_function->num_residuals());

    std::vector<int> block_sizes = cost_function->parameter_block_sizes();
    raw_jacobians = new double *[block_sizes.size()];
    jacobians.resize(block_sizes.size());

    for (int i = 0; i < static_cast<int>(block_sizes.size()); i++) {
        jacobians[i].resize(cost_function->num_residuals(), block_sizes[i]);
        raw_jacobians[i] = jacobians[i].data();
    }
    cost_function->Evaluate(parameter_blocks.data(), residuals.data(), raw_jacobians);


    if (loss_function) {
        double residual_scaling_, alpha_sq_norm_;

        double sq_norm, rho[3];

        sq_norm = residuals.squaredNorm();
        loss_function->Evaluate(sq_norm, rho);
        double sqrt_rho1_ = sqrt(rho[1]);
        if ((sq_norm == 0.0) || (rho[2] <= 0.0)) {
            residual_scaling_ = sqrt_rho1_;
            alpha_sq_norm_ = 0.0;
        } else {
            const double D = 1.0 + 2.0 * sq_norm * rho[2] / rho[1];
            const double alpha = 1.0 - sqrt(D);
            residual_scaling_ = sqrt_rho1_ / (1 - alpha);
            alpha_sq_norm_ = alpha / sq_norm;
        }
        for (int i = 0; i < static_cast<int>(parameter_blocks.size()); i++) {
            jacobians[i] = sqrt_rho1_ * (jacobians[i] - alpha_sq_norm_ * residuals * (residuals.transpose() * jacobians[i]));
        }

        residuals *= residual_scaling_;
    }
}


MarginalizationInfo::~MarginalizationInfo() {

    for (auto it = parameter_block_data.begin(); it != parameter_block_data.end(); ++it)
        delete it->second;
    if (m == 0 || n == 0)return;
    for (int i = 0; i < (int)factors.size(); i++) {

        delete[] factors[i]->raw_jacobians;

        delete factors[i]->cost_function;

        delete factors[i];
    }
}

void MarginalizationInfo::addResidualBlockInfo(ResidualBlockInfo* residual_block_info) {
    factors.emplace_back(residual_block_info);

    std::vector<double*>& parameter_blocks = residual_block_info->parameter_blocks;
    std::vector<int> parameter_block_sizes = residual_block_info->cost_function->parameter_block_sizes();

    for (int i = 0; i < static_cast<int>(residual_block_info->parameter_blocks.size()); i++) {
        double* addr = parameter_blocks[i];
        int size = parameter_block_sizes[i];
        parameter_block_size[reinterpret_cast<long>(addr)] = size;
    }

    for (int i = 0; i < static_cast<int>(residual_block_info->drop_set.size()); i++) {
        double* addr = parameter_blocks[residual_block_info->drop_set[i]];
        parameter_block_drop_idx[reinterpret_cast<long>(addr)] = 0;
    }
    for (int i = 0; i < static_cast<int>(residual_block_info->sparse_set.size()); i++) {
        double* addr = parameter_blocks[residual_block_info->sparse_set[i]];
        parameter_block_sparse_idx[reinterpret_cast<long>(addr)] = 0;
    }

}



int MarginalizationInfo::localSize(int size) const {
    return size == 7 ? 6 : size;
}

int MarginalizationInfo::globalSize(int size) const {
    return size == 6 ? 7 : size;
}

void* ThreadsConstructA(void* threadsstruct) {
    ThreadsStruct* p = ((ThreadsStruct*)threadsstruct);
    for (auto it : p->sub_factors) {
        it->Evaluate();
        for (int i = 0; i < static_cast<int>(it->parameter_blocks.size()); i++) {
            int idx_i = p->parameter_block_idx[reinterpret_cast<long>(it->parameter_blocks[i])];
            int size_i = p->parameter_block_size[reinterpret_cast<long>(it->parameter_blocks[i])];
            if (size_i == 7)size_i = 6;
            Eigen::MatrixXd jacobian_i = it->jacobians[i].leftCols(size_i);
            for (int j = i; j < static_cast<int>(it->parameter_blocks.size()); j++) {
                int idx_j = p->parameter_block_idx[reinterpret_cast<long>(it->parameter_blocks[j])];
                int size_j = p->parameter_block_size[reinterpret_cast<long>(it->parameter_blocks[j])];
                if (size_j == 7)size_j = 6;
                Eigen::MatrixXd jacobian_j = it->jacobians[j].leftCols(size_j);
                if (i == j)p->A.block(idx_i, idx_j, size_i, size_j) += jacobian_i.transpose() * jacobian_j;
                else {
                    p->A.block(idx_i, idx_j, size_i, size_j) += jacobian_i.transpose() * jacobian_j;
                    p->A.block(idx_j, idx_i, size_j, size_i) = p->A.block(idx_i, idx_j, size_i, size_j).transpose();
                }
            }
            p->b.segment(idx_i, size_i) += jacobian_i.transpose() * it->residuals;
        }
    }
    return threadsstruct;
}

void* ThreadsConstructAwithoutb(void* threadsstruct) {
    ThreadsStruct* p = ((ThreadsStruct*)threadsstruct);
    for (auto it : p->sub_factors) {
        it->Evaluate();
        for (int i = 0; i < static_cast<int>(it->parameter_blocks.size()); i++) {
            int idx_i = p->parameter_block_idx[reinterpret_cast<long>(it->parameter_blocks[i])];
            int size_i = p->parameter_block_size[reinterpret_cast<long>(it->parameter_blocks[i])];
            if (size_i == 7)
                size_i = 6;
            Eigen::MatrixXd jacobian_i = it->jacobians[i].leftCols(size_i);
            for (int j = i; j < static_cast<int>(it->parameter_blocks.size()); j++) {
                int idx_j = p->parameter_block_idx[reinterpret_cast<long>(it->parameter_blocks[j])];
                int size_j = p->parameter_block_size[reinterpret_cast<long>(it->parameter_blocks[j])];
                if (size_j == 7)
                    size_j = 6;
                Eigen::MatrixXd jacobian_j = it->jacobians[j].leftCols(size_j);
                if (i == j)
                    p->A.block(idx_i, idx_j, size_i, size_j) += jacobian_i.transpose() * jacobian_j;
                else {
                    p->A.block(idx_i, idx_j, size_i, size_j) += jacobian_i.transpose() * jacobian_j;
                    p->A.block(idx_j, idx_i, size_j, size_i) = p->A.block(idx_i, idx_j, size_i, size_j).transpose();
                }
            }
        }
    }
    return threadsstruct;
}



void MarginalizationInfo::marginalize_exclude() {
    TicToc t_marg;

    int pos = 0;
    int m2 = 0;
    for (auto& it : parameter_block_sparse_idx) {
        if (parameter_block_idx.find(it.first) == parameter_block_idx.end()) {
            parameter_block_idx[it.first] = pos;
            pos += localSize(parameter_block_size[it.first]);
        }
    }
    m2 = pos;
    for (auto& it : parameter_block_drop_idx) {
        if (parameter_block_idx.find(it.first) == parameter_block_idx.end()) {
            parameter_block_idx[it.first] = pos;
            pos += localSize(parameter_block_size[it.first]);
        }
    }
    m = pos - m2;

    for (auto& it : parameter_block_size) {
        if (parameter_block_idx.find(it.first) == parameter_block_idx.end()) {
            parameter_block_idx[it.first] = pos;
            pos += localSize(it.second);
        }
    }
    n = pos - m2 - m;
    assert(m && n);
    A = Eigen::MatrixXd(pos, pos);
    A.setZero();
    TicToc t_thread_summing;
    pthread_t tids[NUM_THREADS];
    ThreadsStruct threadsstruct[NUM_THREADS];
    int i = 0;
    for (auto it : factors) {
        threadsstruct[i].sub_factors.push_back(it);
        i++;
        i = i % NUM_THREADS;
    }
    for (int i = 0; i < NUM_THREADS; i++) {
        TicToc zero_matrix;
        threadsstruct[i].A = Eigen::MatrixXd::Zero(pos, pos);
        threadsstruct[i].parameter_block_size = parameter_block_size;
        threadsstruct[i].parameter_block_idx = parameter_block_idx;
        int ret = pthread_create( &tids[i], NULL, ThreadsConstructAwithoutb, (void*) & (threadsstruct[i]));
        if (ret != 0) {
            printf("pthread_create error");
            abort();
        }
    }
    for ( int i = NUM_THREADS - 1; i >= 0; i--) {
        pthread_join( tids[i], NULL );
        A += threadsstruct[i].A;

    }

    Eigen::MatrixXd Am2m2 = A.block(0, 0, m2, m2) ;
    Eigen::MatrixXd Am2m2_inv = Eigen::VectorXd(Am2m2.diagonal().array().inverse()).asDiagonal();
    Eigen::MatrixXd Am2mn = A.block(0, m2, m2, n + m);
    Eigen::MatrixXd Amnm2 = A.block(m2, 0, n + m, m2);
    Eigen::MatrixXd Amnmn = A.block(m2, m2, n + m, n + m);
    A = Amnmn - Amnm2 * Am2m2_inv * Am2mn;
    Eigen::MatrixXd Amm = A.block(0, 0, m, m);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes2(Amm);
    Eigen::MatrixXd Amm_inv = saes2.eigenvectors() * Eigen::VectorXd((saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array().inverse(), 0)).asDiagonal() * saes2.eigenvectors().transpose();
    Eigen::MatrixXd Amr = A.block(0, m, m, n);
    Eigen::MatrixXd Arm = A.block(m, 0, n, m);
    Eigen::MatrixXd Arr = A.block(m, m, n, n);

    A = Arr - Arm * Amm_inv * Amr;


    // std::cout<<"marginalization time:"<<t_marg.toc()<<",";

    m += m2;
}



void MarginalizationInfo::ResetLinearizationPoint(std::vector<double*>parameters) {
    Eigen::VectorXd dx(n);
    for (int i = 0; i < static_cast<int>(keep_block_size.size()); i++) {
        int size = keep_block_size[i];
        int idx = keep_block_idx[i] - m;
        Eigen::VectorXd x = Eigen::Map<const Eigen::VectorXd>(parameters[i], size);
        Eigen::VectorXd x0 = Eigen::Map<const Eigen::VectorXd>(keep_block_data[i], size);
        if (size != 7)
            dx.segment(idx, size) = x - x0;
        else {
            dx.segment<3>(idx + 0) = x.head<3>() - x0.head<3>();
            dx.segment<3>(idx + 3) = 2.0 * (Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
            if (!((Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).w() >= 0)) {
                dx.segment<3>(idx + 3) = 2.0 * -(Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
            }
        }
        if (size == 1)assert(parameters[i][0] == 0);
        Eigen::Map<Eigen::VectorXd>(keep_block_data[i], size) = Eigen::Map<const Eigen::VectorXd>(parameters[i], size);
    }
    // std::cout<<( linearized_jacobians * dx).transpose()<<std::endl;
    linearized_residuals += linearized_jacobians * dx;
    b += A * dx;



}


void MarginalizationInfo::marginalize(bool initialinformation, bool multithread) {
    multithread = false;
    int pos = 0;

    for (auto& it : parameter_block_drop_idx) {
        if (parameter_block_idx.find(it.first) == parameter_block_idx.end()) {
            parameter_block_idx[it.first] = pos;
            pos += localSize(parameter_block_size[it.first]);
        }
    }
    m = pos;

    for (auto& it : parameter_block_size) {
        if (parameter_block_idx.find(it.first) == parameter_block_idx.end()) {
            parameter_block_idx[it.first] = pos;
            pos += localSize(it.second);
        }
    }

    n = pos - m;
    if ((pos == 0) && !initialinformation) {
        printf("unstable tracking...\n");
        return;
    }
    if (n == 0)return;

    TicToc t_summing;
    A = Eigen::MatrixXd(pos, pos);
    b = Eigen::VectorXd(pos);
    A.setZero();
    b.setZero();
    if (multithread) {
        pthread_t tids[NUM_THREADS];
        ThreadsStruct threadsstruct[NUM_THREADS];
        int i = 0;
        for (auto it : factors) {
            threadsstruct[i].sub_factors.push_back(it);
            i++;
            i = i % NUM_THREADS;
        }
        for (int i = 0; i < NUM_THREADS; i++) {
            TicToc zero_matrix;
            threadsstruct[i].A = Eigen::MatrixXd::Zero(pos, pos);
            threadsstruct[i].b = Eigen::VectorXd::Zero(pos);
            threadsstruct[i].parameter_block_size = parameter_block_size;
            threadsstruct[i].parameter_block_idx = parameter_block_idx;
            int ret = pthread_create( &tids[i], NULL, ThreadsConstructA, (void*) & (threadsstruct[i]));
            if (ret != 0) {
                printf("pthread_create error");
                abort();
            }


        }
        for ( int i = NUM_THREADS - 1; i >= 0; i--) {
            pthread_join( tids[i], NULL );
            A += threadsstruct[i].A;
            b += threadsstruct[i].b;

        }
    } else {
        ThreadsStruct threadsstruct;
        for (auto it : factors)threadsstruct.sub_factors.push_back(it);
        threadsstruct.A = Eigen::MatrixXd::Zero(pos, pos);
        threadsstruct.b = Eigen::VectorXd::Zero(pos);
        threadsstruct.parameter_block_size = parameter_block_size;
        threadsstruct.parameter_block_idx = parameter_block_idx;
        ThreadsConstructA((void*) & (threadsstruct));
        A = threadsstruct.A;
        b = threadsstruct.b;
    }
    // std::cout<<"A0:"<<A<<std::endl;
    // std::cout<<"b0:"<<b<<std::endl;


    if (m != 0) {
        Eigen::MatrixXd Amm = 0.5 * (A.block(0, 0, m, m) + A.block(0, 0, m, m).transpose());
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(Amm);
        Eigen::MatrixXd Amm_inv = saes.eigenvectors() * Eigen::VectorXd((saes.eigenvalues().array() > eps).select(saes.eigenvalues().array().inverse(), 0)).asDiagonal() * saes.eigenvectors().transpose();
        Eigen::VectorXd bmm = b.segment(0, m);
        Eigen::MatrixXd Amr = A.block(0, m, m, n);
        Eigen::MatrixXd Arm = A.block(m, 0, n, m);
        Eigen::MatrixXd Arr = A.block(m, m, n, n);
        Eigen::VectorXd brr = b.segment(m, n);

        A = Arr - Arm * Amm_inv * Amr;
        b = brr - Arm * Amm_inv * bmm;


    }
    // std::cout<<"A:"<<A<<std::endl;
    // std::cout<<"b:"<<b<<std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes2(A);
    Eigen::VectorXd S = Eigen::VectorXd((saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array(), 0));
    Eigen::VectorXd S_sqrt = S.cwiseSqrt();
    linearized_jacobians = S_sqrt.asDiagonal() * saes2.eigenvectors().transpose();
    Eigen::VectorXd delta;
    {
        Eigen::VectorXd S_inv = Eigen::VectorXd((saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array().inverse(), 0));
        Eigen::VectorXd S_inv_sqrt = S_inv.cwiseSqrt();
        linearized_residuals = S_inv_sqrt.asDiagonal() * saes2.eigenvectors().transpose() * b;
    }

    for (auto it : factors) {
        std::vector<int> block_sizes = it->cost_function->parameter_block_sizes();
        for (int i = 0; i < static_cast<int>(block_sizes.size()); i++) {
            long addr = reinterpret_cast<long>(it->parameter_blocks[i]);
            int size = block_sizes[i];
            if (parameter_block_data.find(addr) == parameter_block_data.end()) {
                double* data = new double[size];
                memcpy(data, it->parameter_blocks[i], sizeof(double) * size);
                parameter_block_data[addr] = data;
            }
        }
    }


}



std::vector<double*> MarginalizationInfo::getParameterBlocks() {
    std::vector<double*> keep_block_addr2;
    keep_block_size.clear();
    keep_block_idx.clear();
    keep_block_data.clear();

    for (const auto& it : parameter_block_idx) {
        if (it.second >= m) {
            keep_block_size.push_back(parameter_block_size[it.first]);
            keep_block_idx.push_back(parameter_block_idx[it.first]);
            keep_block_data.push_back(parameter_block_data[it.first]);
            keep_block_addr2.push_back(reinterpret_cast<double*>(it.first));
        }
    }
    sum_block_size = std::accumulate(std::begin(keep_block_size), std::end(keep_block_size), 0);
    keep_block_addr = keep_block_addr2;
    return keep_block_addr;
}


MarginalizationFactor::MarginalizationFactor(MarginalizationInfo* _marginalization_info): marginalization_info(_marginalization_info) {
    int cnt = 0;
    for (auto it : marginalization_info->keep_block_size) {
        mutable_parameter_block_sizes()->push_back(it);
        cnt += it;
    }
    set_num_residuals(marginalization_info->n);
};

bool MarginalizationFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {

    int n = marginalization_info->n;
    int m = marginalization_info->m;
    Eigen::VectorXd dx(n);
    for (int i = 0; i < static_cast<int>(marginalization_info->keep_block_size.size()); i++) {
        int size = marginalization_info->keep_block_size[i];
        int idx = marginalization_info->keep_block_idx[i] - m;
        Eigen::VectorXd x = Eigen::Map<const Eigen::VectorXd>(parameters[i], size);
        Eigen::VectorXd x0 = Eigen::Map<const Eigen::VectorXd>(marginalization_info->keep_block_data[i], size);
        if (size != 7)
            dx.segment(idx, size) = x - x0;
        else {
            dx.segment<3>(idx + 0) = x.head<3>() - x0.head<3>();
            dx.segment<3>(idx + 3) = 2.0 * (Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
            if (!((Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).w() >= 0)) {
                dx.segment<3>(idx + 3) = 2.0 * -(Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
            }
        }
    }
    Eigen::Map<Eigen::VectorXd>(residuals, n) = marginalization_info->linearized_residuals + marginalization_info->linearized_jacobians * dx;


    if (jacobians) {

        for (int i = 0; i < static_cast<int>(marginalization_info->keep_block_size.size()); i++) {
            if (jacobians[i]) {
                int size = marginalization_info->keep_block_size[i], local_size = marginalization_info->localSize(size);
                int idx = marginalization_info->keep_block_idx[i] - m;
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> jacobian(jacobians[i], n, size);
                jacobian.setZero();
                jacobian.leftCols(local_size) = marginalization_info->linearized_jacobians.middleCols(idx, local_size);
            }
        }
    }
    return true;
}


void MarginalizationInfo::setmarginalizeinfo(std::vector<const double*>& parameter_block_addr_,
        std::vector<int>& parameter_block_global_size_, Eigen::MatrixXd& A_, Eigen::VectorXd& b_,
        bool Sqrt) {
    n = m = 1;
    A = A_;
    b = b_;
    Eigen::VectorXd delta;

    if (Sqrt) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes2(A);
        Eigen::VectorXd S = Eigen::VectorXd((saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array(), 0));
        Eigen::VectorXd S_inv = Eigen::VectorXd((saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array().inverse(), 0));
        Eigen::VectorXd S_sqrt = S.cwiseSqrt();
        Eigen::VectorXd S_inv_sqrt = S_inv.cwiseSqrt();
        linearized_jacobians = S_sqrt.asDiagonal() * saes2.eigenvectors().transpose();

        {
            Eigen::VectorXd S_inv = Eigen::VectorXd((saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array().inverse(), 0));
            Eigen::VectorXd S_inv_sqrt = S_inv.cwiseSqrt();
            linearized_residuals = S_inv_sqrt.asDiagonal() * saes2.eigenvectors().transpose() * b;
        }
    }

    for (int i = 0; i < (int)parameter_block_addr_.size(); i++) {
        const double* addr = parameter_block_addr_[i];
        int globalsize = parameter_block_global_size_[i];
        int localsize = localSize(globalsize);
        parameter_block_size[reinterpret_cast<long>(addr)] = globalsize;
        parameter_block_idx[reinterpret_cast<long>(addr)] = n;
        n += localsize;
        double* data = new double[globalsize];
        memcpy(data, addr, sizeof(double) * globalsize);
        parameter_block_data[reinterpret_cast<long>(addr)] = data;
    }
    n -= 1;
    assert(n == A.rows());



}


