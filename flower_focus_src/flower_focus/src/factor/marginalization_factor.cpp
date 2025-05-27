

#include "marginalization_factor.h"
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
        for (int i = 0; i < static_cast<int>(parameter_blocks.size()); i++)
            jacobians[i] = sqrt_rho1_ * (jacobians[i] - alpha_sq_norm_ * residuals * (residuals.transpose() * jacobians[i]));

        residuals *= residual_scaling_;
    }
}


MarginalizationInfo::~MarginalizationInfo() {

    for (auto it = parameter_block_data.begin(); it != parameter_block_data.end(); ++it)
        delete it->second;
    if (m == 0 || n == 0)
        return;
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
        if (parameter_block_size.find(addr) == parameter_block_size.end()) {
            parameter_block_size[addr] = size;
            parameters.push_back(addr);
        }

    }

    for (int i = 0; i < static_cast<int>(residual_block_info->drop_set.size()); i++) {
        double* addr = parameter_blocks[residual_block_info->drop_set[i]];
        if (parameter_block_drop_idx.find(addr) == parameter_block_drop_idx.end()) {
            parameter_block_drop_idx[addr] = 0;
            drop_parameters.push_back(addr);
        }
    }


}




void* ThreadsConstructA(void* threadsstruct) {
    ThreadsStruct* p = ((ThreadsStruct*)threadsstruct);
    for (auto it : p->sub_factors) {
        it->Evaluate();
        for (int i = 0; i < static_cast<int>(it->parameter_blocks.size()); i++) {
            int idx_i = p->parameter_block_idx[it->parameter_blocks[i]];
            int size_i = p->parameter_block_size[it->parameter_blocks[i]];
            if (size_i == 7)
                size_i = 6;
            Eigen::MatrixXd jacobian_i = it->jacobians[i].leftCols(size_i);
            for (int j = i; j < static_cast<int>(it->parameter_blocks.size()); j++) {
                int idx_j = p->parameter_block_idx[it->parameter_blocks[j]];
                int size_j = p->parameter_block_size[it->parameter_blocks[j]];
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
            p->b.segment(idx_i, size_i) += jacobian_i.transpose() * it->residuals;
        }
    }
    return threadsstruct;
}


void MarginalizationInfo::marginalize(bool initialinformation) {

    int pos = 0;
    for (int i = 0; i < (int)drop_parameters.size(); i++) {
        if (parameter_block_idx.find(drop_parameters[i]) == parameter_block_idx.end()) {
            parameter_block_idx[drop_parameters[i]] = pos;
            pos += localSize(parameter_block_size[drop_parameters[i]]);
        }
    }
    m = pos;

    for (int i = 0; i < (int)parameters.size(); i++) {
        if (parameter_block_idx.find(parameters[i]) == parameter_block_idx.end()) {
            parameter_block_idx[parameters[i]] = pos;
            pos += localSize(parameter_block_size[parameters[i]]);
        }
    }

    n = pos - m;
    if ((pos == 0) && !initialinformation) {
        printf("unstable tracking...\n");
        return;
    }
    if (n == 0)
        return;

    TicToc t_summing;
    A = Eigen::MatrixXd(pos, pos);
    b = Eigen::VectorXd(pos);
    A.setZero();
    b.setZero();

    ThreadsStruct threadsstruct;
    for (auto it : factors) threadsstruct.sub_factors.push_back(it);
    threadsstruct.A = Eigen::MatrixXd::Zero(pos, pos);
    threadsstruct.b = Eigen::VectorXd::Zero(pos);
    threadsstruct.parameter_block_size = parameter_block_size;
    threadsstruct.parameter_block_idx = parameter_block_idx;
    ThreadsConstructA((void*) & (threadsstruct));
    A = threadsstruct.A;
    b = threadsstruct.b;


    if (m != 0) {
        Eigen::MatrixXd Amm = 0.5 * (A.block(0, 0, m, m) + A.block(0, 0, m, m).transpose());
        Eigen::MatrixXd Amm_inv;
        if (0) {
            Eigen::VectorXd diag = Eigen::VectorXd(Amm.cols());
            for (int i = 0; i < diag.size(); i++)
                diag(i) = EPSS;
            Amm.diagonal() += diag;
            Eigen::LDLT<Eigen::MatrixXd, Eigen::Upper>ldlt = Amm.selfadjointView<Eigen::Upper>().ldlt();
            Amm_inv = ldlt.solve(Eigen::MatrixXd::Identity(m, m));
#if USE_ASSERT
            Eigen::Index a1, a2;
            ASSERT(ldlt.vectorD().minCoeff(&a1, &a2) > -0.1);
#endif
        } else {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(Amm);
            Amm_inv = saes.eigenvectors()
                      * Eigen::VectorXd((saes.eigenvalues().array() > EPSS).select(saes.eigenvalues().array().inverse(), 0)).asDiagonal()
                      * saes.eigenvectors().transpose();
        }
        Eigen::VectorXd bmm = b.segment(0, m);
        Eigen::MatrixXd Amr = A.block(0, m, m, n);
        Eigen::MatrixXd Arm = A.block(m, 0, n, m);
        Eigen::MatrixXd Arr = A.block(m, m, n, n);
        Eigen::VectorXd brr = b.segment(m, n);

        A = Arr - Arm * Amm_inv * Amr;
        b = brr - Arm * Amm_inv * bmm;
    }
    for (auto it : factors) {
        std::vector<int> block_sizes = it->cost_function->parameter_block_sizes();
        for (int i = 0; i < static_cast<int>(block_sizes.size()); i++) {
            double* addr = it->parameter_blocks[i];
            int size = block_sizes[i];
            if (parameter_block_data.find(addr) == parameter_block_data.end()) {
                double* data = new double[size];
                memcpy(data, it->parameter_blocks[i], sizeof(double) * size);
                parameter_block_data[addr] = data;
            }
        }
    }

    for (const auto& it : parameter_block_idx)
        parameter_block_idx[it.first] -= m;
    // n -= m;
    m = 0;
}


void MarginalizationInfo::marginalize_pointers(std::set<double*>feature_marge_pointer, bool need_jacobian) {



    std::vector<double*> new_addr;
    std::vector<int> new_size;
    std::vector<int>new_idx;

    Eigen::MatrixXd newA;
    Eigen::VectorXd newb;

    Eigen::VectorXd b_plus = b;

    std::unordered_map<const double*, int>addr2idx_old;

    for (int i = 0; i < (int)keep_block_addr.size(); i++)
        addr2idx_old[keep_block_addr[i]] = keep_block_idx[i];


    n = 0; m = 0;
    int pos = 0;

    std::map<int, double*, less<int>>keep_block_addr_mapping;
    std::map<int, int, less<int>>keep_block_size_mapping;

    for (int i = 0; i < (int)keep_block_addr.size(); i++) {
        int idx = keep_block_idx[i];
        keep_block_addr_mapping[idx] = keep_block_addr[i];
        keep_block_size_mapping[idx] = keep_block_size[i];
    }

    for (auto it = keep_block_addr_mapping.begin(); it != keep_block_addr_mapping.end(); it++) {
        double* addr = it->second;
        int size = keep_block_size_mapping[it->first];
        if (feature_marge_pointer.find(addr) != feature_marge_pointer.end()) {
            new_addr.push_back(addr); new_idx.push_back(pos); new_size.push_back(size);
            pos += localSize(size);
        }
    }
    m = pos;
    n = A.cols() - m;
    for (auto it = keep_block_addr_mapping.begin(); it != keep_block_addr_mapping.end(); it++) {
        double* addr = it->second;
        int size = keep_block_size_mapping[it->first];
        if (feature_marge_pointer.find(addr) == feature_marge_pointer.end()) {
            new_addr.push_back(addr); new_idx.push_back(pos); new_size.push_back(size);
            pos += localSize(size);
        }
    }
    ASSERT(pos == A.cols());

    newA = Eigen::MatrixXd::Zero(A.rows(), A.cols());
    newb = Eigen::VectorXd::Zero(b_plus.size());

    for (int i = 0; i < (int)new_addr.size(); i++) {
        const double* addri = new_addr[i];
        int new_idxi = new_idx[i];
        int old_idxi = addr2idx_old[addri];
        int sizei = localSize(new_size[i]);
        newb.segment(new_idxi, sizei) = b_plus.segment(old_idxi, sizei);
        for (int j = i; j < (int)new_addr.size(); j++) {
            const double* addrj = new_addr[j];
            int new_idxj = new_idx[j];
            int old_idxj = addr2idx_old[addrj];
            int sizej = localSize(new_size[j]);
            newA.block(new_idxi, new_idxj, sizei, sizej) = A.block(old_idxi, old_idxj, sizei, sizej);
        }
    }
    newA = newA.selfadjointView<Eigen::Upper>();

    Eigen::VectorXd bmm = newb.segment(0, m);
    Eigen::MatrixXd Amm = newA.block(0, 0, m, m);
    Eigen::MatrixXd Anm = newA.block(m, 0, n, m);
    Eigen::MatrixXd Ann = newA.block(m, m, n, n);
    Eigen::VectorXd bnn = newb.segment(m, n);
    Amm.diagonal().array() += kMinMu;
    Eigen::MatrixXd Amm_inv = Amm.template selfadjointView<Eigen::Upper>().llt().solve(Eigen::MatrixXd::Identity(Amm.rows(), Amm.rows()));


    A = Ann - Anm * Amm_inv * Anm.transpose();
    b = bnn - Anm * Amm_inv * bmm;

    new_addr.clear(); new_size.clear(); new_idx.clear(); pos = 0;
    for (auto it = keep_block_addr_mapping.begin(); it != keep_block_addr_mapping.end(); it++) {
        double* addr = it->second;
        int size = keep_block_size_mapping[it->first];
        if (feature_marge_pointer.find(addr) == feature_marge_pointer.end()) {
            new_addr.push_back(addr); new_idx.push_back(pos); new_size.push_back(size);
            pos += localSize(size);
        }
    }
    ASSERT(pos == n);

    keep_block_size = new_size;
    keep_block_idx = new_idx;
    keep_block_addr = new_addr;
    keep_block_addr_set.clear();
    keep_block_data.clear();

    for (int i = 0; i < (int)keep_block_addr.size(); i++) {
        ASSERT(parameter_block_data.find(keep_block_addr[i]) != parameter_block_data.end());
        keep_block_data.push_back(parameter_block_data[keep_block_addr[i]]);
        keep_block_addr_set[keep_block_addr[i]] = i;
        ASSERT(keep_block_data[keep_block_data.size() - 1] != 0);
    }
    if (need_jacobian)
        GetJacobianAndResidual();

    m = 0;

}


void MarginalizationInfo::setmarginalizeinfo(std::vector<double*>& parameter_block_addr_,
                                             std::vector<int>& parameter_block_global_size_, Eigen::MatrixXd A_, Eigen::VectorXd& b_,
                                             bool Sqrt) {
    n = m = 0;
    A = A_;
    b = b_;


    for (int i = 0; i < (int)parameter_block_addr_.size(); i++) {
        double* addr = parameter_block_addr_[i];
        int globalsize = parameter_block_global_size_[i];
        int localsize = localSize(globalsize);
        parameter_block_size[addr] = globalsize;
        parameter_block_idx[addr] = n;
        n += localsize;
        double* data = new double[globalsize];
        memcpy(data, addr, sizeof(double) * globalsize);
        parameter_block_data[addr] = data;
    }
    ASSERT(n == A.rows());
    ASSERT(m == 0);



}




void MarginalizationInfo::BuildMatrixF(std::vector<double*>pointers) {


    double* x0 = 0;
    double scale_factor = 1;
    if (fix_scale)scale_factor = pointers[scale_factor_i][0];

    x0 = pointers[v_i];
    Eigen::Vector3d v_w_b0 = Eigen::Vector3d(x0[0], x0[1], x0[2]);

    Eigen::Vector3d Bg0 = Eigen::Vector3d(x0[6], x0[7], x0[8]);

    Eigen::Matrix3d R_WI_WC;
    if (ESTIMATE_R_WI_WC) {
        x0 = pointers[R_WI_WC_i];
        R_WI_WC = Eigen::Quaterniond(x0[6], x0[3], x0[4], x0[5]).normalized().toRotationMatrix();
    } else
        R_WI_WC.setIdentity();


    std::vector<Eigen::Vector3d> P_w_ck;
    std::vector<Eigen::Matrix3d> R_w_ck;
    P_w_ck.resize(p_is.size());
    R_w_ck.resize(p_is.size());
    for (int i = 0; i < (int)p_is.size(); i++) {
        x0 = pointers[p_is[i]];
        P_w_ck[i] = Eigen::Vector3d(x0[0], x0[1], x0[2]);
        R_w_ck[i] = Eigen::Quaterniond(x0[6], x0[3], x0[4], x0[5]).toRotationMatrix();
    }

    double sum_distance = 0;
    for (int i = 1; i < (int)p_is.size(); i++)
        sum_distance +=  (P_w_ck[i] - P_w_ck[0]).norm();

    sum_distance *= scale_factor;
    if (!fix)sum_distance = 1;


    std::vector<double>idepth0;
    idepth0.resize(lambda_is.size());
    for (int i = 0; i < (int)lambda_is.size(); i++) {
        x0 = pointers[lambda_is[i]];
        idepth0[i] = x0[0];
    }

    Eigen::Matrix3d R_b_c;
    Eigen::Vector3d P_b_c;
    if (ESTIMATE_EXTRINSIC) {
        x0 = pointers[extrinsic_i];
        P_b_c = Eigen::Vector3d(x0[0], x0[1], x0[2]);
        R_b_c = Eigen::Quaterniond(x0[6], x0[3], x0[4], x0[5]).toRotationMatrix();
    } else {
        R_b_c = RIC[0];
        P_b_c = TIC[0];
    }





    // Eigen::Vector3d v_c0_c0 = (R_w_ck[0].transpose() * R_WI_WC.transpose() *  v_w_b0 + R_b_c.transpose() * Utility::skewSymmetric(gyr0 - Bg0) * P_b_c) / sum_distance;


    {
        matrix_infoF.clear();

        matrix_infoF.push_back(MatirxInfo(keep_block_idx[p_is[0]], keep_block_idx[p_is[0]], R_WI_WC.transpose() / scale_factor - Eigen::Matrix3d::Identity(), false, 3, 3));
        if (fix_scale)
            matrix_infoF.push_back(MatirxInfo(keep_block_idx[p_is[0]], keep_block_idx[scale_factor_i], -P_w_ck[0] / scale_factor, false, 3, 1));
        if (ESTIMATE_R_WI_WC)
            matrix_infoF.push_back(MatirxInfo(keep_block_idx[p_is[0]], keep_block_idx[R_WI_WC_i] + 3, Utility::skewSymmetric(P_w_ck[0])*R_WI_WC.transpose(), false, 3, 3));

        matrix_infoF.push_back(MatirxInfo(keep_block_idx[p_is[0]] + 3, keep_block_idx[p_is[0]] + 3, R_w_ck[0].transpose()*R_WI_WC.transpose() - Eigen::Matrix3d::Identity(), false, 3, 3));
        if (ESTIMATE_R_WI_WC)
            matrix_infoF.push_back(MatirxInfo(keep_block_idx[p_is[0]] + 3, keep_block_idx[R_WI_WC_i] + 3, -R_w_ck[0].transpose()*R_WI_WC.transpose(), false, 3, 3));

        matrix_infoF.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[p_is[0]] + 3, -Utility::skewSymmetric(v_w_b0), false, 3, 3));
        matrix_infoF.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[v_i], R_WI_WC * R_w_ck[0] - Eigen::Matrix3d::Identity(), false, 3, 3));

        if (ESTIMATE_EXTRINSIC) {
            matrix_infoF.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[extrinsic_i], -R_WI_WC * R_w_ck[0] * R_b_c.transpose() * Utility::skewSymmetric(gyr0 - Bg0), false, 3, 3));
            matrix_infoF.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[extrinsic_i] + 3, -R_WI_WC * R_w_ck[0] * Utility::skewSymmetric(R_b_c.transpose() * Utility::skewSymmetric(gyr0 - Bg0) * P_b_c), false, 3, 3));
        }

        // matrix_infoF.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[p_is[0]] + 3, -Utility::skewSymmetric(v_w_b0), false, 3, 3));
        // matrix_infoF.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[v_i], R_WI_WC * R_w_ck[0] * sum_distance - Eigen::Matrix3d::Identity(), false, 3, 3));

        // if (fix)
        //     matrix_infoF.push_back(MatirxInfo(keep_block_idx[v_i], A.cols(), R_WI_WC * R_w_ck[0] * v_c0_c0, false, 3, 1));

        // if (ESTIMATE_EXTRINSIC) {
        //     matrix_infoF.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[extrinsic_i], -R_WI_WC * R_w_ck[0] * R_b_c.transpose() * Utility::skewSymmetric(gyr0 - Bg0), false, 3, 3));
        //     matrix_infoF.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[extrinsic_i] + 3, -R_WI_WC * R_w_ck[0] * Utility::skewSymmetric(R_b_c.transpose() * Utility::skewSymmetric(gyr0 - Bg0) * P_b_c), false, 3, 3));
        // }
        // matrix_infoF.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[v_i] + 6, -R_WI_WC * R_w_ck[0] * R_b_c.transpose() * Utility::skewSymmetric(P_b_c), false, 3, 3));

        if (ESTIMATE_R_WI_WC)
            matrix_infoF.push_back(MatirxInfo(keep_block_idx[R_WI_WC_i] + 3, keep_block_idx[R_WI_WC_i] + 3, R_WI_WC.transpose() - Eigen::Matrix3d::Identity(), false, 3, 3));

        // lambda_l
        for (int j = 0; j < (int)lambda_is.size(); j++) {
            matrix_infoF.push_back(MatirxInfo(keep_block_idx[lambda_is[j]], keep_block_idx[lambda_is[j]], flower_focus::Vector1d(scale_factor / sum_distance - 1), false, 1, 1));
            if (fix_scale)
                matrix_infoF.push_back(MatirxInfo(keep_block_idx[lambda_is[j]], keep_block_idx[scale_factor_i], flower_focus::Vector1d(idepth0[j] / scale_factor), false, 1, 1));
            if (fix)
                matrix_infoF.push_back(MatirxInfo(keep_block_idx[lambda_is[j]], A.cols(), flower_focus::Vector1d(-idepth0[j] / sum_distance), false, 1, 1));
        }

    }




}

void MarginalizationInfo::Getestimate(Eigen::Vector3d& v_c0_c0, Eigen::Vector3d& P_I_ck0, Eigen::Matrix3d& R_I_ck0, double& sum_distance, std::vector<double>& idepth0, Eigen::Matrix3d& R_WI_WC, std::vector<double*>pointers) {




    double* x0 = 0;


    if (ESTIMATE_R_WI_WC) {
        x0 = pointers[R_WI_WC_i];
        R_WI_WC = Eigen::Quaterniond(x0[6], x0[3], x0[4], x0[5]).normalized().toRotationMatrix();
    } else
        R_WI_WC.setIdentity();

    double scale_factor = 1;
    if (fix_scale)scale_factor = pointers[scale_factor_i][0];

    std::vector<Eigen::Vector3d> P_I_ck;
    std::vector<Eigen::Matrix3d> R_I_ck;

    P_I_ck.resize(p_is.size());
    R_I_ck.resize(p_is.size());
    for (int i = 0; i < (int)p_is.size(); i++) {
        x0 = pointers[p_is[i]];
        P_I_ck[i] = R_WI_WC * Eigen::Vector3d(x0[0], x0[1], x0[2]) * scale_factor;
        R_I_ck[i] = R_WI_WC * Eigen::Quaterniond(x0[6], x0[3], x0[4], x0[5]).toRotationMatrix();
    }
    P_I_ck0 = P_I_ck[0];
    R_I_ck0 = R_I_ck[0];


    sum_distance = 0;
    for (int i = 1; i < (int)p_is.size(); i++)
        sum_distance += (P_I_ck[i] - P_I_ck[0]).norm();
    sum_distance *= scale_factor;
    if (!fix)sum_distance = 1;

    idepth0.resize(lambda_is.size());
    for (int i = 0; i < (int)lambda_is.size(); i++) {
        x0 = pointers[lambda_is[i]];
        idepth0[i] = x0[0] * sum_distance / scale_factor;
    }

    x0 = pointers[v_i];
    Eigen::Vector3d v_w_b0 = Eigen::Vector3d(x0[0], x0[1], x0[2]);
    // Eigen::Vector3d Bg0 = Eigen::Vector3d(x0[6], x0[7], x0[8]);

    Eigen::Matrix3d R_b_c;
    Eigen::Vector3d P_b_c;
    if (ESTIMATE_EXTRINSIC) {
        x0 = pointers[extrinsic_i];
        P_b_c = Eigen::Vector3d(x0[0], x0[1], x0[2]);
        R_b_c = Eigen::Quaterniond(x0[6], x0[3], x0[4], x0[5]).toRotationMatrix();
    } else {
        R_b_c = RIC[0];
        P_b_c = TIC[0];
    }

    // v_c0_c0 = (R_I_ck[0].transpose() * v_w_b0 + R_b_c.transpose() * Utility::skewSymmetric(gyr0 - Bg0) * P_b_c) / sum_distance;
    v_c0_c0 = R_I_ck[0].transpose() * v_w_b0 ;

}

void MarginalizationInfo::BuildMatrixG(std::vector<double*>pointers, Eigen::VectorXd dx, Eigen::VectorXd& dy) {


    double* x0 = 0;

    double scale_factor = 1;
    if (fix_scale)scale_factor = pointers[scale_factor_i][0];

    x0 = pointers[v_i];
    Eigen::Vector3d v_w_b0 = Eigen::Vector3d(x0[0], x0[1], x0[2]);

    Eigen::Vector3d Bg0 = Eigen::Vector3d(x0[6], x0[7], x0[8]);

    Eigen::Matrix3d R_WI_WC;
    if (ESTIMATE_R_WI_WC) {
        x0 = pointers[R_WI_WC_i];
        R_WI_WC = Eigen::Quaterniond(x0[6], x0[3], x0[4], x0[5]).normalized().toRotationMatrix();
    } else
        R_WI_WC.setIdentity();

    std::vector<Eigen::Vector3d> P_w_ck;
    std::vector<Eigen::Matrix3d> R_w_ck;
    P_w_ck.resize(p_is.size());
    R_w_ck.resize(p_is.size());
    for (int i = 0; i < (int)p_is.size(); i++) {
        x0 = pointers[p_is[i]];
        P_w_ck[i] = Eigen::Vector3d(x0[0], x0[1], x0[2]);
        R_w_ck[i] = Eigen::Quaterniond(x0[6], x0[3], x0[4], x0[5]).toRotationMatrix();
    }

    double sum_distance = 0;
    for (int i = 1; i < (int)p_is.size(); i++)
        sum_distance += (P_w_ck[i] - P_w_ck[0]).norm();
    sum_distance *= scale_factor;
    if (!fix)sum_distance = 1;


    std::vector<double>idepth0;
    idepth0.resize(lambda_is.size());
    for (int i = 0; i < (int)lambda_is.size(); i++) {
        x0 = pointers[lambda_is[i]];
        idepth0[i] = x0[0];
    }

    Eigen::Matrix3d R_b_c;
    Eigen::Vector3d P_b_c;
    if (ESTIMATE_EXTRINSIC) {
        x0 = pointers[extrinsic_i];
        P_b_c = Eigen::Vector3d(x0[0], x0[1], x0[2]);
        R_b_c = Eigen::Quaterniond(x0[6], x0[3], x0[4], x0[5]).toRotationMatrix();
    } else {
        R_b_c = RIC[0];
        P_b_c = TIC[0];
    }
    dy = Eigen::VectorXd(dx.size() + 1);
    dy.segment(0, dx.size()) = dx;
    for (int i = 0; i < (int)p_is.size(); i++)
        dy.segment(keep_block_idx[p_is[i]], 6).setZero();

    {
        Eigen::Vector3d v_c0_c0_old;
        Eigen::Vector3d P_I_ck0_old;
        Eigen::Matrix3d R_I_ck0_old;
        double sum_distance_old;
        std::vector<double> idepth0_old;
        Eigen::Matrix3d R_WI_WC_old;

        Getestimate(v_c0_c0_old, P_I_ck0_old, R_I_ck0_old, sum_distance_old, idepth0_old, R_WI_WC_old, keep_block_data);

        Eigen::Vector3d v_c0_c0_new;
        Eigen::Vector3d P_I_ck0_new;
        Eigen::Matrix3d R_I_ck0_new;
        double sum_distance_new;
        std::vector<double> idepth0_new;
        Eigen::Matrix3d R_WI_WC_new;
        Getestimate(v_c0_c0_new, P_I_ck0_new, R_I_ck0_new, sum_distance_new, idepth0_new, R_WI_WC_new, keep_block_addr);

        dy.segment(keep_block_idx[v_i], 3) = v_c0_c0_new - v_c0_c0_old;
        dy.segment(keep_block_idx[p_is[0]], 3) = P_I_ck0_new - P_I_ck0_old;


        {
#if 1
            Eigen::Quaterniond tmp = Eigen::Quaterniond(R_I_ck0_new) * Eigen::Quaterniond(R_I_ck0_old).inverse();
#else
            Eigen::Quaterniond tmp = Eigen::Quaterniond(R_I_ck0_old).inverse() * Eigen::Quaterniond(R_I_ck0_new);
#endif

            dy.segment(keep_block_idx[p_is[0]] + 3, 3) = 2.0 * tmp.vec();
            if (!(tmp.w() >= 0))
                dy.segment(keep_block_idx[p_is[0]] + 3, 3) = 2.0 * -tmp.vec();
        }

        dy(dx.size()) = sum_distance_new - sum_distance_old;
        for (int i = 0; i < (int)idepth0_old.size(); i++)
            dy(keep_block_idx[lambda_is[i]]) = idepth0_new[i] - idepth0_old[i];

        if (ESTIMATE_R_WI_WC) {
#if 1
            Eigen::Quaterniond tmp = Eigen::Quaterniond(R_WI_WC_new) * Eigen::Quaterniond(R_WI_WC_old).inverse();
#else
            Eigen::Quaterniond tmp = Eigen::Quaterniond(R_WI_WC_old).inverse() * Eigen::Quaterniond(R_WI_WC_new);
#endif


            dy.segment(keep_block_idx[R_WI_WC_i] + 3, 3) = 2.0 * tmp.vec();
            if (!(tmp.w() >= 0))
                dy.segment(keep_block_idx[R_WI_WC_i] + 3, 3) = 2.0 * -tmp.vec();
        }



    }



    // Eigen::Vector3d v_c0_c0 = (R_w_ck[0].transpose() * R_WI_WC.transpose() * v_w_b0 + R_b_c.transpose() * Utility::skewSymmetric(gyr0 - Bg0) * P_b_c) / sum_distance;


    {
        matrix_infoG.clear();
        matrix_infoG.push_back(MatirxInfo(keep_block_idx[p_is[0]], keep_block_idx[p_is[0]], scale_factor * R_WI_WC - Eigen::Matrix3d::Identity(), false, 3, 3));
        if (fix_scale)
            matrix_infoG.push_back(MatirxInfo(keep_block_idx[p_is[0]], keep_block_idx[scale_factor_i], R_WI_WC * P_w_ck[0], false, 3, 1));
        if (ESTIMATE_R_WI_WC)
            matrix_infoG.push_back(MatirxInfo(keep_block_idx[p_is[0]], keep_block_idx[R_WI_WC_i] + 3, -scale_factor * R_WI_WC * Utility::skewSymmetric(P_w_ck[0]), false, 3, 3));

        matrix_infoG.push_back(MatirxInfo(keep_block_idx[p_is[0]] + 3, keep_block_idx[p_is[0]] + 3, R_WI_WC * R_w_ck[0] - Eigen::Matrix3d::Identity(), false, 3, 3));
        if (ESTIMATE_R_WI_WC)
            matrix_infoG.push_back(MatirxInfo(keep_block_idx[p_is[0]] + 3, keep_block_idx[R_WI_WC_i] + 3, R_WI_WC, false, 3, 3));


        matrix_infoG.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[p_is[0]] + 3, Utility::skewSymmetric(R_w_ck[0].transpose()*R_WI_WC.transpose()*v_w_b0), false, 3, 3));
        matrix_infoG.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[v_i], R_w_ck[0].transpose() *R_WI_WC.transpose() - Eigen::Matrix3d::Identity(), false, 3, 3));


        if (ESTIMATE_EXTRINSIC) {
            matrix_infoG.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[extrinsic_i], R_b_c.transpose() * Utility::skewSymmetric(gyr0 - Bg0), false, 3, 3));
            matrix_infoG.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[extrinsic_i] + 3, Utility::skewSymmetric(R_b_c .transpose()*Utility::skewSymmetric(gyr0 - Bg0) * P_b_c), false, 3, 3));
        }
        if (ESTIMATE_R_WI_WC)
            matrix_infoG.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[R_WI_WC_i] + 3, R_w_ck[0].transpose()*Utility::skewSymmetric(R_WI_WC.transpose()*v_w_b0), false, 3, 3));


        // matrix_infoG.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[p_is[0]] + 3, Utility::skewSymmetric(R_w_ck[0].transpose()*R_WI_WC.transpose()*v_w_b0) / sum_distance, false, 3, 3));
        // matrix_infoG.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[v_i], R_w_ck[0].transpose() *R_WI_WC.transpose() / sum_distance - Eigen::Matrix3d::Identity(), false, 3, 3));
        // if (fix) {
        //     Eigen::Matrix3d tmp = Eigen::Matrix3d::Zero();
        //     for (int i = 1; i < (int)p_is.size(); i++) {
        //         tmp += -v_c0_c0 / sum_distance * scale_factor * (P_w_ck[0] - P_w_ck[i]).transpose() / (P_w_ck[0] - P_w_ck[i]).norm(); //P_w_c0
        //         matrix_infoG.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[p_is[i]], -v_c0_c0 / sum_distance * scale_factor * (P_w_ck[i] - P_w_ck[0]).transpose() / (P_w_ck[i] - P_w_ck[0]).norm(), false, 3, 3));
        //     }
        //     matrix_infoG.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[p_is[0]], tmp, false, 3, 3));
        //     if (fix_scale)
        //         matrix_infoG.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[scale_factor_i], -v_c0_c0 / scale_factor, false, 3, 1));
        // }

        // if (ESTIMATE_EXTRINSIC) {
        //     matrix_infoG.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[extrinsic_i], R_b_c.transpose() * Utility::skewSymmetric(gyr0 - Bg0) / sum_distance, false, 3, 3));
        //     matrix_infoG.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[extrinsic_i] + 3, Utility::skewSymmetric(R_b_c .transpose()*Utility::skewSymmetric(gyr0 - Bg0) * P_b_c) / sum_distance, false, 3, 3));
        // }
        // matrix_infoG.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[v_i] + 6, R_b_c.transpose() * Utility::skewSymmetric(P_b_c) / sum_distance, false, 3, 3));
        // if (ESTIMATE_R_WI_WC)
        //     matrix_infoG.push_back(MatirxInfo(keep_block_idx[v_i], keep_block_idx[R_WI_WC_i] + 3, R_w_ck[0].transpose()*Utility::skewSymmetric(R_WI_WC.transpose()*v_w_b0 / sum_distance), false, 3, 3));

        if (ESTIMATE_R_WI_WC)
            matrix_infoG.push_back(MatirxInfo(keep_block_idx[R_WI_WC_i] + 3, keep_block_idx[R_WI_WC_i] + 3, R_WI_WC - Eigen::Matrix3d::Identity(), false, 3, 3));



        // lambda_l

        for (int j = 0; j < (int)lambda_is.size(); j++) {
            matrix_infoG.push_back(MatirxInfo(keep_block_idx[lambda_is[j]], keep_block_idx[lambda_is[j]], flower_focus::Vector1d(sum_distance / scale_factor - 1), false, 1, 1));
            if (fix) {
                flower_focus::Matrix1_3d tmp = flower_focus::Matrix1_3d::Zero();
                for (int i = 1; i < (int)p_is.size(); i++) {
                    tmp += idepth0[j] * (P_w_ck[0] - P_w_ck[i]).transpose() / (P_w_ck[0] - P_w_ck[i]).norm(); //P_w_c0
                    idepth_pos[i].block(keep_block_idx[lambda_is[j]] - min_idepth_idx, 0, 1, 3) = idepth0[j] * (P_w_ck[i] - P_w_ck[0]).transpose() / (P_w_ck[i] - P_w_ck[0]).norm();
                }
                idepth_pos[0].block(keep_block_idx[lambda_is[j]] - min_idepth_idx, 0, 1, 3) = tmp;

            } else {
                if (fix_scale)matrix_infoG.push_back(MatirxInfo(keep_block_idx[lambda_is[j]], keep_block_idx[scale_factor_i], flower_focus::Vector1d(-idepth0[j] / scale_factor / scale_factor), false, 1, 1));
            }

        }
        if (fix && lambda_is.size()) {
            for (int i = 0; i < (int)idepth_pos.size(); i++)
                matrix_infoG.push_back(MatirxInfo(min_idepth_idx, keep_block_idx[p_is[i]], idepth_pos[i], false, lambda_is.size(), 3));
        }

        if (fix) {
            flower_focus::Matrix1_3d tmp = flower_focus::Matrix1_3d::Zero();
            for (int i = 1; i < (int)p_is.size(); i++) {
                tmp += scale_factor * (P_w_ck[0] - P_w_ck[i]).transpose() / (P_w_ck[0] - P_w_ck[i]).norm(); //P_w_c0
                matrix_infoG.push_back(MatirxInfo(A.cols(), keep_block_idx[p_is[i]], scale_factor * (P_w_ck[i] - P_w_ck[0]).transpose() / (P_w_ck[i] - P_w_ck[0]).norm(), false, 1, 3));
            }
            matrix_infoG.push_back(MatirxInfo(A.cols(), keep_block_idx[p_is[0]], tmp, false, 1, 3));
            if (fix_scale)
                matrix_infoG.push_back(MatirxInfo(A.cols(), keep_block_idx[scale_factor_i], flower_focus::Vector1d(sum_distance / scale_factor), false, 1, 1));
        }

    }

}




void MarginalizationInfo::UpdateMatix(const Eigen::MatrixXd& A, Eigen::MatrixXd& A2) {

    Eigen::MatrixXd A3 = A;
    Eigen::MatrixXd A3_tmp = A3;

    for (int i = 0; i < (int)matrix_infoG.size(); ++i) {
        auto& info = matrix_infoG[i];
        matrix_update(A3_tmp, info.m, A3, info.idx1, info.idx2,  LEFT, info.is_identity, info.b_row, info.b_col);
    }
    A3_tmp = A3;
    for (int i = 0; i < (int)matrix_infoG.size(); ++i) {
        auto& info = matrix_infoG[i];
        matrix_update(A3_tmp, info.m, A3, info.idx1, info.idx2,  RIGHT, info.is_identity, info.b_row, info.b_col);
    }
    A3 = A3.selfadjointView<Eigen::Upper>();
    A2 = A3.block(0, 0, A.cols() - 1, A.cols() - 1);


}



void MarginalizationInfo::UpdateVectorTranspose(const Eigen::VectorXd& b, Eigen::VectorXd& b2) {



    Eigen::VectorXd b3 = b;
    Eigen::VectorXd b3_tmp = b;

    for (int i = 0; i < (int)matrix_infoG.size(); ++i) {
        auto& info = matrix_infoG[i];
        vector_update(b3_tmp, info.m, b3, info.idx1, info.idx2,  RIGHT, info.is_identity, info.b_row, info.b_col);
    }

    b2 = b3.segment(0, b.rows() - 1);



}


void MarginalizationInfo::UpdateVector(const Eigen::VectorXd& b, Eigen::VectorXd& b2, bool onlyF) {

    Eigen::VectorXd b3 = Eigen::VectorXd::Zero(b.rows() + 1);
    b3.segment(0, b.rows()) = b;
    Eigen::VectorXd b3_tmp;

    b3_tmp = b3;
    for (int i = 0; i < (int)matrix_infoG.size(); ++i) {
        auto& info = matrix_infoG[i];
        vector_update(b3_tmp, info.m, b3, info.idx1, info.idx2,  LEFT, info.is_identity, info.b_row, info.b_col);
    }

    b2 = b3;

}



void MarginalizationInfo::compute_residual_rhs(Eigen::VectorXd& residual, Eigen::VectorXd& rhs) {
    ASSERT(m == 0);

    if (have_hist) {

        Eigen::VectorXd dx(n);
        for (int i = 0; i < static_cast<int>(keep_block_size.size()); i++) {
            int size = keep_block_size[i];
            int idx = keep_block_idx[i] - m;
            Eigen::VectorXd x = Eigen::Map<const Eigen::VectorXd>(keep_block_addr[i], size);
            Eigen::VectorXd x0 = Eigen::Map<const Eigen::VectorXd>(keep_block_data[i], size);
            if (size != 7)
                dx.segment(idx, size) = x - x0;
            else {
                dx.segment<3>(idx + 0) = x.head<3>() - x0.head<3>();
                Eigen::Quaterniond tmp = Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5));
                dx.segment<3>(idx + 3) = 2.0 * tmp.vec();
                if (!(tmp.w() >= 0))
                    dx.segment<3>(idx + 3) = 2.0 * -tmp.vec();
            }
        }
        Eigen::VectorXd dy;
        BuildMatrixG(keep_block_addr, dx, dy);


        residual = linearized_jacobiansF * dy + linearized_residuals;
        Eigen::VectorXd tmp;
        UpdateVectorTranspose(FTb, tmp);

        UpdateVectorTranspose(FTAF * dy + FTb, rhs);


    } else {

        Eigen::VectorXd dx(n);
        for (int i = 0; i < static_cast<int>(keep_block_size.size()); i++) {
            int size = keep_block_size[i];
            int idx = keep_block_idx[i] - m;
            Eigen::VectorXd x = Eigen::Map<const Eigen::VectorXd>(keep_block_addr[i], size);
            Eigen::VectorXd x0 = Eigen::Map<const Eigen::VectorXd>(keep_block_data[i], size);
            if (size != 7)
                dx.segment(idx, size) = x - x0;
            else {
                dx.segment<3>(idx + 0) = x.head<3>() - x0.head<3>();

                Eigen::Quaterniond tmp = Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5));

                dx.segment<3>(idx + 3) = 2.0 * tmp.vec();
                if (!(tmp.w() >= 0))
                    dx.segment<3>(idx + 3) = 2.0 * -tmp.vec();
            }
        }


        residual = linearized_jacobians * dx + linearized_residuals;
        rhs = A * dx + b;
    }



}


void MarginalizationInfo::UpdateF() {
    if (have_hist) {

        double sum_distance = 0;
        for (int i = 1; i < (int)p_is.size(); i++) {
            Eigen::Vector3d P_w_ck_i = Eigen::Vector3d(keep_block_data[p_is[i]][0], keep_block_data[p_is[i]][1], keep_block_data[p_is[i]][2]);
            Eigen::Vector3d P_w_ck_0 = Eigen::Vector3d(keep_block_data[p_is[0]][0], keep_block_data[p_is[0]][1], keep_block_data[p_is[0]][2]);
            sum_distance += (P_w_ck_i - P_w_ck_0).norm();
        }
        if (sum_distance > 0.5)fix = true;
        else fix = false;

        BuildMatrixF(keep_block_data);


        max_idepth_idx = -1;
        min_idepth_idx = 100000;

        if (lambda_is.size()) {
            for (int j = 0; j < (int)lambda_is.size(); j++) {
                if (keep_block_idx[lambda_is[j]] > max_idepth_idx)max_idepth_idx = keep_block_idx[lambda_is[j]];
                if (keep_block_idx[lambda_is[j]] < min_idepth_idx)min_idepth_idx = keep_block_idx[lambda_is[j]];
            }
            ASSERT(max_idepth_idx - min_idepth_idx + 1 == (int)lambda_is.size());
        }

        idepth_pos.resize(p_is.size());
        if (lambda_is.size()) {
            for (int i = 0; i < (int)p_is.size(); i++)
                idepth_pos[i] = Eigen::MatrixXd(max_idepth_idx - min_idepth_idx + 1, 3);
        }


    }

    {
        FTAF = Eigen::MatrixXd::Zero(A.rows() + 1, A.rows() + 1);
        FTAF.block(0, 0, A.rows(), A.rows()) = A;
        Eigen::MatrixXd A3_tmp = FTAF;

        for (int i = 0; i < (int)matrix_infoF.size(); ++i) {
            auto& info = matrix_infoF[i];
            matrix_update(A3_tmp, info.m, FTAF, info.idx1, info.idx2,  LEFT, info.is_identity, info.b_row, info.b_col);
        }
        A3_tmp = FTAF;
        for (int i = 0; i < (int)matrix_infoF.size(); ++i) {
            auto& info = matrix_infoF[i];
            matrix_update(A3_tmp, info.m, FTAF, info.idx1, info.idx2,  RIGHT, info.is_identity, info.b_row, info.b_col);
        }
        FTAF = FTAF.selfadjointView<Eigen::Upper>();
    }

    {

        FTb = Eigen::VectorXd::Zero(b.rows() + 1);;
        FTb.segment(0, b.rows()) = b;
        Eigen::VectorXd b3_tmp;

        b3_tmp = FTb;
        for (int i = 0; i < (int)matrix_infoF.size(); ++i) {
            auto& info = matrix_infoF[i];
            vector_update(b3_tmp, info.m, FTb, info.idx1, info.idx2,  RIGHT, info.is_identity, info.b_row, info.b_col);
        }

    }


    {
        Eigen::MatrixXd tmp = Eigen::MatrixXd::Zero(linearized_jacobians.rows() + 1, linearized_jacobians.rows() + 1);
        tmp.block(0, 0, linearized_jacobians.rows(), linearized_jacobians.rows()) = linearized_jacobians;
        Eigen::MatrixXd A3_tmp = tmp;

        for (int i = 0; i < (int)matrix_infoF.size(); ++i) {
            auto& info = matrix_infoF[i];
            matrix_update(A3_tmp, info.m, tmp, info.idx1, info.idx2,  LEFT, info.is_identity, info.b_row, info.b_col);
        }
        ASSERT(tmp.block(tmp.rows() - 1, 0, 1, tmp.rows()).norm() == 0);
        linearized_jacobiansF = tmp.block(0, 0, tmp.rows() - 1, tmp.rows());
    }


}


void MarginalizationInfo::resetLinerizationPoint() {
    if (USE_REPARAM)
        ASSERT(0);
    else {
        Eigen::VectorXd dx(n);
        for (int i = 0; i < static_cast<int>(keep_block_size.size()); i++) {
            int size = keep_block_size[i];
            int idx = keep_block_idx[i] - m;
            Eigen::VectorXd x = Eigen::Map<const Eigen::VectorXd>(keep_block_addr[i], size);
            Eigen::VectorXd x0 = Eigen::Map<const Eigen::VectorXd>(keep_block_data[i], size);
            if (size != 7)
                dx.segment(idx, size) = x - x0;
            else {
                dx.segment<3>(idx + 0) = x.head<3>() - x0.head<3>();
                dx.segment<3>(idx + 3) = 2.0 * (Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
                if (!((Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).w() >= 0))
                    dx.segment<3>(idx + 3) = 2.0 * -(Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
            }
        }
        b = A * dx + b;
        linearized_residuals = linearized_jacobians * dx + linearized_residuals;

        for (int i = 0; i < static_cast<int>(keep_block_size.size()); i++) {
            int size = keep_block_size[i];
            Eigen::Map<Eigen::VectorXd>(keep_block_data[i], size) = Eigen::Map<Eigen::VectorXd>(keep_block_addr[i], size);

        }
    }
}

void MarginalizationInfo::GetJacobianAndResidual() {
    A.diagonal().array() += EPSS;

    Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> llt = A.selfadjointView<Eigen::Upper>().llt();
    linearized_jacobians = llt.matrixL().transpose();
    Eigen::VectorXd inc = llt.solve(b);
    linearized_residuals = linearized_jacobians * inc;
    init_cost = 0.5 * linearized_residuals.squaredNorm();
    LOG_OUT << "gnss cost:" << init_cost << std::endl;
    if (llt.info() != Eigen::Success)
        assert(0);

    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes2(A);
    // Eigen::VectorXd S = Eigen::VectorXd((saes2.eigenvalues().array() > EPSS).select(saes2.eigenvalues().array(), 0));
    // Eigen::VectorXd S_sqrt = S.cwiseSqrt();
    // linearized_jacobians = S_sqrt.asDiagonal() * saes2.eigenvectors().transpose();

    // Eigen::VectorXd S_inv = Eigen::VectorXd((saes2.eigenvalues().array() > EPSS).select(saes2.eigenvalues().array().inverse(), 0));
    // Eigen::VectorXd S_inv_sqrt = S_inv.cwiseSqrt();
    // Eigen::MatrixXd tmp = S_inv_sqrt.asDiagonal() * saes2.eigenvectors().transpose();
    // linearized_residuals = tmp * b;

    // LOG_OUT << "0.5*linearized_residuals.squaredNorm():" << 0.5 * linearized_residuals.squaredNorm() << std::endl;
    // cost0 = 0.5 * linearized_residuals.squaredNorm();
    // A = linearized_jacobians.transpose() * linearized_jacobians;
    // b = linearized_jacobians.transpose() * linearized_residuals;
}

void MarginalizationInfo::getParameterBlocks(bool needJacobian) {
    keep_block_size.clear();
    keep_block_idx.clear();
    keep_block_data.clear();
    keep_block_addr_set.clear();
    keep_block_addr.clear();

    for (const auto& it : parameter_block_idx) {
        if (it.second >= m) {
            keep_block_size.push_back(parameter_block_size[it.first]);
            keep_block_idx.push_back(parameter_block_idx[it.first]);
            keep_block_data.push_back(parameter_block_data[it.first]);
            keep_block_addr.push_back(reinterpret_cast<double*>(it.first));
        }
    }
    for (int i = 0; i < (int)keep_block_addr.size(); i++)
        keep_block_addr_set[keep_block_addr[i]] = i;
    if (needJacobian)
        GetJacobianAndResidual();
    ASSERT(m == 0);


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
                int size = marginalization_info->keep_block_size[i], local_size = localSize(size);
                int idx = marginalization_info->keep_block_idx[i] - m;
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> jacobian(jacobians[i], n, size);
                jacobian.setZero();
                jacobian.leftCols(local_size) = marginalization_info->linearized_jacobians.middleCols(idx, local_size);
            }
        }
    }
    return true;
}
