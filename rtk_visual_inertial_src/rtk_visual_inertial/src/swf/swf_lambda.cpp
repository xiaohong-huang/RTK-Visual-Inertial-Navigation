#include"swf.h"
#include "lambda.h"
#include "../factor/gnss_factor.h"
#include "common_function.h"



void FindReferenceSatellites(mea_t* rover, MarginalizationInfo* marginalization_info, ObsMea** Rsatellites) {
    std::vector<double>biase_floats[6];
    std::vector<ObsMea*>biase_obs[6];
    for (uint8_t i = 0; i < rover->obs_count; i++) {
        ObsMea* d = rover->obs_data + i;
        for (uint8_t f = 0; f < NFREQ; f++) {
            if (!d->RTK_Npoint[f] || d->RTK_Npoint[f]->use )continue;
            long valuepoint = reinterpret_cast<long>(&(d->RTK_Npoint[f]->value));
            if (marginalization_info->parameter_block_idx[valuepoint]) {
                int a = marginalization_info->parameter_block_idx[valuepoint] - marginalization_info->m;
                if (a < 0)continue;
                biase_floats[d->sys * 2 + f].push_back(d->RTK_Npoint[f]->value);
                biase_obs[d->sys * 2 + f].push_back(d);
            }
        }
    }
    for (int sys_f = 0; sys_f < 6; sys_f++) {
        if (biase_floats[sys_f].size()) {
            std::vector<double >biase_floats_cost(biase_floats[sys_f].size(), 0);
            for (int j = 0; j < (int)biase_floats[sys_f].size(); j++) {
                double s = biase_floats[sys_f][j];
                assert(biase_floats_cost[j] == 0);
                for (int i = 0; i < (int)biase_floats[sys_f].size(); i++) {
                    double s2 = biase_floats[sys_f][i] - s;
                    s2 -= round(s2);
                    biase_floats_cost[j] += fabs(s2);
                }
            }

            std::vector<double>biase_float = biase_floats_cost;
            std::sort(biase_float.begin(), biase_float.end());
            double min = biase_float[0];
            bool found = false;
            for (int i = 0; i < (int)biase_floats_cost.size(); i++) {
                if (biase_floats_cost[i] == min) {
                    Rsatellites[sys_f] = biase_obs[sys_f][i];
                    found = true;
                }

            }
            assert(found);

            assert(Rsatellites[sys_f]);
        }
    }
}


Eigen::MatrixXd matrix_reorder(Eigen::MatrixXd A, std::vector<int>ordering) {
    assert(A.cols() == A.rows());
    assert(A.cols() >= (int)ordering.size());
    Eigen::MatrixXd B(A.cols(), A.cols());
    Eigen::MatrixXd C(A.cols(), A.cols());
    for (int i = 0; i < (int)ordering.size(); i++) {
        B.block(i, 0, 1, A.cols()) = A.block(ordering[i], 0, 1, A.cols());
    }
    for (int i = 0; i < (int)ordering.size(); i++) {
        C.block(0, i, A.cols(), 1) = B.block(0, ordering[i], A.cols(), 1);
    }
    return C.block(0, 0, ordering.size(), ordering.size());
}

Eigen::VectorXd vector_reorder(Eigen::VectorXd a, std::vector<int>ordering) {
    // std::cout<<"0"<<std::endl;
    assert(a.size() >= (int)ordering.size());
    // std::cout<<"1"<<std::endl;
    Eigen::VectorXd b(ordering.size());
    // std::cout<<"2"<<std::endl;
    for (int i = 0; i < (int)ordering.size(); i++) {
        b(i) = a(ordering[i]);
    }
    return b;
}

void SWFOptimization::LambdaSearch() {
    static int fixcount;
    static double gnss_last_updatetime = 0;
    if (rovers[rover_count - 1]->ros_time - gnss_last_updatetime > 2) {
        fixcount = 0;
        last_fix = false;
        not_fix_count = 0;
        std::cout << "reset FIX information" << std::endl;
    }
    gnss_last_updatetime = rovers[rover_count - 1]->ros_time;

    MarginalizationInfo* marginalization_info = new MarginalizationInfo();
    marginalization_info->setmarginalizeinfo(parameter_block_addr, parameter_block_global_size, A, b, false);

    if (marginalization_info->n < 6) {
        delete marginalization_info;
        return;
    }

    Eigen::MatrixXd Qy = marginalization_info->A.inverse();
    Eigen::MatrixXd D2(marginalization_info->n, marginalization_info->n);
    Eigen::MatrixXd y(marginalization_info->n, 1);

    for (uint8_t i = 0; i < MAXSATNUM * 2; i++) {
        for (auto Nv_it = rtk_phase_bias_variables[i].begin(); Nv_it != rtk_phase_bias_variables[i].end(); Nv_it++) {
            long valuepoint = reinterpret_cast<long>(&(Nv_it->value));
            if (marginalization_info->parameter_block_idx[valuepoint] && marginalization_info->parameter_block_idx[valuepoint] - marginalization_info->m >= 0) {
                y(marginalization_info->parameter_block_idx[valuepoint] - marginalization_info->m) = Nv_it->value;
            }
        }
    }

    int countD = 0;
    Eigen::MatrixXd D3(marginalization_info->n * rover_count, marginalization_info->n);
    D3.setZero();

    for (uint8_t i = 0; i < MAXSATNUM * 2; i++) {
        for (auto it1 = rtk_phase_bias_variables[i].begin(); it1 != rtk_phase_bias_variables[i].end(); it1++) {
            it1->use = false;
        }
    }

    int last_count = 0;
    int last_ref_count = 0;
    for (int ir = rover_count - 1; ir >= 0; ir--) {
        mea_t* rover = rovers[ir];
        // if(ir!=rover_count-1)continue;
        if (USE_RTK) {

            ObsMea* Rsatellites[6] = {0};
            FindReferenceSatellites(rover, marginalization_info, Rsatellites);
            for (int j = 0; j < 6; j++) {
                if (Rsatellites[j] && ir == rover_count - 1)last_ref_count++;
            }
            for (uint8_t i = 0; i < rover->obs_count; i++) {
                ObsMea* d = rover->obs_data + i;
                uint8_t sys = d->sys;
                for (uint8_t f = 0; f < NFREQ; f++) {
                    if (!d->RTK_Npoint[f] )continue;
                    if (!Rsatellites[sys * 2 + f]) {
                        Rsatellites[sys * 2 + f] = d;
                        continue;
                    }
                    if (d->RTK_Npoint[f]->use || d == Rsatellites[sys * 2 + f])continue;
                    d->RTK_Npoint[f]->use = true;
                    long valuepoint = reinterpret_cast<long>(&(d->RTK_Npoint[f]->value));
                    long bestvaluepoint = reinterpret_cast<long>(&(Rsatellites[sys * 2 + f]->RTK_Npoint[f]->value));

                    if (marginalization_info->parameter_block_idx[valuepoint]) {
                        int a = marginalization_info->parameter_block_idx[valuepoint] - marginalization_info->m;
                        int b = marginalization_info->parameter_block_idx[bestvaluepoint] - marginalization_info->m;

                        if (a < 0 || b < 0)continue;
                        assert(y(a) == d->RTK_Npoint[f]->value);
                        assert(y(b) == Rsatellites[sys * 2 + f]->RTK_Npoint[f]->value);
                        assert(b >= 0);
                        assert(a != b);
                        assert(d->half_flag[f] == d->RTK_Npoint[f]->half_flag);
                        assert(Rsatellites[sys * 2 + f]->half_flag[f] == Rsatellites[sys * 2 + f]->RTK_Npoint[f]->half_flag);
                        assert(((d->half_flag[f]) & 8) && ((Rsatellites[sys * 2 + f]->half_flag[f]) & 8) && ((d->half_flag[f]) & 2) && ((Rsatellites[sys * 2 + f]->half_flag[f]) & 2));

                        if (abs((y(a) - y(b)) - round(y(a) - y(b))) < (last_fix ? 0.2 : 1.4)) { //&&d->sys!=SYS_CMP
                            D3(countD, a) = 1;
                            D3(countD, b) = -1;
                            countD++;
                            //std::cout<<(int)Rsatellites[sys*2+f]->sat<<","<<(int)d->sat<<","<<(int)d->sys<<std::endl;
                            if (ir == rover_count - 1)last_count++;
                        } else {
                            std::cout << "setNtozero" << (int)(d->sat) << "," << (int)(d->sys) << "," << (int)f << "," << y(a) - y(b) << "," << d->el << "," <<
                                      (int)(d->half_flag[f]) << "," << (int)(d->RTK_Npoint[f]->half_flag) << "," << std::endl;
                        }
                    }
                }
            }
        }
    }
    if (last_count + last_ref_count < 6 || last_count < 4)return;

    last_fix = false;
    assert(countD <= marginalization_info->n * 1);
    Eigen::MatrixXd D = D3.block(0, 0, countD, D3.cols());

    if (D.rows() < 4) {
        delete marginalization_info;
        std::cout << "return" << std::endl;
        return;
    }
    Eigen::MatrixXd Qb = D * Qy * D.transpose();
    Eigen::MatrixXd b = D * y;
    std::cout << b.transpose() << std::endl;

    static int count_gns = 0;
    count_gns++;
    double s[2];
    double* F = mat(D.rows(), 2);
    bool search_ok = false;
    int info;


    if (!(info = lambda(Qb.rows(), 2, b.data(), Qb.data(), F, s))) { //nb:5//nb有正负1的行数
        std::cout << s[1] / s[0] << "," << s[0] << "," << s[1] << std::endl;

        Eigen::VectorXd F1 = Eigen::Map<Eigen::VectorXd>(F, D.rows());
        Eigen::VectorXd F2 = Eigen::Map<Eigen::VectorXd>(F + D.rows(), D.rows());

        Eigen::VectorXd e1 = F1 - b;
        Eigen::VectorXd e2 = F2 - b;
        std::vector<int>same_index;
        std::vector<int>different_index;

        for (int i = 0; i < D.rows(); i++) {
            if (fabs(F1(i) - F2(i)) < 1e-2) {
                same_index.push_back(i);
            } else {
                different_index.push_back(i);
            }
        }
        Eigen::MatrixXd Qb2 = Qb;
        for (int i = 0; i < (int)different_index.size(); i++) {
            int i0 = different_index[i];
            e1(i0) = e2(i0) = 0;
            for (int j0 = 0; j0 < D.rows(); j0++) {
                if (i0 == j0)Qb2(i0, j0) = 1;
                else Qb2(i0, j0) = Qb2(j0, i0) = 0;
            }
        }

        double same_cost = e1.transpose() * Qb2.inverse() * e1;
        double s1 = s[1] - same_cost;
        double s0 = s[0] - same_cost;
        if (fabs(s0) < 1e-3)s0 = 1e-3;
        std::cout << s1 / s0 << "," << s0 << "," << s1 << "," << different_index.size() << std::endl;
        assert(s1 > 0 && s0 > 0);


        if (s[0] <= 0.0 || s[1] / s[0] >= 2 || s1 / s0 >= 2) {
            not_fix_count = 0;
            search_ok = true;
            gnss_fix_solution_count++;
            last_fix_time = rovers[rover_count - 1]->ros_time;
        } else {
            not_fix_count += 1;
        }
    }

    double tf[6] = {0};
    bool tfb[6] = {false};
    if (search_ok) {
        fixcount++;
        if (fixcount > FIX_CONTINUE_THRESHOLD) {
            last_fix = true;
            TicToc t_p;
            MarginalizationInfo* marginalization_info2 = new MarginalizationInfo();
            MarginalizationFactor* factormarge;
            ResidualBlockInfo* residual_block_info;

            factormarge = new MarginalizationFactor(last_marg_info);
            residual_block_info = new ResidualBlockInfo(factormarge, NULL, last_marg_info->keep_block_addr,
                    std::vector<int> {}, std::vector<int> {});


            marginalization_info2->addResidualBlockInfo(residual_block_info);

            double* pindex[1000] = {0};
            PBtype* pindexNtype[1000] = {0};

            int numf1 = 0;
            int numf2 = 0;
            for (uint8_t i = 0; i < MAXSATNUM * 2; i++) {

                for (auto Nv_it = rtk_phase_bias_variables[i].begin(); Nv_it != rtk_phase_bias_variables[i].end(); Nv_it++) {
                    long valuepoint = reinterpret_cast<long>(&(Nv_it->value));
                    if (marginalization_info->parameter_block_idx[valuepoint] && marginalization_info->parameter_block_idx[valuepoint] - marginalization_info->m >= 0) {
                        pindex[marginalization_info->parameter_block_idx[valuepoint] - marginalization_info->m] = &(Nv_it->value);
                        pindexNtype[marginalization_info->parameter_block_idx[valuepoint] - marginalization_info->m] = &(*Nv_it);
                        if (i % 2)numf2++;
                        else numf1++;
                    }
                }
            }
            int numcolos = D.cols();


            for (int i = 0; i < last_count; i++) {
                double Fvalue;
                Eigen::MatrixXd diffrowresult;
                {
                    diffrowresult = D.block(i, 0, 1, numcolos);
                    Fvalue = round(F[i]);
                }

                int nump = 0;
                int numn = 0;
                int indexp = -1;
                int indexn = -1;
                for (int i2 = 0; i2 < numcolos; i2++) {
                    if (diffrowresult(0, i2) == -1) {
                        indexn = i2;
                        nump++;
                    }
                    if (diffrowresult(0, i2) == 1) {
                        indexp = i2;
                        numn++;
                    }
                }
                assert(numn == 1 && nump == 1);
                if (numn == 1 && nump == 1) {
                    double* ppoint = pindex[indexp];
                    double* npoint = pindex[indexn];
                    PBtype* ppointNtype = pindexNtype[indexp];
                    PBtype* npointNtype = pindexNtype[indexn];
                    assert(ppointNtype && npointNtype && npoint && ppoint);
                    uint8_t sys = ppointNtype->sys;
                    uint8_t f = ppointNtype->f;

                    if (tfb[sys * 2 + f] == false) {
                        FixedIntegerFactor* factormarge = new FixedIntegerFactor(0, 1 / 0.03);
                        ResidualBlockInfo* residual_block_info = new ResidualBlockInfo(factormarge, NULL,
                                std::vector<double*> {&tf[sys * 2 + f], npoint}, std::vector<int> {0}, std::vector<int> {});
                        marginalization_info2->addResidualBlockInfo(residual_block_info);
                        tfb[sys * 2 + f] = true;
                    }

                    FixedIntegerFactor* factormarge = new FixedIntegerFactor(Fvalue, 1 / 0.03);
                    ResidualBlockInfo* residual_block_info = new ResidualBlockInfo(factormarge, NULL,
                            std::vector<double*> {&tf[sys * 2 + f], ppoint}, std::vector<int> {0}, std::vector<int> {});
                    marginalization_info2->addResidualBlockInfo(residual_block_info);

                }
            }

            fix = true;
            mea_t* rovernew = rovers[rover_count - 1];
            rtk_fix = true;

            printf("                after lambda new pos:%.5f,%.5f,%.5f,%.5f,%d,%d\r\n", rovernew->ros_time, para_pose[image_count + rover_count - 1][0], para_pose[image_count + rover_count - 1][1], para_pose[image_count + rover_count - 1][2], gnss_fix_solution_count, count_gns);

            PhaseBiasSaveAndReset();
            marginalization_info2->marginalize(false, true);
            marginalization_info2->getParameterBlocks();
            PhaseBiasRestore();

            if (last_marg_info) {
                delete last_marg_info;
                if (USE_GLOBAL_OPTIMIZATION) {
                    my_problem.RemoveResidualBlock(marg_residual_block_id);
                }
            }
            last_marg_info = marginalization_info2;
            if (USE_GLOBAL_OPTIMIZATION) {
                MarginalizationFactor* factor = new MarginalizationFactor(last_marg_info);
                marg_residual_block_id = my_problem.AddResidualBlock(factor, 0, last_marg_info->keep_block_addr);
            }
        }

    } else {
        fixcount = 0;
    }


    delete marginalization_info;
    if (F)
        free(F);
}



