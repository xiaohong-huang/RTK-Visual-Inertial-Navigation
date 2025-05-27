#include "gnss_imu_info.h"
#include "../gnss/gnss_ceres_factor.h"
#include "../gnss/gnss_factor.h"
#include "../gnss/lambda.h"

#define EstimatePcorrectionPerio 500
#define Phase_ALL_RESET_COUNT 100
#define FIX_CONTINUE_THRESHOLD 40
#define DISABLE_LONG_BASE_TIME 10


int rover_count_accumulate = 0;
double para_gnss_dt[13] = {0};
enum Optimization_mode {
    MargeExcludeMode,//marginalize the states except the selecting states.
    GNSSMargIncludeMode
};

std::list<PBtype> rtk_phase_bias_variables[MAXSATNUM * 2];
//storing the rover-only phase biases
std::list<PBtype> spp_phase_bias_variables[MAXSATNUM * 2];
//storing the pseudorange correction
std::list<PBtype> pseudorange_correction_variables[MAXSATNUM * 2];

int not_fix_count = 0;
#define ADDCERESRESIDUAL(parameter_block,residual_define,arg2,sparseset){\
    {\
        residual_define\
        problem.AddResidualBlock(factor, 0, parameter_block); \
    }\
}



#define ADDRESIDUAL(parameter_block,residual_define,arg2,sparseset){\
    {\
        bool keep=false;\
        std::vector<double*>parameter_block_vector{parameter_block};\
        std::vector<int>dropset;std::vector<int>keepset;\
        for(int vi=0;vi<(int)parameter_block_vector.size();vi++){\
            if(mode==MargeExcludeMode)keep=true;\
            if(MargePoint.find(parameter_block_vector[vi])!=MargePoint.end()){\
                dropset.push_back(vi);\
                if(mode==GNSSMargIncludeMode)keep=true;\
            }else{\
                keepset.push_back(vi);\
            }\
        }\
        if(keep||arg2){\
            {\
                if(mode==MargeExcludeMode)dropset=keepset;\
                residual_define\
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(factor, 0,parameter_block_vector,dropset,sparseset);\
                marginalization_info->addResidualBlockInfo(residual_block_info);\
            }\
        }\
    }\
}

void GNSSIMUInfo::AddCeresGnssResidual(std::set<double*>MargePoint, MarginalizationInfo* marginalization_info,
                                       ceres::Problem& problem,
                                       ceres::Solver::Options& options, int ir) {

    //gnss residual
    if (USE_GNSS) {
        bool have_base = false;
        mea_t* rover = Rovers[ir];
        double globalxyz[3];
        globalxyz[0] = para_pose[ir][0] + rover->base_xyz[0];
        globalxyz[1] = para_pose[ir][1] + rover->base_xyz[1];
        globalxyz[2] = para_pose[ir][2] + rover->base_xyz[2];
        update_azel(globalxyz, rover);

#define A4s &blackvalue
#define B4s InitialBlackCeresFactor*factor=new InitialBlackCeresFactor(1);
        ADDCERESRESIDUAL(A4s, B4s, false, std::vector<int> {})
#undef A4s
#undef B4s
        if (USE_RTK) {
            have_base = true;
            for (int i = 0; i < rover->obs_count; i++) {
                ObsMea* d = rover->obs_data + i;
                int sys = d->sys;
                for (int f = 0; f < NFREQ; f++) {
                    if (d->RTK_L[f] != 0 && rover_count_accumulate > 1)ASSERT(d->RTK_Npoint[f]);
                    if (!d->RTK_Npoint[f] )continue;
                    if (d->el < AZELMIN)continue;
                    ASSERT(d->el >= AZELMIN);
#define A44 para_pose[ir],&(d->RTK_Npoint[f]->value),para_gnss_dt+sys*2+f
#define B44 RTKCarrierPhaseCeresFactor *factor = new RTKCarrierPhaseCeresFactor(d->satellite_pos,d->RTK_L[f]*lams[d->sys][f],lams[d->sys][f],d->el,rover->br_time_diff,pow(d->RTK_Lstd[f]*lams[d->sys][f],2),rover->base_xyz,true,d->sys,f);
                    ADDCERESRESIDUAL(A44, B44, false, std::vector<int> {})
#undef A44
#undef B44

                }
            }
        }
        if (USE_RTD) {
            have_base = true;
            for (int i = 0; i < rover->obs_count; i++) {
                ObsMea* d = rover->obs_data + i;
                int sys = d->sys;
                ASSERT(sys == 0 || sys == 1 || sys == 2);
                for (int f = 0; f < NFREQ; f++) {
                    if (d->RTK_P[f] == 0.0 || d->SVH != 0 )continue;
                    if (para_gnss_dt[sys * 2 + f] != 0 && d->RTK_Pstd[f] > 2)continue;
                    if (d->el < AZELMIN)continue;
                    ASSERT(d->el >= AZELMIN);
                    have_base = true;
#define A54 para_pose[ir],para_gnss_dt+sys*2+f
#define B54 RTKPseudorangeCeresFactor *factor = new RTKPseudorangeCeresFactor(d->satellite_pos,d->RTK_P[f],  d->el,rover->br_time_diff,pow(d->RTK_Pstd[f],2), rover->base_xyz);
                    ADDCERESRESIDUAL(A54, B54, false, std::vector<int> {})
#undef A54
#undef B54
                }
            }
        }


        for (int i = 0; i < rover->obs_count; i++) {
            ObsMea* d = rover->obs_data + i;
            if (d->SVH != 0)continue;
            if (d->el < AZELMIN)continue;
            if (d->SPP_P[0] != 0)ASSERT(d->el >= AZELMIN);
            if (d->SPP_P[0] != 0.0 && d->SPP_Pstd[0] < 2 && !have_base) {
                ASSERT(!USE_RTD);
                double sin_el = sin(d->el);
                //0.125,0.7,0.35,1 is described in Global Positioning System：Theory and Applications(Volume I),pp481




                double istd = sin_el * sin_el / sqrt(SQR(d->SPP_Pstd[0])+d->ion_var * 0.125 * 0.125 + d->trop_var* 0.7 * 0.7  + d->sat_var* 0.35 * 0.35 + 1) / spp_psr_weight;

                ASSERT(istd > 0);
#define A8 para_pose[ir],para_gnss_dt+6+d->sys*2+0
#define B8 SppPseudorangeCeresFactor *factor = new SppPseudorangeCeresFactor(d->satellite_pos,d->SPP_P[0],istd,rover->base_xyz);
                ADDCERESRESIDUAL(A8, B8, false, std::vector<int> {})
#undef A8
#undef B8

            }
            if (USE_SPP_PHASE && d->SPP_L[0] != 0.0 && d->SPP_Npoint[0]) {
                double lam = lams[d->sys][0];
                double istd = 0;
                double sin_el = sin(d->el);

                istd = sin_el * sin_el / sqrt(SQR(d->SPP_Lstd[0] * lam) + (d->ion_var * 0.125 * 0.125 + d->trop_var * 0.7 * 0.7 + d->sat_var * 0.35 * 0.35));

                ASSERT(istd > 0);
#define A8s para_pose[ir],para_gnss_dt+6+d->sys*2+0,&(d->SPP_Npoint[0]->value)
#define B8s SppCarrierPhaseCeresFactor *factor = new SppCarrierPhaseCeresFactor(d->satellite_pos,d->SPP_L[0]*lam ,istd,rover->base_xyz,lam);
                ADDCERESRESIDUAL(A8s, B8s, false, std::vector<int> {})
#undef A8s
#undef B8s
            }

            if (USE_SPP_CORRECTION && d->SPP_P0[0] != 0.0 && d->SPP_Npoint_PCottections[0]) {
                double lam = lams[d->sys][0];
                double sin_el = sin(d->el);

                ASSERT(d->SPP_P[0] != 0);
                double istd = sin_el * sin_el / sqrt(SQR(d->SPP_Pstd[0])+d->ion_var * 0.125 * 0.125 + d->trop_var* 0.7 * 0.7  + d->sat_var* 0.35 * 0.35) / spp_psr_weight;


                ASSERT(istd > 0);
#define A8ss para_pose[ir],para_gnss_dt+6+d->sys*2+0,&(d->SPP_Npoint_PCottections[0]->value)
#define B8ss SppCarrierPhaseCeresFactor *factor = new SppCarrierPhaseCeresFactor(d->satellite_pos,d->SPP_P0[0] ,istd,rover->base_xyz,lam);
                ADDCERESRESIDUAL(A8ss, B8ss, false, std::vector<int> {})
#undef A8ss
#undef B8ss
            }


        }
        if (USE_DOPPLER) {
            for (int i = 0; i < rover->obs_count; i++) {
                ObsMea* d = rover->obs_data + i;
                if (d->SPP_D[0] == 0.0 || d->SVH != 0)continue;
                if (d->SPP_Dstd[0] > 2)continue;
                if (d->el < AZELMIN)continue;
                double istd = sin(d->el) * sin(d->el) / (d->SPP_Dstd[0] * lams[d->sys][0]);
                ASSERT(istd > 0);

#define A92 para_speed_bias[ir],para_gnss_dt+12,para_pose[ir]
#define B92 SppDopplerCeresFactor*factor=new SppDopplerCeresFactor(d->satellite_vel,d->satellite_pos,para_pose[0],d->SPP_D[0]*lams[d->sys][0],istd,rover->base_xyz);
                ADDCERESRESIDUAL(A92, B92, false, std::vector<int> {})
#undef A92
#undef B92
            }
        }
    }
}


//getting the pointer of the states.
void GNSSIMUInfo::Vector2Double() {


    for (int i = 0; i < (int)Rovers.size(); i++) {
        para_pose[i][0] = Pgs[i].x();
        para_pose[i][1] = Pgs[i].y();
        para_pose[i][2] = Pgs[i].z();

        Eigen::Quaterniond q{Rs[i]};
        para_pose[i][3] = q.x();
        para_pose[i][4] = q.y();
        para_pose[i][5] = q.z();
        para_pose[i][6] = q.w();

        para_speed_bias[i][0] = Vgs[i].x();
        para_speed_bias[i][1] = Vgs[i].y();
        para_speed_bias[i][2] = Vgs[i].z();

        para_speed_bias[i][3] = Bas[i].x();
        para_speed_bias[i][4] = Bas[i].y();
        para_speed_bias[i][5] = Bas[i].z();

        para_speed_bias[i][6] = Bgs[i].x();
        para_speed_bias[i][7] = Bgs[i].y();
        para_speed_bias[i][8] = Bgs[i].z();

    }

}


//saving the states from pointer.
void GNSSIMUInfo::Double2Vector() {
    for (int i = 0; i < (int)Rovers.size(); i++) {
        Rs[i] = Eigen::Quaterniond(para_pose[i][6], para_pose[i][3], para_pose[i][4], para_pose[i][5]).normalized().toRotationMatrix();
        Pgs[i] = Eigen::Vector3d(para_pose[i][0], para_pose[i][1], para_pose[i][2] ) ;
        Vgs[i] = Eigen::Vector3d(para_speed_bias[i][0], para_speed_bias[i][1], para_speed_bias[i][2]);
        Bas[i] = Eigen::Vector3d(para_speed_bias[i][3], para_speed_bias[i][4], para_speed_bias[i][5]);
        Bgs[i] = Eigen::Vector3d(para_speed_bias[i][6], para_speed_bias[i][7], para_speed_bias[i][8]);
    }
}



void GNSSIMUInfo::GnssProcess(mea_t* rover) {

    double globalxyz[3];
    // ASSERT(Pgs[Rovers.size() - 1].norm()>1000&&Pgs[Rovers.size() - 1].norm()<10000);
    globalxyz[0] = Pgs[Rovers.size() - 1].x() + rover->base_xyz[0];
    globalxyz[1] = Pgs[Rovers.size() - 1].y() + rover->base_xyz[1];
    globalxyz[2] = Pgs[Rovers.size() - 1].z() + rover->base_xyz[2];
    update_azel(globalxyz, rover);

    rover_count_accumulate++;
    rtk_fix = false;

    if (!init_gnss) {

        double pos[3];
        Eigen::Matrix3d Rwwg;
        if (ANCHOR_POINT.norm()) {
            LOG_OUT << "ANCHOR_POINT:" << ANCHOR_POINT.transpose() << std::endl;
            ecef2pos(ANCHOR_POINT.data(), pos);
            xyz2enu(pos, Rwwg.data());//transform vector from ecef to enu
        } else {
            ecef2pos(rover->base_xyz, pos);
            xyz2enu(pos, Rwwg.data());//transform vector from ecef to enu
        }

        InitRwgw = Rwwg.transpose();

        base_pos = Eigen::Map<Eigen::Vector3d> (rover->base_xyz);

        Vector2Double();
        ceres::Problem problem;
        ceres::Solver::Options options;
        AddCeresGnssResidual(std::set<double*> {}, 0, problem, options, Rovers.size() - 1);
        options.linear_solver_type = ceres::DENSE_SCHUR;
        // options.trust_region_strategy_type = ceres::DOGLEG;
        options.initial_trust_region_radius = options.max_trust_region_radius = 1e15;
        options.max_num_iterations = 20;
        options.num_threads = 1;
        TicToc t_marg;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        Double2Vector();

    }

    GnssPreprocess(Rovers[Rovers.size() - 1]);

}


void GNSSIMUInfo::CycleSlipDetection(mea_t* rover) {

    if (USE_RTK || USE_RTK) {
#if DISABLE_LONG_BASE_TIME!=0
        if (rover->br_time_diff > DISABLE_LONG_BASE_TIME) {
            for (int i = 0; i < rover->obs_count; i++) {
                ObsMea* d = rover->obs_data + i;
                for (int f = 0; f < NFREQ; f++) {
                    d -> RTK_P[f] = 0;
                    d -> RTK_L[f] = 0;
                }
            }
        }

#endif
    }
    if (USE_SPP_CORRECTION) {
        for (int i = 0; i < rover->obs_count; i++) {
            ObsMea* d = rover->obs_data + i;
            if (d->SPP_P[0] != 0) {
                d->SPP_P0[0] = d->SPP_P[0];
                if (pseudorange_correction_variables[d->sat * 2 + 0].size()) {
                    {
                        auto it = pseudorange_correction_variables[d->sat * 2 + 0].end();
                        it--;
                        it->last_update_time = rover->ros_time;
                        if (it->continue_count > 10)
                            ASSERT(it->value != 0);
                        ASSERT(pseudorange_correction_variables[d->sat * 2 + 0].size() == 1);
                        if (it->continue_count > EstimatePcorrectionPerio) {
                            d->SPP_P0[0] = 0;
                            d->SPP_P[0] += it->value * lams[d->sys][0];
                        } else
                            LOG_OUT << (int)(d->sat) << " is in MyOptimization for SPP corrections" << it->continue_count << "," << it->value << "\r\n";
                    }
                }
            } else
                d->SPP_P0[0] = 0;
        }
    }


    for (int i = 0; i < rover->obs_count; i++) {
        ObsMea* d = rover->obs_data + i;
        if (d->SVH)continue;
        for (uint f = 0; f < NFREQ; f++) {
            d->RTK_Npoint[f] = nullptr;
            d->SPP_Npoint[f] = nullptr;
            d->SPP_Npoint_PCottections[f] = nullptr;

            if (d->RTK_L[f] != 0) {
                if (rtk_phase_bias_variables[d->sat * 2 + f].size()) {//no cycle slip, using the lase phase bias
                    auto it = rtk_phase_bias_variables[d->sat * 2 + f].end();
                    it--;
                    if (rover->ros_time - it->last_update_time < 10)
                        d->RTK_Npoint[f] = &(*it);
                }
            }
            if (d->SPP_L[f] != 0) {
                if (spp_phase_bias_variables[d->sat * 2 + f].size()) {
                    auto it = spp_phase_bias_variables[d->sat * 2 + f].end();
                    it--;
                    if (rover->ros_time - it->last_update_time < 10)
                        d->SPP_Npoint[f] = &(*it);
                }
            }
            if (d-> SPP_P0[f] != 0) {
                if (pseudorange_correction_variables[d->sat * 2 + f].size()) {
                    auto it = pseudorange_correction_variables[d->sat * 2 + f].end();
                    it--;
                    if (rover->ros_time - it->last_update_time < 10)
                        d->SPP_Npoint_PCottections[f] = &(*it);
                }
            }
        }
    }

    Vector2Double();
    int ir = Rovers.size() - 1;

    std::vector<double>error1_rtk(rover->obs_count * 2, 0);
    std::vector<double>error2_rtk[6];
    double median_error_rtk[6];

    std::vector<double>error1_spp(rover->obs_count * 2, 0);
    std::vector<double>error2_spp[6];
    double median_error_spp[6];
    //detecting cycle slip
    for (int i = 0; i < rover->obs_count; i++) {
        ObsMea* d = rover->obs_data + i;
        if (d->SVH)continue;
        double* lam = lams[d->sys];
        int sys = d->sys;
        for (uint f = 0; f < NFREQ; f++) {
            if (d->el < AZELMIN)
                d->RTK_L[f] = d->SPP_L[f] = d->SPP_P0[f] = 0;
            if (d->RTK_Npoint[f]) {
                RTKCarrierPhaseCeresFactor* factor = new RTKCarrierPhaseCeresFactor(d->satellite_pos, d->RTK_L[f]*lam[f], lam[f], d->el, 0, 0, rover->base_xyz, false, d->sys, f);
                double residuals;
                std::vector<double*>parameter_blocks = std::vector<double*> {para_pose[ir], &(d->RTK_Npoint[f]->value), para_gnss_dt + sys * 2 + f};
                factor->Evaluate(parameter_blocks.data(), &residuals, 0);
                error1_rtk[i * 2 + f] = residuals;
                if (d->RTK_Npoint[f]->SLIP_COUNT == d->RTK_SLIP_COUNT[f])
                    error2_rtk[d->sys * 2 + f].push_back(residuals);
                delete factor;
            }

            if (d->SPP_Npoint[f]) {
                RTKCarrierPhaseCeresFactor* factor = new RTKCarrierPhaseCeresFactor(d->satellite_pos, d->SPP_L[f]*lam[f], lam[f], d->el, 0, 0, rover->base_xyz, false, d->sys, f);
                double residuals;
                std::vector<double*>parameter_blocks = std::vector<double*> {
                    para_pose[ir], &(d->SPP_Npoint[f]->value), para_gnss_dt + 6 + d->sys * 2 + 0
                };
                factor->Evaluate(parameter_blocks.data(), &residuals, 0);
                error1_spp[i * 2 + f] = residuals;
                if (d->SPP_Npoint[f]->SLIP_COUNT == d->SPP_SLIP_COUNT[f])
                    error2_spp[d->sys * 2 + f].push_back(residuals);
                delete factor;
            }
        }
    }
    for (int sys = 0; sys < 3; sys++) {
        for (int f = 0; f < NFREQ; f++) {
            if (error2_rtk[sys * 2 + f].size()) {
                std::sort(error2_rtk[sys * 2 + f].begin(), error2_rtk[sys * 2 + f].end());
                median_error_rtk[sys * 2 + f] = error2_rtk[sys * 2 + f][error2_rtk[sys * 2 + f].size() / 2];
                // if(error2_rtk[sys*2+f].size()<=4)median_error_rtk[sys*2+f]+=10;
            }
            if (error2_spp[sys * 2 + f].size()) {
                std::sort(error2_spp[sys * 2 + f].begin(), error2_spp[sys * 2 + f].end());
                median_error_spp[sys * 2 + f] = error2_spp[sys * 2 + f][error2_spp[sys * 2 + f].size() / 2];
                // if(error2_spp[sys*2+f].size()<=3)median_error_spp[sys*2+f]+=10;
            }
        }
    }
    if (not_fix_count > Phase_ALL_RESET_COUNT)LOG_OUT << "                                                                               not_fix_count reach threshold Phase_ALL_RESET_COUNT\r\n";
    //detecting cycle slip
    for (int i = 0; i < rover->obs_count; i++) {
        ObsMea* d = rover->obs_data + i;
        if (d->SVH)continue;
        double* lam = lams[d->sys];
        int sys = d->sys;
        for (uint f = 0; f < NFREQ; f++) {
            if (f == 1)ASSERT(d->RTK_L[f] == 0 && d->SPP_L[f] == 0);
            bool condition3 = false;
            bool condition4 = false;
            if (d->RTK_L[f] != 0) {
                if (USE_RTK && solver_flag == NonLinear && rover_count_accumulate > 1 && d->RTK_Npoint[f] && d->RTK_Npoint[f]->SLIP_COUNT == d->RTK_SLIP_COUNT[f]) {
                    double residuals = error1_rtk[i * 2 + f];
                    if (fabs(residuals - median_error_rtk[sys * 2 + f]) > lam[f] / 2) {
                        condition3 = true;
                        LOG_OUT << "residual_rtk:" << residuals << "," << median_error_rtk[sys * 2 + f] << "," << d->RTK_Npoint[f]->value << "," << (int)(d->sat) << "," << (d->RTK_Lstd[f]) << "," << (int)(d->RTK_Npoint[f]->continue_count) << "," << (int)(d->el / PI * 180) << std::endl;
                    }
                }
            }

            if (d->SPP_L[f] != 0) {
                if (USE_SPP_PHASE && solver_flag == NonLinear && rover_count_accumulate > 1 && d->SPP_Npoint[f] && d->SPP_Npoint[f]->SLIP_COUNT == d->SPP_SLIP_COUNT[f]) {
                    double residuals = error1_spp[i * 2 + f];
                    if (abs((d->SPP_L[f] + d->SPP_Npoint[f]->value)*lam[f] - d->SPP_P[f])*sin(d->el)*sin(d->el) > 10) {
                        condition4 = true;
                        LOG_OUT << "residual_spp2:" << (d->SPP_L[f] + d->SPP_Npoint[f]->value)*lam[f] << "," << d->SPP_P[f] << ","
                                << (d->SPP_L[f] + d->SPP_Npoint[f]->value)*lam[f] - d->SPP_P[f] << "," << d->SPP_Pstd[0] << ","
                                << (int)d->sat << "," << d->SPP_Npoint[f]->continue_count << "," << std::endl;
                    }
                    if (fabs(residuals - median_error_spp[sys * 2 + f]) > lam[f]) {
                        condition4 = true;
                        LOG_OUT << "residual_spp:" << residuals << "," << median_error_spp[sys * 2 + f] << "," << d->SPP_Npoint[f]->value << "," << (int)(d->sat) << "," << (d->SPP_Lstd[f]*lams[d->sys][f]) << "," << (int)(d->SPP_Npoint[f]->continue_count) << "," << (int)(d->el / PI * 180) << std::endl;
                    }
                }
            }

            if (d->RTK_L[f] != 0) {
                if ((!d->RTK_Npoint[f]) || (d->RTK_Npoint[f]->SLIP_COUNT != d->RTK_SLIP_COUNT[f]) || condition3 || not_fix_count > Phase_ALL_RESET_COUNT) {//cycle slip occur, creating new phase bias
                    PBtype n;
                    n.sys = d->sys;
                    n.f = f;
                    n.value = 0;
                    n.continue_count = 0;
                    rtk_phase_bias_variables[d->sat * 2 + f].push_back(n);
                    auto it = rtk_phase_bias_variables[d->sat * 2 + f].end();
                    it--;
                    d->RTK_Npoint[f] = &(*it);
                    d->RTK_Npoint[f]->SLIP_COUNT = d->RTK_SLIP_COUNT[f];
                    d->RTK_Npoint[f]->half_flag = d->half_flag[f];

                }
                if (d->RTK_Npoint[f]) {
                    d->RTK_Npoint[f]->last_update_time = rover->ros_time;
                    ASSERT(d->half_flag[f] == d->RTK_Npoint[f]->half_flag);
                }
            }


            if (d->SPP_L[f] != 0) {
                if ((!d->SPP_Npoint[f]) || (d->SPP_Npoint[f]->SLIP_COUNT != d->SPP_SLIP_COUNT[f]) || condition3 || condition4) {
                    PBtype n;
                    n.sys = d->sys;
                    n.f = f;
                    n.value = 0;
                    n.continue_count = 0;
                    spp_phase_bias_variables[d->sat * 2 + f].push_back(n);
                    auto it2 = spp_phase_bias_variables[d->sat * 2 + f].end();
                    it2--;
                    d->SPP_Npoint[f] = &(*it2);
                    d->SPP_Npoint[f]->SLIP_COUNT = d->SPP_SLIP_COUNT[f];
                    d->SPP_Npoint[f]->half_flag = d->half_flag[f];
                    d->SPP_Npoint[f]->first_observe_time = rover->ros_time;

                }

                if (d->SPP_Npoint[f])
                    d->SPP_Npoint[f]->last_update_time = rover->ros_time;

            }

            if (d->SPP_P0[f] != 0) {
                if (!d->SPP_Npoint_PCottections[f]) {
                    PBtype n;
                    n.sys = d->sys;
                    n.f = f;
                    n.value = 0;
                    n.continue_count = 0;
                    pseudorange_correction_variables[d->sat * 2 + f].push_back(n);
                    auto it2 = pseudorange_correction_variables[d->sat * 2 + f].end();
                    it2--;
                    d->SPP_Npoint_PCottections[f] = &(*it2);

                }

                if (d->SPP_Npoint_PCottections[f])
                    d->SPP_Npoint_PCottections[f]->last_update_time = rover->ros_time;
            }

            if (d->RTK_Npoint[f])d->RTK_Npoint[f]->continue_count++;
            if (d->SPP_Npoint[f])d->SPP_Npoint[f]->continue_count++;
            if (d->SPP_Npoint_PCottections[f])d->SPP_Npoint_PCottections[f]->continue_count++;

        }

    }

}

void GNSSIMUInfo::GenerateGNSSPriorFactor(mea_t* rover) {

    rover->marg_info_gnss = new MarginalizationInfo();

    std::set<double*>RemainPoint{para_pose[Rovers.size() - 1], para_speed_bias[Rovers.size() - 1], &blackvalue}; //posirion,velocity and blackvalue(no use)


    for (int i = 0; i < rover->obs_count; i++) {//phase bias
        ObsMea* d = rover->obs_data + i;
        for (int f = 0; f < NFREQ; f++) {
            if (d->RTK_Npoint[f])
                RemainPoint.insert(&(d->RTK_Npoint[f]->value));
            if (d->SPP_Npoint[f])
                RemainPoint.insert(&(d->SPP_Npoint[f]->value));
            if (d->SPP_Npoint_PCottections[f])
                RemainPoint.insert(&(d->SPP_Npoint_PCottections[f]->value));
        }
    }


    PhaseBiasSaveAndReset();//setting the values of phase biases as 0, for saving computation.
    AddGnssResidual(MargeExcludeMode, RemainPoint, rover->marg_info_gnss, Rovers.size() - 1);
    rover->marg_info_gnss->marginalize(true);
    rover->marg_info_gnss->getParameterBlocks(true);
    PhaseBiasRestore();//restoring the values of the phase biases.

    if (1) {
        TicToc tic;
        //we use the dogleg strategy for global MyOptimization. For stability, we need to initialize the phase biases near the true value.
        //If the levengerg strategy is used for global MyOptimization, this block can be disabled.
        // MarginalizationFactor *factor = new MarginalizationFactor(rover->marg_info_gnss);
        // problem.AddResidualBlock(factor, 0, rover->marg_info_gnss->keep_block_addr);
        ceres::Problem problem;
        ceres::Solver::Options options;
        AddCeresGnssResidual( std::set<double*> {}, 0, problem, options, Rovers.size() - 1);
        // bool needoptimize=false;
        for (int i = 0; i < rover->obs_count; i++) {
            ObsMea* d = rover->obs_data + i;
            for (int f = 0; f < NFREQ; f++) {

                if (d->RTK_Npoint[f] && problem.HasParameterBlock(&(d->RTK_Npoint[f]->value))) {
                    if (d->RTK_Npoint[f]->continue_count > 10)problem.SetParameterBlockConstant(&(d->RTK_Npoint[f]->value));
                    // else needoptimize=true;
                }

                if (d->SPP_Npoint[f] && problem.HasParameterBlock(&(d->SPP_Npoint[f]->value))) {
                    if (d->SPP_Npoint[f]->continue_count > 10)problem.SetParameterBlockConstant(&(d->SPP_Npoint[f]->value));
                    // else needoptimize=true;
                }

                if (d->SPP_Npoint_PCottections[f] && problem.HasParameterBlock(&(d->SPP_Npoint_PCottections[f]->value))) {
                    if (d->SPP_Npoint_PCottections[f]->continue_count > 10)problem.SetParameterBlockConstant(&(d->SPP_Npoint_PCottections[f]->value));
                    // else needoptimize=true;
                }

            }
        }

        problem.SetParameterBlockConstant(para_pose[Rovers.size() - 1]);
        problem.SetParameterBlockConstant(para_speed_bias[Rovers.size() - 1]);
        options.linear_solver_type = ceres::DENSE_SCHUR;
        options.initial_trust_region_radius = options.max_trust_region_radius = 1e15;
        options.max_num_iterations = 2;
        options.num_threads = 1;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        LOG_OUT << summary.BriefReport() << "," << tic.toc() << std::endl;

        for (int i = 0; i < rover->obs_count; i++) {
            ObsMea* d = rover->obs_data + i;
            for (int f = 0; f < NFREQ; f++) {

                if (d->RTK_Npoint[f] && problem.HasParameterBlock(&(d->RTK_Npoint[f]->value)))
                    ASSERT(fabs(d->RTK_Npoint[f]->value) < 1e3);

                if (d->SPP_Npoint[f] && problem.HasParameterBlock(&(d->SPP_Npoint[f]->value)))
                    ASSERT(fabs(d->SPP_Npoint[f]->value) < 1e4);

            }
        }
    }


}
void GNSSIMUInfo::GnssPreprocess(mea_t* rover) {


    TicToc t_initialize;


    CycleSlipDetection(rover);
    GenerateGNSSPriorFactor(rover);

    if (not_fix_count > Phase_ALL_RESET_COUNT)not_fix_count = 0;


}



void GNSSIMUInfo::PhaseBiasSaveAndReset() {
    phase_bias_save.clear();
    for (int i = 0; i < MAXSATNUM * 2; i++) {
        for (auto it1 = rtk_phase_bias_variables[i].begin(); it1 != rtk_phase_bias_variables[i].end(); it1++) {
            phase_bias_save.push_back(it1->value);
            it1->value = 0;
        }
        for (auto it1 = spp_phase_bias_variables[i].begin(); it1 != spp_phase_bias_variables[i].end(); it1++) {
            phase_bias_save.push_back(it1->value);
            it1->value = 0;
        }
        for (auto it1 = pseudorange_correction_variables[i].begin(); it1 != pseudorange_correction_variables[i].end(); it1++) {
            phase_bias_save.push_back(it1->value);
            it1->value = 0;
        }
    }
    phase_bias_save.push_back(blackvalue);
    blackvalue = 0;
}


void GNSSIMUInfo::PhaseBiasRestore() {
    int index = 0;
    for (int i = 0; i < MAXSATNUM * 2; i++) {
        for (auto it1 = rtk_phase_bias_variables[i].begin(); it1 != rtk_phase_bias_variables[i].end(); it1++)
            it1->value = phase_bias_save[index++];
        for (auto it1 = spp_phase_bias_variables[i].begin(); it1 != spp_phase_bias_variables[i].end(); it1++)
            it1->value = phase_bias_save[index++];
        for (auto it1 = pseudorange_correction_variables[i].begin(); it1 != pseudorange_correction_variables[i].end(); it1++)
            it1->value = phase_bias_save[index++];
    }
    blackvalue = phase_bias_save[index++];
    ASSERT(index == (int)phase_bias_save.size());
    phase_bias_save.clear();
}





void GNSSIMUInfo::AddGnssResidual(int mode, std::set<double*>MargePoint, MarginalizationInfo* marginalization_info, int ir) {

    //gnss residual
    if (USE_GNSS) {
        bool have_base = false;
        mea_t* rover = Rovers[ir];
        double globalxyz[3];
        globalxyz[0] = para_pose[ir][0] + rover->base_xyz[0];
        globalxyz[1] = para_pose[ir][1] + rover->base_xyz[1];
        globalxyz[2] = para_pose[ir][2] + rover->base_xyz[2];
        update_azel(globalxyz, rover);
        // ASSERT(Eigen::Vector3d(para_pose[ir]).norm() > 1000 && Eigen::Vector3d(para_pose[ir]).norm() < 10000);

#define A4s &blackvalue
#define B4s InitialBlackFactor2*factor=new InitialBlackFactor2(1);
        ADDRESIDUAL(A4s, B4s, false, std::vector<int> {})
        if (USE_RTK) {
            have_base = true;
            for (int i = 0; i < rover->obs_count; i++) {
                ObsMea* d = rover->obs_data + i;
                int sys = d->sys;
                for (int f = 0; f < NFREQ; f++) {
                    if (d->RTK_L[f] != 0 && rover_count_accumulate > 1)ASSERT(d->RTK_Npoint[f]);
                    if (!d->RTK_Npoint[f] )continue;
                    if (d->el < AZELMIN)continue;
                    ASSERT(d->el >= AZELMIN);
#define A44 para_pose[ir],&(d->RTK_Npoint[f]->value),para_gnss_dt+sys*2+f
#define B44 RTKCarrierPhaseFactor2 *factor = new RTKCarrierPhaseFactor2(d->satellite_pos,d->RTK_L[f]*lams[d->sys][f],lams[d->sys][f],d->el,rover->br_time_diff,pow(d->RTK_Lstd[f]*lams[d->sys][f],2),rover->base_xyz,true,d->sys,f);
                    ADDRESIDUAL(A44, B44, false, std::vector<int> {})

                }
            }
        }
        if (USE_RTD) {
            have_base = true;
            for (int i = 0; i < rover->obs_count; i++) {
                ObsMea* d = rover->obs_data + i;
                int sys = d->sys;
                ASSERT(sys == 0 || sys == 1 || sys == 2);
                for (int f = 0; f < NFREQ; f++) {
                    if (d->RTK_P[f] == 0.0 || d->SVH != 0  || d->RTK_Pstd[f] > 2)continue;
                    if (d->el < AZELMIN)continue;
                    ASSERT(d->el >= AZELMIN);
                    have_base = true;
#define A54 para_pose[ir],para_gnss_dt+sys*2+f
#define B54 RTKPseudorangeFactor2 *factor = new RTKPseudorangeFactor2(d->satellite_pos,d->RTK_P[f],  d->el,rover->br_time_diff,pow(d->RTK_Pstd[f],2), rover->base_xyz);
                    ADDRESIDUAL(A54, B54, false, std::vector<int> {})
                }
            }
        }


        for (int i = 0; i < rover->obs_count; i++) {
            ObsMea* d = rover->obs_data + i;
            if (d->SVH != 0)continue;
            if (d->el < AZELMIN)continue;
            if (d->SPP_P[0] != 0)ASSERT(d->el >= AZELMIN);

            if (d->SPP_P[0] != 0.0 && d->SPP_Pstd[0] < 2 && !have_base) {
                ASSERT(!USE_RTD);
                double sin_el = sin(d->el);

                double istd = sin_el * sin_el / sqrt(SQR(d->SPP_Pstd[0])+d->ion_var * 0.125 * 0.125 + d->trop_var* 0.7 * 0.7  + d->sat_var* 0.35 * 0.35 + 1) / spp_psr_weight;
                LOG_OUT << "gnss weight:" << sin_el << "," << d->SPP_Pstd[0] << "," << d->ion_var + d->trop_var + d->sat_var + 1 << "," << istd << std::endl;



#define A8 para_pose[ir],para_gnss_dt+6+d->sys*2+0
#define B8 SppPseudorangeFactor2 *factor = new SppPseudorangeFactor2(d->satellite_pos,d->SPP_P[0],istd,rover->base_xyz);
                ADDRESIDUAL(A8, B8, false, std::vector<int> {})


            }
            if (USE_SPP_PHASE && d->SPP_L[0] != 0.0 && d->SPP_Npoint[0]) {
                double lam = lams[d->sys][0];
                double istd = 0;
                double sin_el = sin(d->el);

                istd = sin_el * sin_el / sqrt(SQR(d->SPP_Lstd[0] * lam) + (d->ion_var * 0.125 * 0.125 + d->trop_var * 0.7 * 0.7 + d->sat_var * 0.35 * 0.35));


                ASSERT(istd > 0);
#define A8s para_pose[ir],para_gnss_dt+6+d->sys*2+0,&(d->SPP_Npoint[0]->value)
#define B8s SppCarrierPhaseFactor2 *factor = new SppCarrierPhaseFactor2(d->satellite_pos,d->SPP_L[0]*lam ,istd,rover->base_xyz,lam);
                ADDRESIDUAL(A8s, B8s, false, std::vector<int> {})
            }

            if (USE_SPP_CORRECTION && d->SPP_P0[0] != 0.0 && d->SPP_Npoint_PCottections[0]) {
                if (d->SPP_Pstd[0] < 2) {
                    double lam = lams[d->sys][0];
                    double sin_el = sin(d->el);

                    ASSERT(d->SPP_P[0] != 0);

                    double istd = sin_el * sin_el / sqrt(SQR(d->SPP_Pstd[0])+d->ion_var * 0.125 * 0.125 + d->trop_var* 0.7 * 0.7  + d->sat_var* 0.35 * 0.35) / spp_psr_weight;




                    ASSERT(istd > 0);
#define A8ss para_pose[ir],para_gnss_dt+6+d->sys*2+0,&(d->SPP_Npoint_PCottections[0]->value)
#define B8ss SppCarrierPhaseFactor2 *factor = new SppCarrierPhaseFactor2(d->satellite_pos,d->SPP_P0[0] ,istd,rover->base_xyz,lam);
                    ADDRESIDUAL(A8ss, B8ss, false, std::vector<int> {})
                }

            }


        }
        if (USE_DOPPLER) {
            for (int i = 0; i < rover->obs_count; i++) {
                ObsMea* d = rover->obs_data + i;
                if (d->SPP_D[0] == 0.0 || d->SVH != 0)continue;
                if (d->SPP_Dstd[0] > 2)continue;
                if (d->el < AZELMIN)continue;
                double istd = sin(d->el) * sin(d->el) / (d->SPP_Dstd[0] * lams[d->sys][0]);
                ASSERT(istd > 0);

#define A92 para_speed_bias[ir],para_gnss_dt+12,para_pose[ir]
#define B92 SppDopplerFactor2*factor=new SppDopplerFactor2(d->satellite_vel,d->satellite_pos,para_pose[0],d->SPP_D[0]*lams[d->sys][0],istd,rover->base_xyz);
                ADDRESIDUAL(A92, B92, false, std::vector<int> {})
            }
        }
    }


}





void FindReferenceSatellites(mea_t* rover, MarginalizationInfo* marginalization_info, ObsMea** Rsatellites) {

    double max_el[6] = {0};
    for (int i = 0; i < rover->obs_count; i++) {
        ObsMea* d = rover->obs_data + i;
        for (int f = 0; f < NFREQ; f++) {
            if (!d->RTK_Npoint[f] || d->RTK_Npoint[f]->use )continue;
            double* valuepoint = (&(d->RTK_Npoint[f]->value));
            if (marginalization_info->parameter_block_idx.find(valuepoint) != marginalization_info->parameter_block_idx.end()) {
                int a = marginalization_info->parameter_block_idx[valuepoint];
                if (a < 0)continue;

                if (d->el > max_el[d->sys * 2 + f]) {
                    max_el[d->sys * 2 + f] = d->el;
                    Rsatellites[d->sys * 2 + f] = d;
                }
            }
        }
    }
}

int gnss_fix_solution_count = 0;
Eigen::MatrixXd matrix_reorder(Eigen::MatrixXd A, std::vector<int>ordering) {
    assert(A.cols() == A.rows());
    assert(A.cols() >= (int)ordering.size());
    Eigen::MatrixXd B(A.cols(), A.cols());
    Eigen::MatrixXd C(A.cols(), A.cols());
    for (int i = 0; i < (int)ordering.size(); i++)
        B.block(i, 0, 1, A.cols()) = A.block(ordering[i], 0, 1, A.cols());
    for (int i = 0; i < (int)ordering.size(); i++)
        C.block(0, i, A.cols(), 1) = B.block(0, ordering[i], A.cols(), 1);
    return C.block(0, 0, ordering.size(), ordering.size());
}

Eigen::VectorXd vector_reorder(Eigen::VectorXd a, std::vector<int>ordering) {
    // LOG_OUT<<"0"<<std::endl;
    assert(a.size() >= (int)ordering.size());
    // LOG_OUT<<"1"<<std::endl;
    Eigen::VectorXd b(ordering.size());
    // LOG_OUT<<"2"<<std::endl;
    for (int i = 0; i < (int)ordering.size(); i++)
        b(i) = a(ordering[i]);
    return b;
}

void GetPosUsingN(mea_t* current_rover, double* current_pose, bool is_current) {
    ceres::Problem problem; ceres::Solver::Options options;
    double para_gnss_dt_tmp[13];
    memcpy(para_gnss_dt_tmp, para_gnss_dt, sizeof(para_gnss_dt_tmp));


    for (int i = 0; i < current_rover->obs_count; i++) {
        ObsMea* d = current_rover->obs_data + i;
        uint8_t sys = d->sys;
        for (uint8_t f = 0; f < NFREQ; f++) {
            if (!d->RTK_Npoint[f] )continue;
            if (d->el < AZELMIN)continue;
            ASSERT(d->el >= AZELMIN);


            if (is_current && d->RTK_Npoint[f]->last_fix_value) {
                d->RTK_Npoint[f]->tmp = d->RTK_Npoint[f]->last_fix_value;

                RTKCarrierPhaseCeresFactor* factor = new RTKCarrierPhaseCeresFactor(
                    d->satellite_pos,
                    d->RTK_L[f]*lams[d->sys][f],
                    lams[d->sys][f],
                    d->el,
                    current_rover->br_time_diff,
                    pow(d->RTK_Lstd[f]*lams[d->sys][f], 2),
                    current_rover->base_xyz,
                    true,
                    d->sys,
                    f);
                problem.AddResidualBlock(factor, 0, std::vector<double*>({current_pose, &(d->RTK_Npoint[f]->tmp), para_gnss_dt_tmp + sys * 2 + f}));

                InitialBlackCeresFactor* factor2 = new InitialBlackCeresFactor(1e5);
                problem.AddResidualBlock(factor2, 0, std::vector<double*>({ & (d->RTK_Npoint[f]->tmp)}));
            }


            if (!is_current && d->RTK_Npoint[f]->current_fxi_value) {
                d->RTK_Npoint[f]->tmp = d->RTK_Npoint[f]->current_fxi_value;

                RTKCarrierPhaseCeresFactor* factor = new RTKCarrierPhaseCeresFactor(
                    d->satellite_pos,
                    d->RTK_L[f]*lams[d->sys][f],
                    lams[d->sys][f],
                    d->el,
                    current_rover->br_time_diff,
                    pow(d->RTK_Lstd[f]*lams[d->sys][f], 2),
                    current_rover->base_xyz,
                    true,
                    d->sys,
                    f);
                problem.AddResidualBlock(factor, 0, std::vector<double*>({current_pose, &(d->RTK_Npoint[f]->tmp), para_gnss_dt_tmp + sys * 2 + f}));

                InitialBlackCeresFactor* factor2 = new InitialBlackCeresFactor(1e5);
                problem.AddResidualBlock(factor2, 0, std::vector<double*>({ & (d->RTK_Npoint[f]->tmp)}));
            }
        }
    }

    options.initial_trust_region_radius = options.max_trust_region_radius = 1000000000000000.0;
    options.max_num_iterations = 20;
    options.num_threads = 8;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

}
bool last_fix = false;

void GNSSIMUInfo::LambdaSearch(std::vector<mea_t*>& rovers, std::vector<double*>& parameter_block_addr_,
                               std::vector<int>& parameter_block_global_size_, Eigen::MatrixXd A_, Eigen::VectorXd& b_) {
    ASSERT(rovers.size());
    needNreaolve = false;
    // RegenerateGNSSIMUbase();
    static int fixcount;
    static double gnss_last_updatetime = 0;
    if (rovers[rovers.size() - 1]->ros_time - gnss_last_updatetime > 2) {
        fixcount = 0;
        last_fix = false;
        not_fix_count = 0;
        LOG_OUT << "reset FIX information" << std::endl;
    }
    gnss_last_updatetime = rovers[rovers.size() - 1]->ros_time;

    static MarginalizationInfo* marginalization_info = 0;
    marginalization_info = new MarginalizationInfo();
    marginalization_info->setmarginalizeinfo(parameter_block_addr_, parameter_block_global_size_, A_, b_, false);
    int old_size = -1;
    {
        std::set<double*>parameter_head;
        std::set<double*>marge_pointer;
        for (int i = 0; i < MAXSATNUM * 2; i++) {
            for (auto it1 = rtk_phase_bias_variables[i].begin(); it1 != rtk_phase_bias_variables[i].end(); it1++)
                it1->use = false;
        }
        for (int ir = 0; ir < (int)rovers.size(); ir++) {
            mea_t* rover = rovers[ir];
            for (int i = 0; i < rover->obs_count; i++) {
                ObsMea* d = rover->obs_data + i;
                for (int f = 0; f < NFREQ; f++) {
                    if (d->SVH == 0)
                        assert(d->el != 0);
                    if (d->RTK_Lstd[f] <= 0.004 * LTHRE && d->RTK_Npoint[f] &&
                            ((d->half_flag[f]) & 8) && ((d->half_flag[f]) & 2) && d->el >= AZELMIN) { //
                        if (!d->RTK_Npoint[f]->use) {
                            parameter_head.insert(&(d->RTK_Npoint[f]->value));
                            d->RTK_Npoint[f]->use = true;
                        }
                    }
                }
            }
        }
        if (parameter_head.size() == 0)return;
        for (auto it = marginalization_info->parameter_block_idx.begin(); it != marginalization_info->parameter_block_idx.end(); it++) {
            if (parameter_head.find(it->first) == parameter_head.end())
                marge_pointer.insert(it->first);
        }
        ASSERT(marge_pointer.size());
        if (marge_pointer.size() == marginalization_info->parameter_block_idx.size())return;
        marginalization_info->getParameterBlocks(false);
        marginalization_info->marginalize_pointers(marge_pointer, false);
        marginalization_info->parameter_block_idx.clear();
        ASSERT(marginalization_info->keep_block_addr.size() == marginalization_info->keep_block_idx.size());
        for (int i = 0; i < (int)marginalization_info->keep_block_idx.size(); i++)
            marginalization_info->parameter_block_idx[marginalization_info->keep_block_addr[i]] = marginalization_info->keep_block_idx[i];
        ASSERT(marginalization_info->parameter_block_idx.size());
        ASSERT(marginalization_info->parameter_block_idx.size() <= parameter_head.size());


        old_size = marginalization_info->parameter_block_idx.size();


    }

    if (marginalization_info->n < 6) {
        delete marginalization_info;
        return;
    }

    Eigen::MatrixXd Qy = marginalization_info->A.inverse();
    Eigen::MatrixXd y(marginalization_info->n, 1);

    for (int i = 0; i < MAXSATNUM * 2; i++) {
        for (auto Nv_it = rtk_phase_bias_variables[i].begin(); Nv_it != rtk_phase_bias_variables[i].end(); Nv_it++) {
            double* valuepoint = (&(Nv_it->value));
            if (marginalization_info->parameter_block_idx.find(valuepoint) != marginalization_info->parameter_block_idx.end())
                y(marginalization_info->parameter_block_idx[valuepoint]) = Nv_it->value;
        }
    }

    int countD = 0;
    Eigen::MatrixXd D3(marginalization_info->n, marginalization_info->n);
    D3.setZero();

    for (int i = 0; i < MAXSATNUM * 2; i++) {
        for (auto it1 = rtk_phase_bias_variables[i].begin(); it1 != rtk_phase_bias_variables[i].end(); it1++)
            it1->use = false;
    }

    int last_count = 0;
    int last_ref_count = 0;
    for (int ir = rovers.size() - 1; ir >= 0; ir--) {
        mea_t* rover = rovers[ir];
        // if(ir!=rover_count-1)continue;
        if (USE_RTK) {

            ObsMea* Rsatellites[6] = {0};
            FindReferenceSatellites(rover, marginalization_info, Rsatellites);
            for (int j = 0; j < 6; j++) {
                if (Rsatellites[j] && ir == (int)rovers.size() - 1)last_ref_count++;
            }
            for (int i = 0; i < rover->obs_count; i++) {
                ObsMea* d = rover->obs_data + i;
                int sys = d->sys;
                for (int f = 0; f < NFREQ; f++) {
                    if (!d->RTK_Npoint[f] )continue;
                    if (!Rsatellites[sys * 2 + f]) {
                        Rsatellites[sys * 2 + f] = d;
                        continue;
                    }
                    if (d->RTK_Npoint[f]->use || d == Rsatellites[sys * 2 + f])continue;
                    d->RTK_Npoint[f]->use = true;
                    double* valuepoint = (&(d->RTK_Npoint[f]->value));
                    double* bestvaluepoint = (&(Rsatellites[sys * 2 + f]->RTK_Npoint[f]->value));

                    if (marginalization_info->parameter_block_idx.find(valuepoint) != marginalization_info->parameter_block_idx.end()) {
                        int a = marginalization_info->parameter_block_idx[valuepoint];
                        int b = marginalization_info->parameter_block_idx[bestvaluepoint];

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
                            LOG_OUT << (int)Rsatellites[sys * 2 + f]->sat << "," << (int)d->sat << "," << (int)d->sys << std::endl;
                            if (ir == (int)rovers.size() - 1)last_count++;
                        } else {
                            LOG_OUT << "setNtozero" << (int)(d->sat) << "," << (int)(d->sys) << "," << (int)f << "," << y(a) - y(b) << "," << d->el << "," <<
                                    (int)(d->half_flag[f]) << "," << (int)(d->RTK_Npoint[f]->half_flag) << "," << std::endl;
                        }
                    }
                }
            }
        }
    }
    LOG_OUT << last_count << "," << last_ref_count << "," << countD << std::endl;
    if (last_count + last_ref_count < 6 || last_count < 4)return;

    last_fix = false;
    assert(countD <= marginalization_info->n);
    Eigen::MatrixXd D = D3.block(0, 0, countD, D3.cols());

    if (D.rows() < 4) {
        delete marginalization_info;
        LOG_OUT << "return" << std::endl;
        return;
    }
    Eigen::MatrixXd Qb = D * Qy * D.transpose();
    Eigen::MatrixXd b = D * y;
    LOG_OUT << b.transpose() << std::endl;
    ASSERT(old_size == (int)marginalization_info->parameter_block_idx.size());
    static int count_gns = 0;
    count_gns++;
    double s[2];
    double* F = mat(D.rows(), 2);
    bool search_ok = false;
    int info;


    if (!(info = lambda(Qb.rows(), 2, b.data(), Qb.data(), F, s))) { //nb:5//nb有正负1的行数
        LOG_OUT << s[1] / s[0] << "," << s[0] << "," << s[1] << std::endl;

        Eigen::VectorXd F1 = Eigen::Map<Eigen::VectorXd>(F, D.rows());
        Eigen::VectorXd F2 = Eigen::Map<Eigen::VectorXd>(F + D.rows(), D.rows());

        Eigen::VectorXd e1 = F1 - b;
        Eigen::VectorXd e2 = F2 - b;
        std::vector<int>same_index;
        std::vector<int>different_index;

        for (int i = 0; i < D.rows(); i++) {
            if (fabs(F1(i) - F2(i)) < 1e-2)
                same_index.push_back(i);
            else
                different_index.push_back(i);
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
        LOG_OUT << s1 / s0 << "," << s0 << "," << s1 << "," << different_index.size() << std::endl;
        if (s1 <= 0 || s0 <= 0) {
            std::cout << "s1<=0||s0<=0" << std::endl;
            return;
        }
        assert(s1 > 0 && s0 > 0);


        if (s[0] <= 0.0 || s[1] / s[0] >= 2 || s1 / s0 >= 2) {
            not_fix_count = 0;
            search_ok = true;
            gnss_fix_solution_count++;
            last_fix_time = Rovers[Rovers.size() - 1]->ros_time;
        } else
            not_fix_count += 1;
    }

    double tf[6] = {0};
    bool tfb[6] = {false};
    if (search_ok) {

        fixcount++;
        if (fixcount > FIX_CONTINUE_THRESHOLD) {

            last_fix = true;
            TicToc t_p;
            MarginalizationInfo* marginalization_info2 = new MarginalizationInfo();

            ResidualBlockInfo* residual_block_info = new ResidualBlockInfo(new MarginalizationFactor(Rovers[Rovers.size() - 1]->marg_info_gnss), NULL, Rovers[Rovers.size() - 1]->marg_info_gnss->keep_block_addr,
                                                                           std::vector<int> {}, std::vector<int> {});
            marginalization_info2->addResidualBlockInfo(residual_block_info);

            double* pindex[1000] = {0};
            PBtype* pindexNtype[1000] = {0};

            int numf1 = 0;
            int numf2 = 0;
            for (int i = 0; i < MAXSATNUM * 2; i++) {

                for (auto Nv_it = rtk_phase_bias_variables[i].begin(); Nv_it != rtk_phase_bias_variables[i].end(); Nv_it++) {
                    double* valuepoint = (&(Nv_it->value));
                    if (marginalization_info->parameter_block_idx.find(valuepoint) != marginalization_info->parameter_block_idx.end()) {
                        pindex[marginalization_info->parameter_block_idx[valuepoint]] = &(Nv_it->value);
                        pindexNtype[marginalization_info->parameter_block_idx[valuepoint]] = &(*Nv_it);
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
                    int sys = ppointNtype->sys;
                    int f = ppointNtype->f;

                    if (tfb[sys * 2 + f] == false) {
                        FixedIntegerFactor2* factormarge = new FixedIntegerFactor2(0, 1 / 0.03);
                        ResidualBlockInfo* residual_block_info = new ResidualBlockInfo(factormarge, NULL,
                                                                                       std::vector<double*> {&tf[sys * 2 + f], npoint}, std::vector<int> {0}, std::vector<int> {});
                        marginalization_info2->addResidualBlockInfo(residual_block_info);
                        tfb[sys * 2 + f] = true;
                    }

                    FixedIntegerFactor2* factormarge = new FixedIntegerFactor2(Fvalue, 1 / 0.03);
                    ResidualBlockInfo* residual_block_info = new ResidualBlockInfo(factormarge, NULL,
                                                                                   std::vector<double*> {&tf[sys * 2 + f], ppoint}, std::vector<int> {0}, std::vector<int> {});
                    marginalization_info2->addResidualBlockInfo(residual_block_info);

                }
            }

            mea_t* rovernew = rovers[rovers.size() - 1];
            rtk_fix = true;
            cout << "\033[32m";
            std::cout << "------------------------after lambda new pos:" << rovernew->ros_time << ",\t" << Pgs[Pgs.size() - 1].norm() << ",\t" <<  gnss_fix_solution_count << "," << count_gns << "------------------------" << std::endl;
            cout << "\033[0m";
            LOG_OUT << "------------------------after lambda new pos:" << rovernew->ros_time << ",\t" << Pgs[Pgs.size() - 1].norm() << ",\t" <<  gnss_fix_solution_count << "," << count_gns << "------------------------" << std::endl;

            PhaseBiasSaveAndReset();
            marginalization_info2->marginalize(false);
            marginalization_info2->getParameterBlocks(true);
            PhaseBiasRestore();
            delete Rovers[Rovers.size() - 1]->marg_info_gnss;
            Rovers[Rovers.size() - 1]->marg_info_gnss = marginalization_info2;
            RegenerateGNSSIMUbase();

        }

    } else {
        fixcount = 0;
        LOG_OUT << "reset" << std::endl;
    }



    delete marginalization_info;
    if (F)
        free(F);
}


