#include "swf.h"
#include "ceres/residual_block.h"
#include "../gnss/gnss_ceres_factor.h"
#include "../gnss/gnss_factor.h"
// #include "../factor/gnss_ceres_factor.h"




void SWFOptimization::InputGnss(mea_t* rovernew2) {


    if (USE_GNSS) {
        static int count = 0;
        if (count++ < 2)return;
        mea_t* rovernew = (mea_t*)malloc(sizeof(mea_t));
        memcpy(rovernew, rovernew2, sizeof(mea_t));

        rover_buf.push_back(rovernew);
    }
}
void SWFOptimization::InputMag(double t, const Vector3d& Magdata) {
    if (USE_GNSS)
        mag_mean = Magdata;

}


//finding the avaliable gnss observations.
int SWFOptimization::getVariableUseSingleNum(mea_t* obs_data) {
    int i, j = 0;
    ObsMea* datai = obs_data->obs_data;
    for (i = 0; i < obs_data->obs_count; i++) {
        ObsMea* d = datai + i;
        if (d->SVH == 0 && d->SPP_P[0] != 0 && d->SPP_Pstd[0] < 2)
            j += 1;
    }
    return j;
}

