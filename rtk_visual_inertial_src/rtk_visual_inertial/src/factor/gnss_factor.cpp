
#include "gnss_factor.h"
#include "common_function.h"





bool SppPseudorangeFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
    const double* xyz = parameters[0];
    const double* dtgps = parameters[1];
    double xyzglobal[3];
    xyzglobal[0] = xyz[0] + base_pos[0];
    xyzglobal[1] = xyz[1] + base_pos[1];
    xyzglobal[2] = xyz[2] + base_pos[2];
    double sqrt_info = istd;
    double r1, e1[3];
    r1 = distance(xyzglobal, satelite_pos, e1);
    if (residuals)
        residuals[0] = sqrt_info * (r1 + *dtgps - P1);
    if (jacobians) {
        if (jacobians[0]) {
            memset(jacobians[0], 0, sizeof(double) * 7);
            jacobians[0][0] = sqrt_info * e1[0];
            jacobians[0][1] = sqrt_info * e1[1];
            jacobians[0][2] = sqrt_info * e1[2];
        }
        if (jacobians[1]) {
            jacobians[1][0] = sqrt_info * 1;
        }



    }
    // printf("%.2f,%.2f\r\n",r1+*dtgps-P1,sqrt_info);

    return true;

}





bool SppCarrierPhaseFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
    const double* xyz = parameters[0];
    const double* dtgps = parameters[1];
    const double* PB1 = parameters[2];
    double xyzglobal[3];
    xyzglobal[0] = xyz[0] + base_pos[0];
    xyzglobal[1] = xyz[1] + base_pos[1];
    xyzglobal[2] = xyz[2] + base_pos[2];
    double sqrt_info = istd;
    double r1, e1[3];
    r1 = distance(xyzglobal, satelite_pos, e1);
    if (residuals)
        residuals[0] = sqrt_info * (r1 + *dtgps - *PB1 * lam - L1_lam);
    if (jacobians) {

        if (jacobians[0]) {
            memset(jacobians[0], 0, sizeof(double) * 7);
            jacobians[0][0] = sqrt_info * e1[0];
            jacobians[0][1] = sqrt_info * e1[1];
            jacobians[0][2] = sqrt_info * e1[2];
        }
        if (jacobians[1]) {
            jacobians[1][0] = sqrt_info * 1;
        }
        if (jacobians[2]) {
            jacobians[2][0] = -sqrt_info * lam;
        }





    }
    return true;

}




bool FixedIntegerFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
    double sqrt_info = istd;
    if (residuals)
        residuals[0] = sqrt_info * ((parameters[1][0] - parameters[0][0]) - (N21));
    if (jacobians) {
        if (jacobians[0])
            jacobians[0][0] = -sqrt_info;
        if (jacobians[1])
            jacobians[1][0] = sqrt_info;
    }
    return true;
}

double varerr2(double el, double dt, double mea_var) {
    double b = clight * 5e-12 * dt;
    double sinel = sinf(el);
    return (mea_var / sinel / sinel ) + b * b;

}

bool RTKCarrierPhaseFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
    const double* xyz = parameters[0];
    const double* PB1 = parameters[1];
    const double* dtur = parameters[2];

    double r1, e1[3];
    double xyzglobal[3];
    xyzglobal[0] = xyz[0] + base_pos[0];
    xyzglobal[1] = xyz[1] + base_pos[1];
    xyzglobal[2] = xyz[2] + base_pos[2];
    r1 = distance(xyzglobal, satelite_pos, e1);
    double sqrt_info = 1;
    if (use_istd)sqrt_info = 1 / sqrt(varerr2(el, base_rover_time_diff, mea_var));
    if (residuals)
        residuals[0] = sqrt_info * (r1 - *PB1 * lam - L1_lam + *dtur);
    if (jacobians) {
        if (jacobians[0]) {
            memset(jacobians[0], 0, sizeof(double) * 7);
            jacobians[0][0] = sqrt_info * e1[0];
            jacobians[0][1] = sqrt_info * e1[1];
            jacobians[0][2] = sqrt_info * e1[2];
        }
        if (jacobians[1]) {
            jacobians[1][0] = -sqrt_info * lam;
        }
        if (jacobians[2]) {
            jacobians[2][0] = sqrt_info;
        }



    }
    return true;
}

bool RTKPseudorangeFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
    const double* xyz = parameters[0];
    const double* dtur = parameters[1];
    double xyzglobal[3];
    xyzglobal[0] = xyz[0] + base_pos[0];
    xyzglobal[1] = xyz[1] + base_pos[1];
    xyzglobal[2] = xyz[2] + base_pos[2];
    double r1, e1[3];
    r1 = distance(xyzglobal, satelite_pos, e1);
    double sqrt_info = 1 / sqrt(varerr2(el, base_rover_time_diff, mea_var));
    if (residuals) {
        residuals[0] = sqrt_info * (r1 - P1 + *dtur);
    }
    if (jacobians) {
        if (jacobians[0]) {
            memset(jacobians[0], 0, sizeof(double) * 7);
            jacobians[0][0] = sqrt_info * e1[0];
            jacobians[0][1] = sqrt_info * e1[1];
            jacobians[0][2] = sqrt_info * e1[2];
        }
        if (jacobians[1]) {
            jacobians[1][0] = sqrt_info;
        }



    }
    return true;
}





bool SppDopplerFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
    const double* vxyz = parameters[0];
    const double* dt = parameters[1];
    const double* xyz = parameters[2];
    double e[3], xyzglobal[3];
    xyzglobal[0] = xyz[0] + base_pos[0];
    xyzglobal[1] = xyz[1] + base_pos[1];
    xyzglobal[2] = xyz[2] + base_pos[2];
    double rate = velecitydistance(xyzglobal, satelite_pos, vxyz,  satelitev1,  e);

    double sqrt_info = istd;
    if (residuals)
        residuals[0] = sqrt_info * (rate + *dt + D1_lam);
    if (jacobians) {
        if (jacobians[0]) {
            memset(jacobians[0], 0, sizeof(double) * 9);
            jacobians[0][0] = sqrt_info * e[0];
            jacobians[0][1] = sqrt_info * e[1];
            jacobians[0][2] = sqrt_info * e[2];
        }
        if (jacobians[1]) {
            jacobians[1][0] = sqrt_info * 1;
        }
        if (jacobians[2]) {
            memset(jacobians[2], 0, sizeof(double) * 7);
            Eigen::Vector3d e2 = Eigen::Vector3d(xyzglobal) - Eigen::Vector3d(satelite_pos);
            double r = e2.norm();
            e2 /= r;
            Eigen::Vector3d ev = Eigen::Vector3d(vxyz) - Eigen::Vector3d(satelitev1);
            Eigen::MatrixXd J2 = sqrt_info * ev.transpose() * (Eigen::Matrix3d::Identity() - e2 * e2.transpose()) / r;
            jacobians[2][0] = J2(0);
            jacobians[2][1] = J2(1);
            jacobians[2][2] = J2(2);
        }

    }
    return true;

}
