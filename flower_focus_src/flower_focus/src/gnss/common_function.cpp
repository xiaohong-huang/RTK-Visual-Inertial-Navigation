#include "common_function.h"

//wave length
double lams[3][2] = {
    { 0.190293672798364871256993069437, 0.244210213424568250983881512184 },
    { 0.19203948631027648, 0.24834936958430670 },
    { 0.19029367279836487, 0.24834936958430670 }
};



static int ludcmp(double* A, int n, int* indx, double* d) {
    double big, s, tmp, * vv = mat(n, 1);
    int i, imax = 0, j, k;

    *d = 1.0;
    for (i = 0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++) if ((tmp = fabs(A[i + j * n])) > big) big = tmp;
        if (big > 0.0) vv[i] = 1.0 / big;
        else {
            free(vv);
            return -1;
        }
    }
    for (j = 0; j < n; j++) {
        for (i = 0; i < j; i++) {
            s = A[i + j * n];
            for (k = 0; k < i; k++) s -= A[i + k * n] * A[k + j * n];
            A[i + j * n] = s;
        }
        big = 0.0;
        for (i = j; i < n; i++) {
            s = A[i + j * n];
            for (k = 0; k < j; k++) s -= A[i + k * n] * A[k + j * n];
            A[i + j * n] = s;
            if ((tmp = vv[i] * fabs(s)) >= big) {
                big = tmp;
                imax = i;
            }
        }
        if (j != imax) {
            for (k = 0; k < n; k++) {
                tmp = A[imax + k * n];
                A[imax + k * n] = A[j + k * n];
                A[j + k * n] = tmp;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (A[j + j * n] == 0.0) {
            free(vv);
            return -1;
        }
        if (j != n - 1) {
            tmp = 1.0 / A[j + j * n];
            for (i = j + 1; i < n; i++) A[i + j * n] *= tmp;
        }
    }
    free(vv);
    return 0;
}


static void lubksb(const double* A, int n, const int* indx, double* b) {
    double s;
    int i, ii = -1, ip, j;

    for (i = 0; i < n; i++) {
        ip = indx[i];
        s = b[ip];
        b[ip] = b[i];
        if (ii >= 0) for (j = ii; j < i; j++) s -= A[i + j * n] * b[j];
        else if (s) ii = i;
        b[i] = s;
    }
    for (i = n - 1; i >= 0; i--) {
        s = b[i];
        for (j = i + 1; j < n; j++) s -= A[i + j * n] * b[j];
        b[i] = s / A[i + i * n];
    }
}


double satazel(const double* pos, const double* e, double* azel) {
    double az = 0.0, el = PI / 2.0, enu[3];

    if (pos[2] > -RE_WGS84) {
        ecef2enu(pos, e, enu);
        az = dot(enu, enu, 2) < 1E-12 ? 0.0 : atan2(enu[0], enu[1]);
        if (az < 0.0) az += 2 * PI;
        el = asin(enu[2]);
    }
    if (azel) {
        azel[0] = az;
        azel[1] = el;
    }
    return el;
}


double dot(const double* a, const double* b, int n) {
    double c = 0.0;

    while (--n >= 0) c += a[n] * b[n];
    return c;
}


void ecef2pos(const double* r, double* pos) {
    double e2 = FE_WGS84 * (2.0 - FE_WGS84), r2 = dot(r, r, 2), z, zk, v = RE_WGS84, sinp;

    for (z = r[2], zk = 0.0; fabs(z - zk) >= 1E-4;) {
        zk = z;
        sinp = z / sqrt(r2 + z * z);
        v = RE_WGS84 / sqrt(1.0 - e2 * sinp * sinp);
        z = r[2] + v * e2 * sinp;
    }
    pos[0] = r2 > 1E-12 ? atan(z / sqrt(r2)) : (r[2] > 0.0 ? PI / 2.0 : -PI / 2.0);
    pos[1] = r2 > 1E-12 ? atan2(r[1], r[0]) : 0.0;
    pos[2] = sqrt(r2 + z * z) - v;
}


double distance(const double* rr, const double* rs, double* e) {
    double r;
    uint8_t i;

    for (i = 0; i < 3; i++) e[i] = rr[i] - rs[i];
    r = norm(e, 3);
    for (i = 0; i < 3; i++) e[i] /= r;
    return r + OMGE * (rs[0] * rr[1] - rs[1] * rr[0]) / clight;
}


double norm(const double* a, int n) {
    return sqrt(dot(a, a, n));
}


void ecef2enu(const double* pos, const double* r, double* e) {
    double E[9];

    xyz2enu(pos, E);
    matmul("NN", 3, 1, 3, 1.0, E, r, 0.0, e);
}


void xyz2enu(const double* pos, double* E) {
    double sinp = sin(pos[0]), cosp = cos(pos[0]), sinl = sin(pos[1]), cosl = cos(pos[1]);

    E[0] = -sinl;
    E[3] = cosl;
    E[6] = 0.0;
    E[1] = -sinp * cosl;
    E[4] = -sinp * sinl;
    E[7] = cosp;
    E[2] = cosp * cosl;
    E[5] = cosp * sinl;
    E[8] = sinp;
}


void matmul(const char* tr, int n, int k, int m, double alpha,
            const double* A, const double* B, double beta, double* C) {
    double d;
    int i, j, x, f = tr[0] == 'N' ? (tr[1] == 'N' ? 1 : 2) : (tr[1] == 'N' ? 3 : 4);

    for (i = 0; i < n; i++) for (j = 0; j < k; j++) {
            d = 0.0;
            switch (f) {
            case 1:
                for (x = 0; x < m; x++) d += A[i + x * n] * B[x + j * m];
                break;
            case 2:
                for (x = 0; x < m; x++) d += A[i + x * n] * B[j + x * k];
                break;
            case 3:
                for (x = 0; x < m; x++) d += A[x + i * m] * B[x + j * m];
                break;
            case 4:
                for (x = 0; x < m; x++) d += A[x + i * m] * B[j + x * k];
                break;
            }
            if (beta == 0.0) C[i + j * n] = alpha * d;
            else C[i + j * n] = alpha * d + beta * C[i + j * n];
        }
}


void matmul6xx6(double* A, double* B, uint8_t n) {
    uint8_t i, k, l;
    uint16_t j1, j2;
    memset(B, 0, (6 * 6)*sizeof(double));
    for (l = 0; l < 6; l++) {
        for (k = l; k < 6; k++) {
            j1 = k;
            j2 = l;
            for (i = 0; i < n; i++) {
                B[l * 6 + k] += A[j1] * A[j2];
                j1 += 6;
                j2 += 6;
            }
        }
    }
    for (l = 0; l < 6; l++) {
        for (k = l + 1; k < 6; k++) {
            B[k * 6 + l] = B[l * 6 + k];
        }
    }
}


void matmul61(double* A, double* v, double* B, uint8_t n) {
    uint8_t i, l;
    uint16_t j1, j2;
    memset(B, 0, 6 * sizeof(double));
    for (l = 0; l < 6; l++) {
        j1 = l;
        j2 = 0;
        for (i = 0; i < n; i++) {
            B[l] += A[j1] * v[j2];
            j1 += 6;
            j2 += 1;
        }
    }
}


void matmul3xx3(double* A, double* B, uint8_t n) {
    uint8_t i, k, l;
    uint16_t j1, j2;
    memset(B, 0, (3 * 3) * sizeof(double));
    for (l = 0; l < 3; l++) {
        for (k = l; k < 3; k++) {
            j1 = k;
            j2 = l;
            for (i = 0; i < n; i++) {
                B[l * 3 + k] += A[j1] * A[j2];
                j1 += 3;
                j2 += 3;
            }
        }
    }
    for (l = 0; l < 3; l++) {
        for (k = l + 1; k < 3; k++) {
            B[k * 3 + l] = B[l * 3 + k];
        }
    }
}


void matmul31(double* A, double* v, double* B, uint8_t n) {
    uint8_t i, l;
    uint16_t j1, j2;
    memset(B, 0, 3 * sizeof(double));
    for (l = 0; l < 3; l++) {
        j1 = l;
        j2 = 0;
        for (i = 0; i < n; i++) {
            B[l] += A[j1] * v[j2];
            j1 += 3;
            j2 += 1;
        }
    }
}


uint8_t lsq6( double* A,  double* y, int n, int m, double* x, double* Q) {

    uint8_t info;
    if (m < n) return 0;
    double Ay[6];
    matmul61(A, y, Ay, m);
    matmul6xx6(A, Q, m);
    info = matinv(Q, n);
    if (info)
        matmul61(Q, Ay, x, 6);
    return info;
}


uint8_t lsq3(double* A, double* y, int n, int m, double* x, double* Q) {
    uint8_t info;
    if (m < n) return 0;
    double Ay[3];
    double Qinv[9];
    matmul31(A, y, Ay, m);
    matmul3xx3(A, Q, m);
    info = matinv33(Q, Qinv);
    if (info)
        matmul31(Qinv, Ay, x, 3);
    return info;
}


int lsq(double* A, double* y, int n, int m, double* x, double* Q) {
    double* Ay;
    int info;

    if (m < n) return 0;
    Ay = mat(n, 1);
    matmul("NN", n, 1, m, 1.0, A, y, 0.0, Ay); /* Ay=A*y */
    matmul("NT", n, n, m, 1.0, A, A, 0.0, Q);  /* Q=A*A' */
    info = matinv(Q, n);
    if (info) matmul("NN", n, 1, n, 1.0, Q, Ay, 0.0, x); /* x=Q^-1*Ay */
    free(Ay);
    return info;
}


uint8_t matinv33(double* a, double* b) {
    double det = a[0] * a[4] * a[8] + a[1] * a[5] * a[6] + a[3] * a[7] * a[2] - a[2] * a[4] * a[6] - a[5] * a[7] * a[0] - a[1] * a[3] * a[8];
    if (ABS(det) < 1e-2)return 0;
    if (det < 0)printf("wrong!!\r\n");
    b[0] = (a[4] * a[8] - a[7] * a[5]) / det;
    b[1] = (a[2] * a[7] - a[1] * a[8]) / det;
    b[2] = (a[1] * a[5] - a[2] * a[4]) / det;
    b[3] = (a[6] * a[5] - a[3] * a[8]) / det;
    b[4] = (a[0] * a[8] - a[2] * a[6]) / det;
    b[5] = (a[2] * a[3] - a[0] * a[5]) / det;
    b[6] = (a[3] * a[7] - a[4] * a[6]) / det;
    b[7] = (a[1] * a[6] - a[0] * a[7]) / det;
    b[8] = (a[0] * a[4] - a[1] * a[3]) / det;
    return 1;
}


uint32_t ca, cb;
uint8_t matinv33f(float* a, float* b) {
    float det = a[0] * a[4] * a[8] + a[1] * a[5] * a[6] + a[3] * a[7] * a[2] - a[2] * a[4] * a[6] - a[5] * a[7] * a[0] - a[1] * a[3] * a[8];
    if (ABS(det) < 1e-1f)
        return 0;
    b[0] = (a[4] * a[8] - a[7] * a[5]) / det;
    b[1] = (a[2] * a[7] - a[1] * a[8]) / det;
    b[2] = (a[1] * a[5] - a[2] * a[4]) / det;
    b[3] = (a[6] * a[5] - a[3] * a[8]) / det;
    b[4] = (a[0] * a[8] - a[2] * a[6]) / det;
    b[5] = (a[2] * a[3] - a[0] * a[5]) / det;
    b[6] = (a[3] * a[7] - a[4] * a[6]) / det;
    b[7] = (a[1] * a[6] - a[0] * a[7]) / det;
    b[8] = (a[0] * a[4] - a[1] * a[3]) / det;
    return 1;
}


uint8_t matinv(double* A, int n) {
    double d, * B;
    int i, j, * indx;

    indx = imat(n, 1);
    B = mat(n, n);
    matcpy(B, A, n, n);
    if (ludcmp(B, n, indx, &d)) {
        free(indx);
        free(B);
        return 0;
    }
    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) A[i + j * n] = 0.0;
        A[j + j * n] = 1.0;
        lubksb(B, n, indx, A + j * n);
    }
    free(indx);
    free(B);
    return 1;
}


int* imat(int n, int m) {
    int* p;

    if (n <= 0 || m <= 0) return NULL;
    p = (int*)malloc(sizeof(int) * n * m);
    return p;
}


double* mat(int n, int m) {
    double* p;

    if (n <= 0 || m <= 0) return NULL;
    p = (double*)malloc(sizeof(double) * n * m);
    return p;
}


void matcpy(double* A, const double* B, int n, int m) {
    memcpy(A, B, sizeof(double) * n * m);
}


void update_azel(double globalxyz[3], mea_t* rover) {
    for (uint8_t i = 0; i < rover->obs_count; i++) {
        ObsMea* d = rover->obs_data + i;
        if (d->SVH != 0)continue;
        double pos[3], e2[3], azel[2], e[3];
        ecef2pos(globalxyz, pos);
        distance(globalxyz, d->satellite_pos, e);
        e2[0] = -e[0];
        e2[1] = -e[1];
        e2[2] = -e[2];
        satazel(pos, e2, azel);
        d->el = azel[1];

    }
}


double velecitydistance(const double* rr, const double* rs, const double* vr, const double* vs, double* e) {
    double r;
    uint8_t i;
    double ev[3];
    for (i = 0; i < 3; i++) e[i] = rr[i] - rs[i];
    r = norm(e, 3);
    for (i = 0; i < 3; i++) e[i] /= r;
    for (i = 0; i < 3; i++)ev[i] = vr[i] - vs[i];

    return dot(ev, e, 3) + OMGE / clight * (vs[1] * rr[0] + rs[1] * vr[0] - vs[0] * rr[1] - rs[0] * vr[1]);
}

