/******************************
File name: tinympc_task.cpp
Description: 原 H7 任务提取版；保留模型拟合、求解和输出逻辑。
整车状态、计时与临界区依赖改为显式输入和 tinympc_port 接口。
*****************************/

#include "tinympc_task.h"
#include "admm.hpp"
#include <cstring>                  // memset
#include <cstdint>
#include <cmath>                    // std::isfinite

#if defined(__clang__)
#define TINYMPC_RUN_CRASH_WORKAROUND __attribute__((optnone, noinline))
#else
#define TINYMPC_RUN_CRASH_WORKAROUND
#endif



// ============================================================================
// Global task pointer - created in board.cpp via tinympc_test_task_create()
// ============================================================================
TinyMPC_Task* tinympc_ptr = nullptr;

// ============================================================================
// MPC 实时性观测量
//   g_mpc_iter   : 本拍 ADMM 迭代次数（work.iter）
//   g_mpc_solved : 本拍求解状态，1=收敛 / 11=撞 max_iter 未收敛（work.status）
// ============================================================================
volatile float g_mpc_iter   = 0.0f;
volatile float g_mpc_solved = 0.0f;
volatile uint32_t g_mpc_total_us = 0;
volatile uint32_t g_mpc_model_us = 0;
volatile uint32_t g_mpc_solve_us = 0;
volatile uint32_t g_mpc_max_total_us = 0;
volatile uint32_t g_mpc_over_2ms_count = 0;
volatile uint32_t g_mpc_over_4ms_count = 0;

// ============================================================================
// 原轮腿启动参数；正常 run 会使用腿长拟合结果更新。
// ============================================================================

// State transition matrix A (10x10)
tinytype Adyn_data[NSTATES * NSTATES] = {
    // Row 0
    1.0000f, 0.0040f, 0.0f,    0.0f,    -0.0004f, 0.0f,    -0.0004f, 0.0f,    0.0f,    0.0f,
    // Row 1
    0.0f,    1.0000f, 0.0f,    0.0f,    -0.1949f, -0.0004f, -0.1949f, -0.0004f, 0.0f,    0.0f,
    // Row 2
    0.0f,    0.0f,    1.0000f, 0.0040f, 0.0f,     0.0f,     0.0f,     0.0f,     0.0f,    0.0f,
    // Row 3
    0.0f,    0.0f,    0.0f,    1.0000f, -0.0163f, 0.0f,     0.0163f,  0.0f,     0.0f,    0.0f,
    // Row 4
    0.0f,    0.0f,    0.0f,    0.0f,    1.0022f,  0.0040f,  0.0001f,  0.0f,     0.0f,    0.0f,
    // Row 5
    0.0f,    0.0f,    0.0f,    0.0f,    1.1101f,  1.0022f,  0.0580f,  0.0001f,  0.0f,    0.0f,
    // Row 6
    0.0f,    0.0f,    0.0f,    0.0f,    0.0001f,  0.0f,     1.0022f,  0.0040f,  0.0f,    0.0f,
    // Row 7
    0.0f,    0.0f,    0.0f,    0.0f,    0.0580f,  0.0001f,  1.1101f,  1.0022f,  0.0f,    0.0f,
    // Row 8
    0.0f,    0.0f,    0.0f,    0.0f,    0.0f,     0.0f,     0.0f,     0.0f,     1.0001f, 0.0040f,
    // Row 9
    0.0f,    0.0f,    0.0f,    0.0f,    -0.0161f, 0.0f,     -0.0161f, 0.0f,     0.0277f, 1.0001f
};

// Control matrix B (10x4)
tinytype Bdyn_data[NSTATES * NINPUTS] = {
    // Row 0
    0.0001f,  0.0001f,  0.0f,     0.0f,
    // Row 1
    0.0322f,  0.0322f, -0.0048f, -0.0048f,
    // Row 2
    0.0f,     0.0f,     0.0f,     0.0f,
    // Row 3
   -0.0142f,  0.0142f, -0.0004f,  0.0004f,
    // Row 4
   -0.0003f,  0.0f,     0.0001f,  0.0f,
    // Row 5
   -0.1641f, -0.0137f,  0.0275f,  0.0014f,
    // Row 6
    0.0f,    -0.0003f,  0.0f,     0.0001f,
    // Row 7
   -0.0137f, -0.1641f,  0.0014f,  0.0275f,
    // Row 8
    0.0f,     0.0f,     0.0f,     0.0f,
    // Row 9
    0.0007f,  0.0007f, -0.0075f, -0.0075f
};

// Affine vector f (no gravity term for this system)
tinytype fdyn_data[NSTATES] = {
    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
};


// Pre-computed infinite horizon LQR gain Kinf (4x10)
// 初始化: 无穷时域LQR增益矩阵 (使用简单的方阵分解)
tinytype Kinf_data[NINPUTS * NSTATES] = {
    // 4 rows x 10 columns
  -6.4233 , -11.3516 , -3.0662 , -1.2487 , -73.6683 , -4.0772 , -3.8077 , -2.1889 , -11.5755 , -1.6879 , 
  -6.4233 , -11.3516 , 3.0662 , 1.2487 , -3.8077 , -2.1889 , -73.6683 , -4.0772 , -11.5755 , -1.6879 , 
  5.7574 , 9.1717 , -19.4468 , -8.1358 , 131.9737 , 5.4709 , -16.8440 , 0.5308 , -66.8076 , -5.7546 , 
  5.7574 , 9.1717 , 19.4468 , 8.1358 , -16.8440 , 0.5308 , 131.9737 , 5.4709 , -66.8076 , -5.7546 
};

// Pre-computed infinite horizon cost-to-go Pinf (10x10)
// 初始化: 无穷时域代价矩阵 Pinf (单位已×100)
tinytype Pinf_data[NSTATES * NSTATES] = {
    // Row 0
    70057.0f, 59986.0f, 0.0f,     0.0f,     19124.0f, 10916.0f, 19124.0f, 10916.0f, 5880.0f,  3364.0f,
    // Row 1
    59986.0f, 99770.0f, 0.0f,     0.0f,     22532.0f, 18159.0f, 22532.0f, 18159.0f, 7293.0f,  5617.0f,
    // Row 2
    0.0f,     0.0f,     10404.0f, 895.0f,   -995.0f,  -63.0f,   995.0f,   63.0f,    0.0f,     0.0f,
    // Row 3
    0.0f,     0.0f,     895.0f,   375.0f,   -377.0f,  -26.0f,   377.0f,   26.0f,    0.0f,     0.0f,
    // Row 4
    19124.0f, 22532.0f, -995.0f,  -377.0f,  82946.0f, 4728.0f,  5454.0f,  4035.0f,  3424.0f,  1358.0f,
    // Row 5
    10916.0f, 18159.0f, -63.0f,   -26.0f,   4728.0f,  3325.0f,  4035.0f,  3303.0f,  1399.0f,  1029.0f,
    // Row 6
    19124.0f, 22532.0f, 995.0f,   377.0f,   5454.0f,  4035.0f,  82946.0f, 4728.0f,  3424.0f,  1358.0f,
    // Row 7
    10916.0f, 18159.0f, 63.0f,    26.0f,    4035.0f,  3303.0f,  4728.0f,  3325.0f,  1399.0f,  1029.0f,
    // Row 8
    5880.0f,  7293.0f,  0.0f,     0.0f,     3424.0f,  1399.0f,  3424.0f,  1399.0f,  29567.0f, 1622.0f,
    // Row 9
    3364.0f,  5617.0f,  0.0f,     0.0f,     1358.0f,  1029.0f,  1358.0f,  1029.0f,  1622.0f,  428.0f
};

// Pre-computed Quu inverse (4x4)
tinytype Quu_inv_data[NINPUTS * NINPUTS] = {
    // Row 0
    0.3831 ,  -0.0100 ,   0.0281,    0.0002,
   -0.0100 ,   0.3831 ,   0.0002 ,   0.0281,
    0.0281 ,   0.0002 ,   0.8939 ,  -0.0020,
    0.0002 ,   0.0281 ,  -0.0020 ,   0.8939
};

// Pre-computed (A - B*Kinf)^T (10x10)
// 初始化: 闭环状态转移矩阵 (A-B*Kinf)^T (单位已缩放)
tinytype AmBKt_data[NSTATES * NSTATES] = {
    // Row 0
    1.0009f, 0.4689f, 0.0f,    0.0f,   -0.0026f, -1.3090f, -0.0026f, -1.3090f, 0.0002f, 0.0958f,
    // Row 1
    0.0056f, 1.8190f, 0.0f,    0.0f,   -0.0046f, -2.2842f, -0.0046f, -2.2842f, 0.0003f, 0.1541f,
    // Row 2
    0.0f,    0.0f,    0.9998f, -0.1026f, 0.0001f, 0.0465f, -0.0001f, -0.0465f, 0.0f,    0.0f,
    // Row 3
    0.0f,    0.0f,    0.0039f, 0.9580f,  0.0f,    0.0246f, 0.0f,     -0.0246f, 0.0f,    0.0f,
    // Row 4
    0.0057f, 2.8540f, -0.0019f, -0.9467f, 0.9707f, -14.6446f, -0.0026f, -1.3002f, 0.0018f, 0.9059f,
    // Row 5
    0.0005f, 0.2302f, 0.0f,    -0.0248f, 0.0023f, 0.1516f, -0.0009f, -0.4374f, 0.0001f, 0.0496f,
    // Row 6
    0.0057f, 2.8540f, 0.0019f, 0.9467f, -0.0026f, -1.3002f, 0.9707f, -14.6446f, 0.0018f, 0.9059f,
    // Row 7
    0.0005f, 0.2302f, 0.0f,    0.0248f, -0.0009f, -0.4374f, 0.0023f, 0.1516f, 0.0001f, 0.0496f,
    // Row 8
    0.0002f, 0.0984f, 0.0f,    0.0f,   -0.0002f, -0.1216f, -0.0002f, -0.1216f, 0.9981f, -0.9642f,
    // Row 9
    0.0001f, 0.0529f, 0.0f,    0.0f,   -0.0003f, -0.1333f, -0.0003f, -0.1333f, 0.0038f, 0.9156f
};

// Pre-computed affine terms
tinytype APf_data[NSTATES] = {
    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
};

tinytype BPf_data[NINPUTS] = {
    0.0f, 0.0f, 0.0f, 0.0f
};

// Fitted static feed-forward input reference at q=0, qd=0, qdd=0
tinytype Ueq_ref_data[NINPUTS] = {0.0f, 0.0f, 0.0f, 0.0f};

// Initial state for simulation
tinytype xinit_data[NSTATES] = {
    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
};

// Goal state
tinytype xgoal_data[NSTATES] = {
    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
};

// Input constraints
static const tinytype u_min_data[NINPUTS] = {
    -6.0f, -6.0f, -45.0f, -45.0f
};

static const tinytype u_max_data[NINPUTS] = {
    6.0f, 6.0f, 45.0f, 45.0f
};

// -----------------------------------------------------------------------------
// Helpers: evaluate 2D polynomial on (l_left, l_right) using exponent table.
// Coef order follows poly_exp_p/poly_exp_q from leg_fit_coeffs_sparse.h
// -----------------------------------------------------------------------------
static inline void poly2d_prepare(tinytype x1, tinytype x2, tinytype* pow1, tinytype* pow2) {
    pow1[0] = 1.0f;
    pow2[0] = 1.0f;
    for (int i = 1; i <= POLY_MAX_DEG; i++) {
        pow1[i] = pow1[i - 1] * x1;
        pow2[i] = pow2[i - 1] * x2;
    }
}

static inline tinytype poly2d_eval(const tinytype* coef, const tinytype* pow1, const tinytype* pow2) {
    tinytype sum = 0.0f;
    for (size_t i = 0; i < POLY_TERM_COUNT; i++) {
        sum += coef[i] * pow1[poly_exp_p[i]] * pow2[poly_exp_q[i]];
    }
    return sum;
}

static inline bool all_finite(const tinytype* a, int n) {
    for (int i = 0; i < n; i++) {
        if (!std::isfinite(a[i])) {
            return false;
        }
    }
    return true;
}

static inline void set_by_matlab_index(tinytype* M_row_major, int R, int C, uint16_t idx1, tinytype v) {
    // MATLAB idx: 1..R*C, column-major
    const int k = (int)idx1 - 1;
    const int i = k % R;
    const int j = k / R;
    M_row_major[i * C + j] = v; // row-major storage
}

static inline void set_by_lowertri_index(
    tinytype* M_row_major, int R, int C,
    const uint16_t* lower_lin_idx,   // Lp_lower_lin_idx / Lq_lower_lin_idx
    uint16_t idx_lower_1based,       // 1..m (lower-tri packed)
    tinytype v)
{
    const uint16_t matlab_idx = lower_lin_idx[idx_lower_1based - 1];  // 1..R*C
    set_by_matlab_index(M_row_major, R, C, matlab_idx, v);
}

static void build_sparse_mat(tinytype* M, int R, int C,
                             const uint16_t* idx_const, const tinytype* val_const, size_t n_const,
                             const uint16_t* idx_var, const tinytype* coef_var, size_t n_var,
                             const tinytype* pow1, const tinytype* pow2)
{
    // M must be pre-zeroed or pre-initialized by caller.
    for (size_t n = 0; n < n_const; n++) {
        set_by_matlab_index(M, R, C, idx_const[n], val_const[n]);
    }
    for (size_t n = 0; n < n_var; n++) {
        const tinytype* coef = coef_var + n * POLY_TERM_COUNT;
        const tinytype v = poly2d_eval(coef, pow1, pow2);
        set_by_matlab_index(M, R, C, idx_var[n], v);
    }
}

static void build_sparse_mat_lowertri(
    tinytype* M, int R, int C,
    const uint16_t* lower_lin_idx,   // Lp_lower_lin_idx / Lq_lower_lin_idx
    const uint16_t* idx_const, const tinytype* val_const, size_t n_const,
    const uint16_t* idx_var, const tinytype* coef_var, size_t n_var,
    const tinytype* pow1, const tinytype* pow2)
{
    // M must be pre-zeroed by caller.
    for (size_t n = 0; n < n_const; n++) {
        set_by_lowertri_index(M, R, C, lower_lin_idx, idx_const[n], val_const[n]);
    }
    for (size_t n = 0; n < n_var; n++) {
        const tinytype* coef = coef_var + n * POLY_TERM_COUNT;
        const tinytype v = poly2d_eval(coef, pow1, pow2);
        set_by_lowertri_index(M, R, C, lower_lin_idx, idx_var[n], v);
    }
}

static void mat_mul_LT(tinytype* Out, const tinytype* L, int N) {
    // Out = L * L^T, where L is lower-triangular stored full (row-major).
    // This is O(N^3) but N is small (10 or 4).
    std::memset(Out, 0, sizeof(tinytype) * N * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            tinytype sum = 0;
            const int kmax = (j < i) ? j : i;
            for (int k = 0; k <= kmax; k++) {
                sum += L[i * N + k] * L[j * N + k];
            }
            Out[i * N + j] = sum;
            Out[j * N + i] = sum; // sym
        }
    }
}

// -----------------------------------------------------------------------------
// Build fitted matrices for given left/right leg lengths from leg_fit_coeffs_sparse.h
// NOTE: This assumes the header defines the following symbols (as exported by MATLAB):
//   - K_coef_40xT[40][POLY_TERM_COUNT]
//   - Ueq_coef_4xT[4][POLY_TERM_COUNT]
//   - Ad_idx_const[], Ad_val_const[], Ad_idx_var[], Ad_coef_var[][POLY_TERM_COUNT]
//   - Bd_idx_const[], Bd_val_const[], Bd_idx_var[], Bd_coef_var[][POLY_TERM_COUNT]
//   - G_idx_const[],  G_val_const[],  G_idx_var[],  G_coef_var[][POLY_TERM_COUNT]
//   - Am_idx_const[], Am_val_const[], Am_idx_var[], Am_coef_var[][POLY_TERM_COUNT]   (optional)
//   - Lp_idx_const/var + Lp_coef_var (Cholesky lower of P)
//   - Lq_idx_const/var + Lq_coef_var (Cholesky lower of Quu_inv)
// If your header uses different names, rename them here accordingly.
// -----------------------------------------------------------------------------
static void build_leg_fit_all(tinytype l_left, tinytype l_right)
{
    tinytype Kinf_tmp[NINPUTS * NSTATES];
    tinytype Ueq_tmp[NINPUTS];
    tinytype Ad_tmp[NSTATES * NSTATES];
    tinytype Bd_tmp[NSTATES * NINPUTS];
    tinytype Fdyn_tmp[NSTATES];
    tinytype AmBKt_tmp[NSTATES * NSTATES];
    tinytype Pinf_tmp[NSTATES * NSTATES];
    tinytype Quu_tmp[NINPUTS * NINPUTS];

    tinytype powL[POLY_MAX_DEG + 1];
    tinytype powR[POLY_MAX_DEG + 1];
    poly2d_prepare(l_left, l_right, powL, powR);

    // 1) Kinf (4x10) from 40 coefficients
    for (int idx = 0; idx < 40; idx++) {
        const tinytype v = poly2d_eval(&K_coef_40xT[idx][0], powL, powR);
        int row, col;
        if (idx < 10) { row = 0; col = idx; }
        else if (idx < 20) { row = 2; col = idx - 10; }
        else if (idx < 30) { row = 1; col = idx - 20; }
        else { row = 3; col = idx - 30; }
        Kinf_tmp[row * NSTATES + col] = v;
    }

    // 1.5) Ueq (4x1) from fitted static underactuated solve
    for (int i = 0; i < NINPUTS; i++) {
        Ueq_tmp[i] = poly2d_eval(&Ueq_coef_4xT[i][0], powL, powR);
    }

    // 2) Ad (10x10), Bd (10x4), fdyn/G (10x1)
    std::memset(Ad_tmp, 0, sizeof(Ad_tmp));
    std::memset(Bd_tmp, 0, sizeof(Bd_tmp));
    std::memset(Fdyn_tmp, 0, sizeof(Fdyn_tmp));

    build_sparse_mat(Ad_tmp, NSTATES, NSTATES,
                     Ad_idx_const, Ad_val_const, sizeof(Ad_idx_const)/sizeof(Ad_idx_const[0]),
                     Ad_idx_var,   &Ad_coef_var[0][0],  sizeof(Ad_idx_var)/sizeof(Ad_idx_var[0]),
                     powL, powR);

    build_sparse_mat(Bd_tmp, NSTATES, NINPUTS,
                     Bd_idx_const.data(), Bd_val_const.data(), Bd_idx_const.size(),
                     Bd_idx_var,   &Bd_coef_var[0][0],  sizeof(Bd_idx_var)/sizeof(Bd_idx_var[0]),
                     powL, powR);

    build_sparse_mat(Fdyn_tmp, NSTATES, 1,
                     G_idx_const.data(), G_val_const.data(), G_idx_const.size(),
                     G_idx_var,   &G_coef_var[0][0],   sizeof(G_idx_var)/sizeof(G_idx_var[0]),
                     powL, powR);

    // 3) AmBKt = (Ad - Bd*Kinf)^T
    for (int i = 0; i < NSTATES; i++) {
        for (int j = 0; j < NSTATES; j++) {
            tinytype v = Ad_tmp[i * NSTATES + j];
            for (int u = 0; u < NINPUTS; u++) {
                v -= Bd_tmp[i * NINPUTS + u] * Kinf_tmp[u * NSTATES + j];
            }
            AmBKt_tmp[j * NSTATES + i] = v;
        }
    }

    // 4) Pinf from Lp (Cholesky lower), Quu_inv from Lq (Cholesky lower)
    tinytype Lp[NSTATES * NSTATES];
    tinytype Lq[NINPUTS * NINPUTS];
    std::memset(Lp, 0, sizeof(Lp));
    std::memset(Lq, 0, sizeof(Lq));

    build_sparse_mat_lowertri(Lp, NSTATES, NSTATES,
                              Lp_lower_lin_idx,
                              Lp_idx_const.data(), Lp_val_const.data(), Lp_idx_const.size(),
                              Lp_idx_var,   &Lp_coef_var[0][0],  sizeof(Lp_idx_var)/sizeof(Lp_idx_var[0]),
                              powL, powR);

    build_sparse_mat_lowertri(Lq, NINPUTS, NINPUTS,
                              Lq_lower_lin_idx,
                              Lq_idx_const.data(), Lq_val_const.data(), Lq_idx_const.size(),
                              Lq_idx_var,   &Lq_coef_var[0][0],  sizeof(Lq_idx_var)/sizeof(Lq_idx_var[0]),
                              powL, powR);

    mat_mul_LT(Pinf_tmp, Lp, NSTATES);
    mat_mul_LT(Quu_tmp, Lq, NINPUTS);

    // Finite check (keep previous globals if invalid)
    if (!all_finite(Kinf_tmp, NINPUTS * NSTATES) ||
        !all_finite(Ueq_tmp, NINPUTS) ||
        !all_finite(Ad_tmp, NSTATES * NSTATES) ||
        !all_finite(Bd_tmp, NSTATES * NINPUTS) ||
        !all_finite(Fdyn_tmp, NSTATES) ||
        !all_finite(AmBKt_tmp, NSTATES * NSTATES) ||
        !all_finite(Pinf_tmp, NSTATES * NSTATES) ||
        !all_finite(Quu_tmp, NINPUTS * NINPUTS)) {
        return;
    }

    std::memcpy(Kinf_data, Kinf_tmp, sizeof(Kinf_tmp));
    std::memcpy(Ueq_ref_data, Ueq_tmp, sizeof(Ueq_tmp));
    std::memcpy(Adyn_data, Ad_tmp, sizeof(Ad_tmp));
    std::memcpy(Bdyn_data, Bd_tmp, sizeof(Bd_tmp));
    std::memcpy(fdyn_data, Fdyn_tmp, sizeof(Fdyn_tmp));
    std::memcpy(AmBKt_data, AmBKt_tmp, sizeof(AmBKt_tmp));
    std::memcpy(Pinf_data, Pinf_tmp, sizeof(Pinf_tmp));
    std::memcpy(Quu_inv_data, Quu_tmp, sizeof(Quu_tmp));
}
// ============================================================================
//  Create task function - called in board.cpp boardALLInit
// ============================================================================
TINYMPC_TASK_STORAGE static TinyMPC_Task tinympc_instance;

void tinympc_task_create(void)
{
   tinympc_ptr = &tinympc_instance;
}

// ============================================================================
// Constructor
// ============================================================================
TinyMPC_Task::TinyMPC_Task()
{
    for (int i = 0; i < NINPUTS; i++) {
        u_cache0[i] = 0.0f;
    }
    u_new = 0;
}

// ============================================================================
// Init - Setup static TinyMPC structures
// ============================================================================
void TinyMPC_Task::init()
{
    // ========== Initialize Cache with pre-computed Riccati parameters ==========
    cache.rho = 1.0f;

    // Direct copy instead of Map to avoid potential alignment issues
    for (int i = 0; i < NINPUTS; i++) {
        for (int j = 0; j < NSTATES; j++) {
            cache.Kinf(i, j) = Kinf_data[i * NSTATES + j];
        }
    }

    for (int i = 0; i < NSTATES; i++) {
        for (int j = 0; j < NSTATES; j++) {
            cache.Pinf(i, j) = Pinf_data[i * NSTATES + j];
            cache.AmBKt(i, j) = AmBKt_data[i * NSTATES + j];
        }
    }

    for (int i = 0; i < NINPUTS; i++) {
        for (int j = 0; j < NINPUTS; j++) {
            cache.Quu_inv(i, j) = Quu_inv_data[i * NINPUTS + j];
        }
    }

    for (int i = 0; i < NSTATES; i++) {
        cache.APf(i) = APf_data[i];
    }
    for (int i = 0; i < NINPUTS; i++) {
        cache.BPf(i) = BPf_data[i];
    }

    // ========== Initialize Settings ==========
    settings.abs_pri_tol = 1e-3f;
    settings.abs_dua_tol = 1e-3f;
    settings.max_iter = 100;
    settings.check_termination = 1;
    settings.en_state_bound = 1;
    settings.en_input_bound = 1;
    settings.en_state_soc = 0;  // No second order cone constraints
    settings.en_input_soc = 0;
    // ========== Initialize Bounds ==========
    // Set constraint bounds for each timestep
    for (int j = 0; j < NHORIZON; j++) {
        for (int i = 0; i < NSTATES; i++) {
            bounds.x_min(i, j) = -x_set_data[i];
            bounds.x_max(i, j) = x_set_data[i];
        }
    }
    for (int j = 0; j < NHORIZON - 1; j++) {
        for (int i = 0; i < NINPUTS; i++) {
            bounds.u_min(i, j) = u_min_data[i];
            bounds.u_max(i, j) = u_max_data[i];
        }
    }

    // Initialize slack and dual variables to zero
    bounds.z = tiny_MatrixNuNhm1::Zero();
    bounds.znew = tiny_MatrixNuNhm1::Zero();
    bounds.v = tiny_MatrixNxNh::Zero();
    bounds.vnew = tiny_MatrixNxNh::Zero();
    bounds.y = tiny_MatrixNuNhm1::Zero();
    bounds.g = tiny_MatrixNxNh::Zero();

    // ========== Initialize SOCs (not used but need to be set) ==========
    socs.cu[0] = 0.0f;
    socs.cx[0] = 0.0f;
    socs.Acu[0] = 0;
    socs.Acx[0] = 0;
    socs.qcu[0] = 0;
    socs.qcx[0] = 0;

    // ========== Initialize Workspace ==========
    // Direct copy dynamics matrices to avoid alignment issues
    for (int i = 0; i < NSTATES; i++) {
        for (int j = 0; j < NSTATES; j++) {
            work.Adyn(i, j) = Adyn_data[i * NSTATES + j];
        }
        for (int j = 0; j < NINPUTS; j++) {
            work.Bdyn(i, j) = Bdyn_data[i * NINPUTS + j];
        }
        work.fdyn(i) = fdyn_data[i];
        work.Q(i) = Q_data[i];
    }
    for (int i = 0; i < NINPUTS; i++) {
        work.R(i) = R_data[i];
    }

    // Initialize state and input trajectories to zero
    work.x = tiny_MatrixNxNh::Zero();
    work.u = tiny_MatrixNuNhm1::Zero();

    // Initialize linear cost terms to zero
    work.q = tiny_MatrixNxNh::Zero();
    work.r = tiny_MatrixNuNhm1::Zero();

    // Initialize Riccati backward pass terms
    work.p = tiny_MatrixNxNh::Zero();
    work.d = tiny_MatrixNuNhm1::Zero();

    // Initialize reference trajectories
    work.Xref = tiny_MatrixNxNh::Zero();
    work.Uref = tiny_MatrixNuNhm1::Zero();

    for (int i = 0; i < NSTATES; i++) 
       x_meas_data[i] = xinit_data[i];
    for (int i = 0; i < NSTATES; i++) {
       x_ref_data[i] = xgoal_data[i];
    }
		    for (int i = 0; i < NSTATES; i++) {
        work.x(i, 0) = x_meas_data[i];
    }
    // Set reference trajectory (interpolate from initial to goal)
    for (int j = 0; j < NHORIZON; j++) {
        for (int i = 0; i < NSTATES; i++) {
            work.Xref(i, j)=x_ref_data[i];
        }
    }
			
    // Set input reference (thrust to counter gravity)
    for (int j = 0; j < NHORIZON - 1; j++) {
        work.Uref(0, j) = 0.0f;
        work.Uref(1, j) = 0.0f;
        work.Uref(2, j) = 0.0f;  // Upward thrust
				work.Uref(3, j) = 0.0f;  // Upward thrust
    }

    // Initialize terminal cost
    work.p.col(NHORIZON - 1) = -cache.Pinf * work.Xref.col(NHORIZON - 1);

    // Link workspace to bounds and socs
    work.bounds = &bounds;
    work.socs = &socs;

    // Initialize status
    work.status = 0;
    work.iter = 0;
    work.primal_residual_state = 0.0f;
    work.primal_residual_input = 0.0f;
    work.dual_residual_state = 0.0f;
    work.dual_residual_input = 0.0f;

    // ========== Initialize Solver ==========
    solver.settings = &settings;
    solver.cache = &cache;
    solver.work = &work;
}


// ============================================================================
// Run - Execute MPC solve loop
// ============================================================================
TINYMPC_RUN_CRASH_WORKAROUND void TinyMPC_Task::run()
{
    const uint64_t t_start = tinympc_time_us();
	    tinytype l_left = leg_length_left;
	    tinytype l_right = leg_length_right;
	    build_leg_fit_all(l_left, l_right);
	 //==== 加载预计算矩阵到TinyMPC缓存 ====
		//Kinf矩阵
    for (int i = 0; i < NINPUTS; i++) {
        for (int j = 0; j < NSTATES; j++) {
            cache.Kinf(i, j) = Kinf_data[i * NSTATES + j];
        }
    }
    for (int i = 0; i < NSTATES; i++) {
        for (int j = 0; j < NSTATES; j++) {
            cache.Pinf(i, j) = Pinf_data[i * NSTATES + j];
            cache.AmBKt(i, j) = AmBKt_data[i * NSTATES + j];
        }
    }

		//Qinv
    for (int i = 0; i < NINPUTS; i++) {
        for (int j = 0; j < NINPUTS; j++) {
            cache.Quu_inv(i, j) = Quu_inv_data[i * NINPUTS + j];
        }
    }
		//APf向量
    for (int i = 0; i < NSTATES; i++) {
        cache.APf(i) = APf_data[i];
    }
		//BPf向量
    for (int i = 0; i < NINPUTS; i++) {
        cache.BPf(i) = BPf_data[i];
    }
    for (int i = 0; i < NSTATES; i++) {
        for (int j = 0; j < NSTATES; j++) {
            work.Adyn(i, j) = Adyn_data[i * NSTATES + j];
        }
        for (int j = 0; j < NINPUTS; j++) {
            work.Bdyn(i, j) = Bdyn_data[i * NINPUTS + j];
        }
        work.fdyn(i) = fdyn_data[i];
        work.Q(i) = Q_data[i];
    }
    for (int i = 0; i < NINPUTS; i++) {
        work.R(i) = R_data[i];
    }
		const tinytype *x_offest_set = x_bound_center_data;
	
		for (int j = 0; j < NHORIZON; j++) {
				for (int i = 0; i < NSTATES; i++) {
						bounds.x_min(i, j) = x_offest_set[i]-x_set_data[i];
						bounds.x_max(i, j) = x_offest_set[i]+x_set_data[i];
				}
		}
				for (int i = 0; i < NSTATES; i++) {
        work.x(i, 0) = x_meas_data[i];
    }
    // Update reference trajectory
		for (int j = 0; j < NHORIZON; j++) {
				for (int i = 0; i < NSTATES; i++) {
						work.Xref(i, j) = x_ref_data[i];
				}
		}
    // Update input reference from fitted static u_eq(l_left, l_right)
    for (int j = 0; j < NHORIZON - 1; j++) {
        for (int i = 0; i < NINPUTS; i++) {
            tinytype uref = Ueq_ref_data[i];
            if (uref < u_min_data[i]) { uref = u_min_data[i]; }
            if (uref > u_max_data[i]) { uref = u_max_data[i]; }
            work.Uref(i, j) = uref;
			//			 work.Uref(i, j) = 0;
	        }
	    }
    const uint64_t t_model_done = tinympc_time_us();
	    tiny_solve(&solver);
    const uint64_t t_solve_done = tinympc_time_us();
	    g_mpc_iter   = (float)work.iter;     // [RM AWARD 采集] 本拍 ADMM 迭代次数
	    g_mpc_solved = (float)work.status;   // [RM AWARD 采集] 1=收敛, 11=未收敛
	    const uint32_t previous_irq = tinympc_enter_critical();
    for (int i = 0; i < NINPUTS; i++) {
        u_cache0[i] = (float)work.u(i, 0);
    }
	    u_new = 0;
	    tinympc_exit_critical(previous_irq);

    const uint32_t model_us = (uint32_t)(t_model_done - t_start);
    const uint32_t solve_us = (uint32_t)(t_solve_done - t_model_done);
    const uint32_t total_us = (uint32_t)(tinympc_time_us() - t_start);
    g_mpc_model_us = model_us;
    g_mpc_solve_us = solve_us;
    g_mpc_total_us = total_us;
    if (total_us > g_mpc_max_total_us) {
        g_mpc_max_total_us = total_us;
    }
    if (total_us > 2000U) {
        g_mpc_over_2ms_count++;
    }
    if (total_us > 4000U) {
        g_mpc_over_4ms_count++;
    }
}


void TinyMPC_Task::tinympc_set_state(float *x,float *x_set)
{
   for (int i = 0; i < NSTATES; i++) {
				x_meas_data[i] = (tinytype)x[i];
				x_set_data[i] = (tinytype)x_set[i];
   }
}

void TinyMPC_Task::tinympc_set_target(float *xd)
{
   for (int i = 0; i < NSTATES; i++) {
       x_ref_data[i] = (tinytype)xd[i];
   }
}

void TinyMPC_Task::tinympc_set_model_input(float l_left, float l_right,
                                         const float *x_bound_center)
{
    leg_length_left = l_left;
    leg_length_right = l_right;
    for (int i = 0; i < NSTATES; ++i) {
        x_bound_center_data[i] = x_bound_center[i];
    }
}

// 方案3：用 LQR rollout 为 MPC warm start 播种 + 清对偶
// 目的：入 MPC 那一拍让 solver 从“LQR 假设下的可行轨迹”起解，通常 10 次迭代内收敛，
// 避免陈旧 work.x/u 导致第一拍大振荡。
void TinyMPC_Task::seed_mpc_from_lqr(const float *x_meas, const float (*K_lqr)[10])
{
    reset_dual(&solver);  // y, g 归零，避免上次留下的对偶污染

    // 从当前状态出发用 LQR rollout 走 NHORIZON 步
    tinytype x_k[NSTATES];
    for (int i = 0; i < NSTATES; i++) x_k[i] = (tinytype)x_meas[i];

    for (int k = 0; k < NHORIZON - 1; k++) {
        // 写入本步状态
        for (int i = 0; i < NSTATES; i++) work.x(i, k) = x_k[i];

        // u = -K_lqr · x（LQR 对 (X - Xd) 的标准形式，此处 Xd 近似 0 即 regulator）
        tinytype u_k[NINPUTS];
        for (int a = 0; a < NINPUTS; a++) {
            float s = 0;
            for (int b = 0; b < NSTATES; b++) s += K_lqr[a][b] * (float)x_k[b];
            u_k[a] = (tinytype)(-s);
        }
        // 投影到输入边界，确保 seed 本身可行
        if (u_k[0] < -7.0f)  u_k[0] = -7.0f;  else if (u_k[0] > 7.0f)  u_k[0] = 7.0f;
        if (u_k[1] < -7.0f)  u_k[1] = -7.0f;  else if (u_k[1] > 7.0f)  u_k[1] = 7.0f;
        if (u_k[2] < -45.0f) u_k[2] = -45.0f; else if (u_k[2] > 45.0f) u_k[2] = 45.0f;
        if (u_k[3] < -45.0f) u_k[3] = -45.0f; else if (u_k[3] > 45.0f) u_k[3] = 45.0f;
        for (int a = 0; a < NINPUTS; a++) work.u(a, k) = u_k[a];

        // 前向动力学：x_{k+1} = A x + B u + f
        tinytype x_next[NSTATES];
        for (int i = 0; i < NSTATES; i++) {
            tinytype s = fdyn_data[i];
            for (int j = 0; j < NSTATES; j++)  s += Adyn_data[i*NSTATES + j] * x_k[j];
            for (int j = 0; j < NINPUTS; j++)  s += Bdyn_data[i*NINPUTS + j] * u_k[j];
            x_next[i] = s;
        }
        for (int i = 0; i < NSTATES; i++) x_k[i] = x_next[i];
    }
    // 末端状态
    for (int i = 0; i < NSTATES; i++) work.x(i, NHORIZON - 1) = x_k[i];
}


void TinyMPC_Task::get_u(float *u_out)
{
    if(u_out == nullptr) 
		{
       return;
    }
    const uint32_t previous_irq = tinympc_enter_critical();
    if (u_new==1) 
		{
      for (int i = 0; i < NINPUTS; i++) 
			{
         u_out[i] = (float)work.u(i, 1);
      }
      u_new = 0;
    } 
//		else if(u_new==2) 
//		{
//       for (int i = 0; i < NINPUTS; i++) {
//         u_out[i] = (float)work.u(i, 2);
//       }
//			 u_new = 3;
//    }
//		else if(u_new==3) 
//		{
//			for (int i = 0; i < NINPUTS; i++) {
//         u_out[i] =(float)work.u(i, 3);
//       }
//			 u_new = 0;
//		}
    else
    {
        for (int i = 0; i < NINPUTS; i++) {
            u_out[i] = u_cache0[i];
        }
        u_new=1;
    }
    tinympc_exit_critical(previous_irq);
}



extern "C" void tinympc_run_once(void)
{
  tinympc_ptr->run();
}















