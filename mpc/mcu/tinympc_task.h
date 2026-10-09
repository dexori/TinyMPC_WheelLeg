/******************************
File name: tinympc_task.h
Description: 原 H7 TinyMPC 任务提取版，整车输入由调用方显式传入。
*****************************/

#ifndef __TINYMPC_TEST_TASK_H__
#define __TINYMPC_TEST_TASK_H__

#include <cstdint>
#include "tinympc_port.h"

// Use static TinyMPC types from tinympc_stm32
#include "glob_opts.hpp"
#include "types.hpp"
#include "leg_fit_coeffs_sparse.h"
class TinyMPC_Task 
{
public:
    TinyMPC_Task();
    // Static TinyMPC structures (no dynamic allocation)
    TinyCache cache;
    TinySettings settings;
    TinyWorkspace work;
    TinyBounds bounds;
    TinySocs socs;
    TinySolver solver;
    // Cached MPC outputs
		tinytype x_set_data[NSTATES]={10000000000.0f, 10.0f, 5.0f, 10.0f, 1.8f, 4.0f,1.8f, 4.0f, 0.8f, 5.0f};
		tinytype x_meas_data[NSTATES] = {0};
		tinytype x_ref_data[NSTATES] = {0};
    float u_cache0[NINPUTS];
    volatile uint8_t u_new;
    // 原 perception 腿长与 controller 误差坐标边界中心。
    // 所有 setter 和 run 必须由同一执行上下文调用，或由上层同步。
    tinytype leg_length_left = 0.25f;
    tinytype leg_length_right = 0.25f;
    tinytype x_bound_center_data[NSTATES] = {0};

    void init();
    void run();
    void tinympc_set_state(float *x,float *x_set);
    void tinympc_set_target(float *xd);
    void get_u(float *u_out);
    void tinympc_set_model_input(float l_left, float l_right,
                                const float *x_bound_center);

    // 方案3：入 MPC 前用 LQR rollout 为 warm start 播种 + 清对偶
    // K_lqr 是 4x10 行主序增益（与 Controller::K 一致）；x_meas 是当前 10 维状态
    void seed_mpc_from_lqr(const float *x_meas, const float (*K_lqr)[10]);
};

// Global task pointer
extern TinyMPC_Task* tinympc_ptr;
extern tinytype Adyn_data[NSTATES * NSTATES];
extern tinytype Bdyn_data[NSTATES * NINPUTS];
extern tinytype Ueq_ref_data[NINPUTS];
extern volatile float g_mpc_iter;
extern volatile float g_mpc_solved;
extern volatile uint32_t g_mpc_total_us;
extern volatile uint32_t g_mpc_model_us;
extern volatile uint32_t g_mpc_solve_us;
extern volatile uint32_t g_mpc_max_total_us;
extern volatile uint32_t g_mpc_over_2ms_count;
extern volatile uint32_t g_mpc_over_4ms_count;

// Create task function
void tinympc_task_create(void);



#endif
