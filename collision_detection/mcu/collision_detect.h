/* SPDX-License-Identifier: MIT */
#ifndef COLLISION_DETECT_H
#define COLLISION_DETECT_H

#include <stdint.h>
#include "wheelleg_mcg.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    MO_OK = 0, MO_DISABLED, MO_BAD_INPUT, MO_BAD_PARAMS, MO_NUMERICAL_ERROR
} MomentumObsStatus;

typedef struct {
    float K1[MO_DOF], K2[MO_DOF], K3[MO_DOF], K4[MO_DOF];
    float r_threshold[MO_DOF]; /* s: N，其余轴: N·m，按车重新整定 */
    float clamp_max;          /* 原实现同一数值限幅，并非统一物理量 */
    uint8_t legacy_startup_compensation; /* 1 复现 r-Mii；量纲待复核 */
} MomentumObsParams;

typedef struct {
    float qd[MO_DOF];         /* [ds,dyaw,dleft_leg,dright_leg,dpitch] */
    float tau[MO_DOF];        /* 广义力，先按 wheelleg_mcg.h 映射 */
    float dt;                /* s；原工程 0.002 */
    float slip_speed;        /* m/s：轮速与融合速度之差，用于 yaw 门控 */
    uint8_t enabled;         /* 0 复位，包括观测器和判定计数 */
} CollisionInput;

typedef struct {
    float r[MO_DOF], r_out[MO_DOF], int_r[MO_DOF];
    float int_dyn[MO_DOF], int_fe[MO_DOF];
    float debug_norm, debug_max_abs; /* 混合量纲，只作旧实现对拍 */
    uint8_t collision, axis_mask;
    uint8_t startup_cnt, initialized, prev_still;
} MomentumObsState;

typedef struct {
    float yaw_stall_threshold, yaw_max_slip;
    uint32_t yaw_stall_debounce, yaw_stall_cooldown;
    float yaw_deadzone, yaw_release;
    uint32_t yaw_debounce;
} CollisionGateParams;

typedef struct {
    uint32_t stall_count, cooldown_left, hysteresis_count;
    uint8_t yaw_latch;
    uint8_t yaw_stall_event; /* 单拍事件，不是原工程的应用层保持标志 */
    int8_t yaw_stall_dir;
    float yaw_gated_residual; /* 等同原 yawAdmTorque 信号；不执行参考控制 */
} CollisionGateState;

typedef struct {
    MomentumObsParams params;
    MomentumObsState observer;
    CollisionGateParams gate_params;
    CollisionGateState gate;
} CollisionDetector;

/* 一般接入只需要以下三个接口。 */
void collision_init(CollisionDetector *detector);
void collision_reset(CollisionDetector *detector); /* 保留参数 */
MomentumObsStatus collision_update(CollisionDetector *detector,
    const CollisionInput *input, const MomentumObsModel *model);

/* 只使用动量残差、不需要 yaw 判定时，调用以下接口。 */
void momentum_obs_params_default(MomentumObsParams *params);
void momentum_obs_reset(MomentumObsState *state);
/* 单独使用观测器时 slip_speed 不参与计算。模型与速度须来自同一拍。
 * 首次初始化后仍执行原算法的一次积分更新，不自动跳过首拍。 */
MomentumObsStatus momentum_obs_update(MomentumObsState *state,
    const MomentumObsParams *params, const MomentumObsModel *model,
    const CollisionInput *input);


#ifdef __cplusplus
}
#endif
#endif
