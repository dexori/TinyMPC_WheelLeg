/* SPDX-License-Identifier: MIT
 * Extracted from perception.cpp::MomentumObserver().
 * 只包含数值积分、残差和分轴阈值，不访问电机、HAL 或控制器。 */
#include "collision_detect.h"
#include <math.h>
#include <stddef.h>
#include <string.h>

static float clampf(float x, float lo, float hi)
{
    return x < lo ? lo : (x > hi ? hi : x);
}

/* 默认增益与阈值，换车时从这里检查。 */
void momentum_obs_params_default(MomentumObsParams *p)
{
    const MomentumObsParams defaults = {
        {8, 12, 12, 12, 6}, {2, .5f, 2.4f, 2.4f, 1.6f},
        {.15f, .3f, .2f, .2f, .1f}, {1, 1, 1, 1, 1},
        {1, 4, 1, 1, 2}, 300, 1
    };
    if (p != NULL) *p = defaults;
}

void momentum_obs_reset(MomentumObsState *s)
{
    if (s == NULL) return;
    memset(s, 0, sizeof(*s));
    s->prev_still = 1;
}

static int valid_params(const MomentumObsParams *p)
{
    int i;
    if (!isfinite(p->clamp_max) || p->clamp_max <= 0 ||
        p->legacy_startup_compensation > 1) return 0;
    for (i = 0; i < MO_DOF; ++i) {
        if (!isfinite(p->K1[i]) || p->K1[i] <= 0 ||
            !isfinite(p->K2[i]) || p->K2[i] < 0 ||
            !isfinite(p->K3[i]) || p->K3[i] < 0 ||
            !isfinite(p->K4[i]) || p->K4[i] <= 0 ||
            !isfinite(p->r_threshold[i]) || p->r_threshold[i] <= 0) return 0;
    }
    return 1;
}

MomentumObsStatus momentum_obs_update(MomentumObsState *s,
    const MomentumObsParams *cfg, const MomentumObsModel *model,
    const CollisionInput *in)
{
    float p[MO_DOF] = {0}, dyn[MO_DOF] = {0};
    int i, j, still;
    if (s == NULL) return MO_BAD_INPUT;
    if (in == NULL) { momentum_obs_reset(s); return MO_BAD_INPUT; }
    if (!in->enabled) { momentum_obs_reset(s); return MO_DISABLED; }
    if (cfg == NULL || model == NULL || !isfinite(in->dt) ||
        in->dt < .0005f || in->dt > .02f) {
        momentum_obs_reset(s); return MO_BAD_INPUT;
    }
    if (!valid_params(cfg)) { momentum_obs_reset(s); return MO_BAD_PARAMS; }
    for (i = 0; i < MO_DOF; ++i) {
        if (!isfinite(in->qd[i]) || !isfinite(in->tau[i]) ||
            !isfinite(model->G[i])) {
            momentum_obs_reset(s); return MO_BAD_INPUT;
        }
        for (j = 0; j < MO_DOF; ++j) {
            if (!isfinite(model->M[i][j]) || !isfinite(model->C[i][j])) {
                momentum_obs_reset(s); return MO_BAD_INPUT;
            }
        }
    }
    /* p=M*qd，dyn=tau+C^T*qd-G。模型一致性由调用方负责。 */
    for (i = 0; i < MO_DOF; ++i) {
        dyn[i] = in->tau[i] - model->G[i];
        for (j = 0; j < MO_DOF; ++j) {
            p[i] += model->M[i][j] * in->qd[j];
            dyn[i] += model->C[j][i] * in->qd[j];
        }
        if (!isfinite(p[i]) || !isfinite(dyn[i])) {
            momentum_obs_reset(s); return MO_NUMERICAL_ERROR;
        }
    }
    /* 首次调用以当前动量初始化积分状态。 */
    if (!s->initialized) {
        for (i = 0; i < MO_DOF; ++i) {
            s->int_dyn[i] = p[i];
            s->int_fe[i] = p[i];
            s->int_r[i] = s->r[i] = s->r_out[i] = 0;
        }
        s->initialized = 1;
    }
    /* 原动量观测器积分公式，计算本拍残差。 */
    for (i = 0; i < MO_DOF; ++i) {
        float e, fe, candidate;
        s->int_dyn[i] += (dyn[i] + s->r[i]) * in->dt;
        s->int_fe[i] += dyn[i] * in->dt;
        e = p[i] - s->int_dyn[i];
        fe = p[i] - s->int_fe[i];
        s->int_r[i] += (-cfg->K2[i] * s->r[i] + e) * in->dt;
        candidate = cfg->K1[i] * (s->int_r[i] + cfg->K3[i] * fe);
        if (!isfinite(candidate) || !isfinite(s->int_dyn[i]) ||
            !isfinite(s->int_fe[i]) || !isfinite(s->int_r[i])) {
            momentum_obs_reset(s); return MO_NUMERICAL_ERROR;
        }
        s->r[i] = clampf(candidate, -cfg->clamp_max, cfg->clamp_max);
    }
    /* 原启动经验项及分轴阈值；不改变原默认行为。 */
    still = fabsf(in->qd[0]) < .03f && fabsf(in->qd[1]) < .03f &&
            fabsf(in->qd[2]) < .06f && fabsf(in->qd[3]) < .06f &&
            fabsf(in->qd[4]) < .06f;
    if (s->prev_still && !still) s->startup_cnt = 15;
    s->prev_still = (uint8_t)still;
    {
        const float k5 = s->startup_cnt > 0 &&
                         cfg->legacy_startup_compensation ? 1.0f : 0.0f;
        if (s->startup_cnt > 0) --s->startup_cnt;
        s->collision = s->axis_mask = 0;
        s->debug_norm = s->debug_max_abs = 0;
        for (i = 0; i < MO_DOF; ++i) {
            float a;
            s->r_out[i] = cfg->K4[i] * (s->r[i] - k5 * model->M[i][i]);
            a = fabsf(s->r_out[i]);
            s->debug_norm += a * a;
            if (a > s->debug_max_abs) s->debug_max_abs = a;
            if (a > cfg->r_threshold[i]) {
                s->axis_mask |= (uint8_t)(1u << i);
                s->collision = 1;
            }
            if (!isfinite(s->r_out[i]) || !isfinite(s->debug_norm)) {
                momentum_obs_reset(s); return MO_NUMERICAL_ERROR;
            }
        }
        s->debug_norm = sqrtf(s->debug_norm);
    }
    return MO_OK;
}

/* yaw 判定：连续计数、冷却与滞回。 */
void collision_reset(CollisionDetector *d)
{
    if (d == NULL) return;
    momentum_obs_reset(&d->observer);
    memset(&d->gate, 0, sizeof(d->gate));
}

void collision_init(CollisionDetector *d)
{
    const CollisionGateParams defaults = {15, .5f, 5, 2000, 8, 5, 40};
    if (d == NULL) return;
    momentum_obs_params_default(&d->params);
    d->gate_params = defaults;
    collision_reset(d);
}

static int valid_gate(const CollisionGateParams *p)
{
    return isfinite(p->yaw_stall_threshold) && p->yaw_stall_threshold > 0 &&
        isfinite(p->yaw_max_slip) && p->yaw_max_slip > 0 &&
        isfinite(p->yaw_deadzone) && isfinite(p->yaw_release) &&
        p->yaw_deadzone > p->yaw_release && p->yaw_release >= 0 &&
        p->yaw_stall_debounce > 0 && p->yaw_debounce > 0;
}

MomentumObsStatus collision_update(CollisionDetector *d,
    const CollisionInput *in, const MomentumObsModel *model)
{
    MomentumObsStatus status;
    float signal, mag;
    int slip_ok;
    CollisionGateState *g;
    const CollisionGateParams *p;
    if (d == NULL) return MO_BAD_INPUT;
    if (in == NULL) { collision_reset(d); return MO_BAD_INPUT; }
    if (!in->enabled) { collision_reset(d); return MO_DISABLED; }
    if (!isfinite(in->slip_speed)) { collision_reset(d); return MO_BAD_INPUT; }
    if (!valid_gate(&d->gate_params)) {
        collision_reset(d); return MO_BAD_PARAMS;
    }
    status = momentum_obs_update(&d->observer, &d->params, model, in);
    if (status != MO_OK) { collision_reset(d); return status; }
    g = &d->gate; p = &d->gate_params;
    signal = d->observer.r_out[1];
    mag = fabsf(signal);
    slip_ok = fabsf(in->slip_speed) < p->yaw_max_slip;
    /* 原应用标志由上层消费；这里将其改为明确的单拍事件。 */
    g->yaw_stall_event = 0;
    g->yaw_stall_dir = 0;
    if (g->cooldown_left > 0) {
        --g->cooldown_left;
        g->stall_count = 0;
    } else if (mag > p->yaw_stall_threshold && slip_ok) {
        if (g->stall_count < p->yaw_stall_debounce) ++g->stall_count;
        if (g->stall_count >= p->yaw_stall_debounce) {
            g->yaw_stall_event = 1;
            g->yaw_stall_dir = signal > 0 ? 1 : -1;
            g->cooldown_left = p->yaw_stall_cooldown;
            g->stall_count = 0;
        }
    } else {
        g->stall_count = 0;
    }
    /* 原 AdmHysteresis；不添加“偏航速度低”的额外条件。 */
    if (!g->yaw_latch) {
        if (mag > p->yaw_deadzone && slip_ok) {
            if (g->hysteresis_count < p->yaw_debounce) ++g->hysteresis_count;
        } else {
            g->hysteresis_count = 0;
        }
        if (g->hysteresis_count >= p->yaw_debounce) g->yaw_latch = 1;
    } else if (mag < p->yaw_release || !slip_ok) {
        g->yaw_latch = 0;
        g->hysteresis_count = 0;
    }
    g->yaw_gated_residual = g->yaw_latch ?
        signal - (signal > 0 ? p->yaw_release : -p->yaw_release) : 0;
    return MO_OK;
}
