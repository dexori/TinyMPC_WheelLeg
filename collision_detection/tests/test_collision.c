/* SPDX-License-Identifier: MIT */
#include "collision_detect.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(x) do { if (!(x)) { \
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #x); \
    exit(EXIT_FAILURE); } } while (0)

static MomentumObsModel identity_model(void)
{
    MomentumObsModel m = {0};
    int i;
    for (i = 0; i < MO_DOF; ++i) m.M[i][i] = 1;
    return m;
}

static CollisionInput input(void)
{
    CollisionInput in = {0};
    in.enabled = 1;
    in.dt = .002f;
    return in;
}

static void test_zero_and_instances(void)
{
    CollisionDetector a, b;
    CollisionInput in = input();
    MomentumObsModel m = identity_model();
    int n;
    collision_init(&a); collision_init(&b);
    for (n = 0; n < 10000; ++n) {
        CHECK(collision_update(&a, &in, &m) == MO_OK);
        CHECK(a.observer.debug_norm == 0 && !a.observer.collision);
    }
    m.G[1] = 5000; /* 仅构造门控测试信号，不是外力检测精度测试。 */
    CHECK(collision_update(&a, &in, &m) == MO_OK);
    m = identity_model();
    CHECK(collision_update(&b, &in, &m) == MO_OK);
    CHECK(b.observer.r_out[1] == 0 && a.observer.r_out[1] > 15);
}

static void test_first_step_and_startup(void)
{
    CollisionDetector d;
    CollisionInput in = input();
    MomentumObsModel m = identity_model();
    float expected;
    collision_init(&d);
    in.tau[1] = 10;
    CHECK(collision_update(&d, &in, &m) == MO_OK);
    expected = 12 * (-10 * in.dt * in.dt + .3f * (-10 * in.dt));
    CHECK(fabsf(d.observer.r[1] - expected) < 1e-6f);
    in = input(); in.qd[0] = .1f;
    collision_reset(&d);
    CHECK(collision_update(&d, &in, &m) == MO_OK);
    CHECK(d.observer.startup_cnt == 14);
    CHECK(d.observer.r_out[0] == -1); /* 旧 r-Mii 行为，非物理补偿验证 */
    CHECK(d.observer.debug_norm > 0);
    collision_reset(&d);
    d.params.legacy_startup_compensation = 0;
    CHECK(collision_update(&d, &in, &m) == MO_OK);
    CHECK(d.observer.r_out[0] == 0);
}

static void test_gates(void)
{
    CollisionDetector d;
    CollisionInput in = input();
    MomentumObsModel m = identity_model();
    unsigned n;
    m.G[1] = 5000;
    collision_init(&d);
    for (n = 1; n <= 45; ++n) {
        CHECK(collision_update(&d, &in, &m) == MO_OK);
        CHECK(d.observer.axis_mask == (1u << 1));
        CHECK(d.gate.yaw_stall_event == (n == 5));
        if (n == 5) CHECK(d.gate.yaw_stall_dir == 1);
        CHECK(d.gate.yaw_latch == (n >= 40));
    }
    CHECK(d.gate.cooldown_left == 1960);
    CHECK(d.gate.yaw_gated_residual > 0);
    in.slip_speed = .5f; /* 严格小于门槛才通过 */
    CHECK(collision_update(&d, &in, &m) == MO_OK);
    CHECK(!d.gate.yaw_latch && d.gate.yaw_gated_residual == 0);
    in.enabled = 0;
    CHECK(collision_update(&d, &in, NULL) == MO_DISABLED);
    CHECK(!d.observer.initialized && !d.gate.cooldown_left);
    in.enabled = 1; in.slip_speed = 0; m.G[1] = -5000;
    for (n = 0; n < 5; ++n) CHECK(collision_update(&d, &in, &m) == MO_OK);
    CHECK(d.gate.yaw_stall_event && d.gate.yaw_stall_dir == -1);
}

static void test_invalid(void)
{
    CollisionDetector d;
    CollisionInput in = input();
    MomentumObsModel m = identity_model();
    collision_init(&d);
    CHECK(collision_update(NULL, &in, &m) == MO_BAD_INPUT);
    CHECK(collision_update(&d, NULL, &m) == MO_BAD_INPUT);
    CHECK(collision_update(&d, &in, NULL) == MO_BAD_INPUT);
    in.dt = NAN; CHECK(collision_update(&d, &in, &m) == MO_BAD_INPUT);
    in = input(); in.tau[2] = INFINITY;
    CHECK(collision_update(&d, &in, &m) == MO_BAD_INPUT);
    in = input(); in.slip_speed = NAN;
    CHECK(collision_update(&d, &in, &m) == MO_BAD_INPUT);
    in = input(); d.params.K1[0] = 0;
    CHECK(collision_update(&d, &in, &m) == MO_BAD_PARAMS);
    collision_init(&d); d.gate_params.yaw_debounce = 0;
    CHECK(collision_update(&d, &in, &m) == MO_BAD_PARAMS);
    collision_init(&d); m.M[0][0] = FLT_MAX; in.qd[0] = FLT_MAX;
    CHECK(collision_update(&d, &in, &m) == MO_NUMERICAL_ERROR);
    CHECK(!d.observer.initialized && d.observer.r_out[0] == 0);
}

static void test_legacy_model(void)
{
    float q[5] = {0,0,.2f,-.1f,.05f}, qd[5] = {0,0,.6f,-.4f,.2f};
    float plus[5], minus[5], u[4] = {1,2,3,4}, tau[5];
    MomentumObsModel m, mp, mm;
    float defect;
    int i,j;
    CHECK(wheelleg_mcg_legacy(q,qd,&m));
    CHECK(fabsf(m.M[0][0]-24.6191124260355f) < 1e-6f);
    for (i = 0; i < 5; ++i) {
        for (j = 0; j < 5; ++j) CHECK(m.M[i][j] == m.M[j][i]);
        plus[i] = q[i]+.001f*qd[i]; minus[i] = q[i]-.001f*qd[i];
    }
    CHECK(wheelleg_mcg_legacy(plus,qd,&mp));
    CHECK(wheelleg_mcg_legacy(minus,qd,&mm));
    defect = (mp.M[0][2]-mm.M[0][2])/.002f-m.C[0][2]-m.C[2][0];
    CHECK(fabsf(defect) > .1f); /* 已知旧 C 不一致，不能误标为“通过” */
    printf("legacy Mdot-C-C^T [0,2] defect: %g (known issue)\n",(double)defect);
    CHECK(wheelleg_generalized_force(u,.065f,.245f,tau));
    CHECK(fabsf(tau[0]-3/.065f) < 1e-5f);
    CHECK(tau[2] == 2 && tau[3] == 2 && tau[4] == -7);
    CHECK(!wheelleg_generalized_force(u,0,.245f,tau));
    q[0] = NAN; CHECK(!wheelleg_mcg_legacy(q,qd,&m));
    CHECK(m.M[0][0] == 0);
}

int main(void)
{
    test_zero_and_instances(); test_first_step_and_startup();
    test_gates(); test_invalid(); test_legacy_model();
    puts("Collision observer tests passed.");
    return 0;
}
