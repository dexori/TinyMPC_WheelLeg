#include "tinympc_task.h"
#include "tinympc_run.h"
#include "admm.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>

#define CHECK(x) do { if (!(x)) { \
    std::fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #x); \
    std::exit(EXIT_FAILURE); } } while (0)

// PC 替身，仅检查接入逻辑，不验证真实中断/计时或 MCU 最坏耗时。
static std::uint32_t irq_state = 1;
std::uint64_t tinympc_time_us()
{
    static std::uint64_t now;
    now += 100;
    return now;
}
std::uint32_t tinympc_enter_critical()
{
    const std::uint32_t previous = irq_state;
    irq_state = 1;
    return previous;
}
void tinympc_exit_critical(std::uint32_t previous)
{
    irq_state = previous;
}

int main()
{
    tinympc_task_create();
    CHECK(tinympc_ptr != nullptr);
    TinyMPC_Task &task = *tinympc_ptr;
    task.init();
    CHECK(task.cache.rho == 1.0f);
    float x[NSTATES] = {0,0,.10f,0,.08f,0,-.06f,0,.05f,0};
    float x_ref[NSTATES] = {0};
    float half_width[NSTATES] = {1000000,10,5,10,1.8f,4,1.8f,4,.8f,5};
    float center[NSTATES] = {0};
    task.tinympc_set_state(x, half_width);
    task.tinympc_set_target(x_ref);
    task.tinympc_set_model_input(.25f,.27f,center);
    tinympc_run_once();
    CHECK(task.work.status == 1 || task.work.status == 11);
    CHECK(task.work.iter > 0 && task.work.iter <= task.settings.max_iter);
    CHECK(task.work.Adyn.allFinite() && task.work.Bdyn.allFinite());
    CHECK(task.cache.Pinf.isApprox(task.cache.Pinf.transpose()));
    CHECK(task.cache.Quu_inv.isApprox(task.cache.Quu_inv.transpose()));
    CHECK(g_mpc_model_us == 100 && g_mpc_solve_us == 100);
    CHECK(g_mpc_total_us == 300);

    float u[NINPUTS];
    for (int i = 0; i < NINPUTS; ++i) {
        CHECK(std::isfinite(task.u_cache0[i]));
    }
    // 原 get_u 交替返回第 0 / 1 列，此行为保留而非作为正确调度结论。
    task.get_u(u);
    for (int i = 0; i < NINPUTS; ++i) CHECK(u[i] == task.u_cache0[i]);
    task.get_u(u);
    for (int i = 0; i < NINPUTS; ++i) CHECK(u[i] == task.work.u(i,1));
    CHECK(irq_state == 1);  // 进入前关中断，返回后仍然关中断。
    irq_state = 0;
    task.get_u(u);
    CHECK(irq_state == 0);

    // 在线一次项使用原始 Q/R；rho 只乘 slack/dual，不重复加到跟踪权重。
    task.bounds.vnew.setConstant(.2f);
    task.bounds.g.setConstant(.03f);
    task.bounds.znew.setConstant(.4f);
    task.bounds.y.setConstant(.05f);
    task.work.Xref.setConstant(.1f);
    task.work.Uref.setConstant(.3f);
    const float values[] = {.5f,1.0f,2.0f};
    for (float rho : values) {
        task.cache.rho = rho;
        update_linear_cost(&task.solver);
        CHECK(std::fabs(task.work.q(0,0) - (-Q_data[0]*.1f-rho*.17f)) < 2e-6f);
        CHECK(std::fabs(task.work.r(0,0) - (-R_data[0]*.3f-rho*.35f)) < 2e-6f);
    }
    task.cache.rho = 1.0f;  // 这里只测一次项，未用不同 rho 的旧缓存求解。
    std::puts("MCU task input/output and online rho-term tests passed.");
    return 0;
}
