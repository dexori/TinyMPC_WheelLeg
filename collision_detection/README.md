# 碰撞检测

代码从原 H7 工程提取，只计算动量残差和碰撞判据，不执行电机、导纳或脱困动作。

## 只看这两份代码

- [mcu/collision_detect.c](mcu/collision_detect.c)：默认参数 → 动量积分 → 启动修正 → 分轴阈值 → yaw 判定。调用接口在同目录的 `collision_detect.h`。
- [mcu/wheelleg_mcg.c](mcu/wheelleg_mcg.c)：原固定 M/C/G 数值模型和四路力矩映射。

`matlab/` 放 M/C/G 推导和机构参数；`tests/` 放 PC 测试。移植时只加入 `mcu/` 的两个 C 文件，并把该目录加入头文件路径，不需要 HAL、Eigen 或 MPC。

## 怎么接

先包含 `collision_detect.h` 和 `string.h`，静态保存一个 `CollisionDetector`。启动时调用一次 `collision_init(&detector)`，不要每拍重新初始化。

状态和力矩顺序：

```text
q  = [s, yaw, left_leg, right_leg, pitch]
qd = [ds, dyaw, dleft_leg, dright_leg, dpitch]
u  = [左轮力矩, 右轮力矩, 左机体—腿力矩, 右机体—腿力矩]
```

距离用 m，角度用 rad，速度用 m/s 或 rad/s，力矩用 N·m。q、qd 和 u 必须来自同一拍，电流换算和正负号由上层处理。

每拍函数内按这个顺序调用，q/qd/u/slip_speed/dt 由上层提供：

```c
MomentumObsModel model;
CollisionInput input = {0};
if (!wheelleg_mcg_legacy(q, qd, &model) ||
    !wheelleg_generalized_force(u, wheel_radius, half_track, input.tau)) {
    collision_reset(&detector);
    return;  // 模型或输入无效，本拍停止
}
memcpy(input.qd, qd, sizeof(input.qd));
input.dt = dt;
input.slip_speed = slip_speed;
input.enabled = 1;
MomentumObsStatus status = collision_update(&detector, &input, &model);
```

轮半径 `wheel_radius` 和半轮距 `half_track` 单位为 m。换用自己的模型时，只替换 M/C/G 求值，填入 `MomentumObsModel`；只需要残差可调用 `momentum_obs_update()`。

计算有效时 `status == MO_OK`。残差读 `detector.observer.r_out`，分轴超阈值读 `detector.observer.axis_mask`，yaw 单拍事件读 `detector.gate.yaw_stall_event`。纵向残差单位为 N，其余轴为 N·m；这几个量不是已经确认的物理碰撞或外力测量值。

增益和阈值在 `momentum_obs_params_default()`，yaw 计数在 `collision_init()`。默认 500 Hz；计数按调用次数推进，改频率要一起调整。停用时复位；`collision_reset()` 保留参数，`collision_init()` 恢复默认参数。

## MATLAB 入口

需要 MATLAB Symbolic Math Toolbox。在仓库根目录打开 MATLAB：

```matlab
cd collision_detection/matlab
model = MCG(mcg_parameters());
test_mcg();
```

先改 `mcg_parameters.m`；它通过 `balance_canshujisuan.m → balance_Leg.m` 计算机构参数。输出在 `matlab/generated/`，接口为 `M_fun(q)`、`C_fun(q,qd)`、`G_fun(q)`、`Q_fun(u)`，输入为列向量。不会自动覆盖 C 模型，也不含实时腿长变化项。

## 保留的原实现问题

- 原 C 矩阵不满足 `Mdot=C+C'`；MATLAB 已改用 Christoffel 构造，但没有同步替换 C 模型。
- 起动前 15 拍的 `r_out=K4*(r-Mii)` 是量纲待复核的经验项，可用 `legacy_startup_compensation=0` 关闭，之后需要重新整定。
- 分轴阈值没有通用去抖，yaw 判定也没有检查速度是否接近零。模型、增益和阈值需要按车验证，`MO_OK` 不表示可以执行安全控制动作。

编译方法见[根 README](../README.md)。C 测试和原观测器核心 2000 拍对拍通过，只说明实现可复现；MATLAB 完整符号流程与实车效果尚未重新验证。来源见 [NOTICE](NOTICE)，许可证见 [LICENSE](LICENSE)。
