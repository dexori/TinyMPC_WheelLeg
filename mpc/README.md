# MPC 代码说明

本目录只保留原 H7 的 C++ MPC 实现和 MATLAB 参数生成。先看 [mcu/tinympc_task.cpp](mcu/tinympc_task.cpp)，调用接口在 [mcu/tinympc_task.h](mcu/tinympc_task.h)。

## 原 MCU 代码

`mcu/` 里的关键文件：

- `tinympc_task.cpp/.h`：静态工作区初始化、腿长拟合求值、缓存更新、求解、LQR 初值播种和输出缓存。
- `admm.cpp/.hpp`、`types.hpp`、`glob_opts.hpp`：原 C++ ADMM 及固定维度类型。
- `leg_fit_coeffs_sparse.h`：原 H7 使用的四位小数稀疏拟合参数，系数没有重算。
- `tinympc_port.h`：需要底层实现的微秒计时和临界区接口。

移植时把两个 `.cpp` 加入工程，本目录和 Eigen 加入头文件路径。这里没有带完整 Eigen、HAL、controller/perception 或板级工程。PC 检查需要 C++11 和 Eigen 3.4+，在仓库根目录执行：

```bash
cmake -S . -B build-mcu -DCMAKE_BUILD_TYPE=Release
cmake --build build-mcu --parallel
ctest --test-dir build-mcu --output-on-failure
```

状态和输入顺序：

```text
x = [s, ds, yaw, dyaw, left_leg, dleft_leg,
     right_leg, dright_leg, pitch, dpitch]
u = [左轮力矩, 右轮力矩, 左机体—腿力矩, 右机体—腿力矩]
```

距离用 m，角度用 rad，速度用 m/s 或 rad/s，力矩用 N·m。模型周期 4 ms，8 个状态节点、7 段输入；机体—腿力矩到关节电机的映射由调用方完成。

任务调用顺序（变量由调用方提供）：

```cpp
tinympc_task_create();
tinympc_ptr->init();  // 启动时只执行一次

// 每次求解前，传入同一拍的数据。
tinympc_ptr->tinympc_set_state(x_mpc, x_half_width);
tinympc_ptr->tinympc_set_target(x_ref);
tinympc_ptr->tinympc_set_model_input(l_left, l_right, x_bound_center);
tinympc_run_once();
tinympc_ptr->get_u(u_out);
```

原 `controller` 的误差状态需要上层自行构造。`x_bound_center` 对应原任务中负的参考/前馈偏移，不是目标状态的另一个名字；`x_half_width` 是各轴边界半宽。输入更新和 `run()` 需要同一执行上下文或上层同步，不能在中断里边写状态边求解。

底层实现 `tinympc_time_us()`、`tinympc_enter_critical()`、`tinympc_exit_critical()`。STM32 临界区需保存原 PRIMASK，退出时恢复原值。若需要原 DTCM 放置，定义 `TINYMPC_TASK_STORAGE` 并配置链接段；默认不指定段。

保留原逻辑也意味着保留其边界：`get_u()` 交替取预测第 0/1 列，种子轮力矩限幅为 ±7 而求解约束为 ±6；未收敛时仍发布结果，错误模型会沿用旧参数，输出没有额外物理限幅。还有误差坐标仿射项、APf/BPf 与非零 f 的一致性需要核对。这里没有把这些行为宣称为安全或正确的调度方案，上层必须检查状态、同步、时效并做好降级。

提取版在 PC 通过 Debug/Release 接口测试，并与原任务在相同 Eigen 下做了 20 拍不同腿长和参考输入的对照，A/K、首拍输出、状态和迭代数一致。计时与中断在测试中是替身，不代表 MCU 时序、并发或实车效果已经验证。

## rho 到底加在哪里

整理后的 [MATLAB 生成代码](matlab/generate_mpc_parameters.m) 使用 `Q_rho = Q + rho*I`、`R_rho = R + rho*I` 算 K/P/Quu 缓存。Q/R 本身仍是原始跟踪权重，在线一次项使用它们，再单独加 slack/dual 的 rho 项，不重复加到 Q/R。

原工程的 `LQR.m` 则先用原 Q/R 算 K/P，只在 Quu 加 rho。补回的 MCU 任务保留原参数表与 `cache.rho=1`，没有宣称原表已按新方式重新生成。改 rho 时必须重算并验证整套缓存，不能只改 `cache.rho`。新的 MATLAB 输出和原稀疏表接口不同，目前不能直接替换。

## MATLAB 从哪里开始

入口是 [matlab/generate_mpc_fit.m](matlab/generate_mpc_fit.m)，日常参数在 [matlab/wheelleg_parameters.m](matlab/wheelleg_parameters.m)。需要 MATLAB Symbolic Math Toolbox。

在仓库根目录打开 MATLAB：

```matlab
cd mpc/matlab
p = wheelleg_parameters();
p.angle_range_deg = 0:5:50;  % 先用少量点检查流程
result = generate_mpc_fit(p);
```

| 文件 | 作用 |
|---|---|
| `wheelleg_parameters.m` | 机构参数、Q/R、约束、采样与输出配置 |
| `prepare_leg_data.m` | 读取 CAD 表，计算腿长、质心与惯量 |
| `derive_reduced_dynamics.m` | 消去内力，推导动力学 |
| `derive_ab_model.m`、`derive_gravity_model.m` | 连续 A/B 和重力项 |
| `sample_ab_models.m` | 采样左右腿工作点 |
| `generate_mpc_parameters.m` | 离散化、计算缓存、拟合与导出 |

CAD 示例在 `matlab/data/leg_mass_properties.csv`，四列为 `dist_mm,com_y_mm,com_x_mm,Ixx_kgmm2`：腿长、CAD 质心 Y/X、惯量，单位依次为 mm、mm、mm、kg·mm²。腿长必须递增并覆盖采样范围；换机构时同时修改 CAD 表、尺寸、质量和惯量。

结果写入 `matlab/generated/wheelleg_mpc_fit_data.h` 和 `wheelleg_mpc_fit.mat`。这是新的拟合格式，不是原任务使用的 `mcu/leg_fit_coeffs_sparse.h`。不能改名覆盖；接入时需要适配归一化、系数布局和缓存更新，保证一次求解使用同一工作点的一整套数据。

PC 测试只验证数值复现及接口。原任务沿用稳态 Riccati 缓存，终端处理存在近似，不应称为任意终端权重下的精确有限时域 QP 实现。完整 MATLAB 符号流程和 H7 实车性能尚未重新验证。
