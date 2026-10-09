# TinyMPC-WheelLeg

从自己的轮腿 STM32H7 工程中整理了两部分代码：MPC 任务与参数生成、碰撞动量观测器与 M/C/G 推导。两块可以独立使用，不包含整车固件、电机驱动和通信。

## 从哪里看

```text
mpc/                         MPC 模块
├── mcu/                     原 H7 任务、C++ ADMM 和腿长参数表
├── matlab/                  建模、CAD 参数和变腿长拟合
└── tests/                   任务接口测试，PC 使用底层替身

collision_detection/         碰撞检测模块
├── mcu/                     观测器与判定、原 M/C/G 数值模型
├── matlab/                  M/C/G 符号推导及检查
└── tests/                   C 接口测试
```

看 MPC，先打开 [tinympc_task.cpp](mpc/mcu/tinympc_task.cpp)，调用方法见 [MPC README](mpc/README.md)。

看碰撞检测，先打开 [collision_detect.c](collision_detection/mcu/collision_detect.c)，调用接口在 [collision_detect.h](collision_detection/mcu/collision_detect.h)，接入方法见 [碰撞检测 README](collision_detection/README.md)。

原理和实验说明放在论坛，这里只说明代码怎么用。

## 编译测试

需要 C99/C++11 编译器、CMake 3.16+ 和 Eigen 3.4+。Eigen 是 MPC 的外部依赖，没有整份复制进仓库。在仓库根目录执行：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

默认编译 MPC 和碰撞检测。只编译碰撞检测时加 `-DBUILD_MPC=OFF`，不需要 Eigen；只编译 MPC 时加 `-DBUILD_COLLISION=OFF`。也可以用 `-S mpc` 或 `-S collision_detection` 单独构建，使用不同的 build 目录。

上 MCU 时，MPC 需要 `mpc/mcu/` 的两个 `.cpp`、配套头文件、Eigen 和底层适配；碰撞检测只需 `collision_detection/mcu/` 的两个 `.c` 和两个 `.h`。PC 测试和 CMake 不需要带入固件。

## 使用前注意

- MPC 在线使用原 `leg_fit_coeffs_sparse.h`。整理后的 MATLAB 输出采用新格式，rho 缓存生成方式也有区别，不能直接改名覆盖原表，见模块 README。
- 碰撞检测保留原数值模型和启动经验项，两者都有已知问题。残差超阈值不等于已经确认碰撞。
- 两块使用不同的模型，参数不能混用。PC 测试验证实现和接口，不代表 MCU 时序及实车验证通过；完整 MATLAB 符号流程尚未重新运行。

## 许可证

[MIT License](LICENSE)。MPC 算法来源及第三方说明见 [NOTICE](NOTICE)，碰撞代码提取说明见 [collision_detection/NOTICE](collision_detection/NOTICE)。
