/* SPDX-License-Identifier: MIT */
#ifndef WHEELLEG_MCG_H
#define WHEELLEG_MCG_H

#ifdef __cplusplus
extern "C" {
#endif

#define MO_DOF 5

typedef struct {
    float M[MO_DOF][MO_DOF];
    float C[MO_DOF][MO_DOF];
    float G[MO_DOF];
} MomentumObsModel;

/* 原 H7 固定数值系数。C 未通过 Mdot=C+C^T 检查，仅用于复现旧实现。
 * q=[s, yaw, left_leg, right_leg, pitch]，单位 m/rad。
 * 返回 1 表示输入及输出有限，0 表示无效；不代表模型物理验证通过。 */
int wheelleg_mcg_legacy(const float q[MO_DOF], const float qd[MO_DOF],
                        MomentumObsModel *model);

/* 四路力矩均以 N·m 输入，调用方负责电流换算、减速比及正负号。
 * u=[Twl,Twr,Tbl_left,Tbl_right]，Rw 为轮半径，Rl 为半轮距（m）。
 * tau[0] 单位 N，其余 N·m。返回 1 成功，0 输入无效。 */
int wheelleg_generalized_force(const float u[4], float Rw, float Rl,
                               float tau[MO_DOF]);

#ifdef __cplusplus
}
#endif
#endif
