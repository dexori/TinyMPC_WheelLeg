/* SPDX-License-Identifier: MIT
 * 数值函数逐项取自原 perception.cpp。legacy_C 已知不满足标准动量恒等式。
 * 不声称该固定模型覆盖变腿长，也不声称它与重新推导的 C 一致。 */
#include "wheelleg_mcg.h"
#include <math.h>
#include <stddef.h>
#include <string.h>

static void legacy_M(const float q[5], float M[5][5])
{
	(void)q[0];
	(void)q[1];
	const float thll = q[2];
	const float thlr = q[3];
	const float thb = q[4];
	const float cthb = cosf(thb);
	const float cthll = cosf(thll);
	const float cthlr = cosf(thlr);
	const float sthb = sinf(thb);
	const float sthll = sinf(thll);
	const float sthlr = sinf(thlr);
	const float t9 = -0.15824f * cthb;
	const float t10 = 0.4766251173525303f * cthll;
	const float t11 = 0.4766251173525303f * cthlr;
	const float t17 = -0.03234001207326082f * sthb * sthll;
	const float t18 = -0.03234001207326082f * sthb * sthlr;
	const float t22 = 0.2423012752720828f - 0.1167731537513699f * cthll;
	const float t23 = 0.1167731537513699f * cthlr - 0.2423012752720828f;
	const float t24 = 0.09740939111223574f * cthll + 0.09740939111223574f * cthlr +
	                  0.8869211185970738f * sthll * sthlr - 0.2021219683783261f;
	const float m33 = 1.134217461435397f - 0.1948187822244715f * cthll -
	                  0.7654366621452406f * cthll * cthll;
	const float m44 = 1.134217461435397f - 0.1948187822244715f * cthlr -
	                  0.7654366621452406f * cthlr * cthlr;
	M[0][0] = 24.6191124260355f;
	M[0][1] = 0.0f;
	M[0][2] = t10;
	M[0][3] = t11;
	M[0][4] = t9;
	M[1][0] = 0.0f;
	M[1][1] = 1.185467723372781f;
	M[1][2] = t22;
	M[1][3] = t23;
	M[1][4] = 0.0f;
	M[2][0] = t10;
	M[2][1] = t22;
	M[2][2] = m33;
	M[2][3] = t24;
	M[2][4] = t17;
	M[3][0] = t11;
	M[3][1] = t23;
	M[3][2] = t24;
	M[3][3] = m44;
	M[3][4] = t18;
	M[4][0] = t9;
	M[4][1] = 0.0f;
	M[4][2] = t17;
	M[4][3] = t18;
	M[4][4] = 0.56126592f;
}
static void legacy_C(const float q[5], const float qd[5], float C[5][5])
{
	(void)q[0];
	(void)q[1];
	(void)qd[0];
	(void)qd[1];
	const float thll = q[2];
	const float thlr = q[3];
	const float thb = q[4];
	const float dthll = qd[2];
	const float dthlr = qd[3];
	const float dthb = qd[4];
	const float cthb = cosf(thb);
	const float cthll = cosf(thll);
	const float cthlr = cosf(thlr);
	const float sthb = sinf(thb);
	const float sthll = sinf(thll);
	const float sthlr = sinf(thlr);
	const float v11 = -1.906500469410121f * dthll * sthll;
	const float v12 = 0.4670926150054797f * dthll * sthll;
	const float v13 = dthll * (0.38963756444894293f * sthll + 1.530873324290481f * sinf(2.0f * thll));
	const float v14 = dthll * (-0.38963756444894293f * sthll + 3.5476844743882956f * cthll * sthlr);
	const float v15 = -0.1293600482930433f * dthll * cthll * sthb;
	const float v16 = -1.906500469410121f * dthlr * sthlr;
	const float v17 = -0.4670926150054797f * dthlr * sthlr;
	const float v18 = dthlr * (-0.38963756444894293f * sthlr + 3.5476844743882956f * cthlr * sthll);
	const float v19 = dthlr * (0.38963756444894293f * sthlr + 1.530873324290481f * sinf(2.0f * thlr));
	const float v20 = -0.1293600482930433f * dthlr * cthlr * sthb;
	const float v21 = 0.63296f * dthb * sthb;
	const float v23 = -0.1293600482930433f * dthb * cthb * sthll;
	const float v24 = -0.1293600482930433f * dthb * cthb * sthlr;
	C[0][0] = 0.0f; C[0][1] = 0.0f; C[0][2] = v11; C[0][3] = v16; C[0][4] = v21;
	C[1][0] = 0.0f; C[1][1] = 0.0f; C[1][2] = v12; C[1][3] = v17; C[1][4] = 0.0f;
	C[2][0] = 0.0f; C[2][1] = 0.0f; C[2][2] = v13; C[2][3] = v18; C[2][4] = v23;
	C[3][0] = 0.0f; C[3][1] = 0.0f; C[3][2] = v14; C[3][3] = v19; C[3][4] = v24;
	C[4][0] = 0.0f; C[4][1] = 0.0f; C[4][2] = v15; C[4][3] = v20; C[4][4] = 0.0f;
}
static void legacy_G(const float q[5], float G[5])
{
	(void)q[0];
	(void)q[1];
	const float thll = q[2];
	const float thlr = q[3];
	const float thb = q[4];
	G[0] = 0.0f;
	G[1] = 0.0f;
	G[2] = -44.3326322060644f * sinf(thll);
	G[3] = -44.3326322060644f * sinf(thlr);
	G[4] = 1.5523344f * sinf(thb);
}
static void legacy_Q(float Twl, float Twr, float Tbl_l, float Tbl_r, float Q[5])
{
	Q[0] = (200.0f / 13.0f) * (Twl + Twr);
	Q[1] = (49.0f / 13.0f) * (Twr - Twl);
	Q[2] = Tbl_l - Twl;
	Q[3] = Tbl_r - Twr;
	Q[4] = -(Tbl_l + Tbl_r);
}


int wheelleg_mcg_legacy(const float q[MO_DOF], const float qd[MO_DOF],
                        MomentumObsModel *model)
{
    int i, j;
    if (model == NULL) return 0;
    memset(model, 0, sizeof(*model));
    if (q == NULL || qd == NULL) return 0;
    for (i = 0; i < MO_DOF; ++i)
        if (!isfinite(q[i]) || !isfinite(qd[i])) return 0;
    legacy_M(q, model->M);
    legacy_C(q, qd, model->C);
    legacy_G(q, model->G);
    for (i = 0; i < MO_DOF; ++i) {
        if (!isfinite(model->G[i])) goto invalid;
        for (j = 0; j < MO_DOF; ++j)
            if (!isfinite(model->M[i][j]) || !isfinite(model->C[i][j]))
                goto invalid;
    }
    return 1;
invalid:
    memset(model, 0, sizeof(*model));
    return 0;
}

int wheelleg_generalized_force(const float u[4], float Rw, float Rl,
                               float tau[MO_DOF])
{
    int i;
    if (tau == NULL) return 0;
    memset(tau, 0, sizeof(float) * MO_DOF);
    if (u == NULL || !isfinite(Rw) || Rw <= 0 ||
        !isfinite(Rl) || Rl <= 0) return 0;
    for (i = 0; i < 4; ++i) if (!isfinite(u[i])) return 0;
    /* 旧默认几何直接沿用原式以便对拍；其他几何使用显式半径。 */
    if (Rw == .065f && Rl == .245f) {
        legacy_Q(u[0], u[1], u[2], u[3], tau);
    } else {
        tau[0] = (u[0] + u[1]) / Rw;
        tau[1] = (Rl / Rw) * (u[1] - u[0]);
        tau[2] = u[2] - u[0];
        tau[3] = u[3] - u[1];
        tau[4] = -(u[2] + u[3]);
    }
    for (i = 0; i < MO_DOF; ++i) {
        if (!isfinite(tau[i])) {
            memset(tau, 0, sizeof(float) * MO_DOF);
            return 0;
        }
    }
    return 1;
}
