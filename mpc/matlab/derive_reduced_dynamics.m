function model = derive_reduced_dynamics()
%DERIVE_REDUCED_DYNAMICS 写出原始受力方程，消去内力，推导 M、B_tau、G。
%
% 本项目统一使用：
%   M(q)*qdd = B_tau(q)*u - G(q)
%
% q = [position, yaw, left_leg_angle, right_leg_angle, pitch]'.

if exist('sym','file') ~= 2
    error(['derive_reduced_dynamics requires MATLAB Symbolic Math Toolbox. ', ...
        'Install it before running generate_mpc_fit.']);
end

syms position yaw theta_l theta_r pitch real
syms velocity yaw_rate theta_l_rate theta_r_rate pitch_rate real
syms position_accel yaw_accel theta_l_accel theta_r_accel pitch_accel real
syms tau_w_l tau_w_r tau_b_l tau_b_r real
syms wheel_l_accel wheel_r_accel real

syms f_l f_r Fws_l Fws_r Fwh_l Fwh_r Fbs_l Fbs_r Fbh_l Fbh_r real
syms Rw Rl lc mw ml mb Iw Ill Ilr Ib Iz gravity real
syms ll lr lwl lwr lbl lbr real
syms body_offset_l body_offset_r wheel_offset_l wheel_offset_r real

q = [position;yaw;theta_l;theta_r;pitch];
qd = [velocity;yaw_rate;theta_l_rate;theta_r_rate;pitch_rate];
qdd = [position_accel;yaw_accel;theta_l_accel;theta_r_accel;pitch_accel];
u = [tau_w_l;tau_w_r;tau_b_l;tau_b_r];

%% 运动学
body_horizontal_accel = Rw*(wheel_l_accel+wheel_r_accel)/2 ...
    +ll*theta_l_accel/2+lr*theta_r_accel/2;
body_vertical_accel = -(ll/2)*theta_l*theta_l_accel ...
    -(lr/2)*theta_r*theta_r_accel;
leg_l_horizontal_accel = Rw*wheel_l_accel+lwl*theta_l_accel;
leg_r_horizontal_accel = Rw*wheel_r_accel+lwr*theta_r_accel;
leg_l_vertical_accel = body_vertical_accel+lbl*theta_l*theta_l_accel;
leg_r_vertical_accel = body_vertical_accel+lbr*theta_r*theta_r_accel;

eq_position = position_accel == Rw*(wheel_l_accel+wheel_r_accel)/2;
eq_yaw = yaw_accel == Rw*(-wheel_l_accel+wheel_r_accel)/(2*Rl) ...
    -ll*theta_l_accel/(2*Rl)+lr*theta_r_accel/(2*Rl);

%% 动力学方程
eq_wheel_l_force = mw*Rw*wheel_l_accel == f_l-Fws_l;
eq_wheel_r_force = mw*Rw*wheel_r_accel == f_r-Fws_r;
eq_wheel_l_torque = Iw*wheel_l_accel == tau_w_l-f_l*Rw;
eq_wheel_r_torque = Iw*wheel_r_accel == tau_w_r-f_r*Rw;
eq_leg_l_horizontal = ml*leg_l_horizontal_accel == Fws_l-Fbs_l;
eq_leg_r_horizontal = ml*leg_r_horizontal_accel == Fws_r-Fbs_r;
eq_leg_l_vertical = ml*leg_l_vertical_accel == Fwh_l-Fbh_l-ml*gravity;
eq_leg_r_vertical = ml*leg_r_vertical_accel == Fwh_r-Fbh_r-ml*gravity;

wheel_l_sin = taylor(sin(theta_l-wheel_offset_l),theta_l,0,'Order',2);
wheel_r_sin = taylor(sin(theta_r-wheel_offset_r),theta_r,0,'Order',2);
wheel_l_cos = taylor(cos(theta_l-wheel_offset_l),theta_l,0,'Order',2);
wheel_r_cos = taylor(cos(theta_r-wheel_offset_r),theta_r,0,'Order',2);
body_l_sin = taylor(sin(theta_l+body_offset_l),theta_l,0,'Order',2);
body_r_sin = taylor(sin(theta_r+body_offset_r),theta_r,0,'Order',2);
body_l_cos = taylor(cos(theta_l+body_offset_l),theta_l,0,'Order',2);
body_r_cos = taylor(cos(theta_r+body_offset_r),theta_r,0,'Order',2);

eq_leg_l_torque = Ill*theta_l_accel == ...
    Fwh_l*lwl*wheel_l_sin+Fbh_l*lbl*body_l_sin ...
    -Fws_l*lwl*wheel_l_cos-Fbs_l*lbl*body_l_cos-tau_w_l+tau_b_l;
eq_leg_r_torque = Ilr*theta_r_accel == ...
    Fwh_r*lwr*wheel_r_sin+Fbh_r*lbr*body_r_sin ...
    -Fws_r*lwr*wheel_r_cos-Fbs_r*lbr*body_r_cos-tau_w_r+tau_b_r;

eq_body_horizontal = mb*body_horizontal_accel == Fbs_l+Fbs_r;
eq_body_vertical = mb*body_vertical_accel == Fbh_l+Fbh_r-mb*gravity;
eq_body_pitch = Ib*pitch_accel == -(tau_b_l+tau_b_r) ...
    -(Fbs_l+Fbs_r)*lc+(Fbh_l+Fbh_r)*lc*pitch;
eq_body_yaw = Iz*yaw_accel == (-f_l+f_r)*Rl;
eq_equal_vertical_force = Fwh_l == Fwh_r;

%% 消去左右车轮角加速度
wheel_solution = solve([eq_position,eq_yaw], ...
    [wheel_l_accel,wheel_r_accel],'ReturnConditions',false);
replace_wheels = @(expr) simplify(subs(expr, ...
    [wheel_l_accel,wheel_r_accel], ...
    [wheel_solution.wheel_l_accel,wheel_solution.wheel_r_accel]));

%% 用十条线性方程消去十个内部力
internal_forces = [f_l;f_r;Fws_l;Fws_r;Fwh_l;Fwh_r; ...
    Fbs_l;Fbs_r;Fbh_l;Fbh_r];
force_equations = [ ...
    lhs(eq_wheel_l_force)-rhs(eq_wheel_l_force); ...
    lhs(eq_wheel_r_force)-rhs(eq_wheel_r_force); ...
    lhs(eq_wheel_l_torque)-rhs(eq_wheel_l_torque); ...
    lhs(eq_wheel_r_torque)-rhs(eq_wheel_r_torque); ...
    lhs(eq_leg_l_horizontal)-rhs(eq_leg_l_horizontal); ...
    lhs(eq_leg_r_horizontal)-rhs(eq_leg_r_horizontal); ...
    lhs(eq_leg_l_vertical)-rhs(eq_leg_l_vertical); ...
    lhs(eq_leg_r_vertical)-rhs(eq_leg_r_vertical); ...
    lhs(eq_body_vertical)-rhs(eq_body_vertical); ...
    lhs(eq_equal_vertical_force)-rhs(eq_equal_vertical_force)];
force_equations = replace_wheels(force_equations);
[force_A,force_b] = equationsToMatrix(force_equations,internal_forces);
force_solution = simplify(force_A\force_b);

reduced_equations = [ ...
    lhs(eq_body_horizontal)-rhs(eq_body_horizontal); ...
    lhs(eq_body_yaw)-rhs(eq_body_yaw); ...
    lhs(eq_leg_l_torque)-rhs(eq_leg_l_torque); ...
    lhs(eq_leg_r_torque)-rhs(eq_leg_r_torque); ...
    lhs(eq_body_pitch)-rhs(eq_body_pitch)];
reduced_equations = replace_wheels(reduced_equations);
reduced_equations = simplify(subs( ...
    reduced_equations,internal_forces,force_solution));

%% 整理出 M、B_tau 和 G
[M,rhs_all] = equationsToMatrix(reduced_equations,qdd);
M = simplify(M);
rhs_all = simplify(rhs_all);
B_tau = simplify(jacobian(rhs_all,u));
unforced = simplify(rhs_all-B_tau*u);
G = simplify(-subs(unforced,qd,zeros(5,1)));
velocity_terms = simplify(unforced+G);

check = isAlways(simplify(rhs_all-(B_tau*u-G+velocity_terms)) == 0, ...
    'Unknown','false');
if ~all(check(:))
    error('The M/B_tau/G sign check failed.');
end

% 符号函数的参数顺序。sample_ab_models.m 必须按完全相同的顺序传值：
% 1~12  整车常量；13~18 左右腿长和质心距离；19~22 左右质心偏角。
parameters = [Rw;Rl;lc;mw;ml;mb;Iw;Ill;Ilr;Ib;Iz;gravity; ...
    ll;lr;lwl;lwr;lbl;lbr;body_offset_l;body_offset_r; ...
    wheel_offset_l;wheel_offset_r];

model.q = q;
model.qd = qd;
model.u = u;
model.parameters = parameters;
model.M = M;
model.B_tau = B_tau;
model.G = G;
model.dGdq = simplify(jacobian(G,q));
model.velocity_terms = velocity_terms;
end
