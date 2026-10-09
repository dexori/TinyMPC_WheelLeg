function p = wheelleg_parameters()
%WHEELLEG_PARAMETERS 日常需要修改的参数全部放在这个文件。
%
% 状态顺序，也是 Q 对角线顺序：
% [s, ds, yaw, dyaw, left_leg, dleft_leg, ...
%  right_leg, dright_leg, pitch, dpitch]
%
% 输入顺序，也是 R 对角线顺序：
% [left_wheel, right_wheel, left_body_leg, right_body_leg]

root = fileparts(mfilename('fullpath'));

%% 1. 变腿长采样和拟合
p.angle_range_deg = 0:50;  % 快速调试可改成 0:5:50
p.fit_degree = 3;
p.fit_ridge = 1e-10;

%% 2. MPC 参数
p.Ts = 0.004;
p.horizon = 8;
p.rho = 1.0;
p.Q = diag([160,10,100,10,500.38082,3,500.38082,3,3000,1]);
p.R = diag([5.25,5.25,0.5,0.5]);
p.x_min = -[1e6;10;5;10;1.8;4;1.8;4;0.8;5];
p.x_max =  [1e6;10;5;10;1.8;4;1.8;4;0.8;5];
p.u_min = [-6;-6;-45;-45];
p.u_max = [ 6; 6; 45; 45];
p.max_iter = 40;
p.abs_pri_tol = 1e-3;
p.abs_dua_tol = 1e-3;

%% 3. 整车动力学参数
p.Rw = 0.065;       % 车轮半径 [m]
p.Rl = 0.490/2;     % 半轮距 [m]
p.lc = 0.016;       % 机体关节到机体质心 [m]
p.mw = 0.66;        % 单个车轮质量 [kg]
p.ml = 1.605;       % 单侧等效腿质量 [kg]
p.mb = 18.41;       % 机体质量 [kg]
p.Iw = 0.000653;    % 单个车轮惯量 [kg*m^2]
p.Ib = 0.52977;     % 机体俯仰惯量 [kg*m^2]
p.Iz = 0.821;       % 整车航向惯量 [kg*m^2]
p.g = 9.81;

%% 4. 五连杆尺寸 [m]
p.AO = 0.1134;
p.OO1 = 0.1134;
p.AB = 0.135;
p.BO1 = 0.135;
p.l1 = 0.210;
p.l2 = 0.250;

%% 5. 外部 CAD 质心表和输出文件
p.mass_property_csv = fullfile(root,'data','leg_mass_properties.csv');
p.output_header = fullfile(root,'generated','wheelleg_mpc_fit_data.h');
p.output_mat = fullfile(root,'generated','wheelleg_mpc_fit.mat');

%% 6. 静态前馈力矩的数值正则项
p.ueq_regularization = 1e-8;
end
