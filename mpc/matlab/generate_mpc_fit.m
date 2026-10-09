function result = generate_mpc_fit(p)
%GENERATE_MPC_FIT 唯一主程序，MATLAB 代码从这里开始看。
%
% 直接生成：
%   result = generate_mpc_fit();
%
% 临时修改参数：
%   p = wheelleg_parameters();
%   p.rho = 1.0;
%   p.angle_range_deg = 0:5:50;
%   result = generate_mpc_fit(p);

if nargin < 1
    p = wheelleg_parameters();
end

fprintf('[1/4] 读取 CAD 表，计算各角度的腿部参数\n');
leg = prepare_leg_data(p);

fprintf('[2/4] 推导公式，计算各左右腿组合的连续 A/B/G\n');
model = derive_reduced_dynamics();
ab = derive_ab_model(model);
gravity = derive_gravity_model(model);
dynamics = sample_ab_models(p,leg,ab,gravity);

fprintf('[3/4] 计算离散模型和 TinyMPC 矩阵\n');
fprintf('[4/4] 拟合最终矩阵，导出 C 头文件\n');
[mpc,fit] = generate_mpc_parameters(p,dynamics);

result.parameters = p;
result.dynamics = dynamics;
result.mpc = mpc;
result.fit = fit;

if ~isempty(p.output_mat)
    output_folder = fileparts(p.output_mat);
    if ~isempty(output_folder) && ~exist(output_folder,'dir')
        mkdir(output_folder);
    end
    save(p.output_mat,'result');
    fprintf('已保存 %s\n',p.output_mat);
end
end
