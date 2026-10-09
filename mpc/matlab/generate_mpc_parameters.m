function [mpc,fit] = generate_mpc_parameters(p,dynamics)
%GENERATE_MPC_PARAMETERS 由连续 A/B 生成、拟合并导出 TinyMPC 参数。
%
% 这个文件完成生成流程的后半段：
%   continuous Ac/Bc/fc -> discrete Ad/Bd/fd
%   discrete model      -> K/P/Quu_inv
%   all leg samples     -> two-variable polynomial coefficients

% Q、R 是原始跟踪权重。计算 TinyMPC Riccati 矩阵时明确使用
% Q+rho*I 和 R+rho*I。

% 文件末尾只是普通数值计算的局部函数。日常调参只改
% wheelleg_parameters.m，不需要修改这些局部函数。


nx = 10;
nu = 4;
sample_count = dynamics.count;

if ~isequal(size(p.Q),[nx,nx]) || ~isequal(size(p.R),[nu,nu])
    error('Q must be 10x10 and R must be 4x4.');
end
if p.Ts <= 0 || p.rho <= 0
    error('Ts and rho must be positive.');
end

mpc.Ad = zeros(nx,nx,sample_count);
mpc.Bd = zeros(nx,nu,sample_count);
mpc.fd = zeros(nx,1,sample_count);
mpc.K = zeros(nu,nx,sample_count);
mpc.P = zeros(nx,nx,sample_count);
mpc.Quu_inv = zeros(nu,nu,sample_count);
mpc.P_chol = zeros(nx,nx,sample_count);
mpc.Quu_inv_chol = zeros(nu,nu,sample_count);
mpc.u_eq = dynamics.u_eq;

Q_rho = p.Q+p.rho*eye(nx);
R_rho = p.R+p.rho*eye(nu);

for k = 1:sample_count
    [Ad,Bd,fd] = discretize_model( ...
        dynamics.Ac(:,:,k),dynamics.Bc(:,:,k),dynamics.fc(:,:,k),p.Ts);
    [K,P,Quu_inv] = calculate_riccati(Ad,Bd,Q_rho,R_rho);

    mpc.Ad(:,:,k) = Ad;
    mpc.Bd(:,:,k) = Bd;
    mpc.fd(:,:,k) = fd;
    mpc.K(:,:,k) = K;
    mpc.P(:,:,k) = P;
    mpc.Quu_inv(:,:,k) = Quu_inv;
    mpc.P_chol(:,:,k) = chol(P,'lower');
    mpc.Quu_inv_chol(:,:,k) = chol(Quu_inv,'lower');
end

mpc.left_length = dynamics.left_length;
mpc.right_length = dynamics.right_length;

% 只拟合控制器最终使用的矩阵，不先拟合 A/B 再造成第二层拟合误差。
x = [dynamics.left_length,dynamics.right_length];
fit.A = fit_matrix(x,mpc.Ad,p.fit_degree,p.fit_ridge);
fit.B = fit_matrix(x,mpc.Bd,p.fit_degree,p.fit_ridge);
fit.f = fit_matrix(x,mpc.fd,p.fit_degree,p.fit_ridge);
fit.K = fit_matrix(x,mpc.K,p.fit_degree,p.fit_ridge);
fit.P_chol = fit_matrix(x,mpc.P_chol,p.fit_degree,p.fit_ridge);
fit.Quu_inv_chol = fit_matrix(x,mpc.Quu_inv_chol, ...
    p.fit_degree,p.fit_ridge);
fit.u_eq = fit_matrix(x,mpc.u_eq,p.fit_degree,p.fit_ridge);
fit.leg_range = dynamics.leg_range;

fprintf('  A              RMSE %.3e, maximum error %.3e\n', ...
    fit.A.rmse,fit.A.max_error);
fprintf('  B              RMSE %.3e, maximum error %.3e\n', ...
    fit.B.rmse,fit.B.max_error);
fprintf('  f              RMSE %.3e, maximum error %.3e\n', ...
    fit.f.rmse,fit.f.max_error);
fprintf('  K              RMSE %.3e, maximum error %.3e\n', ...
    fit.K.rmse,fit.K.max_error);
fprintf('  P_chol         RMSE %.3e, maximum error %.3e\n', ...
    fit.P_chol.rmse,fit.P_chol.max_error);
fprintf('  Quu_inv_chol   RMSE %.3e, maximum error %.3e\n', ...
    fit.Quu_inv_chol.rmse,fit.Quu_inv_chol.max_error);
fprintf('  u_eq           RMSE %.3e, maximum error %.3e\n', ...
    fit.u_eq.rmse,fit.u_eq.max_error);

export_fit_header(p,fit);
end


function [Ad,Bd,fd] = discretize_model(Ac,Bc,fc,Ts)
% 精确零阶保持，同时离散化重力常量项 fc。
nx = size(Ac,1);
nu = size(Bc,2);
augmented = [Ac,Bc,fc;zeros(nu+1,nx+nu+1)];
discrete = expm(augmented*Ts);
Ad = discrete(1:nx,1:nx);
Bd = discrete(1:nx,nx+1:nx+nu);
fd = discrete(1:nx,end);
end


function [K,P,Quu_inv] = calculate_riccati(A,B,Q_rho,R_rho)
% TinyMPC 使用的无穷时域 Riccati 迭代。
P = Q_rho;
for iteration = 1:10000
    Quu = R_rho+B'*P*B;
    K = Quu\(B'*P*A);
    P_next = Q_rho+A'*P*(A-B*K);
    P_next = (P_next+P_next')/2;
    if max(abs(P_next(:)-P(:))) < 1e-11
        P = P_next;
        break;
    end
    P = P_next;
end
if iteration == 10000
    error('Riccati iteration did not converge. Check A/B/Q/R/rho.');
end

Quu = R_rho+B'*P*B;
K = Quu\(B'*P*A);
Quu_inv = Quu\eye(size(Quu));
end


function fitted = fit_matrix(x,values,degree,ridge)
% 每个矩阵元素都使用同一组“左腿长、右腿长”多项式基。
center = (max(x,[],1)+min(x,[],1))/2;
scale = (max(x,[],1)-min(x,[],1))/2;
if any(scale == 0)
    error('Left and right leg lengths must both contain multiple values.');
end
z = (x-center)./scale;

% 总次数多项式：1, zL, zR, zL^2, zL*zR, zR^2, ...
powers = zeros((degree+1)*(degree+2)/2,2);
row = 0;
for total_degree = 0:degree
    for left_power = total_degree:-1:0
        row = row+1;
        powers(row,:) = [left_power,total_degree-left_power];
    end
end

Phi = zeros(size(x,1),size(powers,1));
for i = 1:size(powers,1)
    Phi(:,i) = z(:,1).^powers(i,1).*z(:,2).^powers(i,2);
end

% 拟合前按 C 语言的行优先顺序排列每个矩阵。
matrix_rows = size(values,1);
matrix_cols = size(values,2);
y = zeros(size(x,1),matrix_rows*matrix_cols);
for k = 1:size(x,1)
    y(k,:) = reshape(values(:,:,k).',1,[]);
end

coefficients = (Phi'*Phi+ridge*eye(size(Phi,2)))\(Phi'*y);
error_value = Phi*coefficients-y;

fitted.coefficients = coefficients;
fitted.powers = powers;
fitted.center = center;
fitted.scale = scale;
fitted.rows = matrix_rows;
fitted.cols = matrix_cols;
fitted.rmse = sqrt(mean(error_value(:).^2));
fitted.max_error = max(abs(error_value(:)));
end


function export_fit_header(p,fit)
% 按 C 求解器使用的行优先顺序写出系数。
output_folder = fileparts(p.output_header);
if ~isempty(output_folder) && ~exist(output_folder,'dir')
    mkdir(output_folder);
end

fid = fopen(p.output_header,'w');
if fid < 0
    error('Cannot create output header: %s',p.output_header);
end

terms = size(fit.A.powers,1);
fprintf(fid,'/* Generated by mpc/matlab/generate_mpc_fit.m. */\n');
fprintf(fid,'#ifndef WHEELLEG_MPC_FIT_DATA_H\n');
fprintf(fid,'#define WHEELLEG_MPC_FIT_DATA_H\n\n');
fprintf(fid,'#define WMPC_FIT_TERMS %d\n',terms);
fprintf(fid,'#define WMPC_FIT_DEGREE %d\n',p.fit_degree);
fprintf(fid,'#define WMPC_FIT_HORIZON %d\n',p.horizon);
fprintf(fid,'#define WMPC_FIT_TS %s\n',float_text(p.Ts));
fprintf(fid,'#define WMPC_FIT_RHO %s\n',float_text(p.rho));
fprintf(fid,'#define WMPC_FIT_MAX_ADMM_ITERATIONS %d\n',p.max_iter);
fprintf(fid,'#define WMPC_FIT_ABS_PRI_TOL %s\n',float_text(p.abs_pri_tol));
fprintf(fid,'#define WMPC_FIT_ABS_DUA_TOL %s\n',float_text(p.abs_dua_tol));
fprintf(fid,'#define WMPC_FIT_LEG_MIN %s\n',float_text(fit.leg_range(1)));
fprintf(fid,'#define WMPC_FIT_LEG_MAX %s\n\n',float_text(fit.leg_range(2)));

write_integer_array(fid,'wmpc_fit_power_left',fit.A.powers(:,1));
write_integer_array(fid,'wmpc_fit_power_right',fit.A.powers(:,2));
write_float_array(fid,'wmpc_fit_leg_center',fit.A.center);
write_float_array(fid,'wmpc_fit_leg_scale',fit.A.scale);
write_float_array(fid,'wmpc_fit_A_coeff',fit.A.coefficients);
write_float_array(fid,'wmpc_fit_B_coeff',fit.B.coefficients);
write_float_array(fid,'wmpc_fit_f_coeff',fit.f.coefficients);
write_float_array(fid,'wmpc_fit_K_coeff',fit.K.coefficients);
write_float_array(fid,'wmpc_fit_P_chol_coeff',fit.P_chol.coefficients);
write_float_array(fid,'wmpc_fit_Quu_inv_chol_coeff', ...
    fit.Quu_inv_chol.coefficients);
write_float_array(fid,'wmpc_fit_u_eq_coeff',fit.u_eq.coefficients);
write_float_array(fid,'wmpc_fit_Q_diag',diag(p.Q));
write_float_array(fid,'wmpc_fit_R_diag',diag(p.R));
write_float_array(fid,'wmpc_fit_x_min',p.x_min);
write_float_array(fid,'wmpc_fit_x_max',p.x_max);
write_float_array(fid,'wmpc_fit_u_min',p.u_min);
write_float_array(fid,'wmpc_fit_u_max',p.u_max);

fprintf(fid,'#endif /* WHEELLEG_MPC_FIT_DATA_H */\n');
fclose(fid);
fprintf('已生成 %s\n',p.output_header);
end


function write_integer_array(fid,name,value)
value = reshape(value.',1,[]);
fprintf(fid,'static const unsigned char %s[%d] = {',name,numel(value));
for i = 1:numel(value)
    if i > 1
        fprintf(fid,', ');
    end
    fprintf(fid,'%d',value(i));
end
fprintf(fid,'};\n\n');
end


function write_float_array(fid,name,value)
% 先转置，使 MATLAB 矩阵写成 C 的行优先数组。
value = reshape(value.',1,[]);
fprintf(fid,'static const float %s[%d] = {\n    ',name,numel(value));
for i = 1:numel(value)
    fprintf(fid,'%s',float_text(value(i)));
    if i < numel(value)
        fprintf(fid,', ');
    end
    if mod(i,6) == 0 && i < numel(value)
        fprintf(fid,'\n    ');
    end
end
fprintf(fid,'\n};\n\n');
end


function text = float_text(value)
text = sprintf('%.9g',value);
if isempty(strfind(text,'.')) && isempty(strfind(lower(text),'e'))
    text = [text,'.0'];
end
text = [text,'f'];
end
