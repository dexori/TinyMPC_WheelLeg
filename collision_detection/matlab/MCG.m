function model = MCG(p, output_dir)
% 从原 MCG.m 提取的能量建模，仅推导 M/C/G 和广义力映射。
% 保留原坐标和运动学近似；C 改由 Christoffel 符号构造。
% 需要 MATLAB Symbolic Math Toolbox。output_dir='' 时不写生成文件。
if nargin < 1 || isempty(p), p = mcg_parameters(); end
if nargin < 2
    output_dir = fullfile(fileparts(mfilename('fullpath')), 'generated');
end
fields = {'Rw','Rl','mw','ml','mb','Iw','Ib','Iz','lc','g', ...
    'll','lr','lwl','lwr','lbl','lbr','I_ll','I_lr'};
for i = 1:numel(fields)
    value = p.(fields{i});
    assert(isnumeric(value) && isscalar(value) && isreal(value) && ...
        isfinite(value), 'Invalid parameter: %s', fields{i});
    % 用有理数避免双精度常数在符号表达式中展开为巨大分母。
    a.(fields{i}) = sym(value, 'r');
end
assert(p.Rw > 0 && p.Rl > 0 && p.ll > 0 && p.lr > 0);
assert(all([p.mw,p.ml,p.mb,p.Iw,p.Ib,p.Iz,p.I_ll,p.I_lr,p.g] > 0));
assert(p.lwl >= 0 && p.lwr >= 0 && p.lbl >= 0 && p.lbr >= 0);
assert(abs(p.lwl+p.lbl-p.ll) < 1e-10 && ...
       abs(p.lwr+p.lbr-p.lr) < 1e-10, 'COM distances do not sum to leg length');

syms s phi thll thlr thb sd dphi dthll dthlr dthb real
syms sdd phidd ddthll ddthlr ddthb real
syms Twl Twr Tbl_l Tbl_r real
q = [s; phi; thll; thlr; thb];
qd = [sd; dphi; dthll; dthlr; dthb];
qdd = [sdd; phidd; ddthll; ddthlr; ddthb];
u = [Twl; Twr; Tbl_l; Tbl_r];

% 沿用原 MCG.m 的小角度轮角逆解，不声称完整三维运动学。
thwl = s/a.Rw - a.Rl/a.Rw*phi - a.ll/(2*a.Rw)*thll + a.lr/(2*a.Rw)*thlr;
thwr = s/a.Rw + a.Rl/a.Rw*phi + a.ll/(2*a.Rw)*thll - a.lr/(2*a.Rw)*thlr;
dthwl = jacobian(thwl,q)*qd;
dthwr = jacobian(thwr,q)*qd;
h_hip = a.ll/2*cos(thll) + a.lr/2*cos(thlr);
s_b = s + a.lc*sin(thb);
h_b = h_hip + a.lc*cos(thb);
s_ll = a.Rw*thwl + a.lwl*sin(thll);
s_lr = a.Rw*thwr + a.lwr*sin(thlr);
h_ll = h_hip - a.lbl*cos(thll);
h_lr = h_hip - a.lbr*cos(thlr);

ds_b = jacobian(s_b,q)*qd; dh_b = jacobian(h_b,q)*qd;
ds_ll = jacobian(s_ll,q)*qd; dh_ll = jacobian(h_ll,q)*qd;
ds_lr = jacobian(s_lr,q)*qd; dh_lr = jacobian(h_lr,q)*qd;
T = a.mw/2*(a.Rw*dthwl)^2 + a.mw/2*(a.Rw*dthwr)^2 ...
    + a.ml/2*(ds_ll^2+dh_ll^2+ds_lr^2+dh_lr^2) ...
    + a.mb/2*(ds_b^2+dh_b^2) ...
    + a.Iw/2*(dthwl^2+dthwr^2) ...
    + a.I_ll/2*dthll^2 + a.I_lr/2*dthlr^2 ...
    + a.Ib/2*dthb^2 + a.Iz/2*dphi^2;
V = a.ml*a.g*(h_ll+h_lr) + a.mb*a.g*h_b;
L = T - V;
dLdqd = jacobian(L,qd).';
E = simplify(jacobian(dLdqd,q)*qd + jacobian(dLdqd,qd)*qdd ...
    - jacobian(L,q).');
M = simplify(jacobian(E,qdd));
G = simplify(jacobian(V,q).');

% 不使用旧 jacobian(rest,qd)*qd -> 再取 Jacobian 的方式。
% 对二次速度项，那样会重复引入倍率，也不保证 Mdot=C+C'。
C = sym(zeros(5));
for i = 1:5
    for j = 1:5
        for k = 1:5
            C(i,j) = C(i,j) + (diff(M(i,j),q(k)) ...
                + diff(M(i,k),q(j)) - diff(M(j,k),q(i))) * qd(k)/2;
        end
    end
end
C = simplify(C);
Mdot = sym(zeros(5));
for k = 1:5, Mdot = Mdot + diff(M,q(k))*qd(k); end

% 映射沿用原工程定义，电机正方向和实际力矩换算须在上层对齐。
Tau = [(Twl+Twr)/a.Rw; a.Rl/a.Rw*(Twr-Twl); ...
    Tbl_l-Twl; Tbl_r-Twr; -(Tbl_l+Tbl_r)];
checks.symmetry = simplify(M-M.');
checks.dynamics = simplify(E-M*qdd-C*qd-G);
checks.momentum = simplify(Mdot-C-C.');
names = fieldnames(checks);
for i = 1:numel(names)
    v = checks.(names{i});
    ok = isAlways(v == sym(zeros(size(v))), 'Unknown', 'false');
    assert(all(ok(:)), 'Symbolic check failed: %s', names{i});
end

model.q = q; model.qd = qd; model.qdd = qdd; model.u = u;
model.M = M; model.C = C; model.G = G; model.Tau = Tau;
model.E = E; model.checks = checks; model.parameters = p;
model.functions.M = matlabFunction(M, 'Vars', {q}, 'Optimize', false);
model.functions.C = matlabFunction(C, 'Vars', {q,qd}, 'Optimize', false);
model.functions.G = matlabFunction(G, 'Vars', {q}, 'Optimize', false);
model.functions.Tau = matlabFunction(Tau, 'Vars', {u}, 'Optimize', false);
if ~isempty(output_dir)
    if ~exist(output_dir,'dir'), mkdir(output_dir); end
    matlabFunction(M,'File',fullfile(output_dir,'M_fun'),'Vars',{q});
    matlabFunction(C,'File',fullfile(output_dir,'C_fun'),'Vars',{q,qd});
    matlabFunction(G,'File',fullfile(output_dir,'G_fun'),'Vars',{q});
    matlabFunction(Tau,'File',fullfile(output_dir,'Q_fun'),'Vars',{u});
    save(fullfile(output_dir,'wheelleg_mcg.mat'),'model');
end
end
