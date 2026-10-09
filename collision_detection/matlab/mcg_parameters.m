function p = mcg_parameters()
% 参数取自原 matlab_balance_tandem/MCG.m，SI 单位。
% 固定腿长工作点，不是实时变腿长模型。
p.Rw = 0.065;
p.Rl = 0.490 / 2;
p.mw = 0.66;
p.ml = 1.605;
p.mb = 19.78;
p.Iw = 0.000653;
p.Ib = 0.56;
p.Iz = 0.895;
p.lc = -0.008;
p.g = 9.81;

% MCG 推导唯一需要的机构参数，不包含 VMC/LQR/MPC。
[L0, la, I] = balance_canshujisuan(60, 60, ...
    0.1134, 0.1134, 0.135, 0.135, 0.210, 0.250, ...
    0.345, 0.3069, 0.011357, 0.002353, 1.188, 0.417);
p.ll = L0; p.lr = L0;
p.lwl = la; p.lwr = la;
p.lbl = L0 - la; p.lbr = L0 - la;
p.I_ll = I; p.I_lr = I;
end
