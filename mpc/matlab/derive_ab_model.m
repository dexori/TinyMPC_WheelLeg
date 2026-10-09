function ab = derive_ab_model(model)
%DERIVE_AB_MODEL 从 M、B_tau、G 推导连续 10 状态 A/B 矩阵。
%
% Aqq = -M0^(-1)*dGdq0
% Bu  =  M0^(-1)*B_tau0

zero_point = zeros(14,1);  % q(5), qdot(5), u(4)
from = [model.q;model.qd;model.u];
M0 = simplify(subs(model.M,from,zero_point));
B_tau0 = simplify(subs(model.B_tau,from,zero_point));
dGdq0 = simplify(subs(model.dGdq,from,zero_point));

Aqq = simplify(-(M0\dGdq0));
Bu = simplify(M0\B_tau0);

Ac = sym(zeros(10,10));
Bc = sym(zeros(10,4));
Ac(1,2)=1; Ac(3,4)=1; Ac(5,6)=1; Ac(7,8)=1; Ac(9,10)=1;
Ac([2,4,6,8,10],[1,3,5,7,9]) = Aqq;
Bc([2,4,6,8,10],:) = Bu;

ab.fA = matlabFunction(Ac,'Vars',{model.parameters});
ab.fB = matlabFunction(Bc,'Vars',{model.parameters});
ab.fM = matlabFunction(M0,'Vars',{model.parameters});
ab.fB_tau = matlabFunction(B_tau0,'Vars',{model.parameters});
end
