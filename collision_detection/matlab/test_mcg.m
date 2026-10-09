function test_mcg()
% MATLAB 内部符号/数值检查；不代表实车外力估计精度。
rng(1);
model = MCG(mcg_parameters(), '');
Mfun = model.functions.M; Cfun = model.functions.C;
for n = 1:20
    q = [0;0;0.2*(2*rand(3,1)-1)];
    qd = 0.5*(2*rand(5,1)-1);
    M = Mfun(q); C = Cfun(q,qd);
    h = 1e-6;
    Mdot = (Mfun(q+h*qd)-Mfun(q-h*qd))/(2*h);
    assert(norm(M-M.','fro') < 1e-10);
    assert(norm(Mdot-C-C.','fro') < 1e-7);
    [~,flag] = chol(M); assert(flag == 0);
end
fprintf('MCG symbolic identities and 20 numeric samples passed.\n');
end
