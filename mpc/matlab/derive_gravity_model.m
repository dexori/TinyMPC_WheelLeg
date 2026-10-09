function gravity = derive_gravity_model(model)
%DERIVE_GRAVITY_MODEL 单独生成直立点的重力向量 G0。

zero_point = zeros(14,1);  % q(5), qdot(5), u(4)
G0 = simplify(subs(model.G,[model.q;model.qd;model.u],zero_point));

gravity.G = model.G;
gravity.G0 = G0;
gravity.fG0 = matlabFunction(G0,'Vars',{model.parameters});
end
