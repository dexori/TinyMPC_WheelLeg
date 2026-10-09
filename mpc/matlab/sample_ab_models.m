function samples = sample_ab_models(p,leg,ab,gravity)
%SAMPLE_AB_MODELS 逐个计算所有“左腿长度 × 右腿长度”的连续 A/B/f。

N = numel(leg);
sample_count = N*N;
samples.count = sample_count;
samples.left_length = zeros(sample_count,1);
samples.right_length = zeros(sample_count,1);
samples.Ac = zeros(10,10,sample_count);
samples.Bc = zeros(10,4,sample_count);
samples.fc = zeros(10,1,sample_count);
samples.u_eq = zeros(4,1,sample_count);
samples.equilibrium_residual = zeros(sample_count,1);

k = 0;
for i = 1:N
    fprintf('  left leg %d/%d, angle %.1f deg\n',i,N,leg(i).angle_deg);
    for j = 1:N
        k = k+1;

        % 这 22 项的顺序必须与 derive_reduced_dynamics.m 末尾一致：
        % 1~12 整车常量，13~18 腿长/质心距离，19~22 质心偏角。
        values = [p.Rw; p.Rl; p.lc; p.mw; p.ml; p.mb; p.Iw; ...
            leg(i).inertia; leg(j).inertia; p.Ib; p.Iz; p.g; ...
            leg(i).length; leg(j).length; ...
            leg(i).wheel_to_com; leg(j).wheel_to_com; ...
            leg(i).body_to_com; leg(j).body_to_com; ...
            leg(i).body_offset; leg(j).body_offset; ...
            leg(i).wheel_offset; leg(j).wheel_offset];

        Ac = double(ab.fA(values));
        Bc = double(ab.fB(values));
        M = double(ab.fM(values));
        B_tau = double(ab.fB_tau(values));
        G = double(gravity.fG0(values));

        u_eq = (B_tau'*B_tau+p.ueq_regularization*eye(4)) ...
            \ (B_tau'*G);
        fc = zeros(10,1);
        fc([2,4,6,8,10]) = -(M\G);

        samples.left_length(k) = leg(i).length;
        samples.right_length(k) = leg(j).length;
        samples.Ac(:,:,k) = Ac;
        samples.Bc(:,:,k) = Bc;
        samples.fc(:,:,k) = fc;
        samples.u_eq(:,:,k) = u_eq;
        samples.equilibrium_residual(k) = norm(B_tau*u_eq-G,inf);
    end
end

samples.leg_range = [min([samples.left_length;samples.right_length]), ...
    max([samples.left_length;samples.right_length])];
fprintf('  Maximum equilibrium residual: %.3e\n', ...
    max(samples.equilibrium_residual));
end
