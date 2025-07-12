function setup_osqp()
    clearvars -global prob_osqp;
    global prob_osqp;

    %% === 车辆参数 ===
    m = 1274 + 71 * 2;
    a = 1.015;
    b = 1.895;
    Iz = 1536.7;
    k1 = -148970;
    k2 = -82204;
    T = 0.01;
    vx_nominal = 50 / 3.6;

    %% === 连续模型并离散化 ===
    A2 = [0, vx_nominal, 1, 0;
          0, 0, 0, 1;
          0, 0, (k1 + k2)/(m * vx_nominal), (a*k1 - b*k2)/(m * vx_nominal) - vx_nominal;
          0, 0, (a*k1 - b*k2)/(Iz * vx_nominal), (a^2*k1 + b^2*k2)/(Iz * vx_nominal)];

    B2 = [0; 0; -k1/m; -a*k1/Iz];

    A1 = eye(4) + T * A2;
    B1 = T * B2;

    %% === 增广系统建模 ===
    Nx = 4;
    Nu = 1;
    A_aug = [A1, B1; zeros(Nu, Nx), eye(Nu)];
    B_aug = [B1; eye(Nu)];

    A = A_aug;
    B = B_aug;

    n = Nx + Nu;   % = 5
    m = Nu;        % = 1
    N = 100;

    %% === 代价函数权重 ===
    Q_diag = [1, 1, 1, 1, 0];  % 最后一位是 u_{k-1}
    R_diag = 500;

    Q  = kron(speye(N), diag(Q_diag));
    QN = diag(Q_diag);
    R  = R_diag * eye(N - 1);

    P = blkdiag(Q, QN, R);
    q = zeros(size(P,1),1);

    %% === 状态/输入约束 ===
    x_min = repmat([-100; -100; -100; -100; -0.44], N, 1);
    x_max = repmat([ 100;  100;  100;  100;  0.44], N, 1);
    u_min = repmat(-0.005, N - 1, 1);
    u_max = repmat( 0.005, N - 1, 1);

    l = [x_min; u_min];
    u = [x_max; u_max];

    %% === 构建 A 矩阵（匹配 OSQP 接口要求） ===
    nv = size(P, 1);  % 决策变量数量
    nc = length(l);   % 约束数量

    A_osqp = speye(nc, nv);  % ✅ A ∈ ℝ^{nc × nv}

    %% === 保存至结构体 ===
    prob_osqp.P = sparse((P + P') / 2);  % 强制对称
    prob_osqp.q = q;
    prob_osqp.A = sparse(A_osqp);
    prob_osqp.l = l;
    prob_osqp.u = u;
    prob_osqp.nx = n;
    prob_osqp.nu = m;
    prob_osqp.N  = N;

    %% === 设置 OSQP 求解器 ===
    prob_osqp.solver = osqp;
    prob_osqp.solver.setup(prob_osqp.P, prob_osqp.q, prob_osqp.A, prob_osqp.l, prob_osqp.u, ...
        'max_iter', 10000, ...
        'eps_abs', 1e-4, ...
        'eps_rel', 1e-4, ...
        'verbose', true);

    disp('[OSQP] Augmented MPC problem setup complete.');
end
