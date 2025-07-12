function [sys,x0,str,ts] = MPC_Final(t,x,u,flag)
switch flag
  case 0
    [sys,x0,str,ts]=mdlInitializeSizes;
  case 2
    sys=mdlUpdate(t,x,u);
  case 3
    sys=mdlOutputs(t,x,u);
  case {1,4,9}
    sys=[];
  otherwise
    error(['Unhandled flag=',num2str(flag)]);
end

function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 4;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 5;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0  = [0.0001,0.0001,0.0001,0.0001]';
str = [];
ts  = [0.01 0];

global U solve_times
U = [0];  % Δ前轮转角
solve_times = [];

global A_cons B_cons lb ub x_start
global A1 B1 u_piao U kesi Q PHI THETA H
global Nx Nu Np Nc Row T

% === 参数设置 ===
Nx = 4;
Nu = 1;
Np = 100;
Nc = 80;
Row = 10;
T = 0.01;

m  = 1274 + 71*2;
a  = 1.015;
b  = 1.895;
Iz = 1536.7;
k1 = -148970;
k2 = -82204;
vx = 50 / 3.6;

q = diag([9, 1, 1, 1]);
Q = kron(eye(Np), q);
R = 500 * eye(Nc * Nu);

% 系统离散线性模型
A2 = [0, vx, 1, 0;
      0, 0, 0, 1;
      0, 0, (k1+k2)/m/vx, (a*k1-b*k2)/m/vx - vx;
      0, 0, (a*k1-b*k2)/Iz/vx, (a^2*k1 + b^2*k2)/Iz/vx];
B2 = [0; 0; -k1/m; -a*k1/Iz];

A1 = A2*T + eye(Nx);
B1 = B2*T;

A = [A1, B1;
     zeros(Nu, Nx), eye(Nu)];
B = [B1;
     eye(Nu)];
C = [eye(Nx), zeros(Nx, Nu)];

% === PHI 与 THETA 构造 ===
PHI = [];
THETA = [];
for i = 1:Np
    PHI = [PHI; C * (A^i)];
    row_theta = [];
    for j = 1:Nc
        if i >= j
            row_theta = [row_theta, C * A^(i-j) * B];
        else
            row_theta = [row_theta, zeros(Nx, Nu)];
        end
    end
    THETA = [THETA; row_theta];
end

% === H矩阵构造 ===
H = [THETA' * Q * THETA + R, zeros(Nc*Nu, 1);
     zeros(1, Nc*Nu), Row];

% === 约束 ===
A_t = tril(ones(Nc));
A_I = kron(A_t, eye(Nu));
Ut = kron(ones(Nc,1), U);
umax = 0.44;
umin = -0.44;
umax_dt = 0.05;
umin_dt = -0.05;
Umax = kron(ones(Nc,1), umax);
Umin = kron(ones(Nc,1), umin);
Umax_dt = kron(ones(Nc,1), umax_dt);
Umin_dt = kron(ones(Nc,1), umin_dt);

A_cons = [A_I, zeros(Nc*Nu,1);
         -A_I, zeros(Nc*Nu,1)];
B_cons = [Umax - Ut;
          -Umin + Ut];
lb = [Umin_dt];
ub = [Umax_dt];

function sys = mdlUpdate(~, x, ~)
sys = x;

function sys = mdlOutputs(t,x,u)
global Nx Nu Np Nc Row T U kesi Q PHI THETA H
global A_cons B_cons lb ub x_start solve_times

% === 当前状态 ===
kesi = [u(1:4); U];

% === 构造代价项 ===
error = PHI * kesi;
f = [2 * error' * Q * THETA, 0];

% === OSQP 求解 ===

persistent prob_osqp;
% 求解
tic;
if isempty(prob_osqp)
    prob_osqp = osqp;
    P = sparse((H + H') / 2);
    q = f(:);
    A_osqp = sparse(A_cons);
    l_osqp = -inf(size(B_cons));
    u_osqp = B_cons;

    prob_osqp.setup(P, q, A_osqp, l_osqp, u_osqp, ...
        'verbose', false, 'eps_abs', 1e-3, 'eps_rel', 1e-3, 'max_iter', 1000);
end

% 更新q
prob_osqp.update('q', f(:));


res = prob_osqp.solve();
solve_time = toc;
solve_times(end+1) = solve_time;

if mod(length(solve_times), 100) == 0
    fprintf('[OSQP] Iter: %d, Avg Time: %.6f s\n', ...
        length(solve_times), mean(solve_times));
end

if strcmp(res.info.status, 'solved')
    u_piao = res.x(1);
else
    warning('OSQP failed: %s', res.info.status);
    u_piao = 0;
end

% 输出控制量
U = kesi(5) + u_piao;
sys = U;
