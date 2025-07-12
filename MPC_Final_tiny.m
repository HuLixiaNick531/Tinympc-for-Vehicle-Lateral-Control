
function [sys,x0,str,ts] = MPC_Final(t,x,u,flag)
% MPC_Final S-function adapted for TinyMPC, preserving original parameters
switch flag
  case 0
    [sys,x0,str,ts] = mdlInitializeSizes;
  case 2
    sys = mdlUpdate(t,x,u);
  case 3
    sys = mdlOutputs(t,x,u);
  case {1,4,9}
    sys = [];
  otherwise
    error(['Unhandled flag = ', num2str(flag)]);
end

function [sys,x0,str,ts] = mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 4;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 5;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0  = [0.0001, 0.0001, 0.0001, 0.0001]';
str = [];
ts  = [0.01 0];
global U solve_times
U = [0]; % Initial control input
solve_times = []
function sys = mdlUpdate(t,x,u)
sys = x;

function sys = mdlOutputs(t,x,u)
global prob U
global solve_times
%判断是否运行了setup程序
if isempty(prob)
    error('TinyMPC problem not initialized. Run setup_tinympc.m first.');
end



% Current states and velocity
vx = u(5);
x = u(1:4); % [lateral error, heading error, lateral vel error, yaw rate error]

% 创建新状态空间矩阵Augmented state: kesi = [x; U]
Nx = 4;
Nu = 1;
kesi = zeros(Nx + Nu, 1);
kesi(1) = u(1);
kesi(2) = u(2);
kesi(3) = u(3);
kesi(4) = u(4);
kesi(5) = U;

tic;
% Set initial state in TinyMPC
verbose_int = int32(0);
prob.set_x0(kesi, verbose_int);

% 测量求解时间

prob.solve(verbose_int);
solve_time = toc;

% 记录求解时间
solve_times(end+1) = solve_time;
if mod(length(solve_times), 100) == 0 % 每100次输出一次平均时间
    fprintf('TinyMPC - Iterations: %d, Avg Solve Time: %.6f s\n', length(solve_times), mean(solve_times));
end


% Get control increment sequence
delta_u = prob.get_u(verbose_int);

if isempty(delta_u) || length(delta_u) < 1
    error('TinyMPC returned an empty or invalid delta_u');
end

% Use first control increment (equivalent to u_piao = X(1))
u_piao = double(delta_u(1));

% Update control input (U = kesi(5) + u_piao)
U = kesi(5) + u_piao;

% % Enforce control input constraint
% u_real = min(max(U, -0.44), 0.44);
u_real = U;

% Debugging output
disp(['Time: ', num2str(t), ...
      ', States: [', num2str(x(1)), ', ', num2str(x(2)), ', ', num2str(x(3)), ', ', num2str(x(4)), ']', ...
      ', vx: ', num2str(vx), ...
      ', Control: ', num2str(u_real), ...
      ', Delta U: ', num2str(delta_u(1,10))]);

% Output
sys = [u_real];