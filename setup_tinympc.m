% setup_tinympc.m
clear all;
close all;

tinympc_matlab_dir = 'C:\Users\13706\Documents\tinympc\tinympc-matlab';
tinympc_dir = fullfile(tinympc_matlab_dir, 'tinympc', 'TinyMPC');
addpath(genpath(tinympc_matlab_dir));

if libisloaded('libtinympcShared')
    unloadlibrary('libtinympcShared');
    disp('Unloaded existing libtinympcShared');
end

global prob;
prob = TinyMPC();

prob.compile_lib(tinympc_dir);

lib_dir = fullfile(tinympc_dir, '\build\src\tinympc\libtinympcShared.dll');
generic_header = fullfile(tinympc_matlab_dir, '/matlab_wrapper/codegen.hpp');
prob = prob.load_lib(lib_dir, generic_header);

% Vehicle parameters (same as original MPC_Final)
m = 1274 + 71*2;
a = 1.015;
b = 1.895;
Iz = 1536.7;
k1 = -148970;
k2 = -82204;
T = 0.01;

% Nominal longitudinal velocity (use a reasonable fixed value)
vx_nominal = 50/3.6; % 13.89 m/s, consistent with typical speed

% State-space matrices (computed at nominal vx)
A2 = [0, vx_nominal, 1, 0;
      0, 0, 0, 1;
      0, 0, (k1+k2)/m/vx_nominal, (a*k1-b*k2)/m/vx_nominal-vx_nominal;
      0, 0, (a*k1-b*k2)/Iz/vx_nominal, (a^2*k1+b^2*k2)/Iz/vx_nominal];
B2 = [0; 0; -k1/m; -a*k1/Iz];
A1 = A2*T + eye(4);
B1 = B2*T;

% Augmented state system (kesi = [x; u_prev])
Nx = 4;
Nu = 1;
A_aug = [A1, B1; zeros(Nu, Nx), eye(Nu)];
B_aug = [B1; eye(Nu)];

A=A_aug;
B=B_aug;

% A = reshape(A_aug.', 1, []);
% B = reshape(B_aug.', 1, []);

n = Nx + Nu;
m = Nu;
N = 100; % Matches Np = 100 from original code

% Cost weights (adapt q to diagonal form for augmented state)
Q_diag = [10000, 10000, 10000, 10000, 0]; % Original q diagonal [1,1,1,1], plus 0 for u_prev
R_diag = 500;            % Matches R = 1000 * eye(Nc * Nu)

% Constraints (same as original)
x_min = [-100, -100, -100, -100, -0.44]; % u_prev in [-0.44, 0.44]
x_max = [100, 100, 100, 100, 0.44];
u_min = -0.005; % Delta u in [-0.005, 0.005]
u_max = 0.005;
x_min = repmat(x_min, 1, N);
x_max = repmat(x_max, 1, N);
u_min = repmat(u_min, 1, N-1);
u_max = repmat(u_max, 1, N-1);

% Solver settings (adjusted for stability)
rho = 1.0;
abs_pri_tol = 0.001;
abs_dual_tol = 0.001;
max_iter = 1000; % Increased for better convergence
check_termination = 1;


% A=[1, 0.138888888888889, 0.0100000000000000, 0,	0, 0, 1, 0, 0.0100000000000000, 0, 0, 0, 0.882453898305085, -0.136564127871940, 1.05204802259887, 0, 0, 0.00214216281642481, 0.789782097421748,	0.983956204854558, 0, 0, 0, 0, 1];
% A=[1, 0.138888888888889, 0.0100000000000000, 0, 0];
% B=[0, 0, 1.05204802259887, 0.983956204854558, 1];
prob = prob.setup(n, m, N, A, B, Q_diag, R_diag, x_min, x_max, u_min, u_max, ...
    'rho', rho, 'abs_pri_tol', abs_pri_tol, 'abs_dual_tol', abs_dual_tol, ...
    'max_iter', max_iter, 'check_termination', check_termination);

output_dir = fullfile(tinympc_matlab_dir, 'generated_code');
prob.tiny_codegen(tinympc_dir, output_dir);
prob.compile_lib(output_dir);

prob.unload_lib('libtinympcShared');
compiled_lib_dir = fullfile(output_dir, '/build/tinympc/libtinympcShared.dll');
compiled_header = fullfile(tinympc_matlab_dir, '/matlab_wrapper/wrapper.hpp');
prob = prob.load_lib(compiled_lib_dir, compiled_header);

disp('TinyMPC setup complete. Run the Simulink model with MPC_Final.');