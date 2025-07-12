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
    error.(['Unhandled flag=',num2str(flag)]);

end

function [sys,x0,str,ts]=mdlInitializeSizes

sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 4;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 5;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed
sys = simsizes(sizes);
x0  = [0.0001,0.0001,0.0001,0.0001]';
str = [];
ts  = [0.01 0];
global U solve_times
U=[0]; %Δ前轮转角
solve_times = []; % 全局变量，记录每次求解时间
global A_cons B_cons lb ub x_start
global A1 B1 u_piao U kesi Q PHI THETA H

global Nx Nu Np Nc Row T
Nx=4;%状态量4个，横向误差，航向角误差，侧向速度误差，横摆角速度误差
Nu=1;%控制量1个前轮转角
Np=100;%预测步长50
Nc=80;%控制步长30
Row=10;%松弛因子
T=0.01;%采用时间0.01
%%车辆参数
m=1274+71*2;%整车质量
a=1.015;%质心到前轴的距离
b=1.895;%质心到后轴的距离
Iz=1536.7;%绕Z轴的转动惯量
% k1=-112600;%前轴侧偏刚度,注意，这边的侧偏刚度已经乘以2；
% k2=-89500;%后轴侧偏刚度 
%通过拟合和魔术公式得到的精确版（可以写一下）；
k1=-148970;
k2=-82204;

%%控制器设计
vx=50/3.6;%%%%%这边多设计几个传入值，k1,k2,Mu;



%可以进一步调节这里的q和R
% q=[100,0,0,0;
%    0,1,0,0;
%    0,0,1,0;
%    0,0,0,1];
q=[10,0,0,0;
   0,1,0,0;
   0,0,1,0;
   0,0,0,1];
Q=kron(eye(Np),q);
% R=eye(Nc*Nu);
R=500*eye(Nc*Nu);

%雅可比矩阵获得线性化汽车动力学式子
A2=[0,vx,1,0;
   0,0,0,1;
   0,0,(k1+k2)/m/vx,(a*k1-b*k2)/m/vx-vx;
   0,0,(a*k1-b*k2)/Iz/vx,(a^2*k1+b^2*k2)/Iz/vx];
B2=[0;0;-k1/m;-a*k1/Iz];

%前向欧拉法获得
A1=A2*T+eye(Nx);
B1=B2*T;

%针对变量kesi，重新设计线性系统
A_cell=cell(2,2);
A_cell{1,1}=A1;
A_cell{1,2}=B1;
A_cell{2,1}=zeros(Nu,Nx);
A_cell{2,2}=eye(Nu);
A=cell2mat(A_cell);

B_cell=cell(2,1);
B_cell{1,1}=B1;
B_cell{2,1}=eye(Nu);
B=cell2mat(B_cell);

C_cell=cell(1,2);
C_cell{1,1}=eye(Nx);
C_cell{1,2}=zeros(Nx,Nu);
C=cell2mat(C_cell);


%% 预测时域矩阵升维
PHI_cell=cell(Np,1);
for i=1:1:Np
    PHI_cell{i,1}=C*A^i;
end
PHI=cell2mat(PHI_cell);

THETA_cell=cell(Np,Nc);
for i=1:1:Np
    for j=1:1:Nc
        if  i>=j
            THETA_cell{i,j}=C*A^(i-j)*B;
        else
            THETA_cell{i,j}=zeros(Nx,Nu);
        end
    end
end
THETA=cell2mat(THETA_cell);
% 获得H和f
H_cell=cell(2,2);
H_cell{1,1}=THETA'*Q*THETA+R;
H_cell{1,2}=zeros(Nc*Nu,1);
H_cell{2,1}=zeros(1,Nc*Nu);
H_cell{2,2}=Row; %松弛因子
H=cell2mat(H_cell);

% 约束条件：
% 控制量约束：
A_t=zeros(Nc,Nc);%工具人矩阵
for i=1:1:Nc
    for j=1:1:Nc
        if i>=j
            A_t(i,j)=1;
        else
            A_t(i,j)=0;
        end
    end
end
A_I=kron(A_t,eye(Nu));
% 控制量升维度
Ut=kron(ones(Nc,1),U);
% 转角限制
umax=0.44;
umin=-0.44;
% 转角变化率限制
umax_dt=0.05;
umin_dt=-0.05;
% 两者升维度
Umax=kron(ones(Nc,1),umax);
Umin=kron(ones(Nc,1),umin);
Umax_dt=kron(ones(Nc,1),umax_dt);
Umin_dt=kron(ones(Nc,1),umin_dt);
% 不等式子AX<=b
A_cons_cell=cell(2,2);
A_cons_cell{1,1}=A_I;
A_cons_cell{1,2}=zeros(Nu*Nc,1);
A_cons_cell{2,1}=-A_I;
A_cons_cell{2,2}=zeros(Nu*Nc,1);
A_cons=cell2mat(A_cons_cell);

B_cons_cell=cell(2,1);
B_cons_cell{1,1}=Umax-Ut;
B_cons_cell{2,1}=-Umin+Ut;
B_cons=cell2mat(B_cons_cell);
% 不等式lb<x<Ub
lb=[Umin_dt];
ub=[Umax_dt];


function sys=mdlUpdate(t,x,u)
warning off all
sys = x;

function sys=mdlOutputs(t,x,u)
global Nx Nu Np Nc Row T U kesi PHI Q THETA H
global A_cons B_cons lb ub x_start solve_times
warning off all
%% 定义全局变量
% global A1 B1 u_piao U kesi
% global solve_times
% Nx=4;%状态量4个，横向误差，航向角误差，侧向速度误差，横摆角速度误差
% Nu=1;%控制量1个前轮转角
% Np=100;%预测步长50
% Nc=80;%控制步长30
% Row=10;%松弛因子
% T=0.01;%采用时间0.01
% %%车辆参数
% m=1274+71*2;%整车质量
% a=1.015;%质心到前轴的距离
% b=1.895;%质心到后轴的距离
% Iz=1536.7;%绕Z轴的转动惯量
% % k1=-112600;%前轴侧偏刚度,注意，这边的侧偏刚度已经乘以2；
% % k2=-89500;%后轴侧偏刚度 
% %通过拟合和魔术公式得到的精确版（可以写一下）；
% k1=-148970;
% k2=-82204;
% 
% %%控制器设计
% vx=u(5);%%%%%这边多设计几个传入值，k1,k2,Mu;
% 
% 
% kesi=zeros(Nx+Nu,1);
% kesi(1)=u(1);
% kesi(2)=u(2);
% kesi(3)=u(3);
% kesi(4)=u(4);
% kesi(5)=U;
% %可以进一步调节这里的q和R
% % q=[100,0,0,0;
% %    0,1,0,0;
% %    0,0,1,0;
% %    0,0,0,1];
% q=[10,0,0,0;
%    0,1,0,0;
%    0,0,1,0;
%    0,0,0,1];
% Q=kron(eye(Np),q);
% % R=eye(Nc*Nu);
% R=1000*eye(Nc*Nu);
% 
% %雅可比矩阵获得线性化汽车动力学式子
% A2=[0,vx,1,0;
%    0,0,0,1;
%    0,0,(k1+k2)/m/vx,(a*k1-b*k2)/m/vx-vx;
%    0,0,(a*k1-b*k2)/Iz/vx,(a^2*k1+b^2*k2)/Iz/vx];
% B2=[0;0;-k1/m;-a*k1/Iz];
% 
% %前向欧拉法获得
% A1=A2*T+eye(Nx);
% B1=B2*T;
% 
% %针对变量kesi，重新设计线性系统
% A_cell=cell(2,2);
% A_cell{1,1}=A1;
% A_cell{1,2}=B1;
% A_cell{2,1}=zeros(Nu,Nx);
% A_cell{2,2}=eye(Nu);
% A=cell2mat(A_cell);
% 
% B_cell=cell(2,1);
% B_cell{1,1}=B1;
% B_cell{2,1}=eye(Nu);
% B=cell2mat(B_cell);
% 
% C_cell=cell(1,2);
% C_cell{1,1}=eye(Nx);
% C_cell{1,2}=zeros(Nx,Nu);
% C=cell2mat(C_cell);
% 
% 
% %% 预测时域矩阵升维
% PHI_cell=cell(Np,1);
% for i=1:1:Np
%     PHI_cell{i,1}=C*A^i;
% end
% PHI=cell2mat(PHI_cell);
% 
% THETA_cell=cell(Np,Nc);
% for i=1:1:Np
%     for j=1:1:Nc
%         if  i>=j
%             THETA_cell{i,j}=C*A^(i-j)*B;
%         else
%             THETA_cell{i,j}=zeros(Nx,Nu);
%         end
%     end
% end
% THETA=cell2mat(THETA_cell);
% % 获得H和f
% H_cell=cell(2,2);
% H_cell{1,1}=THETA'*Q*THETA+R;
% H_cell{1,2}=zeros(Nc*Nu,1);
% H_cell{2,1}=zeros(1,Nc*Nu);
% H_cell{2,2}=Row; %松弛因子
% H=cell2mat(H_cell);
% error=PHI*kesi;
% 
% f_cell=cell(2,1);
% f_cell{1,1}=2*error'*Q*THETA;
% f_cell{1,2}=0;
% f=cell2mat(f_cell);
% % 约束条件：
% % 控制量约束：
% A_t=zeros(Nc,Nc);%工具人矩阵
% for i=1:1:Nc
%     for j=1:1:Nc
%         if i>=j
%             A_t(i,j)=1;
%         else
%             A_t(i,j)=0;
%         end
%     end
% end
% A_I=kron(A_t,eye(Nu));
% % 控制量升维度
% Ut=kron(ones(Nc,1),U);
% % 转角限制
% umax=0.44;
% umin=-0.44;
% % 转角变化率限制
% umax_dt=0.05;
% umin_dt=-0.05;
% % 两者升维度
% Umax=kron(ones(Nc,1),umax);
% Umin=kron(ones(Nc,1),umin);
% Umax_dt=kron(ones(Nc,1),umax_dt);
% Umin_dt=kron(ones(Nc,1),umin_dt);
% % 不等式子AX<=b
% A_cons_cell=cell(2,2);
% A_cons_cell{1,1}=A_I;
% A_cons_cell{1,2}=zeros(Nu*Nc,1);
% A_cons_cell{2,1}=-A_I;
% A_cons_cell{2,2}=zeros(Nu*Nc,1);
% A_cons=cell2mat(A_cons_cell);
% 
% B_cons_cell=cell(2,1);
% B_cons_cell{1,1}=Umax-Ut;
% B_cons_cell{2,1}=-Umin+Ut;
% B_cons=cell2mat(B_cons_cell);
% % 不等式lb<x<Ub
% lb=[Umin_dt];
% ub=[Umax_dt];

kesi=zeros(Nx+Nu,1);
kesi(1)=u(1);
kesi(2)=u(2);
kesi(3)=u(3);
kesi(4)=u(4);
kesi(5)=U;



error=PHI*kesi;


f_cell=cell(2,1);


f_cell{1,1}=2*error'*Q*THETA;

f_cell{1,2}=0;
f=cell2mat(f_cell);


% 二次规划问题
x_start=zeros(Nc+1,1);%加入一个起始点
options=optimset('Algorithm','interior-point-convex');

% 测量求解时间
tic;
[X, ~, exitflag] = quadprog(H, f, A_cons, B_cons, [], [], lb, ub, x_start, options);
solve_time = toc;


% 记录求解时间
solve_times(end+1) = solve_time;
if mod(length(solve_times), 100) == 0 % 每100次输出一次平均时间
    fprintf('Quadprog - Iterations: %d, Avg Solve Time: %.6f s\n', length(solve_times), mean(solve_times));
end


% 赋值输出
u_piao=X(1);

% disp(['Time: ', num2str(t), ', y: ', num2str(u1)]);
% disp(['Time: ', num2str(t), ', vx: ', num2str(u_piao)]);

U=kesi(5)+u_piao;


u_real=U;


sys = u_real;

% 
% 
% 
% 
% function [sys,x0,str,ts] = MPC_Final(t,x,u,flag)
% % MPC_Final S-function adapted for TinyMPC, preserving original parameters
% switch flag
%   case 0
%     [sys,x0,str,ts] = mdlInitializeSizes;
%   case 2
%     sys = mdlUpdate(t,x,u);
%   case 3
%     sys = mdlOutputs(t,x,u);
%   case {1,4,9}
%     sys = [];
%   otherwise
%     error(['Unhandled flag = ', num2str(flag)]);
% end
% 
% function [sys,x0,str,ts] = mdlInitializeSizes
% sizes = simsizes;
% sizes.NumContStates  = 0;
% sizes.NumDiscStates  = 4;
% sizes.NumOutputs     = 1;
% sizes.NumInputs      = 5;
% sizes.DirFeedthrough = 1;
% sizes.NumSampleTimes = 1;
% sys = simsizes(sizes);
% x0  = [0.0001, 0.0001, 0.0001, 0.0001]';
% str = [];
% ts  = [0.01 0];
% global U solve_times
% U = [0]; % Initial control input
% solve_times = []
% function sys = mdlUpdate(t,x,u)
% sys = x;
% 
% function sys = mdlOutputs(t,x,u)
% global prob U
% global solve_times
% %判断是否运行了setup程序
% if isempty(prob)
%     error('TinyMPC problem not initialized. Run setup_tinympc.m first.');
% end
% 
% 
% 
% % Current states and velocity
% vx = u(5);
% x = u(1:4); % [lateral error, heading error, lateral vel error, yaw rate error]
% 
% % 创建新状态空间矩阵Augmented state: kesi = [x; U]
% Nx = 4;
% Nu = 1;
% kesi = zeros(Nx + Nu, 1);
% kesi(1) = u(1);
% kesi(2) = u(2);
% kesi(3) = u(3);
% kesi(4) = u(4);
% kesi(5) = U;
% 
% % Set initial state in TinyMPC
% verbose_int = int32(0);
% prob.set_x0(kesi, verbose_int);
% 
% % 测量求解时间
% tic;
% prob.solve(verbose_int);
% solve_time = toc;
% 
% % 记录求解时间
% solve_times(end+1) = solve_time;
% if mod(length(solve_times), 100) == 0 % 每100次输出一次平均时间
%     fprintf('TinyMPC - Iterations: %d, Avg Solve Time: %.6f s\n', length(solve_times), mean(solve_times));
% end
% 
% 
% % Get control increment sequence
% delta_u = prob.get_u(verbose_int);
% 
% if isempty(delta_u) || length(delta_u) < 1
%     error('TinyMPC returned an empty or invalid delta_u');
% end
% 
% % Use first control increment (equivalent to u_piao = X(1))
% u_piao = double(delta_u(1));
% 
% % Update control input (U = kesi(5) + u_piao)
% U = kesi(5) + u_piao;
% 
% % % Enforce control input constraint
% % u_real = min(max(U, -0.44), 0.44);
% u_real = U;
% 
% % Debugging output
% disp(['Time: ', num2str(t), ...
%       ', States: [', num2str(x(1)), ', ', num2str(x(2)), ', ', num2str(x(3)), ', ', num2str(x(4)), ']', ...
%       ', vx: ', num2str(vx), ...
%       ', Control: ', num2str(u_real), ...
%       ', Delta U: ', num2str(delta_u(1,10))]);
% 
% % Output
% sys = [u_real];