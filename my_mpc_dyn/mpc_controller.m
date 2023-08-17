function [sys,x0,str,ts] = MY_MPCController(t,x,u,flag)

switch flag,
 case 0
  [sys,x0,str,ts] = mdlInitializeSizes;
 case 2
  sys = mdlUpdates(t,x,u);
 case 3
  sys = mdlOutputs(t,x,u);
 case {1,4,9}
  sys = [];
 otherwise
  error(['unhandled flag = ',num2str(flag)]); % Error handling
end

function [sys,x0,str,ts] = mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 3; % this parameter doesn't matter
sizes.NumOutputs     = 1;
sizes.NumInputs      = 13;
sizes.DirFeedthrough = 1; % Matrix D is non-empty.
sizes.NumSampleTimes = 1;
sys = simsizes(sizes); 
x0 =[0;0;0];   
% Initialize the discrete states.
str = [];             % Set str to an empty matrix.
global T;
global car_steer;
car_steer=0;
T=0.05;
ts  = [T 0];       % sample time: [period, offset]

function sys = mdlUpdates(t,x,u)
sys = x;

function sys = mdlOutputs(t,x,u)
tic
% 参考变量处理
Np=30;
Nx=6;
Nu=1;
[x_ref,y_ref,phi_ref]=deal(u(1),u(2),u(3));
[x_dot_ref,y_dot_ref,phi_dot_ref]=deal(u(4),u(5),u(6));
delta_ref=u(7);
if abs(x_dot_ref)<0.001
    x_dot_ref=0.001*abs(x_dot_ref+0.01)/(x_dot_ref+0.01);
end
if abs(y_dot_ref)<0.001
    y_dot_ref=0.001*abs(y_dot_ref+0.01)/(y_dot_ref+0.01);
end
if abs(phi_dot_ref)<0.001
    phi_dot_ref=0.001*abs(phi_dot_ref+0.01)/(phi_dot_ref+0.01);
end
% 车辆状态初始化
global T;
[car_x,car_y,car_phi]=deal(u(8),u(9),u(10));
[car_x_dot,car_y_dot,car_phi_dot]=deal(u(11),u(12),u(13));
global car_steer;
% 构造矩阵
Q=10*eye(Nx*Np);
R=100*eye(Nu*Np);
A=zeros(Nx,Nx);
B=zeros(Nx,Nu);
if 1
    a11=1;
    a12=0;
    a13=-1.0*T*(y_dot_ref*cos(phi_ref) + x_dot_ref*sin(phi_ref));
    a14=T*cos(phi_ref);
    a15=-1.0*T*sin(phi_ref);
    a16=0;
    a21=0;
    a22=1;
    a23=T*(x_dot_ref*cos(phi_ref) - 1.0*y_dot_ref*sin(phi_ref));
    a24=T*sin(phi_ref);
    a25=T*cos(phi_ref);
    a26=0;
    a31=0;
    a32=0;
    a33=1;
    a34=0;
    a35=0;
    a36=T;
    a41=0;
    a42=0;
    a43=0;
    a44=(77.65*T*delta_ref*(1.23*phi_dot_ref + y_dot_ref))/x_dot_ref^2 + 1.0;
    a45=T*(phi_dot_ref - (77.65*delta_ref)/x_dot_ref);
    a46=T*(y_dot_ref - (95.67*delta_ref)/x_dot_ref);
    a51=0;
    a52=0;
    a53=0;
    a54=-1.0*T*(phi_dot_ref + (0.00116*(92043.6*phi_dot_ref - 62700.0*y_dot_ref))/x_dot_ref^2 - (77.65*(1.23*phi_dot_ref + y_dot_ref))/x_dot_ref^2);
    a55=1.0 - (150.43*T)/x_dot_ref;
    a56=-1.0*T*(x_dot_ref - 11.16/x_dot_ref);
    a61=0;
    a62=0;
    a63=0;
    a64=T*((39.48*(1.23*phi_dot_ref + y_dot_ref))/x_dot_ref^2 + (0.000239*(270240*phi_dot_ref - 184087*y_dot_ref))/x_dot_ref^2);
    a65=(4.60*T)/x_dot_ref;
    a66=1.0 - (113.37*T)/x_dot_ref;
    A=[[a11,a12,a13,a14,a15,a16]
        [a21,a22,a23,a24,a25,a26]
        [a31,a32,a33,a34,a35,a36]
        [a41,a42,a43,a44,a45,a46]
        [a51,a52,a53,a54,a55,a56]
        [a61,a62,a63,a64,a65,a66]];
end
if 1
    b11=0;
    b21=0;
    b31=0;
    b41=T*((267600*delta_ref)/1723 - (133800*((154*phi_dot_ref)/125 + y_dot_ref))/(1723*x_dot_ref));
    b51=77.65*T;
    b61=39.48*T;
    B=[b11;b21;b31;b41;b51;b61];
end
if 1
    PHI_cell=cell(Np,1);
    THETA_cell=cell(Np,Np);
    for i=1:Np
        PHI_cell{i,1}=A^i;
        for j=1:Np
            if j<=i
                THETA_cell{i,j}=(A^(i-j))*B;
            else
                THETA_cell{i,j}=zeros(Nx,Nu);
            end
        end
    end
    PHI=cell2mat(PHI_cell);
    THETA=cell2mat(THETA_cell);
    H=2*(THETA'*Q*THETA+R);
    [x_e,y_e,phi_e]=deal(car_x-x_ref,car_y-y_ref,car_phi-phi_ref);
    [x_dot_e,y_dot_e,phi_dot_e]=deal(car_x_dot-x_dot_ref,car_y_dot-y_dot_ref,car_phi_dot-phi_dot_ref);
    e=[x_e;y_e;phi_e;x_dot_e;y_dot_e;phi_dot_e];
    E=PHI*e;
    f=2*E'*Q*THETA;
end
if 1
    A_t=zeros(Np,Np);
    for i=1:Np
        for j=1:Np
            if i<=j
                A_t(i,j)=1;
            else
                A_t(i,j)=0;
            end
        end
    end
    A_I=kron(A_t,eye(Nu));
    Ut=kron(ones(Np,1),car_steer);
    umin=[-0.436];
    umax=[0.436];
    delta_umin = [-0.0082];
    delta_umax = [ 0.0082];
    Umin=kron(ones(Np,1),umin);
    Umax=kron(ones(Np,1),umax);
    delta_Umin=kron(ones(Np,1),delta_umin);
    delta_Umax=kron(ones(Np,1),delta_umax);
    A_cons_cell={A_I;-A_I};
    b_cons_cell={Umax-Ut;Ut-Umin};
    A_cons=cell2mat(A_cons_cell);
    b_cons=cell2mat(b_cons_cell);
    lb=delta_Umin;
    ub=delta_Umax;
end

% 滚动优化
options=optimset('Algorithm','interior-point-convex');
Uk=quadprog(H,f,A_cons,b_cons,[],[],lb,ub,[],options);
% fprintf('uk1:6.3%f, delta_ref:6.3%f',Uk(1),delta_ref)
car_steer=Uk(1)+delta_ref;
% car_steer=delta_ref;
sys=car_steer;
toc
 
    
    
    