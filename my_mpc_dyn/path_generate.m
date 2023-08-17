  clc;
  clear;
  N=400;
  T=0.05;

  shape=2.4;%参数名称，用于参考轨迹生成
  dx1=25;dx2=21.95;%没有任何实际意义，只是参数名称
  dy1=4.05;dy2=5.7;%没有任何实际意义，只是参数名称
  Xs1=27.19;Xs2=56.46;%参数名称
  phi_ref=zeros(N,1);%用于保存预测时域内的期望轨迹
  Y_ref=zeros(N,1);%用于保存预测时域内的期望轨迹
  X_phi=0.05:0.5:200;
  X_ref=X_phi';

for p=1:1:N
    %%双移线追踪
      z1=shape/dx1*(X_ref(p,1)-Xs1)-shape/2;
      z2=shape/dx2*(X_ref(p,1)-Xs2)-shape/2;
      Y_ref(p,1)=dy1/2*(1+tanh(z1))-dy2/2*(1+tanh(z2));
      phi_ref(p,1)=atan(dy1*(1/cosh(z1))^2*(1.2/dx1)-dy2*(1/cosh(z2))^2*(1.2/dx2));
      Yita_ref_cell{p,1}=[phi_ref(p,1);Y_ref(p,1)];
end

Phi_ref=phi_ref;

X_dot_ref=zeros(N,1);
Y_dot_ref=zeros(N,1);
Phi_dot_ref=zeros(N,1);
for i=2:N-1
    X_dot_ref(i)=(X_ref(i+1)-X_ref(i))/T;
    Y_dot_ref(i)=(Y_ref(i+1)-Y_ref(i))/T;
    Phi_dot_ref(i)=(Phi_ref(i+1)-Phi_ref(i))/T;
end
X_dot_ref(1)=X_dot_ref(2);
Y_dot_ref(1)=Y_dot_ref(2);
Phi_dot_ref(1)=Phi_dot_ref(2);
X_dot_ref(N)=X_dot_ref(N-1);
Y_dot_ref(N)=Y_dot_ref(N-1);
Phi_dot_ref(N)=Phi_dot_ref(N-1);

Ax_ref=zeros(N,1);
Ay_ref=zeros(N,1);
for i=2:N-1
    Ax_ref(i)=(X_dot_ref(i+1)-X_dot_ref(i))/T;
    Ay_ref(i)=(Y_dot_ref(i+1)-Y_dot_ref(i))/T;
end
Ax_ref(1)=Ax_ref(2);
Ay_ref(1)=Ay_ref(2);
Ax_ref(N)=Ax_ref(N-1);
Ay_ref(N)=Ay_ref(N-1);

L=2.6;
Detal_ref=zeros(N,1);
for i=1:N
    Vr=sqrt(X_dot_ref(i)^2+Y_dot_ref(i)^2);
    Detal_ref(i)=atan(L*Phi_dot_ref(i)/Vr);
end

clearvars -except X_ref Y_ref Phi_ref X_dot_ref Y_dot_ref Phi_dot_ref Detal_ref;
% plot(X_ref, Y_ref,'b--','LineWidth',1);

