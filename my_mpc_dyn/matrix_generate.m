%% �˳����ܣ����ݼ򻯶���ѧģ��(����С�Ƕȼ�����)�����ſ˱Ⱦ���
% �汾��V1.0����дʱ��2013.12.11
% ��������ſ˱Ⱦ������복�����в���������صģ������仯�ˣ��ſ˱Ⱦ���Ҳ����Ӧ�仯��
% ������ֻ��һ������ֵ���������б�Ҫ������������
clc
clear all;
%% ����Ϊ����
%������������ 
syms x_dot y_dot phi phi_dot Y X;%����״̬��
syms delta_f  %ǰ��ƫ��,������
%syms sf sr;%�ֱ�Ϊǰ���ֵĻ�����,��Ҫ�ṩ
Sf=0.2; Sr=0.2;
%syms a b;%ǰ���־��복�����ĵľ��룬�������в���
a=1.232;b=1.468;
%syms C_cf C_cr C_lf C_lr;%�ֱ�Ϊǰ���ֵ��ݺ����ƫ�նȣ��������в���
Ccf=66900;Ccr=62700;Clf=66900;Clr=62700;
%syms m g I;%mΪ����������gΪ�������ٶȣ�IΪ������Z���ת���������������в���
m=1723;g=9.8;I=4175;

% ��������ѧģ��
dy_dot=-x_dot*phi_dot+2*(Ccf*(delta_f-(y_dot+a*phi_dot)/x_dot)+Ccr*(b*phi_dot-y_dot)/x_dot)/m;
dx_dot=y_dot*phi_dot+2*(Clf*Sf+Clr*Sr+Ccf*delta_f*(delta_f-(y_dot+phi_dot*a)/x_dot))/m;
%dphi_dot=dphi_dot;
dphi_dot=(2*a*Ccf*(delta_f-(y_dot+a*phi_dot)/x_dot)-2*b*Ccr*(b*phi_dot-y_dot)/x_dot)/I;
Y_dot=x_dot*sin(phi)+y_dot*cos(phi);
X_dot=x_dot*cos(phi)-y_dot*sin(phi);

% �ſ˱Ⱦ������
% f=[dy_dot;dx_dot;phi_dot;dphi_dot;Y_dot;X_dot];%����ѧģ��
% kesi=[y_dot,x_dot,phi,phi_dot,Y,X];%ϵͳ״̬��
f=[X_dot;Y_dot;phi_dot;dx_dot;dy_dot;dphi_dot];%����ѧģ��
kesi=[X,Y,phi,x_dot,y_dot,phi_dot];%ϵͳ״̬��
v=delta_f;
R=jacobian(f,kesi);%����A(t)-����
R2=jacobian(f,v);%����B(t)-����

% ���ƾ������(����������ת��Ϊ��ɢ���󣬲��ý����㷨  A=I+T*A(t),B=T*B(t))
I=eye(6);
syms T;
A=I+T*R;
B=T*R2;
A1=vpa(A,3);
B1=vpa(B,3);


