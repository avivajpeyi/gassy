
%constants in cgs
G = 6.67e-8;
c = 29979245800;
Rsol = 696342*1e5;
Msol = 1.98855*1e33;
AU = 14959787070000;

M=15*Msol
% parameters
R = 111.3642*Rsol; %UIS twice
% M = 1.4*Msol; %UIS companion
m = 1*Msol;%UIS core
% e = 0.9;
e=0;
Rs = G*(M+m)/c^2;
a = 0.7*R;% Update in Spiralling113... (UIS)
% a = 5e1*Rs;

%time-step
dt = 1e5;
Tend = 1e9;
Tend2 = 5e6;
Tend3=1e3

%solution of the equations of motion
[X,Y,v_x,v_y,t] = spiralling113MESA15(M,m,a,e,dt,Tend2);%e = 0 always
% ---------------------------------------------------------------------------
% [X2,Y2,v_x2,v_y2,t2] = spiralling113MESAsecond(M,m,sqrt(X(end)^2+Y(end)^2),X(end),Y(end),v_x(end),v_y(end),dt,Tend2);
% 
% L1 = length(t);
% L2 = length(t2);
% X_all = zeros(L1 + L2,1);
% Y_all = zeros(L1 + L2,1);
% v_x_all = zeros(L1 + L2,1);
% v_y_all = zeros(L1 + L2,1);
% X_all(1:L1) = X; Y_all(1:L1) = Y; v_x_all(1:L1) = v_x; v_y_all(1:L1) = v_y;
% X_all(L1+1:end) = X2; Y_all(L1+1:end) = Y2; v_x_all(L1+1:end) = v_x2; v_y_all(L1+1:end) = v_y2;
% X = X_all; Y = Y_all; v_x = v_x_all; v_y = v_y_all; t_all = zeros(L1+L2,1); t_all(1:L1) = t; t_all(L1+1:end) = t2; t = t_all;
% ---------------------------------------------------------------------------

