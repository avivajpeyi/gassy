% MAIN ORBIT EVOLVER 


%constants in cgs
G = 6.67e-8;
c = 29979245800;
Rsol = 696342*1e5;
Msol = 1.98855*1e33;
AU = 14959787070000;

M=15*Msol;
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
Tend = 1e9; % full spiral
Tend2 = 5e6; % few orbits
Tend3=1.5e6; % few steps

%solution of the equations of motion
[X,Y,v_x,v_y,t] = spiralling113MESA15(M,m,a,e,dt,Tend);%e = 0 always
out = struct('X', X, 'Y', Y, 'v_x', v_x, 'v_y', v_y, 't', t);
% save('spiral_out', 'out');


plot(X, Y)
title('1M object sprilaing into a 15M object')
xlabel('X [units]') 
ylabel('Y [units]') 
saveas(gcf,'orbit_new','png')




