%% FUNCTION TO EVOLVE 2 BODIES

function [X,Y,v_x,v_y,t] = spiralling113MESA15(M,m,r,e,~,Tend)
%Solves the equations of motion for r two-body system in r common envelope
%with polytropic index p, i.e. P ~ rho^(1+1/p).
G = 6.67e-8;
Rsol = 696342*1e5;
Msol = 1.98855*1e33;

% -------------------------------------------------------------------------

%Dynamical friction law


global rho;
global C;
global M_e;



load profile15Msol.mat rho q c_s;

Rho = rho;
q = q*Rsol; c_s = c_s*100;
R = max(q); q = smooth(q);


%Smoothing

rho = @(r) (r<=R).*real(interp1(q,Rho,r,'linear','extrap')) + (r > R).*0;

C = @(r) abs(interp1(q,c_s,r,'linear','extrap'));

M_e = @(x) (x<=R).*integral(@(r) 4*pi*r.^2.*rho(r),0,x) + (x > R).*integral(@(r) 4*pi*r.^2.*rho(r),0,R);


% -------------------------------------------------------------------------

%initial conditions
T2 = 4*pi^2*r^3/(G*(M+m));
T = sqrt(T2);
L = sqrt(G*(M+m+M_e(r))*r*(1-e^2));
y0 = [(1+e) 0 0 L*T/(r^2*(1+e))];

%ode113 stuff

tspan = [0 Tend]/T;
opts = odeset('RelTol',1e-10,'Stats','on');


%Integrator
[t,y] = ode113(@odefun,tspan,y0,opts);


X = r*y(:,1);
Y = r*y(:,3);
v_x = r*y(:,2)/T;
v_y = r*y(:,4)/T;
t = t*T;
end

