function [X,Y,v_x,v_y,t] = spiralling113MESA15(M,m,a,e,~,Tend)
%Solves the equations of motion for a two-body system in a common envelope
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
T2 = 4*pi^2*a^3/(G*(M+m));
T = sqrt(T2);
L = sqrt(G*(M+m+M_e(a))*a*(1-e^2));
y0 = [(1+e) 0 0 L*T/(a^2*(1+e))];

%ode113 stuff

tspan = [0 Tend]/T;
opts = odeset('RelTol',1e-10,'Stats','on');


%Integrator
[t,y] = ode113(@odefun,tspan,y0,opts);


X = a*y(:,1);
Y = a*y(:,3);
v_x = a*y(:,2)/T;
v_y = a*y(:,4)/T;
t = t*T;
end

function dydt = odefun(t,y)
dydt = zeros(1,4)';
%parameters
G = 6.67e-8;
c = 29979245800;
Msol = 1.98855*1e33;
Rsol = 696342*1e5;
R = 111.3642*Rsol;
M=15*Msol;
% M = 1.4*Msol;
m = 1*Msol;
mu = M*m/(M+m);
a = 0.7*R;
T = sqrt(4*pi^2*a^3/(G*(M+m)));


global rho;
global C;
global M_e;

r = sqrt(y(1)^2 + y(3)^2);
v = sqrt(y(2)^2 + y(4)^2);

b_90 = max([G*M/(v*a/T)^2 1e-1*Rsol]);
N = R/b_90;

c_s = C(a*r);
Mach = (v*a/T)/c_s;
% Note f is cut at 1 so as not to diverge
h1 = max([1/N 1e-2]);

h = (1-(2-h1)*exp(-2+2*h1)/(N^2*h1))^(-1/2) - 1;

%Ostriker
if (Mach < 1-h1)
    I = 0.5*log((1+Mach)/(1-Mach)) - Mach;
else
    if((Mach >= 1-h1) && (Mach < 1 + h))
        I = 0.5*log((2-h1)/h1) - 1 + h1;
%         disp('Speed of sound reached');
    else
        I = 0.5*log(1-1/Mach^2) + log(N);
        if (Mach < 1)
            disp('Wrong condition');
        end
    end
end
f = I*4*pi*G^2*M^2*rho(a*r)/((a*v/T)^3);
if (I<0)
    disp('boo');
end



rdot = (y(1)*y(2) + y(3)*y(4))/r;
rdot = rdot*a/T;
nu = mu/(M+m);

v = a*v/T;
r = a*r;

A = 1/c^2*(-3*rdot^2*nu/2 + v^2 + 3*nu*v^2 - G*(M+m)*(4+2*nu)/r) + ...
    1/c^4*(15*rdot^4*nu/8 - 45*rdot^4*nu^2/8 - 9*rdot^2*nu*v^2/2 + 6*rdot^2*nu^2*v^2 + 3*nu*v^4 - 4*nu^2*v^4 +...
   G*(M+m)/r*(-2*rdot^2 - 25*rdot^2*nu -2*rdot^2*nu^2 - 13*nu*v^2/2 + 2*nu^2*v^2) + ...
   (G*(M+m)/r)^2*(9 + 87*nu/4) ) + ...
    1/c^5*(-24*nu*rdot*v^2*G*(M+m)/(5*r) - (136*rdot*nu/15)*(G*(M+m)/r)^2);

B = 1/c^2*(-4*rdot + 2*rdot*nu) + ...
    1/c^4*(9*rdot^3*nu/2 + 3*rdot^3*nu^2 -15*rdot*nu*v^2/2 - 2*rdot*nu^2*v^2 + G*(M+m)/r*(2*rdot + 41*rdot*nu/2 + 4*rdot*nu^2)) + ...
    1/c^5*((8*nu*v^2/5)*G*(M+m)/r + 24*nu/5*(G*(M+m)/r)^2);

r = r/a;

dydt(1) = y(2);
dydt(2) = -4*pi^2*((1+A)*y(1)/r + B*y(2)*a/T)/r^2 - T*f*y(2)/M -4*pi^2*((M_e(a*r)-m)/(M+m))*y(1)/r^3;
dydt(3) = y(4);
dydt(4) = -4*pi^2*((1+A)*y(3)/r + B*y(4)*a/T)/r^2 - T*f*y(4)/M -4*pi^2*((M_e(a*r)-m)/(M+m))*y(3)/r^3;



Rt = 1e-1*Rsol*max([(M/m)^(1/3) (m/M)^(1/3)]);


if (r*a <= Rt)
    dydt(1) = 0;
    dydt(2) = 0;
    dydt(3) = 0;
    dydt(4) = 0;
    disp(t*T);
end

end
