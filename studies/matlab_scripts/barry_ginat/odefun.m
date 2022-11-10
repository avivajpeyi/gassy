
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
r = 0.7*R;
T = sqrt(4*pi^2*r^3/(G*(M+m)));


global rho;
global C;
global M_e;

r = sqrt(y(1)^2 + y(3)^2);
v = sqrt(y(2)^2 + y(4)^2);

b_90 = max([G*M/(v*r/T)^2 1e-1*Rsol]);
N = R/b_90;


c_s = C(r*r);
Mach = (v*r/T)/c_s;
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
f = I*4*pi*G^2*M^2*rho(r*r)/((r*v/T)^3);
if (I<0)
    disp('boo');
end


rdot = (y(1)*y(2) + y(3)*y(4))/r;
rdot = rdot*r/T;
nu = mu/(M+m);

v = r*v/T;
r = r*r;

A = 1/c^2*(-3*rdot^2*nu/2 + v^2 + 3*nu*v^2 - G*(M+m)*(4+2*nu)/r) + ...
    1/c^4*(15*rdot^4*nu/8 - 45*rdot^4*nu^2/8 - 9*rdot^2*nu*v^2/2 + 6*rdot^2*nu^2*v^2 + 3*nu*v^4 - 4*nu^2*v^4 +...
   G*(M+m)/r*(-2*rdot^2 - 25*rdot^2*nu -2*rdot^2*nu^2 - 13*nu*v^2/2 + 2*nu^2*v^2) + ...
   (G*(M+m)/r)^2*(9 + 87*nu/4) ) + ...
    1/c^5*(-24*nu*rdot*v^2*G*(M+m)/(5*r) - (136*rdot*nu/15)*(G*(M+m)/r)^2);

B = 1/c^2*(-4*rdot + 2*rdot*nu) + ...
    1/c^4*(9*rdot^3*nu/2 + 3*rdot^3*nu^2 -15*rdot*nu*v^2/2 - 2*rdot*nu^2*v^2 + G*(M+m)/r*(2*rdot + 41*rdot*nu/2 + 4*rdot*nu^2)) + ...
    1/c^5*((8*nu*v^2/5)*G*(M+m)/r + 24*nu/5*(G*(M+m)/r)^2);

r = r/r;

dydt(1) = y(2);
dydt(2) = -4*pi^2*((1+A)*y(1)/r + B*y(2)*r/T)/r^2 - T*f*y(2)/M -4*pi^2*((M_e(r*r)-m)/(M+m))*y(1)/r^3;
dydt(3) = y(4);
dydt(4) = -4*pi^2*((1+A)*y(3)/r + B*y(4)*r/T)/r^2 - T*f*y(4)/M -4*pi^2*((M_e(r*r)-m)/(M+m))*y(3)/r^3;



Rt = 1e-1*Rsol*max([(M/m)^(1/3) (m/M)^(1/3)]);


if (r*r <= Rt)
    dydt(1) = 0;
    dydt(2) = 0;
    dydt(3) = 0;
    dydt(4) = 0;
    disp(t*T);
end

end