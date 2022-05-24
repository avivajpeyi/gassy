% Run only after back_reaction_of_gravy_waves. This plots the wave-form as
% a function of time. 


kpc = 3.086e21; 
% D = 50000*kpc;
D = 10*kpc;
c = 29979245800;
G = 6.67e-8;
Me = 15*Msol;
rg = 2*G*Me/c^2;
mu = M*m/(M+m);
phi0=0.7

x = X; y = Y;

% second derivative through equations of motion




r = sqrt(x.^2 + y.^2);
rdot = (x.*v_x + y.*v_y)./r;
v = sqrt(v_x.^2 + v_y.^2);
nu = mu/(M+m);

load profile15Msol.mat rho q c_s;
Rho = rho;
q = q*Rsol; c_s = c_s*100;
R = max(q); q = smooth(q);


%Smoothing

rho = @(r) (r<=R).*real(interp1(q,Rho,r,'linear','extrap')) + (r > R).*0;
C = @(r) abs(interp1(q,c_s,r,'linear','extrap'));

M_e = @(x) (x<=R).*integral(@(r) 4*pi*r.^2.*rho(r),0,x) + (x > R).*integral(@(r) 4*pi*r.^2.*rho(r),0,R);




c_s = real(C(r));
Mach = v./c_s;
b_90 = v;
len = length(v);
for j=(1:len)
    b_90(j) = max([G*M/v(j)^2 1e-1*Rsol]);
end
N = R./b_90;
I = zeros(1,len);
for j=(1:len)
    if (Mach(j) < 1)
        I(j) = 0.5*log((1+Mach(j))/(1-Mach(j))) - Mach(j);   
    else
        if(Mach(j) == 1)
                disp('Exactly speed of sound reached');
                return;
        else
                I(j) = 0.5*log(1-1/Mach(j)^2) + log(N(j));   
        end
    end
end
f = 4*pi*G^2*M^2*I'.*rho(r)./(v.^3);

A = 1/c^2*(-3*rdot.^2*nu/2 + v.^2 + 3*nu*v.^2 - G*(M+m)*(4+2*nu)./r) + ...
    1/c^4*(15*rdot.^4*nu/8 - 45*rdot.^4*nu^2/8 - 9*rdot.^2*nu.*v.^2/2 + 6*rdot.^2*nu^2.*v.^2 + 3*nu*v.^4 - 4*nu^2*v.^4 +...
   G*(M+m)./r.*(-2*rdot.^2 - 25*rdot.^2*nu -2*rdot.^2*nu^2 - 13*nu*v.^2/2 + 2*nu^2*v.^2) + ...
   (G*(M+m)./r).^2*(9 + 87*nu/4) ) + ...
    1/c^5*(-24*nu*rdot.*v.^2*G*(M+m)./(5*r) - (136*rdot*nu/15).*(G*(M+m)./r).^2);

B = 1/c^2*(-4*rdot + 2*rdot*nu) + ...
    1/c^4*(9*rdot.^3*nu/2 + 3*rdot.^3*nu^2 -15*rdot*nu.*v.^2/2 - 2*rdot*nu^2.*v.^2 + G*(M+m)./r.*(2*rdot + 41*rdot*nu/2 + 4*rdot*nu^2)) + ...
    1/c^5*((8*nu*v.^2/5)*G*(M+m)./r + 24*nu/5*(G*(M+m)./r).^2);
M_ee = r;
% potent = r;
for j=(1:len)
    M_ee(j) = M_e(r(j));
%     potent(j) = -4*pi*G*(1/r(j)*integral(@(x) rho(x).*x.^2,0,r(j)) + integral(@(x) rho(x).*x,r(j),R));
end
mu = M*Me/(M+Me);
x_ddot = -G*(M+m)*((1+A).*x./r + B.*v_x)./r.^2 - f.*v_x/M -G*(M_ee-m).*x./r.^3;
y_ddot = -G*(M+m)*((1+A).*y./r + B.*v_y)./r.^2 - f.*v_y/M -G*(M_ee-m).*y./r.^3;

M_ddot11 = 2*mu*(x.*x_ddot + v_x.^2); M_ddot12 = mu*(x.*y_ddot + 2*v_x.*v_y + y.*x_ddot); M_ddot22 = 2*mu*(y.*y_ddot + v_y.^2);

% M_11 = mu*x.^2; M_12 = mu*x.*y; M_22 = mu*y.^2;
phi = phi0;
theta = i;
% M_11dot = my_diff(M_11,order)./my_diff(t,order); M_12dot = my_diff(M_12,order)./my_diff(t,order);M_22dot = my_diff(M_22,order)./my_diff(t,order);
% M_ddot11 = my_diff(M_11dot,order)./my_diff(t(1:end-order),order); M_ddot12 = my_diff(M_12dot,order)./my_diff(t(1:end-order),order);M_ddot22 = my_diff(M_22dot,order)./my_diff(t(1:end-order),order);
% phi(end - 2*order + 1:end) = [];

h_plus = G/(D*c^4)*( M_ddot11.*((cos(phi)).^2 - (sin(phi)).^2*(cos(i))^2) + M_ddot22.*((sin(phi)).^2 - (cos(phi)).^2*(cos(i))^2) - M_ddot12.*(sin(2*phi))*(1+(cos(i))^2) );
h_cross =  G/(D*c^4)*( ( M_ddot11 - M_ddot22).*sin(2*phi)*cos(i) + 2*M_ddot12.*cos(2*phi)*cos(i));

h_plus_av = G/(D*c^4)*( M_ddot11.*(pi - pi^2/2) + M_ddot22.*(pi - pi^2/2));
