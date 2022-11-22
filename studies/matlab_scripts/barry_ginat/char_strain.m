global M;
Msol = 1.98855*1e33;
load lisa.mat f hf;
loglog(f,sqrt(f).*hf,'k');
hold on;
f_decigo_bbo = linspace(1e-3,5e1,1e3);
f_p = 7.36;
hf_decigo = sqrt(7.05e-48*(1+(f_decigo_bbo/f_p).^2) + 4.8e-51*f_decigo_bbo.^(-4).*(1+(f_decigo_bbo/f_p).^2).^(-1) + 5.33e-52*f_decigo_bbo.^(-4));
hf_bbo = sqrt(2e-49*(f_decigo_bbo).^2 +4.58e-49 + 1.26e-51*f_decigo_bbo.^(-4));
loglog(f_decigo_bbo,sqrt(f_decigo_bbo).*hf_decigo,'k');
loglog(f_decigo_bbo,sqrt(f_decigo_bbo).*hf_bbo,'k');
i = 0; %inclination relative to Earth
phi0 = 0;%orientation in x-y plain
indic = 0;

SN = zeros(6,1);
SN_dec = SN;
SN_bbo = SN;

%M = 0.6 White dwarf =====================================================
M = 0.6*Msol;

% 15
Dt = 1;
two_five_PN_15;
wave_form_15;
if (~all(diff(r)))
    j = find(diff(r) == 0,1,'first');
    t(j+1:end) = []; X(j+1:end) = [];Y(j+1:end) = [];v_x(j+1:end) = [];v_y(j+1:end) = [];
    wave_form_15;
end
%
% save 15-06.mat;
% load 15-06.mat;
% Dt = 1;

noise;
loglog(freq,2*freq.*sqrt(Sig)*(t(end)-t(n_start))/(0.4*L));
SN(1) = SNR;
SN_dec(1) = SNR_dec;
SN_bbo(1) = SNR_bbo;

% 8
two_five_PN_Hila;
wave_form_Hila;
if (~all(diff(r)))
    j = find(diff(r) == 0,1,'first');
    t(j+1:end) = []; X(j+1:end) = [];Y(j+1:end) = [];v_x(j+1:end) = [];v_y(j+1:end) = [];
    wave_form_Hila;
end
% save 8-06.mat;
% load 8-06.mat;
noise;
loglog(freq,2*freq.*sqrt(Sig)*(t(end)-t(n_start))/(0.4*L),'r');
SN(2) = SNR;
SN_dec(2) = SNR_dec;
SN_bbo(2) = SNR_bbo;

% 5
% Dt = 1;
two_five_PN_MESA5;
wave_form_5;
if (~all(diff(r)))
    j = find(diff(r) == 0,1,'first');
    t(j+1:end) = []; X(j+1:end) = [];Y(j+1:end) = [];v_x(j+1:end) = [];v_y(j+1:end) = [];
    wave_form_5;
end
% save 5-06.mat;
% load 5-06.mat;
% Dt = 0.1;
noise;
loglog(freq,2*freq.*sqrt(Sig)*(t(end)-t(n_start))/(0.4*L),'g');
SN(3) = SNR;
SN_dec(3) = SNR_dec;
SN_bbo(3) = SNR_bbo;

xlabel('f (Hz)');
ylabel('Characteristic Strain');
title('$m = 0.6~M_\odot$','Interpreter','Latex');
saveas(figure(1),'char_strain_06.fig');
saveas(figure(1),'characteristic_strain_06.eps','epsc');
hold off;
close(figure(1));

% M = 1.4 neutron star
loglog(f,sqrt(f).*hf,'k');
hold on;
loglog(f_decigo_bbo,sqrt(f_decigo_bbo).*hf_decigo,'k');
loglog(f_decigo_bbo,sqrt(f_decigo_bbo).*hf_bbo,'k');
% Dt = 1;
M = 1.4*Msol;
% 15
two_five_PN_15;
wave_form_15;
if (~all(diff(r)))
    j = find(diff(r) == 0,1,'first');
    t(j+1:end) = []; X(j+1:end) = [];Y(j+1:end) = [];v_x(j+1:end) = [];v_y(j+1:end) = [];
    wave_form_15;
end
% save 15-14.mat;
% load 15-14.mat
noise;
loglog(freq,2*freq.*sqrt(Sig)*(t(end)-t(n_start))/(0.4*L));
SN(4) = SNR;
SN_dec(4) = SNR_dec;
SN_bbo(4) = SNR_bbo;

% Dt = 1;
% 8
two_five_PN_Hila;
wave_form_Hila;
if (~all(diff(r)))
    j = find(diff(r) == 0,1,'first');
    t(j+1:end) = []; X(j+1:end) = [];Y(j+1:end) = [];v_x(j+1:end) = [];v_y(j+1:end) = [];
    wave_form_Hila;
end
% save 8-14.mat;
noise;
loglog(freq,2*freq.*sqrt(Sig)*(t(end)-t(n_start))/(0.4*L),'r');
SN(5) = SNR;
SN_dec(5) = SNR_dec;
SN_bbo(5) = SNR_bbo;

% 5
two_five_PN_MESA5;
wave_form_5;
if (~all(diff(r)))
    j = find(diff(r) == 0,1,'first');
    t(j+1:end) = []; X(j+1:end) = [];Y(j+1:end) = [];v_x(j+1:end) = [];v_y(j+1:end) = [];
    wave_form_5;
end
% save 5-14.mat;
% load 5-14.mat;
noise;
loglog(freq,2*freq.*sqrt(Sig)*(t(end)-t(n_start))/(0.4*L),'g');
SN(6) = SNR;
SN_dec(6) = SNR_dec;
SN_bbo(6) = SNR_bbo;

xlabel('f (Hz)');
ylabel('Characteristic Strain');
title('$m = 1.4~M_\odot$','Interpreter','Latex');
saveas(figure(1),'char_strain_14.fig');
saveas(figure(1),'characteristic_strain_14.eps','epsc');
hold off;
close(figure(1));

SN_puff = zeros(6,1);
SN_dec_puff = SN_puff;
SN_bbo_puff = SN_puff;

% Puffing up =============================================================
lambda = 10;
loglog(f,sqrt(f).*hf,'k');
hold on;
loglog(f_decigo_bbo,sqrt(f_decigo_bbo).*hf_decigo,'k');
loglog(f_decigo_bbo,sqrt(f_decigo_bbo).*hf_bbo,'k');
% M = 0.6 White dwarf
M = 0.6*Msol;
% Dt = 1;
% 15
two_five_PN_MESA15_puff;
wave_form_mesa15_puff;
if (~all(diff(r)))
    j = find(diff(r) == 0,1,'first');
    t(j+1:end) = []; X(j+1:end) = [];Y(j+1:end) = [];v_x(j+1:end) = [];v_y(j+1:end) = [];
    wave_form_mesa15_puff;
end
% save 15-06_puff.mat;
% load 15-06_puff.mat;
noise;
loglog(freq,2*freq.*sqrt(Sig)*(t(end)-t(n_start))/(0.4*L));
SN_puff(1) = SNR;
SN_dec_puff(1) = SNR_dec;
SN_bbo_puff(1) = SNR_bbo;

% 8
two_five_PN_Hila_puff;
wave_form_Hila_puff;
if (~all(diff(r)))
    j = find(diff(r) == 0,1,'first');
    t(j+1:end) = []; X(j+1:end) = [];Y(j+1:end) = [];v_x(j+1:end) = [];v_y(j+1:end) = [];
    wave_form_Hila_puff;
end
% save 8-06_puff.mat;
% load 8-06_puff.mat;
noise;
loglog(freq,2*freq.*sqrt(Sig)*(t(end)-t(n_start))/(0.4*L),'r');
SN_puff(2) = SNR;
SN_dec_puff(2) = SNR_dec;
SN_bbo_puff(2) = SNR_bbo;

% Dt = 1;
% 5
two_five_PN_MESA5_puff;
wave_form_5_puff;
if (~all(diff(r)))
    j = find(diff(r) == 0,1,'first');
    t(j+1:end) = []; X(j+1:end) = [];Y(j+1:end) = [];v_x(j+1:end) = [];v_y(j+1:end) = [];
    wave_form_5_puff;
end
% save 5-06_puff.mat;
% load 5-06_puff.mat;
% ind = 1;
noise;
% ind = 0;
loglog(freq,2*freq.*sqrt(Sig)*(t(end)-t(n_start))/(0.4*L),'g');
SN_puff(3) = SNR;
SN_dec_puff(3) = SNR_dec;
SN_bbo_puff(3) = SNR_bbo;


xlabel('f (Hz)');
ylabel('Characteristic Strain');
title('$m = 0.6~M_\odot, \lambda = 10$','Interpreter','Latex');
saveas(figure(1),'char_strain_06_puff.fig');
saveas(figure(1),'characteristic_strain_06_puff.eps','epsc');
hold off;
close(figure(1));

% M = 1.4 neutron star
loglog(f,sqrt(f).*hf,'k');
hold on;
loglog(f_decigo_bbo,sqrt(f_decigo_bbo).*hf_decigo,'k');
loglog(f_decigo_bbo,sqrt(f_decigo_bbo).*hf_bbo,'k');
% Dt = 1;
M = 1.4*Msol;
%15
two_five_PN_MESA15_puff;
wave_form_mesa15_puff;
if (~all(diff(r)))
    j = find(diff(r) == 0,1,'first');
    t(j+1:end) = []; X(j+1:end) = [];Y(j+1:end) = [];v_x(j+1:end) = [];v_y(j+1:end) = [];
    wave_form_mesa15_puff;
end
% save 15-14_puff.mat;
noise;
loglog(freq,2*freq.*sqrt(Sig)*(t(end)-t(n_start))/(0.4*L));
SN_puff(4) = SNR;
SN_dec_puff(4) = SNR_dec;
SN_bbo_puff(4) = SNR_bbo;

% Dt = 1;
% 8
two_five_PN_Hila_puff;
wave_form_Hila_puff;
if (~all(diff(r)))
    j = find(diff(r) == 0,1,'first');
    t(j+1:end) = []; X(j+1:end) = [];Y(j+1:end) = [];v_x(j+1:end) = [];v_y(j+1:end) = [];
    wave_form_Hila_puff;
end
% save 8-14_puff.mat;
noise;
loglog(freq,2*freq.*sqrt(Sig)*(t(end)-t(n_start))/(0.4*L),'r');
SN_puff(5) = SNR;
SN_dec_puff(5) = SNR_dec;
SN_bbo_puff(5) = SNR_bbo;

% 5
two_five_PN_MESA5_puff;
wave_form_5_puff;
if (~all(diff(r)))
    j = find(diff(r) == 0,1,'first');
    t(j+1:end) = []; X(j+1:end) = [];Y(j+1:end) = [];v_x(j+1:end) = [];v_y(j+1:end) = [];
    wave_form_5_puff;
end
% save 5-14_puff.mat;
% ind = 1;
noise;
% ind = 0;
loglog(freq,2*freq.*sqrt(Sig)*(t(end)-t(n_start))/(0.4*L),'g');
SN_puff(6) = SNR;
SN_dec_puff(6) = SNR_dec;
SN_bbo_puff(6) = SNR_bbo;

xlabel('f (Hz)');
ylabel('Characteristic Strain');
title('$m = 1.4~M_\odot, \lambda = 10$','Interpreter','Latex');
saveas(figure(1),'char_strain_14_puff.fig');
saveas(figure(1),'characteristic_strain_14_puff.eps','epsc');
hold off;
close(figure(1));
