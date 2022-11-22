
n_start = 300;
if (indic == 1)
     n_start = find(t<=0.9*t(end),1,'last');
end

signal = 0.5*(1+cos(i)^2)*cos(2*phi0)*h_plus + cos(i)*sin(2*phi0)*h_cross;
% sign = interp1(t,signal,(t(n_start):Dt:t(end)),'linear','extrap');
% h = fft(sign)';
% L = length(h);
% S = abs(h).^2;
% S(floor(L/2+1):end) = [];
% df = 1/(L*Dt);
% freq = df*(0:(L-1));
% freq(floor(L/2+1):end) = [];

[S,Sig,freq,L,df] = fourier_bit(n_start,Dt,signal,t);
load lisa.mat f hf;
S_n = @(x) interp1(f,hf.^2,x,'linear','extrap');

% Integral
% S = @(x) interp1(freq,S,x,'linear','extrap');
func = @(freq) S(freq)./S_n(freq);
f_min = max([df min(f)]);
f_max = min([df*(L-1)/2 max(f)]);
SNR = 2*Dt*sqrt(integral(func,f_min,f_max));
S_n_dec = @(x) interp1(f_decigo_bbo,hf_decigo.^2,x,'linear','extrap');
f_min = max([df min(f_decigo_bbo)]);
f_max = min([df*(L-1)/2 max(f_decigo_bbo)]);
func = @(freq) S(freq)./S_n_dec(freq);
SNR_dec = 2*Dt*sqrt(integral(func,f_min,f_max));
S_n_bbo = @(x) interp1(f_decigo_bbo,hf_bbo.^2,x,'linear','extrap');
func = @(freq) S(freq)./S_n_bbo(freq);
SNR_bbo = 2*Dt*sqrt(integral(func,f_min,f_max));

% Riemann sum
% func = S./S_n(freq);
% SNR = 2*sqrt(sum(func)*Dt^2*df);

%fprintf('SNR = %f \n',double(vpa(SNR)));

% semilogx(min(f):df:max(f),func(min(f):df:max(f)));


