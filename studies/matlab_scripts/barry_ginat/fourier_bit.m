function [Signal,Sig,frequence,L,df] = fourier_bit(n_start,Dt,signal,t)

sign = interp1(t,signal,(t(n_start):Dt:t(end)),'linear','extrap');
h = fft(sign)';
L = length(h);
S = abs(h).^2;
S(floor(L/2+1):end) = [];
df = 1/(L*Dt);
freq = df*(0:(L-1));
freq(floor(L/2+1):end) = [];

frequence = linspace(min(freq),max(freq),1e5);
Sig = interp1(freq,S,frequence,'linear','extrap');
Signal = @(x) interp1(freq,S,x,'linear','extrap');
end

