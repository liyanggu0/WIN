function [ H , h] = estimate_h(y,x,t,fc)
%ESTIMATE_H Summary of this function goes here
%   Detailed explanation goes here
%Get the autocorrelation of the reference signal
c = 3e8; %speed of light
x_real = cos(fc.*t.*2.*pi).*x;
x_sim = exp(-1i*2*pi*fc*t).*x_real; 

fsampleSim = 1/(t(2) - t(1));
fs = 2.4e6;
%[filterCoefb, filterCoefa] = fir1(10,2*fs/fsampleSim);
[filterCoefb, filterCoefa] = fir1(100, 0.01);
Filter = freqz(filterCoefb, filterCoefa, floor(length(t)/2));
if rem(length(t),2)
    Filter = [fliplr(Filter') 1 Filter'];
else
    Filter = [fliplr(Filter') Filter'];
end
X_sim = fft(x_sim);
X_sim = X_sim.*ifftshift(Filter);
x_sim = ifft(X_sim);
x = real(x_sim);
y = real(y);
autoCorRef = xcorr(x);
%autoCorRef = autoCorRef(end/2+1:end);
%% Get the channel

crossCor = xcorr(y,x);
% % h = deconvwnr(crossCor(end/2+0.5:end),autoCorRef(end/2+0.5:end),1);
% % h = fliplr(h);
% % h(1:floor(end/2)) = [];
% % h(length(y)-length(x):end) = [];
autoCorRef(end+1:length(crossCor)) = zeros(1,length(crossCor)-length(autoCorRef));
% Px = fft(autoCorRef);
% Pyx = fft(crossCor);
% Px = Px/max(Px);
% Pyx = Pyx/max(Pyx);
x(end+1:length(y)) = zeros(1,length(y)-length(x));
h = deconvwnr(crossCor(end/2+0.5:end),autoCorRef(end/2+0.5:end),1);
Px = (fft(x)).^2;
Pyx = fft(crossCor(end/2:end));
H = Px./Pyx;
Filter = freqz(filterCoefb, filterCoefa, floor(length(H)/2));
if rem(length(H),2)
    Filter = [fliplr(Filter') 1 Filter'];
else
    Filter = [fliplr(Filter') Filter'];
end
H = H.*ifftshift(Filter) .* ifftshift(Filter);
h = ifft(H);
h = h(end/2:end/2+(5/c)/dt);
% H = fft(h);
end

