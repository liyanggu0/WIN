function [ H ] = estimate_h(y,x,dt)
%ESTIMATE_H Summary of this function goes here
%   Detailed explanation goes here
%Get the autocorrelation of the reference signal
c = 3e8; %speed of light
autoCorRef = xcorr(x);
%autoCorRef = autoCorRef(end/2+1:end);
%% Get the channel
crossCor = xcorr(y,x);
h = deconvwnr(crossCor(end/2+0.5:end),autoCorRef(end/2+0.5:end),1);
h(1:length(autoCorRef)/2) = [];
%h() = [];
%h = diff(h);
% Px = fft(autoCorRef);
% Pyx = fft(crossCor);
% H = Px./Pyx;
% h = ifft(H);
% h = h(end/2:end/2+(5/c)/dt);
H = fft(h);
end

