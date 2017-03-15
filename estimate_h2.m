function [h] = estimate_h2(yr,yc,x,fc,t)
%ESTIMATE_H2 Summary of this function goes here
%   Detailed explanation goes here
y = yr.*cos(2*pi*fc*t)+yc.*sin(2*pi*fc*t);
rxx = xcorr(x);
ryx = xcorr(y,x);
h = deconvwnr(ryx,rxx);
h = h(end/2:end/2+length(y)-length(x));
end

