function [outr, outq, outt] = rtlSim(t, sIn, fcar, amp )
%RTLSIM according to figure 5 of lab 6.pdf

fs = 2.4e6;                %Sample frequency
s = amp * sIn;              %Low noise amp

[~, tmax] = max(s);
phase = t(tmax);
sr = exp(-1i*2*pi*fcar*t+phase).*s;  %Shift the frequency to the left so that f0 =fc


%% Low pass filter
N   = length(t);        
fsampleSim = 1/(t(2) - t(1));

% [filterCoefb, filterCoefa] = fir1(10, 2*fs/fsampleSim);
[filterCoefb, filterCoefa] = fir1(100, 0.01);
Filter = freqz(filterCoefb, filterCoefa, floor(length(t)/2));
if rem(length(t),2)
    Filter = [fliplr(Filter') 1 Filter'];
else
    Filter = [fliplr(Filter') Filter'];
end
Sr = fft(sr);
Sr = Sr.*ifftshift(Filter);
sr = ifft(Sr);

%% Quantizer
%Downsample to fs
D = round(fsampleSim/fs);
sr = downsample(sr,D);
outt = downsample(t,D);
outr = fi(real(sr),1,8);
outq = fi(imag(sr),1,8); %8 bit signed float
