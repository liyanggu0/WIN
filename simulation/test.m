close all
tmin = 0;
dt = 0.0000001;
tmax = 0.005;
fc = 100*10^3; %100 kHz
fbase = 1000; %1000 Hz
Rb = 50000; %100 Hz
filter = [zeros(1,1000) 0.3 zeros(1,25) 0.2 zeros(1,80) 0.4 zeros(1,500) 0.1];

t = tmin:dt:tmax;

st = createTransmitter(t,fc,Rb,fbase);
s = channel(st, filter);


%% plot s to t and plot s in f domain
subplot(2,1,1)
tplot = tmin:dt:tmax+length(filter)*dt;
plot(tplot,s);
subplot(2,1,2)
fmin = -(length(tplot)-1)/2;
fmax = (length(tplot)-1)/2;
df = 1;
f = fmin:df:fmax;
plot(f,fftshift(abs(fft(s)))) 
