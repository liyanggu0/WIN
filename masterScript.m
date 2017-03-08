close all
clearvars
%% Variables
tmin = 0;
dt = 0.0000001;
tmax = 0.001;
fc = 100*10^3; %100 kHz
fb = 1000; %1000 Hz
Rb = 50000; %50 kHz
filter1 = [zeros(1,1000) 0.05 zeros(1,250) 0.02 zeros(1,80) 0];
filter2 = [zeros(1,1025) 0.2 zeros(1,245) -0.03 zeros(1,28) 0.05];
SNR = inf; %5; %dB awgn added at the receiver
RTLAmp = 10;
t = tmin:dt:tmax;

%% Load
%{[y, Fs] = audioread(With Human Between.wav);
%[y, Fs] = audioread(Without Objects.wav);

%% Transmitter
bits = Createbitstream(Rb,t);               %Create bitstream
sBase = createBPSK(t,bits,fb,Rb);           %Create (BPSK) signal
stransmit = transmit(t,fc,sBase);           %Transmit signal over carrier

%% Simulation
[srec1, t1] = channel(stransmit, filter1, SNR, t);
[srec2, t2] = channel(stransmit, filter2, SNR, t);

[si1, sq1, ts1] = rtlSim(t1, srec1, fc, RTLAmp);
[si2, sq2, ts2] = rtlSim(t2, srec2, fc, RTLAmp);

%% initial nulling
y1 = double(si1 + 1i*sq1); %Should be done with a wiener deconv like in SDSP project
y2 = double(si2 + 1i*sq2);

upSamp = round(1/(dt*2.4e6));
si1 = double(si1);
si1 = upsampleZOH(si1,upSamp); %2.4e6 from sampling rate RTL
sq1 = double(sq1);
sq1 = upsampleZOH(sq1,upSamp); %2.4e6 from sampling rate RTL

H1I = estimate_h(si1,sBase,dt);
H1Q = estimate_h(sq1,sBase,dt);
H1 = (H1Q+H1I)/2;
H2 = estimate_h(y2,stransmit,dt);
lh1 = length(H1);
lh2 = length(H2);
if lh1 >= lh2
    H1(lh2+1:lh1) = [];
else
    H2(lh1+1:lh2) = [];
end
P = -H1./H2; %Add extra filter??
p = ifft(P);


%Transmit simultaniously
[srec1, t1] = channel(stransmit, filter1, SNR, t);
tnew = [t t(end)+dt:dt:ts2(end)];
bits = Createbitstream(Rb,tnew);               %Create bitstream
sBase = createBPSK(tnew,bits,fb,Rb);               %Create (BPSK) signal
stransmitLong = transmit(tnew,fc,sBase);           %Transmit signal over carrier
lp = length(p);
lstransmitLong = length(stransmitLong);
if lp >= lstransmitLong
    p(lstransmitLong+1:lp) = [];
else
    stransmitLong(lp+1:lstransmitLong) = [];
end
t2in = [t t(end)+dt:dt:length(p)*dt-dt];
[srec2, t2] = channel(p.*stransmitLong, filter2, SNR, t2in);

if length(srec1) >= length(srec2)
    srec2(end+1:length(srec1)) = 0;
    trec = t1;
else
    srec1(end+1:length(srec2)) = 0;
    trec = t2;
end
srec = srec1 + srec2;
[si, sq, ts] = rtlSim(trec, srec, fc, RTLAmp);
%% plot s to t and plot s in f domain
figure
%subplot(2,1,1)
plot(ts,si);
hold on
plot(ts,sq);
% subplot(2,1,2)
% fmin = -(length(t)-1)/2;
% fmax = (length(t)-1)/2;
% df = 1;
% f = fmin:df:fmax;
% plot(f,fftshift(abs(fft(s))))
