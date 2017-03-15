close all
clearvars
%% Variables
tmin = 0;
dt = 0.0000001;
tmax = 0.01;
fc = 100*10^3; %100 kHz
fb = 1000; %1002 Hz
Rb = 1000; %50 kHz
filter1 = [0 zeros(1,1000) 1 zeros(1,250) 0.4 zeros(1,80) 0.00001];
filter2 = [zeros(1,1025) 0.2 zeros(1,245) -0.03 zeros(1,28) 0.05];
SNR = inf; %dB awgn added at the receiver (not completly correct)
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
%y1 = double(si1 + 1i*sq1); %Should be done with a wiener deconv like in SDSP project
%y2 = double(si2 + 1i*sq2);

%Upsample to simulation sampletime (ZOH-process)
upSamp = round(1/(dt*2.4e6));
si1 = double(si1);
si1 = upsampleZOH(si1,upSamp); %2.4e6 from sampling rate RTL
sq1 = double(sq1);
sq1 = upsampleZOH(sq1,upSamp); %2.4e6 from sampling rate RTL
tup1 = 0:dt:length(si1)*dt-dt;
si2 = double(si2);
si2 = upsampleZOH(si2,upSamp); %2.4e6 from sampling rate RTL
sq2 = double(sq2);
sq2 = upsampleZOH(sq2,upSamp); %2.4e6 from sampling rate RTL
tup2 = 0:dt:length(si2)*dt-dt;

% H1I = estimate_h(si1,sBase,dt);
% H1Q = estimate_h(sq1,sBase,dt);
%H1 = (H1Q+H1I)/2;
% [H1, h1] = estimate_h(si1+1i.*sq1,sBase,t, fc);
h1 = estimate_h2(si1,sq1,stransmit,fc,tup1);
h2 = estimate_h2(si2,sq2,stransmit,fc,tup2);

%Zero pad the shorter signal
lh1 = length(h1);
lh2 = length(h2);
if lh1 >= lh2
    h2(end+1:length(h1)) = 0;
else
    h1(end+1:length(h2)) = 0;
end
H1 = fft(h1);
H2 = fft(h2);
P = -H1./H2;
p = ifft(P);
%p = deconvwnr(h1,h2); %-H1./H2; %Add extra filter??


%Transmit simultaniously
[srec1, t1] = channel(stransmit, filter1, SNR, t);
tnew = [t t(end)+dt:dt:ts2(end)];
bits = Createbitstream(Rb,tnew);               %Create bitstream
sBase = createBPSK(tnew,bits,fb,Rb);               %Create (BPSK) signal
stransmitLong = transmit(tnew,fc,sBase);           %Transmit signal over carrier
lp = length(p);
lstransmitLong = length(stransmitLong);
if lp >= lstransmitLong
    stransmitLong(lstransmitLong+1:lp) = zeros(1,lp-lstransmitLong);
else
    p(lp+1:lstransmitLong) = zeros(1,lstransmitLong-lp);
end
P = fft(p);
StransmitLong = fft(stransmitLong);
Stransmit2 = P.*StransmitLong;
stransmit2 = ifft(Stransmit2);
t2in = [t t(end)+dt:dt:length(p)*dt-dt];
[srec2, t2] = channel(stransmit2, filter2, SNR, t2in);

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
plot(filter1);
hold on
plot(h1);
% subplot(2,1,2)
% fmin = -(length(t)-1)/2;
% fmax = (length(t)-1)/2;
% df = 1;
% f = fmin:df:fmax;
% plot(f,fftshift(abs(fft(s))))
