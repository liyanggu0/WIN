function [out, tout] = channel(in, h, SNR, t)
%TO DO add white gaussian noise to the channel
lin = length(in);
lh = length(h);
out(lin+lh) = 0;
for i = 1:length(h)
    if h(i) ~= 0
        out(i:lin+i-1) = out(i:lin+i-1) + in*h(i);
    end
end

%% Add noise
P = rms(out).^2;
P = 10.*log(P); %[dBW]
out = awgn(out, SNR, P);
dt = t(2) - t(1);
tout = [t t(end)+dt:dt:t(end)+length(h)*dt];
