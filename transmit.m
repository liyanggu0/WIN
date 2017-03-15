function [s] = transmit(t, fc, sBase)
%CREATETRANSMITTER creates a transmitter transmitting a 101010 bitstream
%with bitrate Rb with baseband frequency fb and carrier frequency fc. If fb
%is not given the BPSK will use the carrier for the shift keying.

sCar = createCarrier(t,fc,0);
s = sCar .* sBase;
end

