function [out] = createBPSK(t, bitStream, fbase, Rb)
%createBPSK only works if 1/(Rb*dt) is an integer
dt = t(2) - t(1);
bitLength = round(1/(Rb*dt));

%expand bitstream
bitStream2(1) = 1;
for i = 1:length(bitStream)
    bitStream2 = [bitStream2; ones(bitLength,1).*bitStream(i)];
end
bitStream2 = bitStream2(1:length(t));
%create signal with 1 --> cos and 0 --> -cos
out = bitStream2'.*cos(fbase.*t.*2.*pi) - (1-bitStream2').*cos(fbase.*t.*2.*pi);


