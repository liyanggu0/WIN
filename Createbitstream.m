function [ out ] = Createbitstream(Rb, t)
%Creates bitstream in the form of "01010101" till end of time
tTotal = t(end) - t(1);
streamLength = ceil(tTotal*Rb);
for i = 1:streamLength
    if rem(i,2) == 0
        out(i) = 0;
    else
        out(i) = 1;
    end
end
