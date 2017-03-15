function [s] = upsampleZOH(s,up)
%UPSAMPLE Summary of this function goes here
%   Detailed explanation goes here
s = upsample(s,up); %2.4e6 from sampling rate RTL
for i = 1:length(s) %ZOH process
    delay = rem(i,up)-1;
    if delay == -1
        delay = 3;
    end
    s(i) = s(i-delay);
end

end

