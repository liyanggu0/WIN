function [out] = createCarrier(t, fc, phase)
out = cos(fc.*t.*2.*pi+phase);
end

