function dB = rms2spl(X)
% X is either (N x 1) vector, input signal or a scalar
pa = rms(X);
pref = 2e-5;
ratio = pa/pref;

dB = 20*log10(ratio);
end