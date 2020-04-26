function prms = spl2rms(spl)
% Convert SPL (dB) to Prms (Pa of rms)
% For example, SPL 94 dB is equivalent to Prms 1 Pa

prms = (2e-5)*10^(spl/20);


end