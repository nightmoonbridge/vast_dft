function varout = getVariables(general,array,zone)
% varout = GETVARIABLES(general,array,zone)
% This function returns some variables in order to calculate the following
% frequency information:
%  1. The transfer function matrix between loudspks and ctr pts
%  2. The transfer function vector between vir src and ctr pts
%
% The output 'varout' contains
%  1. lenir            : the length of the impulse response
%  2. multiplierfactor : a multiplication factor for nfft
%  3. nfft             : the length of fft
%  4. dF               : the frequency resolution
%  5. pidx             : the frequency bin index set
%

narginchk(1,3)

if ~isfield(general.setup,'simrir')
    warning('''general.setup.simrir'' is missing.')
    general.setup.rir = false;
end

if general.setup.simrir
    lenir = general.lenIResponse;
else
    try
        lenir = general.lenRir;
    catch
        lenir = general.lenIResponse;
    end
end

multiplierfactor = ceil(lenir/general.lenConFilter/10)*10;
nfft = multiplierfactor*general.lenConFilter;

dF = general.fs/general.lenConFilter;

pickidx = dF/(general.fs/nfft);
pidx = (0:pickidx:nfft/2)+1;  % incl. the DC and the Nyquist frequency
freq = (0:dF:general.fs/2);

if ~general.incl_dcnyq
    pidx = pidx(2:end-1);  % excl. the DC and the Nyquist frequency

    freq = freq(2:end-1);  % excl. the DC and the Nyquist frequency
end

varout.lenir = lenir;
varout.multiplierfactor = multiplierfactor;
varout.nfft = nfft;
varout.dF = dF;

varout.pickidx = pickidx;
varout.pidx = pidx;

varout.freq = freq;

if nargin >= 2
    if ~isempty(array) && isfield(array,'numLoudspk')
        varout.nloudspks = array.numLoudspk;
    end
    if nargin == 3 && ~isempty(zone) && isfield(zone,'number')
        varout.nzones = zone.number;
        varout.nctrpts = zone.numCtrPts;
    end
end

varout.Kbins = length(freq);
varout.LJ = array.numLoudspk*general.lenConFilter;
varout.nmethods = general.nmethods;
varout.fs = general.fs;
varout.sfactor = general.sfactor;

varout = orderfields(varout);
end