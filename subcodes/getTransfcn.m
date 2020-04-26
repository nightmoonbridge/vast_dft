function [Hml, Dm, hml, dm] = getTransfcn(general,array,zone,rirs,drirs,figoption)
% [hml, dm, Hml, Dm] = GETTRANSFCN(general,array,zone,rirs,drirs)
% A function GETTRANSFCN returns the transfer functions from room impulse
% response functions in the time domain
% The first two output variables are for the selected transfer functions
% The last two output variables are for the all transfer functions with a
% proper amount of frequency bins
%
% Hml : Selected transfer function between loudspk l to ctr pts m
% Dm  : Selected transfer function between virtual src z to ctr pts m
% hml : Transfer function between loudspk l to ctr pts m
% dm  : Transfer function between virtual src z to ctr pts m
%

if nargin < 6
    figoption = false;
end

try
    incl_dcnyq = general.incl_dcnyq;
catch
    incl_dcnyq = false;
end


% if general.setup.simrir
%     lenir = general.lenIResponse;
% else
%     lenir = general.lenRir;
% end
lenir = size(drirs{1},1);

multiplierfactor = ceil(lenir/general.lenConFilter/10)*10;
nfft = multiplierfactor*general.lenConFilter;

dF = general.fs/general.lenConFilter;

pickidx = dF/(general.fs/nfft);
pidx = (0:pickidx:nfft/2)+1;  % incl. the DC and the Nyquist frequency

if ~incl_dcnyq
    pidx = pidx(2:end-1);  % excl. the DC and the Nyquist frequency
end

hml = cellfun(@(x) zeros(nfft,zone.numCtrPts,array.numLoudspk), ...
    cell(zone.number,1), 'UniformOutput', false);
dm = cellfun(@(x) zeros(nfft,zone.numCtrPts), ...
    cell(zone.number,1), 'UniformOutput', false);

% Convert RIRs and IR into the frequency domain
% rirs   >> hml
% drirs  >> hz
for ss = 1:zone.number
    for midx = 1:zone.numCtrPts
        for lidx = 1:array.numLoudspk
            inlidx = (1:lenir) + (lidx-1)*lenir;
            hml{ss}(:,midx,lidx) = fft(rirs{ss}(inlidx,midx),nfft);
        end
        dm{ss}(:,midx) = fft(drirs{ss}(:,midx),nfft);
    end
end

Hml = cellfun(@(x) zeros(length(pidx),zone.numCtrPts,array.numLoudspk), ...
    cell(zone.number,1), 'UniformOutput', false);
Dm = cellfun(@(x) zeros(length(pidx),zone.numCtrPts), ...
    cell(zone.number,1), 'UniformOutput', false);
for sridx = 1:zone.number
    Hml{sridx}(:,:,:) = hml{sridx}(pidx,:,:);
    Dm{sridx} = dm{sridx}(pidx,:);
end

if figoption
    freq = (0:dF:general.fs/2)';
    if ~incl_dcnyq
        freq = freq(2:end-1);
    end
    cmic = round(median(1:zone.numCtrPts));
    cloudspk = round(median(1:array.numLoudspk));
    figure('Name','Magnitude spectrum of RIR for ACC in the frequency domain')
    stem(freq,abs(Hml{1}(:,cmic,cloudspk)),'LineStyle','none','DisplayName','Chosen RIR')
    hold on
    plot((0:nfft-1)/nfft*general.fs,abs(hml{1}(:,cmic,cloudspk)), 'DisplayName', 'Measured RIR')
    xlim([0 general.fs/2])
    grid on
    title('Magnitude spectrum')
    xlabel('Frequency (Hz)'),ylabel('Amplitude')
    legendhitcallback
end

end