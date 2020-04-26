function [ctrfilter, perform, rankcheck] = calculatefVAST(general, array, zone, ctrfilter, rirs, drirs, dummy, bandoption, taroption)

if nargin < 8
    bandoption = 'narrow'; % 'narrow', 'broadindi'
end

% The convex optimization problem properties
% opttype    : 'min_sb' or 'min_sd'
% const      : if 'min_sb', then 'sd'
%              if 'min_sd', then 'sb'
% tarval     : the target value, by default, 1e-5
if ~isfield(ctrfilter,'cvxopt_properties')
    ctrfilter.cvxopt_properties.opttype = 'min_sb';
    ctrfilter.cvxopt_properties.const = 'sd';
    ctrfilter.cvxopt_properties.tarval = 1e-5;
    ctrfilter.cvxopt_properties.findopt = true;
    ctrfilter.cvxopt_properties.initpara = 0;
    ctrfilter.cvxopt_properties.optpara = 0;
end

nzones = zone.number;
nloudspks = array.numLoudspk;
nctrpts = zone.numCtrPtsBr;

if ctrfilter.incl_dcnyq
    general.incl_dcnyq = true;
else
    general.incl_dcnyq = false;
end

[Hml, Dm, hml, dm] = getTransfcn(general,array,zone,rirs,drirs);

dF = general.fs/general.lenConFilter;
freq = (0:dF:general.fs/2)';

if ~ctrfilter.incl_dcnyq
    freq = freq(2:end-1);  % excl. the DC and the Nyquist frequency
end

% calculate the performance at a single frequency for VAST-NF
try
    taridx = taroption.taridx;
catch
    taridx = 16; % the index of 1kHz when df = 66.66667 Hz
end
try
    tarfreq = taroption.tarfreq;
catch
    tarfreq = (taridx-1)*dF;
end

Kbins = length(freq);

priorAC = cellfun(@(x) zeros(Kbins,1), ...
    cell(zone.number,1), 'UniformOutput', false);
inputpow = priorAC;
nsde = priorAC;
sde = priorAC;
postAC = priorAC;
iRbi = priorAC;
iRdi = priorAC;
nre = priorAC;
re = priorAC;

q = cellfun(@(x) zeros(nloudspks, Kbins), ...
    cell(nzones,1), 'UniformOutput', false);

q_fvast = cellfun(@(x) zeros(nloudspks*general.lenConFilter,1), ...
    cell(nzones,1), 'UniformOutput', false);

rankcheck = cellfun(@(x) zeros(Kbins,2), ...
    cell(nzones,1), 'UniformOutput', false);

switch bandoption
    case 'narrow'
        taroption.taridx = taridx;
        taroption.tarfreq = tarfreq;
        % Independent Kbins frequency bins
        [q, ctrfilter] = getfVASTnarrow(ctrfilter, Hml, Dm, taroption);

    case 'broadindi'
        [q, ctrfilter] = getfVASTbroadindi(ctrfilter, Hml, Dm, taroption);

    case 'time'
        nfft = general.lenConFilter;
        for zidx = 1:nzones
            for lidx = 1:nloudspks
                inidx = (1:nfft) + (lidx-1)*nfft;
                getfft = fft(ctrfilter.conFilter{zidx}(inidx), nfft);
                
                if ctrfilter.incl_dcnyq
                    q{zidx}(lidx,:) = getfft(1:nfft/2+1);
                else
                    q{zidx}(lidx,:) = getfft(2:nfft/2);
                end
            end
        end
    otherwise
end

if taroption.journal_exp_1
    return
else

    Hbq_Kbins = zeros(nctrpts,Kbins);
    d_Kbins = zeros(nctrpts,Kbins);

    for fbinidx = 1:Kbins
        H = cell(nzones,1);
        ssidx = flipud(perms(1:nzones));
        for sridx = 1:nzones
            H{sridx} = squeeze(Hml{sridx}(fbinidx,:,:));
        end
        for sridx = 1:nzones
            Hb = H{ssidx(sridx,1)};
            Hd = H{ssidx(sridx,2)};

            d = Dm{sridx}(fbinidx,:).';
            hz = d;

            Rb = Hb'*Hb;
            Rd = Hd'*Hd;
            rb = Hb'*hz;

            qf = q{sridx}(:,fbinidx);
            pfm = pfm_mtx_mse;

            Hbq_Kbins(:,fbinidx) = Hb*qf;
            d_Kbins(:,fbinidx) = d;

            % Using an object function
            pfm.getpriorAC(Rb,Rd);
            pfm.getpostAC(Rb,Rd,qf);
            pfm.getnsde(Rb,Hb,hz,qf);
            pfm.getsde(Rb,Hb,hz,qf);
            pfm.getiRi(Rb,true);
            pfm.getiRi(Rd,false);

            pfm.getnre(Rd,qf);
            pfm.getre(Rd,qf);


            priorAC{sridx}(fbinidx) = pfm.priorAC;
            postAC{sridx}(fbinidx) = pfm.postAC;

            nsde{sridx}(fbinidx) = pfm.nsde;
            sde{sridx}(fbinidx) = pfm.sde;

            nre{sridx}(fbinidx) = pfm.nre;
            re{sridx}(fbinidx) = pfm.resiEner;

            iRbi{sridx}(fbinidx) = pfm.iRbi;
            iRdi{sridx}(fbinidx) = pfm.iRdi;
        end
    end

    perform.inputpow = inputpow;
    perform.priorAC = priorAC;
    perform.postAC = postAC;
    perform.nsde = nsde;
    perform.sde = sde;
    perform.iRbi = iRbi;
    perform.iRdi = iRdi;

    perform.nre = nre;
    perform.re = re;

    perform = orderfields(perform);

    for sridx = 1:nzones
        fidx = 1:general.lenConFilter;
        for lidx = 1:nloudspks
            if ctrfilter.incl_dcnyq
               Q = [q{sridx}(lidx,:),flip(conj(q{sridx}(lidx,2:end-1)))].';
            else
                Q = [0,q{sridx}(lidx,:),0,flip(conj(q{sridx}(lidx,:)))].';
            end
            q_fvast{sridx}(fidx) = real(ifft(Q));
            fidx = fidx + general.lenConFilter;
        end
    end

    ctrfilter.conFilter = q_fvast;

    misfig = true;
    mcenter = round(median(1:nctrpts));
    if misfig%figoption
        maxcont = max([10*log10(postAC{1});10*log10(postAC{2})]);
        mincont = min([10*log10(postAC{1});10*log10(postAC{2})]);

        figure('Name', 'AC with respect to frequency','Position', [406 321 880 560])

        semilogx(freq, 10*log10(postAC{1}),'r', freq, 10*log10(postAC{2}),'k:');
        xlabel('Frequency (Hz)');ylabel('AC (dB)');grid on;ylim([mincont/1.1 maxcont*1.1]);xlim([100 general.fs/2]);

        % Reproduced sound pressure field
        minspl = min([10*log10(abs(d_Kbins(mcenter,:))),10*log10(abs(Hbq_Kbins(mcenter,:)))]);
        maxspl = max([10*log10(abs(d_Kbins(mcenter,:))),10*log10(abs(Hbq_Kbins(mcenter,:)))]);

        figure('Name', ['Reprod & desired sound pressure at the center_' bandoption])
        semilogx(freq, 10*log10(abs(d_Kbins(mcenter,:)))+94, 'DisplayName', 'desired')
        hold on
        semilogx(freq, 10*log10(abs(Hbq_Kbins(mcenter,:)))+94, 'DisplayName', 'reprod')
        hold off
        ylabel('SPL (dB) ref. 20 uPa')
        xlabel('Frequency (Hz)')
        ylim([round(minspl/10)*10-20 round(maxspl/10)*10+20]+94)
        xlim([100 general.fs/2])
        grid on
        legendhitcallback
    end
end

end