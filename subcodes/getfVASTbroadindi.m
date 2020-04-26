function [q, ctrfilter] = getfVASTbroadindi(ctrfilter, Hml, Dm, taroption)
% Broadband approach without inter-frequency information

nzones = size(Hml,1);
[Kbins, nctrpts, nloudspks] = size(Hml{nzones});

q = cellfun(@(x) zeros(nloudspks, Kbins), ...
    cell(nzones,1), 'UniformOutput', false);

if strcmpi(ctrfilter.cvxopt_properties.opttype, 'min_sb')
    minSb = true;
else
    minSb = false;
end

% Options for plots
taroption.broadband = true;
taroption.nloudspks = nloudspks;
taroption.tightfigure = true;

if strcmpi(ctrfilter.cvxopt_properties.opttype, 'min_sd')
    if strcmpi(ctrfilter.cvxopt_properties.const, 'sb')
        taroption.nsb = false;
    else
        taroption.nsb = true;
    end
end
if strcmpi(ctrfilter.cvxopt_properties.opttype, 'min_sb')
    if strcmpi(ctrfilter.cvxopt_properties.const, 'sd')
        taroption.nsd = false;
    else
        taroption.nsd = true;
    end
end
taroption.frequency = true;

stfreqidx = 1;%4;
edfreqidx = Kbins;%33;
ssidx = flipud(perms(1:nzones));
for sridx = 1:2%1:nzones
    RRb = cell(Kbins,1);
    RRd = cell(Kbins,1);
    rrb = cell(Kbins,1);
    hhz = cell(Kbins,1);
    
    H = cell(nzones,1);
    vard = zeros(2,1);
    for fbinidx = stfreqidx:edfreqidx%1:Kbins
        for zidx = 1:nzones
            H{zidx} = squeeze(Hml{zidx}(fbinidx,:,:));
        end
        
        Hb = H{ssidx(sridx,1)};
        Hd = H{ssidx(sridx,2)};
        
        d = Dm{sridx}(fbinidx,:).';
        
        RRb{fbinidx} = Hb'*Hb;
        RRd{fbinidx} = Hd'*Hd;
        rrb{fbinidx} = Hb'*d;
        hhz{fbinidx} = d;
        vard(sridx) = d'*d;
    end
    Rb = blkdiag(RRb{:});
    Rd = blkdiag(RRd{:});
    rb = vertcat(rrb{:});
    hz = vertcat(hhz{:});
    
    % Joint diagonalization
    [U, D] = jdiag(Rb, Rd, 'vector', true);
    D = real(D);
    V = ctrfilter.V;
    
    d_part = real(D(1:V));
    U_part = U(:,1:V);
    uvrb = U_part'*rb;

    % Calculate the optimal coefficient a
    if ctrfilter.cvxopt_properties.findopt
        if strcmpi(ctrfilter.cvxopt_properties.opttype, 'min_sd')
            ctrfilter = getOptPara(uvrb,d_part,hz,ctrfilter,[1, sridx]);
        else
            allones = ones(Kbins*nloudspks,1);
            iRdi = abs(allones'*Rd*allones);
            ctrfilter = getOptPara(uvrb,d_part,iRdi,ctrfilter,[1, sridx]);
        end
    else
        ctrfilter.cvxopt_properties.optpara(1,sridx) = ctrfilter.mu;
    end
    
    mu_new = ctrfilter.cvxopt_properties.optpara(1,sridx);
    ctrfilter.mu = mu_new;
    
    qf = U_part*(uvrb./(mu_new + d_part));

    if taroption.journal_exp_1
%         calinfo.mucan = logspace(-15,15,101);
%         calinfo.Vcan = 1:Kbins:nloudspks*Kbins;
%         
%         calinfo.wRd = Rd;
%         calinfo.whz = d;
%         calinfo.wrb = wrb;
%         calinfo.U = U;
%         calinfo.D = D;
%         calinfo.sridx = sridx;
%         
%         pfmmtx = get_mse_pfm(calinfo);
%         q = pfmmtx;
    else
        for fbinidx = 1:Kbins
            q{sridx}(:,fbinidx) = qf((1:nloudspks)+(fbinidx-1)*nloudspks);
        end
    end
    
    % Plot the cost functions
    if sridx == 1
%     costfcnplot(U,D,rb,abs(hz'*hz),minSb,true,taroption)
    end

end

end