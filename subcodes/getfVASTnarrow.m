function [q, ctrfilter] = getfVASTnarrow(ctrfilter, Hml, Dm, taroption)
% Independent Kbins frequency bins
taridx = taroption.taridx;

try
    journal_exp_1 = taroption.journal_exp_1;
catch
    journal_exp_1 = false;
end

nzones = size(Hml,1);
[Kbins, nctrpts, nloudspks] = size(Hml{nzones});

allones = ones(nloudspks,1);

q = cellfun(@(x) zeros(nloudspks, Kbins), ...
    cell(nzones,1), 'UniformOutput', false);

if strcmpi(ctrfilter.cvxopt_properties.opttype, 'min_sb')
    minSb = true;
else
    minSb = false;
end

% Since this is the narrowband approach
taroption.broadband = false;
taroption.tightfigure = true;
taroption.nloudspks = nloudspks;
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

if journal_exp_1
    ff = taridx;
    H = cell(nzones,1);
    ssidx = flipud(perms(1:nzones));
    for sridx = 1:nzones
        H{sridx} = squeeze(Hml{sridx}(ff,:,:));
    end
    
    for sridx = 1:nzones
        Hb = H{ssidx(sridx,1)};
        Hd = H{ssidx(sridx,2)};
        
        d = Dm{sridx}(ff,:).';
        
        Rb = Hb'*Hb;
        Rd = Hd'*Hd;
        rb = Hb'*d;
        
        % Joint diagonalization
        [U, D] = jdiag(Rb, Rd, 'vector', true);
        
        calinfo.mucan = logspace(-15,15,101);
        calinfo.Vcan = 1:16;
        
        calinfo.wRd = Rd;
        calinfo.whz = d;
        calinfo.wrb = rb;
        calinfo.U = U;
        calinfo.D = D;
        calinfo.sridx = sridx;
        calinfo.tarfreq = taroption.tarfreq;
        
        pfmmtx = get_mse_pfm(calinfo);
        q = pfmmtx;
    end

else
    
    for ff = 1:Kbins
        H = cell(nzones,1);
        ssidx = flipud(perms(1:nzones));
        for sridx = 1:nzones
            H{sridx} = squeeze(Hml{sridx}(ff,:,:));
        end

        for sridx = 1:nzones
            Hb = H{ssidx(sridx,1)};
            Hd = H{ssidx(sridx,2)};

            d = Dm{sridx}(ff,:).';

            Rb = Hb'*Hb;
            Rd = Hd'*Hd;
            rb = Hb'*d;

            % Joint diagonalization
            [U, D] = jdiag(Rb, Rd, 'vector');

            V = ctrfilter.V;

            d_part = real(D(1:V));
            U_part = U(:,1:V);
            uvrb = U_part'*rb;

            % Calculate the optimal coefficient a
            if ctrfilter.cvxopt_properties.findopt
                if strcmpi(ctrfilter.cvxopt_properties.opttype, 'min_sd')
                    ctrfilter = getOptPara(uvrb,d_part,d,ctrfilter,[ff, sridx]);
                else
                    iRdi = abs(allones'*Rd*allones);
                    ctrfilter = getOptPara(uvrb,d_part,iRdi,ctrfilter,[ff, sridx]);
                end

                mu_new = ctrfilter.cvxopt_properties.optpara(ff, sridx);
            else
                mu_new = ctrfilter.mu;
            end

            ctrfilter.mu = mu_new;

            qf = U_part*(uvrb./(mu_new + d_part));
            q{sridx}(:,ff) = qf;

%             % Plot the cost functions
%             if sridx == 1 && ff == taridx
%                 costfcnplot(U,D,wrb,abs(whz'*whz),minSb,true,taroption)
%             end

        end

    end
end

end