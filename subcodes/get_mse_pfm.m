function pfmmtx = get_mse_pfm(calinfo)

mucan = calinfo.mucan;
Vcan = calinfo.Vcan;

sridx = calinfo.sridx;

wRd = calinfo.wRd;
whz = calinfo.whz;
wrb = calinfo.wrb;
U = calinfo.U;
D = calinfo.D;

lv = length(Vcan);
lm = length(mucan);
oAC = zeros(length(Vcan),length(mucan));
nsde = zeros(length(Vcan),length(mucan));
nre = zeros(length(Vcan),length(mucan));

sde = nsde;
re = nre;

ovt = ones(size(wRd,1),1);
iRdi = abs(ovt'*wRd*ovt);

desvar = abs(whz'*whz);
for midx = 1:lm
    mu = mucan(midx);
    for vidx = 1:lv
        V = Vcan(vidx);
        d_part = real(D(1:V));
        U_part = U(:,1:V);
        uvrb = U_part'*wrb;
        
        qRbq = d_part.*abs(uvrb).^2 ./ (mu + d_part).^2;
        
        qRdq = abs(uvrb).^2 ./ (mu + d_part).^2;
        
        oAC(vidx,midx) = sum(qRbq) / sum(qRdq);
        
        sde(vidx,midx) = desvar - sum(((2*mu + d_part).*abs(uvrb).^2)./ (mu + d_part).^2);
        
        re(vidx,midx) = sum(qRdq);
        
        nsde(vidx,midx) = sde(vidx,midx) / desvar;
        
        nre(vidx,midx) = re(vidx,midx) / iRdi;
    end
end

pfmmtx.oAC = oAC;
pfmmtx.nsde = nsde;
pfmmtx.nre = nre;
pfmmtx.sde = sde;
pfmmtx.re = re;
pfmmtx.Vcan = Vcan;
pfmmtx.mucan = mucan;


if sridx == 1
    % Note -------------------------------------------------------------- %
    % In this implementation, 
    % the simulated RIRs by rir_generator are used.
    % Therefore, the results are similar but different from those shown in
    % the manuscript.
    % -------------------------------------------------------------- Note %
    
    str_ref = {'*Note*',...
        'In this implementation, the simulated RIRs by rir_generator are used.'};
    str_ref1 = 'Therefore, the results are similar but different from those shown in Fig. ';
    bgc = [1,1,1,0.8];
    % %%
    % Fig. 5 (a)
    str1 = [str_ref(:)', {[str_ref1 '5 (a).']}];
    figname = ['exp1_rank_vs_oAC_zone' num2str(sridx) '_at_' num2str(calinfo.tarfreq/1000) 'kHz'];
    figure('Name', figname)
    plot(Vcan,10*log10(oAC(:,[1,51,lm])))
    xlim([1 16])
    ylim([-3 17])
    legend({'\mu = 0', '\mu = 1', '\mu = \infty'},'Interpreter','tex')
    xlabel('Subspace rank $V$', 'Interpreter', 'latex')
    ylabel('oAC (dB)')
    grid minor
    text(1.2,-1,str1,'BackgroundColor',bgc)
    title(figname)
    
    % %%
    % Fig. 5 (b)
    str2 = [str_ref(:)', {[str_ref1 '5 (b).']}];
    figname = ['exp1_rank_vs_SDE_zone' num2str(sridx) ''];
    figure('Name', figname)
    plot(Vcan,10*log10(sde(:,[1,51,lm])))
    xlim([1 16])
    ylim([-20 -5])
    legend({'\mu = 0', '\mu = 1', '\mu = \infty'},'Interpreter','tex')
    xlabel('Subspace rank $V$', 'Interpreter', 'latex')
    ylabel('SDE (dB)')
    grid minor
    text(1.2, -18, str2, 'BackgroundColor', bgc)
    title(figname)
    

    % %%
    % Fig. 6
    str3 = [str_ref(:)', {[str_ref1 '6.']}];
    figname = ['exp1_rank_vs_RE_zone' num2str(sridx) ''];
    figure('Name',figname)
    plot(Vcan,10*log10(re(:,[1,51,lm])))
    xlim([1 16])
    % ylim([-50 -20])
    legend({'\mu = 0', '\mu = 1', '\mu = \infty'},'Interpreter','tex')
    xlabel('Subspace rank $V$', 'Interpreter', 'latex')
    ylabel('RE (dB)')
    grid minor
    text(1.2, -150, str3, 'BackgroundColor', bgc)
    title(figname)
    
    % %%
    % Fig. 7 (a)
    str4 = [str_ref(:)', {[str_ref1 '7 (a).']}];
    figname = ['exp1_mu_vs_oAC_zone' num2str(sridx) ''];
    figure('Name', figname)
    semilogx(mucan,10*log10(oAC([1,2,4,8,16],:)'))
    xlim(10.^[-7 7])
    ylim([-3 17])
    legend({'V = 1', 'V = 2', 'V = 4', 'V = 8', 'V = 16'})
    xlabel('\mu', 'Interpreter', 'tex')
    ylabel('oAC (dB)')
    grid minor
    text(10^(-6.8), 2, str4, 'BackgroundColor', bgc)
    title(figname)

    % %%
    % Fig. 7 (b)
    str5 = [str_ref(:)', {[str_ref1 '7 (b).']}];
    figname = ['exp1_mu_vs_nSDE_zone' num2str(sridx) ''];
    figure('Name',figname)
    semilogx(mucan,10*log10(nsde([1,2,4,8,16],:)'))
    xlim(10.^[-7 7])
    ylim([-20 3])
    legend({'V = 1', 'V = 2', 'V = 4', 'V = 8', 'V = 16'})
    xlabel('\mu', 'Interpreter', 'tex')
    ylabel('nSDE (dB)')
    grid minor
    text(10^(-6.8), -12, str5, 'BackgroundColor', bgc)
    title(figname)

    % %%
    % Fig. 7 (c)
    str6 = [str_ref(:)', {[str_ref1 '7 (c).']}];
    figname = ['exp1_mu_vs_nRE_zone' num2str(sridx) ''];
    figure('Name',figname)
    semilogx(mucan,10*log10(nre([1,2,4,8,16],:)'))
    xlim(10.^[-7 7])
    ylim([-70 -5])
    legend({'V = 1', 'V = 2', 'V = 4', 'V = 8', 'V = 16'})
    xlabel('\mu', 'Interpreter', 'tex')
    ylabel('nRE (dB)')
    grid minor
    text(10^(-6.8), -55, str6, 'BackgroundColor', bgc)
    title(figname)
    
    
    
    % %%
%     figname = ['exp1_rank_vs_nSDE_zone' num2str(sridx) ''];
%     figure('Name', figname)
%     plot(Vcan,10*log10(nsde(:,[1,51,lm])))
%     xlim([1 16])
%     ylim([-20 0])
%     legend({'\mu = 0', '\mu = 1', '\mu = \infty'},'Interpreter','tex')
%     xlabel('Subspace rank $V$', 'Interpreter', 'latex')
%     ylabel('nSDE (dB)')
%     title(figname)

%     % %%
%     figname = ['exp1_rank_vs_nRE_zone' num2str(sridx) ''];
%     figure('Name',figname)
%     plot(Vcan,10*log10(nre(:,[1,51,lm])))
%     xlim([1 16])
%     % ylim([-50 -20])
%     legend({'\mu = 0', '\mu = 1', '\mu = \infty'},'Interpreter','tex')
%     xlabel('Subspace rank $V$', 'Interpreter', 'latex')
%     ylabel('nRE (dB)')
%     title(figname)




%     figname = ['exp1_mu_vs_SDE_zone' num2str(sridx) ''];
%     figure('Name',figname)
%     semilogx(mucan,10*log10(sde([1,2,4,8,16],:)'))
%     % xlim([1 16])
%     % ylim([2 18])
%     legend({'V = 1', 'V = 2', 'V = 4', 'V = 8', 'V = 16'})
%     xlabel('\mu', 'Interpreter', 'tex')
%     ylabel('SDE (dB)')
%     title(figname)
% 
%     % %%
%     figname = ['exp1_rank_vs_RE_zone' num2str(sridx) ''];
%     figure('Name',figname)
%     plot(Vcan,10*log10(re(:,[1,51,lm])))
%     xlim([1 16])
%     % ylim([-50 -20])
%     legend({'\mu = 0', '\mu = 1', '\mu = \infty'},'Interpreter','tex')
%     xlabel('Subspace rank $V$', 'Interpreter', 'latex')
%     ylabel('RE (dB)')
%     title(figname)
% 
%     figname = ['exp1_mu_vs_RE_zone' num2str(sridx) ''];
%     figure('Name',figname)
%     semilogx(mucan,10*log10(re([1,2,4,8,16],:)'))
%     % xlim([1 16])
%     % ylim([-50 -20])
%     legend({'V = 1', 'V = 2', 'V = 4', 'V = 8', 'V = 16'})
%     xlabel('\mu', 'Interpreter', 'tex')
%     ylabel('RE (dB)')
%     title(figname)
end
end