function [pmsolution, accsolution, accpmsolution, oSpaCorrMtx, oSpaCorrVec, sigDistort, eveMtx, evaVec] = ...
    getVASTcontrolfilter(general, zone, array, irArray, irVirsrc, regnum, vastoption, ctrfilter)
narginchk(5,8)

if nargin < 7
    vastoption = true;
    if nargin < 6
        regnum = 0;
    end
end

if vastoption
    general.lenInput = min(general.lenInput, general.lenConFilter - 1 + ...
        max(general.lenIResponse, general.lenRir));
end

[oSpaCorrMtx, oSpaCorrVec, sigDistort] = getVAST_spaCorrInfo(general, zone, array, irArray, irVirsrc);

accsolution = cell(zone.number,1);
accsolution{1} = (oSpaCorrMtx{1}+oSpaCorrMtx{2})\oSpaCorrVec{1};
accsolution{2} = (oSpaCorrMtx{1}+oSpaCorrMtx{2})\oSpaCorrVec{2};

Vmax = ctrfilter.Vmax;

if strcmpi(ctrfilter.cvxopt_properties.opttype, 'min_sb')
    minSb = true;
else
    minSb = false;
end

pmsolution = cell(zone.number,1);
% accsolution = cell(zone.number,1);
accpmsolution = cell(zone.number,1);

ssidx = flipud(perms(1:zone.number));

vast_pm_cnst = zeros(Vmax,1);
vast_accpm_cnst = zeros(Vmax,1);

taroption.broadband = true;
taroption.nloudspks = array.numLoudspk;
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

taroption.frequency = false;

% Almost the same control filters except for the amplitude of ACC,
% But the VAST framework is faster than traditional since it calculates
% the joint diagonalization once for each scenario
if vastoption
    
    eveMtx = cellfun(@(x) zeros(general.lenConFilter*array.numLoudspk), ...
        cell(zone.number, 1), 'UniformOutput', false);
    evaVec = cellfun(@(x) zeros(general.lenConFilter*array.numLoudspk,1), ...
        cell(zone.number, 1), 'UniformOutput', false);

    % VAST framework
    for sridx = 1:zone.number
        [eveMtx{sridx}, evaVec{sridx}] = jdiag(oSpaCorrMtx{ssidx(sridx,1)}, ...
            oSpaCorrMtx{ssidx(sridx,2)}, 'vector', true);
        
        V = ctrfilter.V;
        if ctrfilter.cvxopt_properties.findopt
            D = evaVec{sridx};
            U = eveMtx{sridx};
            rb = oSpaCorrVec{sridx};
            d_part = real(D(1:V));
            U_part = U(:,1:V);
            uvrb = U_part'*rb;
            hz = sqrt(sigDistort(sridx));
            ctrfilter = getOptPara(uvrb,d_part,hz,ctrfilter,[1, sridx]);
            
            
            mu_new = ctrfilter.cvxopt_properties.optpara(1, sridx);
            ctrfilter.mu = mu_new;
            
%         if minSb %strcmpi(ctrfilter.cvxopt_properties.const,'sd')
            accpmsolution{sridx} = U_part*(uvrb./(mu_new + d_part));
%         else
%             accpmsolution{sridx} = U_part*(mu_new*uvrb./(1 + mu_new*d_part));
%         end

        else
            mu_new = ctrfilter.mu;
            % Calculation of coefficients
            for vidx = 1:Vmax
                vast_pm_cnst(vidx) = (eveMtx{sridx}(:,vidx)'*oSpaCorrVec{sridx})/(1 + evaVec{sridx}(vidx));
                vast_accpm_cnst(vidx) = (eveMtx{sridx}(:,vidx)'*oSpaCorrVec{sridx})/(mu_new + evaVec{sridx}(vidx));
            end
%             accsolution{sridx} = eveMtx{sridx}(:,1);
            pmsolution{sridx} = eveMtx{sridx}*vast_pm_cnst;
            accpmsolution{sridx} = eveMtx{sridx}(:,1:V)*vast_accpm_cnst(1:V);
        end
    end
else
    mu = ctrfilter.mu;
    kappa = mu/(1+mu);
    % traditional approach
    pmsolution{1} = (oSpaCorrMtx{1} + oSpaCorrMtx{2} + regnum*eye(array.numLoudspk*general.lenConFilter))\oSpaCorrVec{1};
    pmsolution{2} = (oSpaCorrMtx{1} + oSpaCorrMtx{2} + regnum*eye(array.numLoudspk*general.lenConFilter))\oSpaCorrVec{2};
    
    [eveMtx,~] = jdiag(oSpaCorrMtx{1}, oSpaCorrMtx{2} + regnum*eye(array.numLoudspk*general.lenConFilter), 'vector');
    accsolution{1} = eveMtx(:,1);
    [eveMtx,~] = jdiag(oSpaCorrMtx{2}, oSpaCorrMtx{1} + regnum*eye(array.numLoudspk*general.lenConFilter), 'vector');
    accsolution{2} = eveMtx(:,1);
    
    accpmsolution{1} = ((1-kappa)*oSpaCorrMtx{1} + kappa*oSpaCorrMtx{2} + regnum*eye(array.numLoudspk*general.lenConFilter))\((1-kappa)*oSpaCorrVec{1});
    accpmsolution{2} = ((1-kappa)*oSpaCorrMtx{2} + kappa*oSpaCorrMtx{1} + regnum*eye(array.numLoudspk*general.lenConFilter))\((1-kappa)*oSpaCorrVec{2});
end

end