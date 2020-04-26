function [oSpaCorrMtx, oSpaCorrVec, sigDistort] = getVAST_spaCorrInfo(general, zone, array, irArray, irVirsrc)
% For ICASSP 2018 implementation, see getVAST_spaCorrInfo_old.m
narginchk(5,5)
nargoutchk(2,3)

% The input signal x[n] is the Kronecker delta function
% xin = [1;zeros(8*general.lenSegment-1,1)];

% xin = [1;zeros(2*general.lenSegment-1,1)];  % original

xin = [1;zeros(general.lenSegment-1,1)];
% rng default
% xin = randn(general.lenSegment,1);
xin = repmat(xin,2,1);

oSpaCorrMtx = cellfun(@(x) zeros(general.lenConFilter*array.numLoudspk),...
                        cell(zone.number,1),'UniformOutput',false);
oSpaCorrVec = cellfun(@(x) zeros(general.lenConFilter*array.numLoudspk,1),...
                        cell(zone.number,1),'UniformOutput',false);
sigDistort = zeros(zone.number,1);

% if general.setup.simrir
%     lenir = general.lenIResponse;
% else
%     lenir = general.lenRir;
% end
lenir = size(irVirsrc{1},1);
                    
for sridx = 1:zone.number
    y = zeros(length(xin)+general.lenPredata-1, zone.numCtrPts,array.numLoudspk);
%     y = zeros(length(xin)+general.lenPredata+general.lenConFilter-1, zone.numCtrPts,array.numLoudspk);
    
    %     y = zeros(2*general.lenSegment+general.lenPredata+general.lenConFilter-1, zone.numCtrPts,array.numLoudspk);
    d = y(:,:,1);
    for midx = 1:zone.numCtrPts
        nidx = 1:lenir;
        for lidx = 1:array.numLoudspk
            % Take each ml impulse response function
            hml = irArray{sridx}(nidx,midx);
            
%             y(1+general.lenPredata+general.lenConFilter-1:end, midx, lidx) = filter(hml,1,xin);
            xx = filter(hml,1,xin);
            y(:, midx, lidx) = xx(1:end-1);
            nidx = nidx + lenir;
        end
        hmz = irVirsrc{sridx}(:,midx);
        
        xy = filter(hmz,1,xin);
        d(:, midx) = xy(1:end-1);
%         d(1+general.lenPredata+general.lenConFilter-1:end, midx) = filter(hmz,1,xin);
    end
    [oSpaCorrMtx{sridx}, oSpaCorrVec{sridx}, sigDistort(sridx)] = getSpaCorrInfo(general, zone, y, d, 'bright', 'vast');
end

end