function [spaCorrMtx, varargout] = getSpaCorrInfo(general, zone, wUctrReprodSoundField, wDesiredSoundField, zonetype, vasttype)
% GETSPACORRINFO returns the spatial correlation matrix and vector.
% It depends on how one defines the variable zonetype:
%   zonetype == 'bright' or 'br'
%     > spaCorrMtx (br)
%     > spaCorrVec (br)
%   zonetype == 'dark' or 'dk'
%     > spaCorrMtx (dk)
%
% Common variables
%   - general
%   - zone
% For the bright zone
%   - wUctrReprodSoundField:   (2N+J-1) x Mb x L
%   - wDesiredSoundField:      (2N+J-1) x Mb
% For the dark zone
%   - wUctrReprodSoundField:   (2N+J-1) x Md x L
%
% The output arguments are as follows:
% For the bright zone
%   - spaCorrMtx:              LJ x LJ
%   - spaCorrVec:              LJ x 1
%   - sigDistort:              scalar
% For the dark zone
%   - spaCorrMtx:              LJ x LJ
%
% This function can be tested by the following syntax:
%   results = run(test_getSpaCorrInfo)
%
% See test_getSpaCorrInfo.m
%

narginchk(5,6)
if nargin < 6
    vasttype = 'apvast';
end
varargout = cell(nargout-1,1);
nargoutchk(1,3)

[ndata, ~, numLoudspk] = size(wUctrReprodSoundField);

if strcmpi(vasttype, 'apvast')
    datalength = 2*general.lenSegment + general.lenPredata;
    ym = zeros(general.lenConFilter*numLoudspk,datalength);
elseif strcmpi(vasttype, 'pvast')
    datalength = general.lenInput;
    ym = zeros(general.lenConFilter*numLoudspk,datalength);
else % vast, need to check -- 4/oct-2019
    datalength = general.lenSegment + general.lenPredata;
%     datalength = 2*general.lenSegment + general.lenPredata;
    ym = zeros(general.lenConFilter*numLoudspk,datalength);
end

% ym = zeros(general.lenConFilter*numLoudspk,2*general.lenSegment);
switch lower(zonetype)
    case {'bright', 'br'}
        % normalization factor
%         nfactor = 1/(ndata*zone.numCtrPtsBr);
%         nfactor = 2/ndata;
        nfactor = 1;
%         nfactor_rb = 1/(datalength*zone.numCtrPtsBr);
        
        % Spatial correlation matrix
        spaCorrMtx = zeros(general.lenConFilter*numLoudspk);
        
        % Spatial correlation vector
        spaCorrVec = zeros(general.lenConFilter*numLoudspk,1);
        
        % Signal distortion energy
        sigDistort = 0;
        
        for mbridx = 1:zone.numCtrPtsBr
            %%% length of data
            %%%   indata:   (2N+J-1) x L
            %%%   dmdata:   (2N+J-1) x 1
            indata = squeeze(wUctrReprodSoundField(:,mbridx,:));
            %             dmdata = wDesiredSoundField(:,mbridx);
            dmdata = wDesiredSoundField(general.lenConFilter:end,mbridx);
            jidx = 1:general.lenConFilter;
            for lidx = 1:numLoudspk
                datain = indata(:,lidx);
                yml = toeplitz(flip(datain(1:general.lenConFilter)),datain(general.lenConFilter:end));
%                 cdata = convmtx(indata(:,lidx)',general.lenConFilter);
%                 yml = cdata(:,general.lenConFilter:2*general.lenSegment+general.lenConFilter-1);
                
                ym(jidx,:) = yml;
                
                jidx = jidx + general.lenConFilter;
                clear datain
            end
            spaCorrMtx = spaCorrMtx + ym*ym';
            spaCorrVec = spaCorrVec + sum(ym.*dmdata',2);
            sigDistort = sigDistort + sum(dmdata.^2);
        end
        spaCorrMtx = nfactor*spaCorrMtx;
        spaCorrVec = nfactor*spaCorrVec;
        sigDistort = nfactor*sigDistort;
        varargout{1} = spaCorrVec;
        varargout{2} = sigDistort;
        
    case {'dark', 'dk'}
        % normalization factor
        nfactor = 1/(ndata*zone.numCtrPtsDk);
        
        % Spatial correlation matrix
        spaCorrMtx = zeros(general.lenConFilter*numLoudspk);
        
        for mdkidx = 1:zone.numCtrPtsDk
            indata = squeeze(wUctrReprodSoundField(:,mdkidx,:));   % (2N+J-1) x L
            jidx = 1:general.lenConFilter;
            for lidx = 1:numLoudspk
                datain = indata(:,lidx);
                yml = toeplitz(flip(datain(1:general.lenConFilter)),datain(general.lenConFilter:end));
%                 cdata = convmtx(indata(:,lidx)',general.lenConFilter);
%                 yml = cdata(:,general.lenConFilter:2*general.lenSegment+general.lenConFilter-1);
                
                ym(jidx,:) = yml;
                
                jidx = jidx + general.lenConFilter;
                clear datain
            end
            spaCorrMtx = spaCorrMtx + ym*ym';
        end
        spaCorrMtx = nfactor*spaCorrMtx;
end
end
