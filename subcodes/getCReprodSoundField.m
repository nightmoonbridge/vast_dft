function controlledReprodSoundField = getCReprodSoundField ...
    (general, array, zone, fxInputSignals, irMeasured, onlycenter)
% GETCREPRODSOUNDFIELD returns the controlled reproduced sound field for
% each input signal. "controlledReprodSoundField" is a structure whose
% fields are ~.Br and ~.Dk. They consist of two cell arrays, respectively.
% Each cell array is a 3D array such that
%    X: the length of input source
%    Y: the number of control points
%
%
numPts = size(irMeasured{1},2);
cenPt = round(median(1:numPts));
nzone = zone.number;
numPtsBr = numPts;
numPtsDk = (nzone-1)*numPts;
controlledReprodSoundField.Br = cellfun(@(x) zeros(general.lenInput, ...
    numPtsBr), cell(zone.number,1), 'UniformOutput', false);
controlledReprodSoundField.Dk = cellfun(@(x) zeros(general.lenInput, ...
    numPtsDk), cell(zone.number,1), 'UniformOutput', false);

lenir = size(irMeasured{1},1)/array.numLoudspk;

tidx = 1:zone.number;

for sridx = 1:zone.number                                    % source index
    bridx = sridx;
    dkidx = tidx(tidx ~= bridx);
    irBr = horzcat(irMeasured{bridx});                 % in the time domain
    irDk = horzcat(irMeasured{dkidx});                 % in the time domain
    zonedata = fxInputSignals{sridx};                  % in the time domain
    
    if onlycenter
        for lidx = 1:array.numLoudspk                   % loudspeaker index
            irbridx = (1:lenir) + (lidx-1)*lenir;
            emittedSignalBr = filter(irBr(irbridx,cenPt),1,zonedata(:,lidx));
            
            controlledReprodSoundField.Br{sridx}(:,cenPt) = ...
                controlledReprodSoundField.Br{sridx}(:,cenPt) + emittedSignalBr;
            
            irdkidx = (1:lenir) + (lidx-1)*lenir;
            emittedSignalDk = filter(irDk(irdkidx,cenPt),1,zonedata(:,lidx));
            
            controlledReprodSoundField.Dk{sridx}(:,cenPt) = ...
                controlledReprodSoundField.Dk{sridx}(:,cenPt) + emittedSignalDk;
        end
    else
        for mbridx = 1:numPtsBr                           % control point index
            for lidx = 1:array.numLoudspk                   % loudspeaker index
                irbridx = (1:lenir) + (lidx-1)*lenir;
                emittedSignalBr = filter(irBr(irbridx,mbridx),1,zonedata(:,lidx));

                controlledReprodSoundField.Br{sridx}(:,mbridx) = ...
                    controlledReprodSoundField.Br{sridx}(:,mbridx) + emittedSignalBr;
            end
        end

        for mdkidx = 1:numPtsDk                           % control point index
            for lidx = 1:array.numLoudspk                   % loudspeaker index
                irdkidx = (1:lenir) + (lidx-1)*lenir;
                emittedSignalDk = filter(irDk(irdkidx,mdkidx),1,zonedata(:,lidx));

                controlledReprodSoundField.Dk{sridx}(:,mdkidx) = ...
                    controlledReprodSoundField.Dk{sridx}(:,mdkidx) + emittedSignalDk;
            end
        end
    end
end

end