function [mergedCtrReprodSoundField, gmax] = getMergedCtrReprodSoundField(ctrReprodSoundField)
% function mergedCtrReprodSoundField = getMergedCtrReprodSoundField(ctrReprodSoundField)
%  This function merges the bright and dark zones at the given zones for
%  the different audio inputs. 
%  For example, there are two zones, labeled zone A and zone B, 
%  the reproduced sound field at zone A is the sum between the followings:
%    the bright zone for input A + the dark zone for input B
%  The reproduced sound field at zone B also can be done in the same manner
%
%  > Input variable
%   ctrReprodSoundField : a structure, there are two fields: 'Br', 'Dk'
%
%  > Output variable
%   mergedCtrReprodSoundField : a cell array with the number of zones
%   gmax                      : the maximum value out of all data points
%
% code: 22nd/July-2019
%

if ~isstruct(ctrReprodSoundField)
    error('Invalid structure for the controlled reproduced sound field')
end
[nzone, ~] = size(ctrReprodSoundField.Br);
[nsample, nctrpts] = size(ctrReprodSoundField.Br{nzone});

ssidx = flipud(perms(1:nzone));
mergedCtrReprodSoundField = cellfun(@(x) zeros(nsample, nctrpts), ...
    cell(nzone,1), 'UniformOutput', false);

maxdatapoint = zeros(nctrpts,nzone);

for sridx = 1:nzone
    for midx = 1:nctrpts
        mergedCtrReprodSoundField{sridx}(:,midx) = ctrReprodSoundField.Br{ssidx(sridx,1)}(:,midx)+ctrReprodSoundField.Dk{ssidx(sridx,2)}(:,midx);
        maxdatapoint(midx,sridx) = max(abs(mergedCtrReprodSoundField{sridx}(:,midx)));
    end
end
gmax = max(max(maxdatapoint));

end