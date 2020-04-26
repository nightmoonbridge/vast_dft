function irVirsrc = generateImpulseResponseVirtualSource(general, array, zone, idx)
% GENERATEIMPULSERESPONSEVIRTUALSOURCE returns the impulse response
% functions for the virtual source. The order of the bright and dark zones 
% are included. The desired sound field should be the input signal itself 
% in the free field condition.

c = general.c;
if general.sfactor == 1
    fs = general.fs;
    K = general.lenIResponse;
else
    fs = general.tfs;
    K = general.lenRir*general.sfactor;
end
M = zone.numCtrPts;

rv = sqrt(sum((zone.geometry(idx).ctrPtsPositions-array.virsrcPosition(idx,:)).^2,2));

vidx = round(fs*rv/c);
% vidx = floor(fs*rv/c);
va = 1./(4*pi*rv);
tVa = zeros(K,M);
for ii = 1:zone.numCtrPts
    tidx = vidx(ii);
    defpt = zeros(K,1);
    defpt(1) = 1;
    
    % This index shifting is changed to tidx from tidx-1.
    % This is to match the index from McRoomSim and RIRgenerator.
    defpt = circshift(defpt,tidx);
    tVa(:,ii) = defpt*va(ii);
end
irVirsrc = tVa;


end