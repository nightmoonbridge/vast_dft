function irArray = generateImpulseResponseArray(general, array, zone, room, fullflag, idx)
% GENERATEIMPULSERESPONSE returns the impulse responses according to the
% given setup. The variable "zone" is updated.
%%%%%% New version from 1st/Feb/2018 %%%%%  Start

if fullflag
    fs = general.tfs;
    if general.setup.simrir
        lenir = general.lenIResponse;
    else
        lenir = general.lenRir;
    end
    lenir = general.sfactor*lenir;
else
    fs = general.fs;
    lenir = general.lenIResponse;
end

irArray = zeros(array.numLoudspk*lenir, zone.numCtrPts);
if ~room.reverberation
    % No reverberation is considered.
    rl = zeros(zone.numCtrPts,array.numLoudspk);
    for ii = 1:array.numLoudspk
        rl(:,ii) = sqrt((zone.geometry(idx).ctrPtsPositions(:,1) - ...
            array.loudspkPositions(ii,1)).^2 ...
            + (zone.geometry(idx).ctrPtsPositions(:,2) - ...
            array.loudspkPositions(ii,2)).^2);
    end
    
    %%% Index allocation
    pidx = round(fs*rl/general.c);
    denominator = 1./(4*pi*rl);
    for ii = 1:zone.numCtrPts
        tidx = pidx(ii,:) + (0:array.numLoudspk-1)*lenir;
        for jj = 1:array.numLoudspk
            defpt = zeros(lenir,1);
            defpt(1) = 1;
            
            % This index shifting is changed to tidx(jj) from tidx(jj)-1.
            % This is to match the index from McRoomSim and RIRgenerator.
            defpt = circshift(defpt,tidx(jj));
%             defpt = circshift(defpt,tidx(jj)-1);
            irArray((jj-1)*lenir+1:jj*lenir,ii) = defpt*denominator(ii,jj);
        end
    end
else
    % Reverberation is considered.
    
    % Since rir_generator is the most reliable implementation for now,
    % rir-generator is used to generate room impulse responses.
    % rir-generator can be found from the following link:
    % https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator
    roominfo = [room.xlength,room.ywidth,room.zheight];
%     nsamples = 0.5*fs;
    tfs = 48000;
    
    for lidx = 1:array.numLoudspk
        mic = [zone.geometry(idx).ctrPtsPositions, ...
            array.height*ones(zone.numCtrPts,1)];
        src = [array.loudspkPositions(lidx,:), array.height];
        
        % Generation of RIR (room impulse response)
        %     h = rir_generator(general.c, fs, mic, src, roominfo, room.T60, nsamples);
        %     h = h';
        
        h48k = rir_generator(general.c, tfs, mic, src, roominfo, room.T60, tfs*0.5);
        h48k = h48k';
        if fullflag
            hrsmp = h48k;
        else
            h48k = h48k(1:tfs*(lenir/fs),:);
            hrsmp = (tfs/fs)*resample(h48k,fs,tfs);
        end
        
        irArray((lidx-1)*lenir+1:lidx*lenir,:) = hrsmp(1:lenir,:);
    end
end


end