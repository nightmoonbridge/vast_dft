function [spaCorrMtx, varargout] = getSpaCorrInfoLS_fast_version6(rir_romm, rir_virtual,numLoudspk,J,K, M,zonetype,Rupper)
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
% 
% narginchk(5,5)
% varargout = cell(nargout-1,1);
% nargoutchk(1,3)

% [data_len, ~, numLoudspk] = size(wUctrReprodSoundField);
% N=length(inputdata);
R_dim=K+J-1;

% ym = zeros(general.lenConFilter*numLoudspk,2*general.lenSegment);
switch lower(zonetype)
    case {'bright', 'br'}
        % normalization factor
        nfactor = 1/M;
        

        R_Y2=0;
        r_final=0;
        

        
        cat_rir_virtual=[rir_virtual;zeros(K+J-1-size(rir_virtual,1),size(rir_virtual,2))];        
        Kappa=0;
        for ii=1:M
            %     tic
            U=zeros(numLoudspk*J,(K+J-1)); r_value=zeros(numLoudspk*J,1); U2=zeros(numLoudspk*J,(K+J-1));
            a_m_set=reshape(rir_romm(:,ii),K,numLoudspk);
                        
            dx=Rupper*cat_rir_virtual(:,ii);
            
            
            Kappa=Kappa+cat_rir_virtual(:,ii)'*dx;
            for ll=1:numLoudspk
                row_inx=(ll-1)*J+1:(ll)*J;
                
 
                hRIR=toeplitz([a_m_set(1,ll);zeros(J-1,1)],[a_m_set(:,ll);zeros(J-1,1)]);
                Y=hRIR*Rupper;

                r_part=hRIR*dx;
                
                
                U(row_inx,:)=Y;
                U2(row_inx,:)=hRIR;
                r_value(row_inx)=r_part;
            end
%             tic
            R_Y2=R_Y2+U*U2';
%             toc
            
            
            r_final=r_final+r_value;
            
        end
         

        spaCorrMtx = nfactor*R_Y2;
        spaCorrVec = nfactor*r_final;
        sigDistort = nfactor*Kappa;
        varargout{1} = spaCorrVec;
        varargout{2} = sigDistort;
%         varargout{3} = filtered_d;
        
    case {'dark', 'dk'}
        % normalization factor
       % normalization factor
        nfactor = 1/M;

        R_Y2=0;
        
        
        for ii=1:M
            %     tic
            U=zeros(numLoudspk*J,(K+J-1));  U2=zeros(numLoudspk*J,(K+J-1));
            a_m_set=reshape(rir_romm(:,ii),K,numLoudspk);
            for ll=1:numLoudspk
                row_inx=(ll-1)*J+1:(ll)*J;
                hRIR=toeplitz([a_m_set(1,ll);zeros(J-1,1)],[a_m_set(:,ll);zeros(J-1,1)]);
                Y=hRIR*Rupper;
                
                U(row_inx,:)=Y;
                U2(row_inx,:)=hRIR;
            end
            
            R_Y2=R_Y2+U*U2';

            
        end
        spaCorrMtx = nfactor*R_Y2;

end