classdef pfm_mtx_mse < handle
% Previously an object named 'PerformanceMetric' was used.
% This is a modified version of the object 'PerformanceMetric' in order to
% the journal paper in preparation.

    properties
        % Prior acoustic contrast
        priorAC = 0;
        
        % Posterior acoustic contrast
        postAC = 0;
        
        % Signal distortion in the bright zone
        % this is not only a performance metric but also a MSE
        sde = 0;
        
        % Normalized signal distortion in the bright zone
        nsde = 0;
        
        % Residual energy in the dark zone
        resiEner = 0;
        
        % Normalized residual energy in the bright zone
        nre = 0;
        
        
        % without any control
        iRbi = 0;
        iRdi = 0;
    end
    methods
        function getpriorAC(obj,Rb,Rd,dboption)
            if nargin < 4
                dboption = false;
            end
            n = size(Rb,1);
            ovt = ones(n,1);
            obj.priorAC = (abs(ovt'*Rb*ovt)+eps)/(abs(ovt'*Rd*ovt)+eps);
            if dboption
                obj.priorAC = 10*log10(obj.priorAC);
            end
        end
        
        function getpostAC(obj,Rb,Rd,q,dboption)
            if nargin < 5
                dboption = false;
            end
            obj.postAC = (abs(q'*Rb*q)+eps)/(abs(q'*Rd*q)+eps);
            if dboption
                obj.postAC = 10*log10(obj.postAC);
            end
        end
        
        function getnsde(obj,Rb,Hb,hz,q)
            obj.nsde = real((q'*Rb*q - 2*q'*Hb'*hz + hz'*hz + eps)/(hz'*hz + eps));
        end

        function getsde(obj,Rb,Hb,hz,q)
            obj.sde = real(q'*Rb*q - 2*q'*Hb'*hz + hz'*hz + eps);
        end
        
        function getre(obj,Rd,q)
            obj.resiEner = real(q'*Rd*q);
        end
        
        function getnre(obj,Rd,q)
            n = size(Rd,1);
            ovt = ones(n,1);
            obj.nre = real((abs(q'*Rd*q) + eps)/(abs(ovt'*Rd*ovt) + eps));
        end

        function getiRi(obj,R,isbright)
            n = size(R,1);
            ovt = ones(n,1);
            iRi = abs(ovt'*R*ovt)+eps;
            if isbright
                obj.iRbi = iRi;
            else
                obj.iRdi = iRi;
            end
        end
    end

end