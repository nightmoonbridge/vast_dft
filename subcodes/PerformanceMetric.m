classdef PerformanceMetric < handle
    properties
        % Prior acoustic contrast
        priorAC = 0;
        
        % Initial field correction
        initFC = 0;
        
        % Posterior field correction
        postFC = 0;
        
        % Initial energy ratio between the desired and uncontrolled
        % sound fields in the bright zone
        initER = 0;        
        
        % Posterior acoustic contrast
        postAC = 0;
        
        % Energy reduction in the dark zone
        egred_D = 0;
        
        % Energy reduction in the bright zone
        egred_B = 0;
        
        % Signal reduction in the bright zone
        sgred_B = 0;
        
        % Normalized signal distortion in the bright zone
        nsde = 0;
        
        % Signal distortion in the bright zone
        sde = 0;
        
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
        
        function getinitFC(obj,Hb,hz)
            n = size(Hb,2);
            ovt = ones(n,1);
            obj.initFC = ((Hb*ovt)'*hz)/(hz'*hz);
        end
        
        function getpostFC(obj,Hb,hz,q)
            obj.postFC = ((Hb*q)'*hz)/(hz'*hz);
        end
        
        function getinitER(obj,Rb,hz)
            n = size(Rb,1);
            ovt = ones(n,1);
            obj.initER = real((hz'*hz)/(ovt'*Rb*ovt));
        end
        
        function getegred(obj,R,q,isbright)
            n = size(R,1);
            ovt = ones(n,1);
            ratio = (abs(ovt'*R*ovt)+eps)/(abs(q'*R*q)+eps);
            if isbright
                obj.egred_B = ratio;
            else
                obj.egred_D = ratio;
            end
        end
        
        function getsgred(obj,Rb,hz,q)
            obj.sgred_B = real((hz'*hz+eps)/(q'*Rb*q+eps));
        end
        
        function getnsde(obj,Rb,Hb,hz,q)
            obj.nsde = real((q'*Rb*q - 2*q'*Hb'*hz + hz'*hz + eps)/(hz'*hz + eps));
        end

        function getsde(obj,Rb,Hb,hz,q)
            obj.sde = real(q'*Rb*q - 2*q'*Hb'*hz + hz'*hz + eps);
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