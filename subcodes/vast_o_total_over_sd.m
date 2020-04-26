function result=vast_o_total_over_sd(RX,RN,r,dsquare,ranks,K_SET)
[U,EIG_VALUE]=jdiag(RX,RN,'vector');
EIG_VALUE_sub_vector=EIG_VALUE(1:ranks);
c=U(:,1:ranks)'*r;
U_part=U(:,1:ranks);

% minimum_value=max(-sum(c.^2.*(1./EIG_VALUE_sub_vector))+dsquare,eps);

normalized_minimum_value=max((-sum(c.^2.*(1./EIG_VALUE_sub_vector))+dsquare)/dsquare,1e-3);

db_min=10*log10(normalized_minimum_value); 

% sd_set=linspace(db_min+1,-3,10);
sd_set=linspace(-10,-3,10);

K_SET=10.^(sd_set/10)*dsquare;

for var_cnt=1:length(K_SET)
    K=K_SET(var_cnt);
    if normalized_minimum_value*dsquare>K   %% the constraints can be satisfied or not
        conFilter_RRCG(:,var_cnt)=nan(size(r));
        mu_RRCG(var_cnt)=nan;
    else
        
        if K>=dsquare  %% contains the minimum or not
            conFilter_RRCG(:,var_cnt)=zeros(size(r));
            mu_RRCG(var_cnt)=nan;
        else
            x0=0;
            f=@(x) sum(EIG_VALUE_sub_vector.*c.^2./(EIG_VALUE_sub_vector+x).^2) ...
            -2*sum(c.^2./(EIG_VALUE_sub_vector+x))+dsquare-K;
            g=@(x) -2*sum(EIG_VALUE_sub_vector.*c.^2./(EIG_VALUE_sub_vector+x).^3) ...
                +2*sum(c.^2./(EIG_VALUE_sub_vector+x).^2);
            
            x=newton_root_sd(f,g,x0);
            mu=x;
            
            
            conFilter = U_part*([1./(mu+EIG_VALUE_sub_vector)].*c);
            
            conFilter_RRCG(:,var_cnt) = conFilter;
            mu_RRCG(var_cnt) = mu;
        end
        
        
    end
    
    
    %     acc_RRCG(sigma_square)=db10(acc_measure(conFilter));
    %     sd_RRCG(sigma_square)=sd_measure(conFilter);
    %     pow_RRCG(sigma_square)=conFilter'*RN*conFilter;
end
result.conFilter_RRCG=conFilter_RRCG;
result.mu_RRCG=mu_RRCG;
% result.acc_RRCG=acc_RRCG;
% result.sd_RRCG=sd_RRCG;
% result.pow_RRCG=pow_RRCG;


end



function x0=newton_root_sd(f,fprime,x0,maxIterations,tolerance)
if nargin<=3
    maxIterations=1e4;
    tolerance=1e-7;
end

for i = 1 : maxIterations
    
    y = f(x0);
%     if (y>=0)
%         break;
%     end
    
    
    yprime = fprime(x0);
    
    if(abs(yprime) < eps) %Don't want to divide by too small of a number
        % denominator is too small
        break; %Leave the loop
    end
    
    x1 = x0 - y/yprime; %Do Newton's computation
    
    if(norm(x1 - x0)/norm(x0) <= tolerance) %If the result is within the desired tolerance
        break; %Done, so leave the loop
    end
    
    x0 = x1; %Update x0 to start the process again
    
end
end
