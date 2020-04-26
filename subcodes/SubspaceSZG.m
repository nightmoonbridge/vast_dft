classdef SubspaceSZG< handle
% https://github.com/LimingShi/Subspace-Sound-Zones-Gneration
    properties (SetAccess=immutable) 
%         N
%         L
%         F
%         pitchBoundsOuter = [0.0 0.5] % The outer bounds that are acceptable
%         
%         % default values
%         epsilon = 0.0
%         dcIsIncluded = false
%         epsilon_ref = 0.0
%         
%         % Precomputed quantities
%         crossCorrelationVectors
%         Gamma1
%         Gamma2
%         fftShiftVector
%          uncontroled_dk_energy
    end
    
    properties (SetAccess=private)
        irArray
        irVirsrc
        fs
        K
        M1
        M2
        L
        q_uncontrol;
        inputpower;
        uncontroled_dk_energy=cell(1,2);
        spaCorrMtxBr=cell(1,2);
        spaCorrVecBr=cell(1,2);
        spaCorrMtxDk=cell(1,2);
        sigDistort=zeros(1,2);
        U_total_gevd=cell(1,2);
        EIG_VALUE_total_gevd=cell(1,2);
        U_total_cg=cell(1,2);
        EIG_VALUE_total_cg=cell(1,2);      
        mu_LB_ac=nan(1,2);
        mu_UB_ac=nan(1,2);
        mu_LB_sd=nan(1,2);
        mu_UB_sd=nan(1,2);
        mu_LB_er=nan(1,2);
        mu_UB_er=nan(1,2);
        min_ac_db=nan(1,2);
        max_ac_db=nan(1,2);
        min_sd_db=nan(1,2);
        max_sd_db=nan(1,2);
        min_er_db=nan(1,2);
        max_er_db=nan(1,2);
        FILTER_COEFFS=cell(1,2);
        measured_ac;
        measured_sd;
        measured_tir;
        measured_er
        ReprodSoundField_FINAL;
        
    end

    methods

        
        function obj = SubspaceSZG(irVirsrc,irArray,FS)
%           
            obj.irArray=irArray;
            obj.irVirsrc=irVirsrc;
            obj.fs=FS;
            obj.K=size(irVirsrc{1},1);
            obj.M1=size(irArray{1},2);
            obj.M2=size(irArray{2},2);
            obj.L=round(size(irArray{1},1)/size(irVirsrc{1},1));                
        end
        
        function computeStatistics(obj,J,indata,method_used)                   
            for sridx = 1:2                                
                %% data_autocorrelation matrix
                if strcmp(method_used,'DATA')
                    clm_values=[indata{sridx}.xin(1);zeros(obj.K+J-2,1)];
                    row_values=[indata{sridx}.xin];
                    [R]=fast_autocorrelation_R(clm_values,row_values);
                    obj.inputpower(sridx)=R(1,1);
                
                end
                if strcmp(method_used,'WHITE')
                    R=eye(obj.K+J-1);              
                    obj.inputpower(sridx)=R(1,1);
                end
                
                
                
                %% spatial-transformation to obtain spatial autocorrelation
                if sridx==1
                    [obj.spaCorrMtxBr{sridx}, obj.spaCorrVecBr{sridx}, obj.sigDistort(sridx)]=getSpaCorrInfoLS_fast_version6...
                        (obj.irArray{sridx},obj.irVirsrc{sridx}, obj.L,J,obj.K,obj.M1,'bright',R);
                    
                    obj.spaCorrMtxDk{sridx} = getSpaCorrInfoLS_fast_version6(obj.irArray{mod(sridx,2)+1}, [],obj.L,J,obj.K,obj.M2, 'dark',R);
                else
                    [obj.spaCorrMtxBr{sridx}, obj.spaCorrVecBr{sridx}, obj.sigDistort(sridx)]=getSpaCorrInfoLS_fast_version6...
                        (obj.irArray{sridx},obj.irVirsrc{sridx}, obj.L,J,obj.K,obj.M2,'bright',R);
                    
                    obj.spaCorrMtxDk{sridx} = getSpaCorrInfoLS_fast_version6(obj.irArray{mod(sridx,2)+1}, [],obj.L,J,obj.K,obj.M1, 'dark',R);
                end
                
                obj.q_uncontrol=kron(ones(obj.L,1),[1;zeros(J-1,1)]);
                
                if strcmp(method_used, 'DATA')
                    if sridx==1
    %                     [NoControlSoundField.Br{sridx}]=get_produced_x_signal(indata{sridx},obj.q_uncontrol,obj.irArray{sridx},obj.M1,obj.L,obj.K);
                        [NoControlSoundField.dk{sridx}]=get_produced_x_signal(indata{sridx},obj.q_uncontrol,obj.irArray{mod(sridx,2)+1},obj.M2,obj.L,obj.K);
                        obj.uncontroled_dk_energy{sridx}=norm(NoControlSoundField.dk{sridx},'fro').^2;
    %                     obj.uncontroled_dk_energy{sridx}=norms(NoControlSoundField.dk{sridx}).^2;
                    else
    %                     [NoControlSoundField.Br{sridx}]=get_produced_x_signal(indata{sridx},obj.q_uncontrol,obj.irArray{sridx},obj.M2,obj.L,obj.K);
                        [NoControlSoundField.dk{sridx}]=get_produced_x_signal(indata{sridx},obj.q_uncontrol,obj.irArray{mod(sridx,2)+1},obj.M1,obj.L,obj.K);
                        obj.uncontroled_dk_energy{sridx}=norm(NoControlSoundField.dk{sridx},'fro').^2;
    %                     obj.uncontroled_dk_energy{sridx}=norms(NoControlSoundField.dk{sridx}).^2;
                    end
                end
                
            end
            
            
        end
        
        
        
        function formSubspace(obj,method1,Iterations)                   
            for sridx = 1:2                                
                %% data_autocorrelation matrix
                if strcmp(method1,'GEVD')
                     [obj.U_total_gevd{sridx},obj.EIG_VALUE_total_gevd{sridx}]=jdiag(obj.spaCorrMtxBr{sridx},obj.spaCorrMtxDk{sridx},'vector');                                    
                end
                
                if strcmp(method1,'CG')
                    x_0=zeros(length(obj.spaCorrVecBr{sridx}),1);
                    R_IN=obj.spaCorrMtxBr{sridx}+obj.spaCorrMtxDk{sridx};
                    r_IN=obj.spaCorrVecBr{sridx};
%                     p_=1./diag(R_IN);
                     p_=1./1;
                    P_R_IN=p_.*R_IN;
                    P_r_IN=p_.*r_IN;
                    
                    [~,R_total_cg]=conjgrad(P_R_IN,P_r_IN,x_0,Iterations);
                    R_total_cg=R_total_cg./norms(R_total_cg);
                    R_total_cg=R_total_cg(:,1:Iterations);
                    [obj.U_total_cg{sridx},obj.EIG_VALUE_total_cg{sridx}]=Subspace_rotate_2_EVD(obj.spaCorrMtxBr{sridx},obj.spaCorrMtxDk{sridx},R_total_cg);

                end                              
            end            
        end
        
        function mu_final=obtainRegularization(obj,subspaceMethod,ranks,constrainMethod,gain)
            if ~strcmp(constrainMethod,'TIR')
                for sridx = 1:2
                    if strcmp(subspaceMethod,'GEVD')
                        U_total=obj.U_total_gevd{sridx}(:,1:ranks);
                        EIG_VALUE_total=obj.EIG_VALUE_total_gevd{sridx}(1:ranks);
                    end
                    if strcmp(subspaceMethod,'CG')
                        U_total=obj.U_total_cg{sridx}(:,1:ranks);
                        EIG_VALUE_total=obj.EIG_VALUE_total_cg{sridx}(1:ranks);
                    end
                    
                    %% data_autocorrelation matrix
                    c_ful=(U_total'*obj.spaCorrVecBr{sridx}).^2;
                    switch constrainMethod
                        case 'AC'
                            acc_min_num=sum(c_ful./EIG_VALUE_total);
                            acc_min_dem=sum(c_ful./EIG_VALUE_total./EIG_VALUE_total);
                            min_ac=acc_min_num./acc_min_dem;
                            obj.min_ac_db(sridx)=db10(min_ac);
                            
                            acc_min_num=c_ful'*EIG_VALUE_total;
                            acc_min_dem=sum(c_ful);
                            max_ac=acc_min_num./acc_min_dem;
                            obj.max_ac_db(sridx)=db10(max_ac);
                            
                            
                            if gain<min_ac
                                mu=0;
                                obj.mu_LB_ac(sridx)=mu;
                                obj.mu_UB_ac(sridx)=inf;
                            else
                                if gain>=max_ac
                                    mu=nan;
                                    obj.mu_LB_ac(sridx)=nan;
                                    obj.mu_UB_ac(sridx)=nan;
                                else
                                    f=@(mu_vaule) sum(c_ful.*EIG_VALUE_total./(EIG_VALUE_total+mu_vaule)./(EIG_VALUE_total+mu_vaule))/...
                                        sum(c_ful./(EIG_VALUE_total+mu_vaule)./(EIG_VALUE_total+mu_vaule))-gain;
                                    term1=(c_ful*c_ful');
                                    term2=(EIG_VALUE_total-EIG_VALUE_total').^2;
                                    term12=term1.*term2;
                                    g=@(x) 2*sum(sum(triu(term12./repmat((EIG_VALUE_total+x).^3,1,length(EIG_VALUE_total))./repmat((EIG_VALUE_total+x)'.^3,length(EIG_VALUE_total),1))))/sum(c_ful./((EIG_VALUE_total+x).^2))^2;
                                    mu_0=0;
                                    mu=newton_root_sd(f,g,mu_0);
                                    obj.mu_LB_ac(sridx)=mu;
                                    obj.mu_UB_ac(sridx)=inf;
                                end
                            end
                            mu_final(sridx)=obj.mu_LB_ac(sridx);
                        case 'SD'
                            
                            minimum_value_full=sum(c_ful./EIG_VALUE_total);
                            min_sd=(obj.sigDistort(sridx)-minimum_value_full)/obj.sigDistort(sridx);
                            obj.min_sd_db(sridx)=db10(min_sd);
                            
                            max_sd=1;
                            obj.max_sd_db(sridx)=0;
                            
                            
                            
                            if gain<min_sd
                                mu=nan;
                                obj.mu_LB_sd(sridx)=nan;
                                obj.mu_UB_sd(sridx)=nan;
                            else
                                if gain>=max_sd
                                    mu=nan;
                                    obj.mu_LB_sd(sridx)=nan;
                                    obj.mu_UB_sd(sridx)=nan;
                                else
                                    f=@(x) (- sum((EIG_VALUE_total+2*x)./((EIG_VALUE_total+x).^2).*c_ful)+obj.sigDistort(sridx))/obj.sigDistort(sridx)-gain;
                                    g=@(x) 2.*x.*sum(c_ful./((EIG_VALUE_total+x).^3))/obj.sigDistort(sridx);
                                    mu_0=1;
                                    mu=newton_root_sd(f,g,mu_0);
                                    obj.mu_LB_sd(sridx)=0;
                                    obj.mu_UB_sd(sridx)=mu;
                                end
                            end
                            mu_final(sridx)=obj.mu_UB_sd(sridx);
                        case 'ER'
                            dem_value_full=sum(c_ful./EIG_VALUE_total./EIG_VALUE_total);
                            uncontrolled_energy=(obj.q_uncontrol'*obj.spaCorrMtxDk{sridx}*obj.q_uncontrol);
                            min_er=uncontrolled_energy/dem_value_full;
                            obj.min_er_db(sridx)=db10(min_er);
                            
                            %                         max_er=inf;
                            obj.max_er_db(sridx)=inf;
                            
                            
                            
                            if gain<min_er
                                mu=0;
                                obj.mu_LB_er(sridx)=0;
                                obj.mu_UB_er(sridx)=inf;
                            else
                                f=@(x) uncontrolled_energy/sum(c_ful./(EIG_VALUE_total+x)./(EIG_VALUE_total+x))-gain;
                                g=@(x) 2*uncontrolled_energy/(sum(c_ful./(EIG_VALUE_total+x)./(EIG_VALUE_total+x)).^2)*sum(c_ful./((EIG_VALUE_total+x).^3));
                                mu_0=0;
                                mu=newton_root_sd(f,g,mu_0);
                                obj.mu_LB_er(sridx)=mu;
                                obj.mu_UB_er(sridx)=inf;
                            end
                            mu_final(sridx)=obj.mu_LB_er(sridx);
                            
                            
                            
                    end
                    
                end
            else
                
                for sridx = 1:2
                    if strcmp(subspaceMethod,'GEVD')
                        U_total{sridx}=obj.U_total_gevd{sridx}(:,1:ranks);
                        EIG_VALUE_total{sridx}=obj.EIG_VALUE_total_gevd{sridx}(1:ranks);
                    end
                    if strcmp(subspaceMethod,'CG')
                        U_total{sridx}=obj.U_total_cg{sridx}(:,1:ranks);
                        EIG_VALUE_total{sridx}=obj.EIG_VALUE_total_cg{sridx}(1:ranks);
                    end
                end
                
                ratio_power=obj.inputpower(2)/obj.inputpower(1);
                c_square_1=(U_total{1}'*obj.spaCorrVecBr{1}).^2; 
                c_square_1=c_square_1*ratio_power;
                c_square_2=(U_total{2}'*obj.spaCorrVecBr{2}).^2;
                mu_final=check_roots_exist(c_square_1,c_square_2,EIG_VALUE_total,gain(1),gain(2));
                
                
                
% iter=1;
% sample_set=0:.1:1;
% % sample_set=1e3;
% finalv=zeros(size(sample_set));
% finalv2=zeros(size(sample_set));
% for values=sample_set
% mu_vaule11_try=values;
% a_value=sum(c_square_1.*EIG_VALUE_total{1}./(EIG_VALUE_total{1}+mu_vaule11_try)./(EIG_VALUE_total{1}+mu_vaule11_try));
% f_try=@(mu_vaule2) a_value/gain(1)-sum(c_square_2./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2));
% f_try_prime=@(mu_vaule2) 2*sum(c_square_2./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2));
% x0=0;
% finalv(iter)=newton_root_sd(f_try,f_try_prime,x0);
% 
% b_value=sum(c_square_1./(EIG_VALUE_total{1}+mu_vaule11_try)./(EIG_VALUE_total{1}+mu_vaule11_try));
% f_try2=@(mu_vaule2) b_value*gain(2)-sum(c_square_2.*EIG_VALUE_total{2}./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2));
% f_try2_prime=@(mu_vaule2) 2*sum(c_square_2.*EIG_VALUE_total{2}./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2));
% x0=0;
% finalv2(iter)=newton_root_sd(f_try2,f_try2_prime,x0);
% 
% iter=iter+1;
% end
% 
%  plot(sample_set,finalv)
% hold on;plot(sample_set,finalv2)
                
                
%                 
%                 if roots_exist
%                     f1=@(mu_vaule1,mu_vaule2) sum(c_square_1.*EIG_VALUE_total{1}./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1))/sum(c_square_2./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2))-gain(1);
%                     f2=@(mu_vaule1,mu_vaule2) sum(c_square_2.*EIG_VALUE_total{2}./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2))/sum(c_square_1./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1))-gain(2);
%                     
%                     
%                     g11=@(mu_vaule1,mu_vaule2) -2*sum(c_square_1.*EIG_VALUE_total{1}./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1))/sum(c_square_2./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2));
%                     g12=@(mu_vaule1,mu_vaule2) 2*sum(c_square_1.*EIG_VALUE_total{1}./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1))...
%                         *sum(c_square_2./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2))/(sum(c_square_2./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2)))^2;
%                     g21=@(mu_vaule1,mu_vaule2) 2*sum(c_square_2.*EIG_VALUE_total{2}./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2))...
%                         *sum(c_square_1./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1))/(sum(c_square_1./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1)))^2;
%                     g22=@(mu_vaule1,mu_vaule2) -2*sum(c_square_2.*EIG_VALUE_total{2}./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2))/sum(c_square_1./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1));
%                     x0=[0;0];
%                     
%                     mu_final=newton_root_sd_dim2(f1,f2,g11,g12,g21,g22,x0);
%                 else
%                     disp('The solution does not exist');
%                     mu_final=[nan,nan];
%                 end
                
                
            end
        end
        
        
        
        
        
        
        function obtain_filter(obj,subspaceMethod,ranks,mu_value)
             for sridx = 1:2
                if strcmp(subspaceMethod,'GEVD')
                    U_total=obj.U_total_gevd{sridx}(:,1:ranks);
                    EIG_VALUE_total=obj.EIG_VALUE_total_gevd{sridx}(1:ranks);
                end
                if strcmp(subspaceMethod,'CG')
                    U_total=obj.U_total_cg{sridx}(:,1:ranks);
                    EIG_VALUE_total=obj.EIG_VALUE_total_cg{sridx}(1:ranks);
                end
                c=U_total'*obj.spaCorrVecBr{sridx};
                obj.FILTER_COEFFS{sridx}=U_total*([1./(mu_value(sridx)+EIG_VALUE_total)].*c);
             end
        end
        function process(obj,indata)
            for sridx = 1:2
                %      filter_vast{sridx}=CG_o_total_over_ac(spaCorrMtxBr1{sridx},spaCorrMtxDk{sridx},spaCorrVecBr1{sridx},sigDistort1(sridx),round((L*J)/8),K_SET);
                if sridx==1
                    [ReprodSoundField.Br{sridx}]=get_produced_x_signal(indata{sridx},obj.FILTER_COEFFS{sridx},obj.irArray{sridx},obj.M1,obj.L,obj.K);
                    [ReprodSoundField.dk{sridx}]=get_produced_x_signal(indata{sridx},obj.FILTER_COEFFS{sridx},obj.irArray{mod(sridx,2)+1},obj.M2,obj.L,obj.K);
                    
%                     [NoControlSoundField.Br{sridx}]=get_produced_x_signal(indata{sridx},obj.q_uncontrol,obj.irArray{sridx},obj.M1,obj.L,obj.K);
%                     [NoControlSoundField.dk{sridx}]=get_produced_x_signal(indata{sridx},obj.q_uncontrol,obj.irArray{mod(sridx,2)+1},obj.M2,obj.L,obj.K);
                else
                     [ReprodSoundField.Br{sridx}]=get_produced_x_signal(indata{sridx},obj.FILTER_COEFFS{sridx},obj.irArray{sridx},obj.M2,obj.L,obj.K);
                    [ReprodSoundField.dk{sridx}]=get_produced_x_signal(indata{sridx},obj.FILTER_COEFFS{sridx},obj.irArray{mod(sridx,2)+1},obj.M1,obj.L,obj.K);
                    
%                      [NoControlSoundField.Br{sridx}]=get_produced_x_signal(indata{sridx},obj.q_uncontrol,obj.irArray{sridx},obj.M2,obj.L,obj.K);
%                      [NoControlSoundField.dk{sridx}]=get_produced_x_signal(indata{sridx},obj.q_uncontrol,obj.irArray{mod(sridx,2)+1},obj.M1,obj.L,obj.K);
                end
                
                
                obj.measured_ac{sridx}=db10(mean(sum(norms(ReprodSoundField.Br{sridx}).^2))/mean(sum(norms(ReprodSoundField.dk{sridx}).^2)));

                inputdata=indata{sridx}.xin;
                rlen = size(obj.irVirsrc{sridx},1)+length(inputdata)-1;
                rlen_p2 = 2^nextpow2(rlen);
                XX = fft(inputdata,rlen_p2);
                filter_coefs=obj.irVirsrc{sridx};
                YY = fft(filter_coefs,rlen_p2);
                rr = ifft(XX.*YY,'symmetric');
                filtered_d = rr(1:length(inputdata),:);    %column vector                
                obj.measured_sd{sridx}=db10(mean((norms(filtered_d-ReprodSoundField.Br{sridx}).^2./(norms(filtered_d).^2))));
                
                obj.measured_er{sridx}=db10(mean((obj.uncontroled_dk_energy{sridx}./(norms(ReprodSoundField.dk{sridx}).^2))));
            end
            
            for sridx = 1:2
                obj.measured_tir{sridx}=db10(mean((norms(ReprodSoundField.Br{sridx}).^2./(norms(ReprodSoundField.dk{mod(sridx,2)+1}).^2))));
                obj.ReprodSoundField_FINAL{sridx}(:,1)=ReprodSoundField.Br{sridx}(:,1)+ReprodSoundField.dk{mod(sridx,2)+1}(:,1);
            end
        end
    end
end
function y=log_sumsum_exp_ls(x)
max_temp=max(max(x));
inf_ind=isinf(max_temp);
y=log(sum(sum(exp(x-max_temp))))+max_temp;
y(inf_ind)=-inf;
end