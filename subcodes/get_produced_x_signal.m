function[final_output]=get_produced_x_signal(indata,filter_con,irArray,M,L,K)

%%
% tic
% for ll=1:L
%     row_inx=(ll-1)*J+1:(ll)*J;
%     CON_filtered_x(:,ll)=filter(filter_con.conFilter_RRCG(row_inx),1,indata.xin(1:end));
% end
% 
% 
% 
% for ii=1:M
%     for ll=1:L
%         row_inx=(ll-1)*K+1:(ll)*K;
%         a=irArray(row_inx,ii);
%         filtered_x(:,ii,ll)=filter(a,1,CON_filtered_x(:,ll));              
%     end
% end
% final_output=squeeze(sum(filtered_x,3));
% toc



% tic
J=round(length(filter_con)/L);
rlen = J+length(indata.xin)-1;
rlen_p2 = 2^nextpow2(rlen);
XX = fft(indata.xin,rlen_p2);
filter_coefs=reshape(filter_con,J,L);
YY = fft(filter_coefs,rlen_p2);
rr = ifft(XX.*YY,'symmetric');
CON_filtered_x1 = rr(1:length(indata.xin),:);    %column vector


rlen = K+length(indata.xin)-1;
rlen_p2 = 2^nextpow2(rlen);
filtered_x=zeros(length(indata.xin),M,L);
for ll=1:L
    XX = fft(CON_filtered_x1(:,ll),rlen_p2);
    row_inx=(ll-1)*K+1:(ll)*K;
    filter_coefs=irArray(row_inx,:);
    YY = fft(filter_coefs,rlen_p2);
    rr = ifft(XX.*YY,'symmetric');
    filtered_x(:,:,ll) = rr(1:length(indata.xin),:);    %column vector
end
final_output=squeeze(sum(filtered_x,3));
% toc





end