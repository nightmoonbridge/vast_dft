function R1=fast_autocorrelation_R(clm_values,row_values)
%% The function compute the outer product of the toeplitize matrix U, i.e., U*U'/M= 1/M\sum_{m=1}^{M} u_m*u_m^T
%  U=toeplitz(clm_values,row_values)
%  R1=U*U'/length(row_values);
%  Using the naive implementation, when U is n \times m, U*U' takes n^2 \times m operations. Also, we
%  need to create U first.
%  However, for this implementation, we only need to performs n*m
%  operations.
%  Written by Liming Shi (ls@create.aau.dk), in Aalborg University, Audio analysis lab,
%  Denmark, 

dim_=length(clm_values);
R1=zeros(dim_,dim_);
row_1=row_values;
for jj=2:dim_
    row_2=[clm_values(jj:-1:2);row_1(1:end-jj+1)];
    R1(1,jj)=sum(row_1.*row_2);
end

for ii=2:dim_
    for jj=ii+1:dim_
        R1(ii,jj)=R1(ii-1,jj-1)+clm_values(ii)*clm_values(jj)-row_1(end-ii+2)*row_1(end-jj+2);
    end
end
R1=R1+R1';
R1(1,1)=sum(row_1.*row_1);
for ii=2:dim_
      R1(ii,ii)=R1(ii-1,ii-1)+clm_values(ii)*clm_values(ii)-row_1(end-ii+2)*row_1(end-ii+2);
end
R1=R1/length(row_values);
end