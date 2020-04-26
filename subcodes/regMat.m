function [l,dif,Al] = regMat(A,b)
% Eigenvalue problem
N = size(A,1);
[ev,el] = eig(A);
eigenVects = cell(1,N);
eigenLams = cell(1,N);

if ~issorted(diag(el)) == 0
    [el,I] = sort(diag(el),'descend');
    ev = ev(:, I);
end
for ii = 1:N
    eigenVects{ii} = ev(:,ii);
    eigenLams{ii} = el(ii);
end

ii = 1;
Al = zeros(N);
while ii <= size(A,1)
    tem = (eigenVects{ii}*eigenVects{ii}')/eigenLams{ii};
    Al = Al + tem;
    dif = norm(b-A*Al*b,2)/norm(b,2);
    if dif < 0.001
%         l = ii;
        break
    end
    ii = ii+1;    
end
l = ii;

end

