function symmat = gsmat(n,density,rc)
% Generate Symmetric Matrix with specific eigenvalues
% This is the same with a function "sprandsym", but in a full matrix form.
% SEE also sprandsym
%
% sprandsym(2,1,[1 0].')
%
% ans =
% 
%    (1,1)       0.8758
%    (2,1)       0.3298
%    (1,2)       0.3298
%    (2,2)       0.1242
%
%
% gsmat(2,1,[1 0].')
%
% ans =
%     0.9992    0.0291
%     0.0291    0.0008
symmat = full(sprandsym(n,density,rc));

end

