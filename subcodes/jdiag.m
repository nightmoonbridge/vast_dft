function [U,D] = jdiag(A, B, evaOption, eigOption)
% JDIAG    Joint diagonalization function
%  [U,D] = JDIAG(A, B, evaOption)
% JDIAG returns the eigenvectors and the eigenvalues from Au = dBu
% where u is an eigenvector and d is an eigenvalue, respectively.
% Both are in a range of 1 <= d, u <= dim(A or B).
% U gives you the joint diagonalization property such that inv(B)*A*U = U*D
%                                U'*A*U = D
%                                U'*B*U = I
% where
%   I is the identity matrix,
%   D is the diagonal matrix whose elements are the eigenvalues
%     (typically in descending order),
%   U is the eigenvector matrix corresponding to D,
%
% and this has a relationship described as follows:
%           diag(U'*A*U) = diag(Ueve'*A*Ueve)./diag(Ueve'*B*Ueve)
%
% Although this gives you a similar solution from [Ueve,Ueva] = eig(B\A),
% the order of the eigenvalues can be different from each other.
%
% JDIAG input arguments:
% A                              - a (semi) positive definite matrix
% B                              - a positive definite matrix
% evaOption                      - 'vector' returns D as a vector, diag(D)
%                                - 'matrix' returns D as a diag. matrix
%
% Latest update   :     21st/October-2019
% Taewoong Lee (tlee at create.aau.dk)
%
% This was modified from the code 'jeig.m' provided in the following book:
%  [1] J. Benesty, M. G. Christensen, and J. R. Jensen,
%    Signal enhancement with variable span linear filters. Springer, 2016.
%
%  DOI: 10.1007/978-981-287-739-0
%
%
% For example,
%  rng default
%  A = full(sprandsym(3,1,[3 4 5]));
%  B = full(sprandsym(3,1,[10 20 30]));
%  [U,D] = JDIAG(A,B);
%
% U'*A*U                                U'*B*U
% ans =                                 ans =
%     0.4313    0.0000   -0.0000            1.0000    0.0000   -0.0000
%     0.0000    0.1662   -0.0000            0.0000    1.0000   -0.0000
%    -0.0000   -0.0000    0.1395           -0.0000   -0.0000    1.0000
%
% [Ueve,Ueva] = eig(B\A);
% Ueve'*A*Ueve                          Ueve'*B*Ueve
% ans =                                 ans =
%     4.4291   -0.0000   -0.0000           10.2682    0.0000   -0.0000
%    -0.0000    3.2703   -0.0000            0.0000   19.6714   -0.0000
%    -0.0000    0.0000    3.9718           -0.0000   -0.0000   28.4816
%
% diag(Ueve'*A*Ueve)./diag(Ueve'*B*Ueve)
% ans =
%     0.4313
%     0.1662
%     0.1395
%
if nargin < 4
    eigOption = false;
    if nargin < 3
        evaOption = 'matrix';
    end
end
if eigOption
    [X, d1] = eig(A, B, 'chol', 'vector');
    
    [D,dind] = sort(d1,'descend');
%     D = diag(dd);
    U = X(:,dind);
else
    [Bc,pd] = chol(B,'lower');  % B = Bc*tranpose(Bc)
    argname = char(inputname(2));
    
    if pd ~= 0
        error(['Matrix ', argname ,' seems NOT a positive definite.']);
    elseif pd == 0
        % Matrix B is a Positive definite.
        C = Bc\A/Bc';           % C = inv(Bc)*A*inv(conjugate transpose(Bc))
        [U,T] = schur(C);
        X = Bc'\U;              % X = inv(conjugate transpose(Bc))*U;
        
        [dd,dind] = sort(diag(T),'descend');
        D = diag(dd);
        U = X(:,dind);
    end
    
    switch lower(evaOption)
        case 'vector'
            D = dd;
        otherwise
    end
end
end