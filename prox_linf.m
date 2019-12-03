function [X, iter] = prox_linf(B,lambda)

% The proximal operator of the l_{inf,inf} norm of a matrix
% l_{inf,inf} norm is the maximum absolute value of all elements of a matrix 
%
% min_X lambda*||X||_{inf,inf}+0.5*||X-B||_2^2
%
% version 1.0 - 1/12/2019
%
% Written by Wenyu HU
%

[m,n] = size(B);
X = zeros(m,n);
[Bse,I] = sort(abs(B(:)),'descend'); 
Bse = [Bse;0];
iter=1;
for iter =1:m*n
    T = sum(Bse(1:iter))-lambda;
    T = T /iter;
    
    if T>=Bse(iter+1) & T<=Bse(iter)
        X(I(1:iter)) = sign(B(I(1:iter))).*T;
        X(I(iter+1:end)) = B(I(iter+1:end));
        break;
    end
    
    iter = iter+1;
end

%     if iter==m*n
%         X = zeros(m,n);
%     end



 