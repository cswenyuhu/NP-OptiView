function S=compute_norm(X,r,s)
% compute \|X\|_r,s^s

[m,n]=size(X);

S=0;
if r==inf & s==inf
    S=max(max(abs(X)));
else
    for k=1:m
        S=S+norm(X(k,:),r)^s;
    end
end
end
