% The nearest polynomial to multiple given polynomials with a given zero: A
% unified optimization approach
% Nov. 27, 2019

clc
clear all
close all
format long

m = 3; % function number
n = m ; % basis number
z = 2.0;   % -1,2 % given zero  

ThrV = 1e-10;
beta = 0.1;
beta_max =1.0e8;
gamma=1.6;
isConverged=0;
iterMax = 1000;
r = 2;
s = 2;

% compute function, basis
% e_j(x), Ve
Ve=zeros(n,1);
for j=1:n
    Ve(j)=z^(j-1);
end
norme=norm(Ve,2);
% Vf0
f0=@(x)-2-x+x^2;
bc=[-2;-1;1];
Vf0=zeros(m,1);
for i=1:m
    Vf0(i)=f0(z);  % vector f0
end
% A
if m<=n 
    % initialize A, Vf, B
    A=zeros(m,n); % 3X3
    A(1,:) = [0.01, 0, 0];
    A(2,:) = [0, -0.02,  0];
    A(3,:) = [0, 0.00, 0.03];
    A = double(A);
    
    Vf=zeros(m,1);
    Vf=Vf0+A*Ve;
    
    B=zeros(m,m);
    B=eye(m)-ones(m,1)*ones(m,1)'/m;
    
    
    X=zeros(m,n);
    Y=zeros(m,n);
    Vc=zeros(n,1);
    Lambda1 = zeros(m,n);
    Lambda2 = zeros(m,1);
    Lambda3 = zeros(m,n);
    iter=0;
    
    tic
    while ~isConverged
        iter=iter+1;
        
        % update X
        pre_X=X;
        temp = B'*(B*A+Lambda1)+(Vf+Lambda2)*Ve'+Y+Lambda3;
        X = lyap(eye(m)+B'*B,Ve*Ve',-temp);
        
        % update Y
        temp = X-Lambda3;
        
        switch lower([num2str(r),',',num2str(s)])
            case 'inf,inf'
                [Y,  ] = prox_linf(temp, 1/beta);
            case '1,1'
                for i=1:m
                    Y(i,:) = prox_l1(temp(i,:),1/beta);
                end
            case '2,1'
                for i=1:m
                    Y(i,:) = prox_l21(temp(i,:)',1/beta);
                end
            case '2,2'
                Y = beta*temp/(2+beta);
            case 'inf,1'
                for i=1:m
                    [Y(i,:), ] = prox_linf(temp(i,:),1/beta);
                end
            case ohterwise
            disp(['Unknown norm'])
        end
        
        % updeate multipliers
        Lambda1 = Lambda1-B*(X-A);
        Lambda2 = Lambda2-(X*Ve-Vf);
        Lambda3 = Lambda3-(X-Y);
        beta = min(gamma*beta, beta_max);
        
        % Stopping check
        if iter>=1
            normX=compute_norm(X,r,s);
            relGap=abs(compute_norm(pre_X,r,s)-compute_norm(X,r,s))/compute_norm(X,r,s);
            Res1 = B*(X-A);  % BX=BA
            Res2 = X*Ve-Vf;  % Xe=f
            Res3 = X-Y  ;    % X=Y
            disp(['#iter ', num2str(iter), ': normX= ', num2str(normX), ', relGap= ', num2str(relGap)]);
            if max([norm(Res1,Inf),norm(Res2,Inf),norm(Res3,Inf),relGap])<=ThrV  || iter>=iterMax
                isConverged = 1;
            end
        end
    end
    TimeX=toc
    
    % Compute c
    Vc = (A-X)'*ones(m,1);
    Vc_val = Vc/m
    bf = bc+Vc_val
    
    % Compute final f(z)
   f_val = f0(z)+Vc_val'*Ve
%    f_val2 = bf'*Ve
   norm_val = compute_norm(X,r,s)
   iter_val = iter

   
   Res1 = B*(X-A);  % BX=BA
   Res2 = X*Ve-Vf;  % Xe=f
   Res3 = X-Y;      % X=Y
   
else
    disp('Error: m> n');
end







