function obj = f_min2(x,beta,Data,s,T,nbasis)

x = reshape(x,T-s,s);
% Break up the data
Y = Data(:,1:T); Z = Data(:,T+1:2*T); XSplines = Data(:,2*T+1:end);
n= length(Z);

% And calculate the M functions
for i= 1 : T
    A(:,i) = Y(:,i) - XSplines(:,1+nbasis*(i-1):nbasis*i) * beta; 
end

M = A(:,1:s) + A(:,s+1:end) * x;

% And the objective function
K1 = zeros(n,n);
K2 = [];
obj = 0;
for i = 1 : s
    K2 = [K2 K1+(repmat(Z(:,i),1,n)-repmat(Z(:,i)',n,1)).^2];
    K1 = K2(:,(i-1)*n+1:i*n);
    obj = obj + ((M(:,i)' * sqrt(K1) * M(:,i))^2)/n^2;
end

