function Bobj = f_boot_min(x, BData,s,T,BCon,nbasis)

beta = x(1:nbasis);
phi = reshape(x(nbasis+1:end),T-s,s);
%phi = [0.5 0.5;1.1 0.2;0.1 0.6]';

% Break up the data as it is coming from bootstrap sample
BY = BData(1:length(BData)/2,1:T); 
BZ = BData(1:length(BData)/2,T+1:2*T); 
BXS = BData(1:length(BData)/2,2*T+1:end);

Y = BData(length(BData)/2+1:end,1:T); 
Z = BData(length(BData)/2+1:end,T+1:2*T); 
XSplines = BData(length(BData)/2+1:end,2*T+1:end);

n= length(Z);

% And calculate the M functions
for i= 1 : T
    BA(:,i) = BY(:,i) - BXS(:,1+nbasis*(i-1):nbasis*i) * beta;
    A(:,i) = Y(:,i) - XSplines(:,1+nbasis*(i-1):nbasis*i) * beta; 
end

BM = BA(:,1:s) + BA(:,s+1:end) * phi;
M = A(:,1:s) + A(:,s+1:end) * phi;

% And the objective function
K1 = zeros(n,n);
K2 = [];
BK1 = zeros(n,n);
BK2 = [];
KK = [];
obj = 0;
Bobj = 0;
for i = 1 : s
    K2 = [K2 K1+(repmat(Z(:,i),1,n)-repmat(Z(:,i)',n,1)).^2];
    BK2 = [BK2 BK1+(repmat(BZ(:,i),1,n)-repmat(BZ(:,i)',n,1)).^2];
    K1 = K2(:,(i-1)*n+1:i*n);
    BK1 = BK2(:,(i-1)*n+1:i*n);
    KK = [KK (M(:,i)' * sqrt(K1) * M(:,i))/n];
    obj = obj + ((M(:,i)' * sqrt(K1) * M(:,i))^2)/n^2;
    Bobj = Bobj + ((BM(:,i)' * sqrt(BK1) * BM(:,i))/n-(n-1)/n*KK(i))^2;
end

Bobj = Bobj + BCon * obj / (n^2);