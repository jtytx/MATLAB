clear
clc

load MCData;


% Set the parameters we will run this under
L = 0;                 % The constraint value         
B = [10^2; 10^2];      % The smoothness constraint, first on beta, then on bootstrap
k_n = 1;               % The number of knots for Bsplines
Ord_X = 4;             % The order of the Bspline
s = T - R;

% And the structures we'll save this in.
quant = NaN(Replications,1); 
xflag = NaN(Replications,1);
fval = NaN(Replications,1);

% Find out how many bases we actually have
A = 0:(k_n+1); A = augknt(A,Ord_X); 
A = spcol(A, Ord_X, 1:10); nbasis = size(A,2);
Save_opt = zeros(size(A,2), Replications);
%Save_opt = zeros(size(A,2)+ R*(T-R), Replications);
clear A;


% Moving on to the smoothness constraint, we need to set up a grid of
% points to evaluate the derivative of the splines
div4g = 0.01; knots = 0:1/(k_n+1):1; knots = norminv(knots(2:end-1));
knots = 2*(normcdf(knots/3) - 0.5); knots = [-1 knots 1]; knots = augknt(knots,Ord_X);
T1 = (-1:div4g:1)'; T1 = repmat(T1,Ord_X,1); T1 = sort(T1);
T2 = spcol(knots, Ord_X, T1);

XGrid = zeros(length(-1:div4g:1), nbasis, Ord_X);
for i = 1:Ord_X
    for j = 1:length(-1:div4g:1)
        XGrid(j,:,i) = T2(Ord_X*(j-1) + i,:);
    end
end

% Now with this grid we'll do numerical integration
LM = zeros(size(XGrid,2));
for d = 1:size(XGrid,3)
    for j = 1:size(XGrid,2)
        A = repmat(XGrid(:,j,d), 1, size(XGrid,2)).*XGrid(:,:,d);
        LM(j,:) = LM(j,:) + trapz(A)*div4g;
    end
end

clear A T1 T2 i j d div4g;


% Set the starting values of the first optimization
%x0 = zeros(nbasis,1)*1;
%x0 = zeros(nbasis + R*(T-R),1)*1;



for k = 1 : Replications
    tic
    x0 = [2 -4 0 4 -2]';
    %x0 = [2 -4 0 4 -2 0.1 0.6]';
    % Evaluate the X at the different splines. Knots put at the percentiles
    XSplines = [];
    for t = 1+T*(k-1) : T*k
        XSplines = [XSplines spcol(knots, Ord_X, X(:,t))];
    end
    
    % And put the data together
    Data = [Y(:,1+T*(k-1):T*k) Z(:,1+T*(k-1):T*k) XSplines];
    
    % Finish setting up the constraints
    Am = []; Bm = []; LB = []; UB = [];
    
    
    % This is the value of the function at the point 0
    Aeq = spcol(knots, Ord_X, 0);
    %Aeq = [spcol(knots, Ord_X, 0) zeros(1,R*(T-R))];
    %BootAeq = [Aeq zeros(1,size(XSplines,2)); zeros(1,size(XSplines,2)) Aeq];
    
    % And optimization options - for full sample problem
    options = optimset('fmincon'); options.GradObj = 'off'; options.Display = 'off';
    options.GradConstr = 'off'; options.DerivativeCheck = 'on';
    options.MaxFunEvals = 10000; options.MaxIter = 10000;
    options.Algorithm = 'active-set';
    options.TolFun = 1e-8;
     
%     opts = optimoptions(@fmincon,'Algorithm','sqp','GradObj','off','Display','off','GradConstr','off','DerivativeCheck','on','MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-8);
%     gs = GlobalSearch;
    % Now look at the bounds
    Beq = L;
    
    % Set up the functions to be called
    fobj = @(x)f_min(x,Data,s,T,nbasis);
    fcon = @(x)f_nonlcon(x,LM,B(1),nbasis);
    
    
    % And finally call the optimizer
    [Save_opt(:,k),fval(k),xflag(k),output] = fmincon(fobj, x0, Am, Bm, Aeq, Beq, LB, UB, fcon,options);
    
%     problem = createOptimProblem('fmincon','x0',x0,'objective',fobj,'Aeq',Aeq,'beq',Beq,'nonlcon',fcon,'options',opts);
%     [Save_opt(:,k),fval(k),xflag(k),output] = run(gs,problem);
    
    %Calculate the bootstrap statistic. Get sample, etc.
    r = 0; BI = int64(1 + rand(n,r)*(n-1));
    BXflag = zeros(r,1);BVal = zeros(r,1);
    
    %Now get the bootstrap critical value - but only if you converged
    if xflag(k) > 0
        
        %Set the starting value for the bootstrap
        %x0 = ones(nbasis,1);
        %x0 = zeros(nbasis + R*(T-R),1)*1;
        %x0 = [2 -4 0 4 -2 0.5 0.5 1.1 0.2 0.1 0.6]';
        x0 = Save_opt(:,k);
        
        BCon = n^(4/3);
        
        
        for b = 1:r
            %Put together the bootstrap sample
            BY = Y(BI(:,b),1+T*(k-1):T*k); BZ = Z(BI(:,b),1+T*(k-1):T*k);
            
            %And the Xsplines
            BXS = XSplines*0;
            for c = 1 : size(XSplines,2)
                A = XSplines(:,c);
                BXS(:,c) = A(BI(:,b));
            end
            
            BData = [BY BZ BXS; Data];
            
            %St up the functions we will call
            fbootobj = @(x)f_boot_min(x,BData,s,T,BCon,nbasis);
            fbootcon = @(x)f_boot_nonlcon(x,LM,B(1),nbasis);
            
            %Call the optimizer
            [x1, BVal(b), BXflag(b)] = fmincon(fbootobj, x0, Am, Bm, Aeq, Beq, LB, UB, fbootcon, options);
            
        end
        %Get the quantile we are looking for
        BVal = BVal(BXflag > 0);
        quant(k) = sum(BVal <= fval(k))/length(BVal);
    
    else                
        quant(k,:) = NaN;
    end
    k
    toc
    % And use the optimum as the starting point for the next optimization
    %x0 = Save_opt(:,k);
end
    

%save MCOutput quant fval Save_opt;








