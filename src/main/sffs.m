function [rParaOut,fx,x,fXbest,Xbest]=sffs(p,k,A,d,c,diagInv,iPara,rPara,IotherPara,IstopCondPara,targetfbest)
% 8/4/2026 new implementation for sffs

% parameters for quad. min.
todebug=0;
whichalgo=IotherPara(23); % selected algo. to fit the data in the reduced dim.
isTF=0; % whether to use target fbest during quad_min in the reduced space, we are not using target fbest inside quad.min. call
isSoftStop=1; % to use soft stop while computing min. in the reduced space or not
firstDigit=floor( IotherPara(21)/(10^(length(num2str(IotherPara(21)))-1)) );  % first digit of IotherPara(21) will tell use which version to use. 
iPara(2)=IotherPara(21)-firstDigit*( 10^(length(num2str(IotherPara(21)))-1) ); % modify iPara for refinement
num_fcalls=0;
stopflag=0; 
timeLimit=abs(IstopCondPara(6)); % hard cputime limit for algo.
itime=cputime; % to save the cputime for the algorithm
iter=0;
psuppk=[];
Xbest=cell(1,k);  %set to hold the best set of ki predictors 
fXbest=inf(1,k); % will store the best f value f(X(1:i)) for i=1:k 
I_curr=[]; % current set of predictors
fXarray=zeros(1,k); % fXarray(i) is the rss for predictors I_curr(1:i)
%Step 0 : Sequentially select the first 2 features using forward stepwise selection
[~,xhat,fXarray(1),I_curr]=fs(p,1,A,d,c,diagInv,iPara,rPara,IotherPara,IstopCondPara,targetfbest);
if todebug==1
fprintf('p=%d, tmax=%d ========================= \n',p,k);    
fprintf('X=');printArray(I_curr,'%d');
fprintf('fX=');printArray(fXarray,'%1.5f');
end

if k==1
   x=xhat;fx=fXarray(1);
   rParaOut.numOfIter=iter;rParaOut.necConMaxVioRefQM=-1;rParaOut.stopflag=stopflag;rParaOut.nfcall=num_fcalls;
   return
end
ki=1; % current features in the set X
Xbest{ki}=I_curr;
fXbest(ki)=fXarray(ki);

    while ki<k
        iter=iter+1;
        if todebug==1
        fprintf('iter=%d; ki=%d ============== \n',iter,ki);
        fprintf('X=');printArray(I_curr,'%d');
        fprintf('fX=');printArray(fXarray(1:ki),'%1.9f');
        end
        %Step 1 : Inclusion
        % istarlocal is the variable selected to be in the support in the reduced space
        Ic=setdiff(1:p,I_curr); % complement of the current support set X(1:ki)
        [num_fcalls,~,fXkplus1,istar]=mostSigniXc(fXarray(ki),I_curr,Ic,p,ki, A,d,c,diagInv,  iPara,rPara,whichalgo,isSoftStop, isTF,targetfbest, num_fcalls);
        I_curr=union(I_curr,istar,'stable');
        ki=ki+1;
        fXarray(ki)=fXkplus1;
        if todebug==1
        fprintf('istar=%d\n',istar);
        fprintf('Xkplus1=');printArray(I_curr,'%d');
        end
        
        % save the support if X has select k features at some point
        if ki==k
           psuppk=Xbest{ki}; 
        end
        
        % update the global best for size ki 
        if fXkplus1<fXbest(ki)
            Xbest{ki}=I_curr;
            fXbest(ki)=fXkplus1;
        end

        %Step 2 : Conditional Exclusion
        cond_exclusion=1;
        while cond_exclusion==1 && ki>2
            % find a predictor with the least significant value
            [num_fcalls,fstar,xr]=leastSigniX(fXkplus1,I_curr,p,ki, A,d,c,diagInv, iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest,num_fcalls);
            
            if fstar<fXbest(ki-1)
                % accept the new set with one less predictor
                I_curr=setdiff(I_curr,xr,'stable');
                ki=ki-1;
                fXarray(ki)=fstar;
                
                % update global best
                Xbest{ki}=I_curr;
                fXbest(ki)=fstar;
            else
                cond_exclusion=0; % exit the while loop
            end
            
            ctime=(cputime-itime)/60; % intermediate check for cputime limit
            if ctime > timeLimit  % if the cputime becomes greater than timeLimit=|IstopCondPara(6)| min, then stop
               stopflag=6;
               rParaOut.numOfIter=iter;rParaOut.necConMaxVioRefQM=-1;rParaOut.stopflag=stopflag;rParaOut.nfcall=num_fcalls;
               break;
            end
        end
        
        %fprintf('ki=%d, |supp|=%d \n',ki,length(X(X~=0)));
    end % while ki=/k
    
    if stopflag==6
       nred=length(Xbest{ki});
       if nred<k && ~isempty(psuppk)  % if selected support is less than k, use previous set with k support if available
          X=psuppk;
          nred=length(X);
       end
    else, nred=k;X=Xbest{ki};   
    end
    % fit the data
    x=zeros(p,1);
    sortX=sort(X);
    [Q,diagD]=eig(A(sortX,sortX),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
    Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
    Sstruct.Q=Q(:, (nred-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
    Sstruct.D=diagD( (nred-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
    Sstruct.diagInv=diagInv(sortX);
    [~,x(sortX),fx]=quad_min_algo_pool(nred,A(sortX,sortX),d(sortX),c,[],[],zeros(nred,1),Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);
    rParaOut.numOfIter=iter;rParaOut.necConMaxVioRefQM=-1;rParaOut.stopflag=stopflag;rParaOut.nfcall=num_fcalls;
end %==========================================================================================================================================================================

%% mostSigniXc
function [num_fcalls,xout,fstar,istar]=mostSigniXc(fX,X,Xc,p,ki, A,d,c,diagInv,  iPara,rPara,whichalgo,isSoftStop, isTF,targetfbest,num_fcalls)
% find the most significant feature w.r.t the set X complement
% i.e. find the feature x* such that S(x*)= max (f(X) - f(X + xj) )  for xj in Xc
% ki=no. of entries in X or predictors already in the model
    pS=-inf; 
    for i=1:(p-ki)
       Xplus1=union(X,Xc(i)); % include one feature from Xc at a time;
       idx_kiplus1=false(1,p);idx_kiplus1(Xplus1)=1; % flag of 0 and 1 indicating reduce dim

       [Q,diagD]=eig(A(idx_kiplus1,idx_kiplus1),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
       Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
       Sstruct.Q=Q(:, (ki+1-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
       Sstruct.D=diagD( (ki+1-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
       Sstruct.diagInv=diagInv(idx_kiplus1);

       [~,xopt,fXplus1]=quad_min_algo_pool(ki+1,A(idx_kiplus1,idx_kiplus1),d(idx_kiplus1),c,[],[],zeros(ki+1,1),Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);
       num_fcalls=num_fcalls+1;
       S=fX-fXplus1;
       if pS<S
          pS=S;istar=Xc(i); fstar=fXplus1; xout=xopt; 
       end
    end

end

%% leastSigniX
function [num_fcalls,fstar,istar]=leastSigniX(fX,X,p,ki, A,d,c,diagInv,  iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest,num_fcalls) 
% find the least significant features w.r.t. the set X
% i.e. find the feature x* such that S(x*)= min ( f(X-xj)-f(X) ) for xj in X

    pS=inf; 
    for i=1:ki
       Xmin1=setdiff(X,X(i),'stable'); % delete one feature from X at a time;
       idx_kimin1=false(1,p);idx_kimin1(Xmin1)=1; % flag of 0 and 1 indicating reduce dim
       
       [Q,diagD]=eig(A(idx_kimin1,idx_kimin1),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
       Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
       Sstruct.Q=Q(:, (ki-1-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
       Sstruct.D=diagD( (ki-1-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
       Sstruct.diagInv=diagInv(idx_kimin1);
       [~,~,fXmin1]=quad_min_algo_pool(ki-1,A(idx_kimin1,idx_kimin1),d(idx_kimin1),c,[],[],zeros(ki-1,1),Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);
       num_fcalls=num_fcalls+1;
       S=fXmin1-fX;
       if S<pS
          pS=S;istar=X(i); fstar=fXmin1; 
       end
    end

end

function [fstar, istar] = leastSigniX_v2(fX, X, p, ki, A, d, c, lb, ub, diagInv, iPara, rPara, whichalgo, isSoftStop, isTF, targetfbest) 
% find the least significant feature w.r.t. the set X
% i.e. find the feature x* such that S(x*) = min (f(X - xj) - f(X)) for xj in X

pS = inf; 
istar = -1;
fstar = fX;

% Pre-allocate initial guess for ki-1 dimension
x0_init = zeros(ki - 1, 1);

for i = 1:ki
    % Fast candidate feature removal without calling setdiff
    Xmin1 = X;
    Xmin1(i) = []; % removes the i-th feature directly in O(ki) time
    
    % Fast logical indexing mask
    idx_kimin1 = false(1, p);
    idx_kimin1(Xmin1) = true;
    
    % Sub-matrix extract
    A_sub = A(idx_kimin1, idx_kimin1);
    
    % Eigenvalue decomposition for reduced size (ki-1) x (ki-1)
    [Q, diagD] = eig(A_sub, 'vector');
    
    % Direct boolean filtering of non-zero eigenvalues
    nonZero_idx = diagD > rPara(1);
    
    Sstruct.nNonZeroEig = sum(nonZero_idx);
    Sstruct.Q = Q(:, nonZero_idx);
    Sstruct.D = diagD(nonZero_idx);
    Sstruct.diagInv = diagInv(idx_kimin1);
    
    % Solve quadratic sub-problem for ki-1 features
    [~, ~, fXmin1] = quad_min_algo_pool(ki - 1, A_sub, d(idx_kimin1), c, ...
        lb(idx_kimin1), ub(idx_kimin1), x0_init, Sstruct, iPara, rPara, ...
        whichalgo, isSoftStop, isTF, targetfbest);
        
    S = fXmin1 - fX;
    if S < pS
        pS = S;
        istar = X(i); 
        fstar = fXmin1; 
    end
end

end

function [xout, fstar, istar] = mostSigniXc_v2(fX, X, Xc, p, ki, A, d, c, lb, ub, diagInv, iPara, rPara, whichalgo, isSoftStop, isTF, targetfbest)
% find the most significant feature w.r.t the set X complement
% i.e. find the feature x* such that S(x*) = max (f(X) - f(X + xj)) for xj in Xc

pS = -inf; 
istar = -1;
fstar = fX;
xout = [];

numCandidates = length(Xc);

% Pre-allocate index logical array mask for base set X
base_mask = false(1, p);
base_mask(X) = true;

% Zero vector initialization for initial guess
x0_init = zeros(ki + 1, 1);

for i = 1:numCandidates
    candidate = Xc(i);
    
    % Construct candidate index mask without calling union/setdiff inside loop
    idx_kiplus1 = base_mask;
    idx_kiplus1(candidate) = true;
    
    % Sub-matrix extract
    A_sub = A(idx_kiplus1, idx_kiplus1);
    
    % Eigenvalue decomposition for the current sub-matrix
    [Q, diagD] = eig(A_sub, 'vector'); 
    
    % Extract non-zero eigenvalues and corresponding eigenvectors
    nonZero_idx = diagD > rPara(1);
    
    Sstruct.nNonZeroEig = sum(nonZero_idx);
    Sstruct.Q = Q(:, nonZero_idx);
    Sstruct.D = diagD(nonZero_idx);
    Sstruct.diagInv = diagInv(idx_kiplus1);
    
    % Solve quadratic sub-problem
    [~, xopt, fXplus1] = quad_min_algo_pool(ki + 1, A_sub, d(idx_kiplus1), c, ...
        lb(idx_kiplus1), ub(idx_kiplus1), x0_init, Sstruct, iPara, rPara, ...
        whichalgo, isSoftStop, isTF, targetfbest);
    
    S = fX - fXplus1;
    if S > pS
        pS = S;
        istar = candidate; 
        fstar = fXplus1; 
        xout = xopt; 
    end
end

end

%% mostSigniX
function [fstar,istar]=mostSigniX(fX,X,p,ki, A,d,c,lb,ub,  iPara,rPara,IotherPara)
% find the most significant feature w.r.t. the set X
% i.e. find the feature x* such that S(x*)= max (f(X-xj) - f(X) )  for xj in X

    pS=-inf; 
    for i=1:ki
       Xmin1=setdiff(X,X(i)); % delete one feature from X at a time;
       idx_kimin1=false(1,p);idx_kimin1(Xmin1)=1; % flag of 0 and 1 indicating reduce dim
       [~,~,fXmin1]=quad_min_algo_pool(ki-1,A(idx_kimin1,idx_kimin1),d(idx_kimin1),c,lb(idx_kimin1),ub(idx_kimin1),zeros(ki-1,1),iPara,rPara,IotherPara(23),1,0,[]);
       S=fXmin1-fX;
       if pS<S
          pS=S;istar=X(i); fstar=fXmin1; 
       end
    end

end

%% leastSigniXc
function [fstar,istar]=leastSigniXc(fX,X,p,ki, A,d,c,lb,ub,  iPara,rPara,IotherPara)
% find the least significant feature in X complement
% i.e. the feature x* such that S(x*)= min ( f(X)-f(X+xj) ) for xj in Xc

    pS=inf; 
    Xc=setdiff(1:p,X); 
    for i=1:(p-ki)
       Xplus1=union(X,Xc(i)); % delete one feature from X at a time;
       idx_kiplus1=false(1,p);idx_kiplus1(Xplus1)=1; % flag of 0 and 1 indicating reduce dim
       [~,~,fXplus1]=quad_min_algo_pool(ki+1,A(idx_kiplus1,idx_kiplus1),d(idx_kiplus1),c,lb(idx_kiplus1),ub(idx_kiplus1),zeros(ki+1,1),iPara,rPara,IotherPara(23),1,0,[]);
       S=fX-fXplus1;
       if S<pS
          pS=S;istar=Xc(i); fstar=fXplus1; 
       end
    end

end
















