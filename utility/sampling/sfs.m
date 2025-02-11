function [rParaOut,X,xtilde,fxtilde,x,fx]=sfs(p,lb,ub,A,d,c,k,xRelaxedOpt,ixstar,diagInv,targetfbest ,IotherPara,IstopCondPara, iPara,rPara,textfileName,todebug)
% 21 May24

% parameters for quad. min.
whichalgo=IotherPara(23); % selected algo. to fit the data in the reduced dim.
if ismember(whichalgo,12:20)
   isBoxed=0;  % are we using unconstrained or constrained version, if isBoxed=1 constrained, isBoxed=0 unconstrained 
else
   isBoxed=1;  % are we using unconstrained or constrained version, if isBoxed=1 constrained, isBoxed=0 unconstrained  
end
isTF=0; % whether to use target fbest during quad_min in the reduced space, we are not using target fbest inside quad.min. call
isSoftStop=1; % to use soft stop while computing min. in the reduced space or not
firstDigit=floor( IotherPara(21)/(10^(length(num2str(IotherPara(21)))-1)) );  % first digit of IotherPara(21) will tell use which version to use. 
iPara(2)=IotherPara(21)-firstDigit*( 10^(length(num2str(IotherPara(21)))-1) ); % modify iPara for refinement
stopflag=0; 
timeLimit=abs(IstopCondPara(6)); % hard cputime limit for algo.
itime=cputime; % to save the cputime for the algorithm
stop=0;
numOfIter=0;
[~,X]=maxk(abs(xRelaxedOpt),k);  % find the largest Tm components of xRelaxedOpt
 % the initial point has been chosen as the truncated xRelaxedopt with
X=sort(X');      
     
[Q,diagD]=eig(A(X,X),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
Sstruct.Q=Q(:, (k-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
Sstruct.D=diagD( (k-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
Sstruct.diagInv=diagInv(X);
if isBoxed==0
   [~,xhat,fX]=quad_min_algo_pool(k,A(X,X),d(X),c,[],[],xRelaxedOpt(X),Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest); 
else
   [~,xhat,fX]=quad_min_algo_pool(k,A(X,X),d(X),c,lb(X),ub(X),xRelaxedOpt(X),Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);  
end
xtilde=zeros(p,1);xtilde(X)=xhat; fxtilde=fX;

if p==k
   x=xtilde;fx=fxtilde;
   stop=1;rParaOut.numOfIter=numOfIter;rParaOut.necConMaxVioRefQM=-1;rParaOut.stopflag=stopflag;
   return;
end
     
Xc=setdiff(1:p,X);  % possible set of variables to be added from
    while stop==0
        numOfIter=numOfIter+1;
        % minimize the gain by dropping one variable
        [fXminusi,istar]=leastSigniX(fX,X,p,k, A,d,c,lb,ub,diagInv, iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);
        % i is index to drop,  fXminusi=f(X - istar)

        % maximize the reduction by adding one variable to Xminusi
        Xminusi=setdiff(X,istar);
        [~,fXplusj,jstar]=mostSigniXc(fXminusi,Xminusi,Xc,p,k, A,d,c,lb,ub,diagInv,  iPara,rPara,whichalgo,isSoftStop, isTF,targetfbest);
        % jstar index to add, fXplusj=f(X + jstar)

        % Swap if possible
        if fX>fXplusj
           X=union(Xminusi,jstar);  % new support
           Xc=union( setdiff(Xc,jstar) , istar); 
           fX=fXplusj;
        else
           stop=1;rParaOut.numOfIter=numOfIter;rParaOut.necConMaxVioRefQM=-1;rParaOut.stopflag=stopflag;
           break;
        end

        % hard stop of CPU time limit
        ctime=(cputime-itime)/60;
         if ctime >= timeLimit
            stop=1;stopflag=6;
            rParaOut.numOfIter=numOfIter;rParaOut.necConMaxVioRefQM=-1;rParaOut.stopflag=stopflag;
            break
         end
        
    end % while stop=0
    
    
    % fit the data
    x=zeros(p,1);
    sortX=X; % X is already sorted
    [Q,diagD]=eig(A(sortX,sortX),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
    Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
    Sstruct.Q=Q(:, (k-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
    Sstruct.D=diagD( (k-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
    Sstruct.diagInv=diagInv(sortX);
    [~,x(sortX),fx]=quad_min_algo_pool(k,A(sortX,sortX),d(sortX),c,lb(sortX),ub(sortX),xRelaxedOpt(sortX),Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);

end %==========================================================================================================================================================================

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

%% mostSigniXc
function [xout,fstar,istar]=mostSigniXc(fX,X,Xc,p,ki, A,d,c,lb,ub,diagInv,  iPara,rPara,whichalgo,isSoftStop, isTF,targetfbest)
% find the most significant feature w.r.t the set X complement
% i.e. find the feature x* such that S(x*)= max (f(X) - f(X + xj) )  for xj in Xc
% ki=no. of entries in X or predictors already in the model
    pS=-inf; 
    for i=1:(p-ki)
       Xplus1=union(X,Xc(i)); % include one feature from Xc at a time;
       idx_kiplus1=false(1,p);idx_kiplus1(Xplus1)=1; % flag of 0 and 1 indicating reduce dim

       [Q,diagD]=eig(A(idx_kiplus1,idx_kiplus1),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
       Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
       Sstruct.Q=Q(:, (ki-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
       Sstruct.D=diagD( (ki-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
       Sstruct.diagInv=diagInv(idx_kiplus1);

       [~,xopt,fXplus1]=quad_min_algo_pool(ki,A(idx_kiplus1,idx_kiplus1),d(idx_kiplus1),c,lb(idx_kiplus1),ub(idx_kiplus1),zeros(ki,1),Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);
       S=fX-fXplus1;
       if pS<S
          pS=S;istar=Xc(i); fstar=fXplus1; xout=xopt; 
       end
    end

end

%% leastSigniX
function [fstar,istar]=leastSigniX(fX,X,p,ki, A,d,c,lb,ub,diagInv,  iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest) 
% find the least significant features w.r.t. the set X
% i.e. find the feature x* such that S(x*)= min ( f(X-xj)-f(X) ) for xj in X
    
    if ki==1
       istar=X;fstar=c;return;
    end
    pS=inf; 
    for i=1:ki
       Xmin1=setdiff(X,X(i),'stable'); % delete one feature from X at a time;
       idx_kimin1=false(1,p);idx_kimin1(Xmin1)=1; % flag of 0 and 1 indicating reduce dim
       
       [Q,diagD]=eig(A(idx_kimin1,idx_kimin1),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
       Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
       Sstruct.Q=Q(:, (ki-1-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
       Sstruct.D=diagD( (ki-1-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
       Sstruct.diagInv=diagInv(idx_kimin1);
       [~,~,fXmin1]=quad_min_algo_pool(ki-1,A(idx_kimin1,idx_kimin1),d(idx_kimin1),c,lb(idx_kimin1),ub(idx_kimin1),zeros(ki-1,1),Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);
       S=fXmin1-fX;
       if S<pS
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
















