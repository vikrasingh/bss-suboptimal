function [stopflag,xhat,fhat,X]=fs1(p,k,A,d,c,lb,ub,diagInv,iPara,rPara,IotherPara,IstopCondPara,targetfbest)
% 19May24
%adding those features one at a time that maximizes the reduction in RSS

whichalgo=IotherPara(23); % selected algo. to fit the data in the reduced dim.
isTF=0; % whether to use target fbest during quad_min in the reduced space, we are not using target fbest inside quad.min. call
isSoftStop=1; % to use soft stop while computing min. in the reduced space or not
firstDigit=floor( IotherPara(21)/(10^(length(num2str(IotherPara(21)))-1)) );  % first digit of IotherPara(21) will tell use which version to use. 
iPara(2)=IotherPara(21)-firstDigit*( 10^(length(num2str(IotherPara(21)))-1) ); % modify iPara for refinement
stopflag=0; 
timeLimit=abs(IstopCondPara(6)); % hard cputime limit for algo.
itime=cputime; % to save the cputime for the algorithm
fXarray=zeros(1,k); % will store the f value f(X(1:i)) for i=1:k
X=[];
Xc=1:p;
ki=1;
[xhat,fXarray(ki),istar]=mostSigniXc(c,X,Xc,p,0, A,d,c,lb,ub,diagInv,  iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);
X=istar; Xc=setdiff(Xc,istar);
if k==1
   xout=zeros(p,1);xout(istar)=xhat;xhat=xout; fhat=fXarray; return; 
end

    while ki<k 
       [xhat,fXarray(ki+1),istar]=mostSigniXc(fXarray(ki),X,Xc,p,ki, A,d,c,lb,ub,diagInv,  iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);
       Xc=setdiff(Xc,istar); 
       X=union(X,istar,'stable');  % stable for not sort the entries after union 
       ki=ki+1;
    
       % check stopping condition
       ctime=(cputime-itime)/60; % intermediate check for cputime limit
       if ctime > timeLimit  % if the cputime becomes greater than timeLimit=|IstopCondPara(6)| min, then stop
           stopflag=6;
           break;
       end
    end
%xhat is the final x     
fhat=fXarray(ki); % final fx
X,fXarray
X=sort(X);
xout=zeros(p,1);xout(X)=xhat;xhat=xout;


end %========================================================================================================


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
       Sstruct.Q=Q(:, (ki+1-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
       Sstruct.D=diagD( (ki+1-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
       Sstruct.diagInv=diagInv(idx_kiplus1);

       [~,xopt,fXplus1]=quad_min_algo_pool(ki+1,A(idx_kiplus1,idx_kiplus1),d(idx_kiplus1),c,lb(idx_kiplus1),ub(idx_kiplus1),zeros(ki+1,1),Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);
       S=fX-fXplus1;
       if pS<S
          pS=S;istar=Xc(i); fstar=fXplus1; xout=xopt; 
       end
    end

end