function [rParaOut,X,xtilde,fxtilde,x,fx]=AG92(p,A,d,c,k,xRelaxedOpt,diagInv,targetfbest ,IotherPara,IstopCondPara, iPara,rPara)
% 21 May24
%fprintf('tmax=%d \n',k);
% parameters for quad. min.
whichalgo=IotherPara(23); % selected algo. to fit the data in the reduced dim.
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
[~,xhat,fX]=quad_min_algo_pool(k,A(X,X),d(X),c,[],[],xRelaxedOpt(X),Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest); 
xtilde=zeros(p,1);xtilde(X)=xhat; fxtilde=fX;

if p==k
   x=xtilde;fx=fxtilde;
   stop=1;rParaOut.numOfIter=numOfIter;rParaOut.necConMaxVioRefQM=-1;rParaOut.stopflag=stopflag;
   return;
end

if k==1 || (p-k)==1
   [rParaOut,X,xtilde,fxtilde,x,fx]=AG9(p,[],[],A,d,c,k,xRelaxedOpt,diagInv,targetfbest ,IotherPara,IstopCondPara, iPara,rPara);
   return;
end
ds=2;  % for this subroutine, implementation is for a general ds value
if k<=ds || (p-k)<ds
   ds=1; 
   %ds=min(k,p-k); % protection 
end
Xc=setdiff(1:p,X);  % possible set of variables to be added from
ndsCombSupp=nchoosek(k,ds); % no. of possible 2 combinations from supp set
ndsCombSuppC=nchoosek(p-k,ds);  % no. of possible 2 combinations from suppC set
alldsCombSupp=nchoosek(X,ds);  % get all the possible combinations of ds variables from supp set
alldsCombSuppC=nchoosek(Xc,ds);  % get all the possible combinations of ds variables from suppC set
     
    while stop==0
        numOfIter=numOfIter+1;
        % minimize the gain by dropping ds variable
        [fXminusds,~,J]=leastSigniX(fX,X,k,ds,ndsCombSupp,alldsCombSupp, A,d,c,diagInv, iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);
        % istar are the indices to drop,  fXminusi=f(X - istar)

        % maximize the reduction by adding ds variable to Xminusi
        Xminusds=setdiff(X,J);
        [~,fXplusj,~,Q]=mostSigniXc(fXminusds,Xminusds,k-ds,ds,ndsCombSuppC,alldsCombSuppC, A,d,c,diagInv,  iPara,rPara,whichalgo,isSoftStop, isTF,targetfbest);
        % jstar index to add, fXplusj=f(X + jstar)

        % Swap if possible
        if fX>fXplusj
           X=union(Xminusds,Q);  % new support
           Xc=union( setdiff(Xc,Q) , J); 
           for i=1:ds
              alldsCombSupp(alldsCombSupp==J(i))=Q(i); 
              alldsCombSuppC(alldsCombSuppC==Q(i))=J(i);
           end
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
    [~,x(sortX),fx]=quad_min_algo_pool(k,A(sortX,sortX),d(sortX),c,[],[],xRelaxedOpt(sortX),Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);

end %==========================================================================================================================================================================

%% mostSigniXc
function [xout,fstar,qstar,Qsupp]=mostSigniXc(fX,X,ki,ds,ndsCombSuppC,alldsCombSuppC, A,d,c,diagInv,  iPara,rPara,whichalgo,isSoftStop, isTF,targetfbest)
% card(X)=ki, 
% find the most significant feature w.r.t the set X complement
% i.e. find the feature x* such that S(x*)= max (f(X) - f(X + xj) )  for xj in Xc
% ki=no. of entries in X or predictors already in the model
    pS=-inf; 
    x0=zeros(ki+ds,1);
    for i=1:ndsCombSuppC
       Xplusds=union(X,alldsCombSuppC(i,:)); % include one feature from Xc at a time;

       [Q,diagD]=eig(A(Xplusds,Xplusds),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
       Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
       Sstruct.Q=Q(:, (ki+ds-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
       Sstruct.D=diagD( (ki+ds-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
       Sstruct.diagInv=diagInv(Xplusds);

       [~,xopt,fXplusds]=quad_min_algo_pool(ki+ds,A(Xplusds,Xplusds),d(Xplusds),c,[],[],x0,Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);
       S=fX-fXplusds;
       if pS<S
          pS=S; Qsupp=alldsCombSuppC(i,:); qstar=i; fstar=fXplusds; xout=xopt; 
       end
    end

end % end mostSigniXc===================================================================

%% leastSigniX
function [fstar,jstar,J]=leastSigniX(fX,X,ki,ds,ndsCombSupp,alldsCombSupp, A,d,c,diagInv,  iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest) 
% find the two least significant features w.r.t. the set X
% i.e. find the feature x* such that S(J)= min ( f(X-J)-f(X) ) for xj in X
    
    if ki==1
       J=X;jstar=1;fstar=c;return;
    end
    pS=inf; 
    
    x0=zeros(ki-ds,1);
    for i=1:ndsCombSupp
       Xminds=setdiff(X,alldsCombSupp(i,:),'stable'); % delete ds feature from X at a time;
      
       [Q,diagD]=eig(A(Xminds,Xminds),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
       Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
       Sstruct.Q=Q(:, (ki-ds-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
       Sstruct.D=diagD( (ki-ds-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
       Sstruct.diagInv=diagInv(Xminds);
       [~,~,fXminds]=quad_min_algo_pool(ki-ds,A(Xminds,Xminds),d(Xminds),c,[],[],x0,Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);
       S=fXminds-fX;
       if S<pS
          pS=S;J=alldsCombSupp(i,:);jstar=i; fstar=fXminds; 
       end
    end

end % end leastSigniX=======================================================================================

















