function [rParaOut,bestVal,betaOut]=dfo(p,n,k,Q,q,c,diagInv,normxRelaxedOpt,L,iPara,rPara,IotherPara)
% VS:23 June24, to use multiple initial points and choose the results corresponding the best solution
   
   num_runs=min(IotherPara(24),p-k+1); % no. of runs with different initial support for GA
   fbest_list=zeros(1,num_runs);
   pre_fbest=inf;
   [~,initial_supp]=sort(abs(normxRelaxedOpt),'descend');
   for irun=1:num_runs
       cpuStart=cputime;
       pbeta=normxRelaxedOpt;
       if irun>1
          pbeta(initial_supp(irun-1))=0;
       end
       [stopflag,bestVal,betaOut]=dfo1(p,n,k,Q,q,c,diagInv,pbeta,L,iPara,rPara,IotherPara);
       fbest_list(irun)=bestVal;
       cpuEnd=cputime-cpuStart;
       if bestVal<pre_fbest
          pre_fbest=bestVal;prebetaOut=betaOut;prestopflag=stopflag;precpu=cpuEnd; 
       end
   end
   bestVal=pre_fbest;betaOut=prebetaOut;stopflag=prestopflag;
   rParaOut.stopflag=stopflag; rParaOut.fbestList=fbest_list; rParaOut.cpusec=precpu;
   

end
%=======================================================================
function  [stopflag,bestVal,betaOut]=dfo1(p,n,k,Q,q,c,diagInv,pbeta,L,iPara,rPara,IotherPara)
% Sec3. Algo.1 page 830 from Bertsimas 2016 BSS via a modern optimization lens
% based on first order discrete necessary conditions.

% L>=l where || f(x)-f(y) ||<= l || x-y ||, in the ref. p833 l=max. eigenvalue of (X'X)
% L=max(eig(X'*X)) where Q=X'X for us
%pbeta is the previous solution if available, else a vector of zeros with first tmax no. of 1

% 4 Sep 2023, adding the structure below to be passed to quadratic minimization package
%  S is a structure S.Q = Q matrix from spectral decomposition such that A=QDQ'  where Q is an orthonormal matrix and D is diagonal with eig. value in the decreasing order 
%                  S.nNonZeroEig = no. of non zero eigenvalues  
%                  S.D = diag vector  of order nNonZeroEig x 1 of the D matrix 
%                  S.diagInv= 1/diag(A)  nx1 vector

    maxiter=1000; % as used in their implementation
    iter=0; % initialization
    nruns=20; % 
    epstol=1e-4; % this tolerance is suggested in the reference p833 
    stopflag=0;  % =0 means full convergence, =5 means, maxiter has been reached for atleast 1 of the random runs
    %1. initial solution b with ||b||_0 <=k
    if isempty(pbeta)
       if p<n
           sol=Q\(-q);
       else
          sol=q./(-diag(Q));
       end
       [~,isupp]=maxk(abs(sol),k);  % find the largest k components 
       pbeta=zeros(p,1);pbeta(isupp)=sol(isupp);
    else 
       [~,isupp]=maxk(abs(pbeta),k);  % find the largest k components 
       sol=zeros(p,1);sol(isupp)=pbeta(isupp);
       pbeta=sol;
    end
    beta0=pbeta;
    bestVal=Inf; % initialization
    betaOut=beta0; % final solution to output
    for i=1:nruns
        %2. run the loop until converges 
        while iter<=maxiter
            
            gradvec=0.5*( Q*pbeta+q );
            v=pbeta-(1/L)*gradvec;   % v in Hk(v)
            [~,idx]=sort(abs(v),'descend');
            
            beta=zeros(p,1);beta(idx(1:k))=v(idx(1:k));  % new beta
            currentSupp=idx(1:k);
           
            [Qorthonormal,diagD]=eig(Q(currentSupp,currentSupp),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
            Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
            Sstruct.Q=Qorthonormal(:, (k-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
            Sstruct.D=diagD( (k-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
            Sstruct.diagInv=diagInv(currentSupp);

            [~,beta(currentSupp),~]=quad_min_algo_pool(k,Q(currentSupp,currentSupp),q(currentSupp),c,[],[],beta(currentSupp),Sstruct,iPara,rPara,IotherPara(23),1,0,[]);  % isSoftStop=1, isTF=0 means do not use targetfbest
            
            df=norm(beta-pbeta,1)/max(norm(beta,1)); % using L_1 norm
            
            if df<epstol
               break; 
            end
            pbeta=beta;
            iter=iter+1;
        end
        if iter==maxiter, stopflag=5;end
        currentVal=fx(beta,p,Q,q,c);
        
        if currentVal<bestVal
           bestVal=currentVal;
           betaOut=beta;
        end
        pbeta=beta0+2*rand(p,1)*max(abs(beta0));  % start the next run from a new random point
        iter=0; % reset the counter
    end

end
%============================================================
