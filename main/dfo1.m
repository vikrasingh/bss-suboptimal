function  [stopflag,bestVal,betaOut]=dfo1(p,n,X,y,k,Q,q,c,lb,ub,diagInv,XtX,pbeta,L,iPara,rPara,IotherPara)
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

    maxiter=iPara(2); % max. iter. limit
    iter=0; % initialization
    nruns=50; 
    epstol=1e-4; % this tolerance is suggested in the reference p833 
    stopflag=0;  % =0 means full convergence, =5 means, maxiter has been reached for atleast 1 of the random runs
    %1. initial solution b with ||b||_0 <=k
    if isempty(pbeta)
       if p<n
%           sol=XtX\(X'*y);
           sol=Q\(-q);
       else
%           rowsum=zeros(p,1);
%           for j=1:p
%               rowsum(j)=sum(X(:,j).^2);
%           end
%           sol=(X'*y)./rowsum;  
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
%     bestSupp=find(beta0); % current best support of the solution
    for i=1:nruns
        %2. run the loop until converges 
        while iter<=maxiter
            gradvec=0.5*( Q*pbeta+q );
%             gradvec=-X'*(y-X*pbeta);
            v=pbeta-(1/L)*gradvec;   % v in Hk(v)
            [~,idx]=sort(abs(v),'descend');
            
            beta=zeros(p,1);beta(idx(1:k))=v(idx(1:k));  % new beta
            currentSupp=idx(1:k);
            % polish the coefficients
%             beta(currentSupp)=XtX(currentSupp,currentSupp)\(X(:,currentSupp)'*y);
            
            % 5 Sep 2023
            [Qorthonormal,diagD]=eig(Q(currentSupp,currentSupp),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
            Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
            Sstruct.Q=Qorthonormal(:, (k-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
            Sstruct.D=diagD( (k-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
            Sstruct.diagInv=diagInv(currentSupp);

            if isempty(lb)  % in this case ub will also be empty, flag for unconstrained algo.
               [~,beta(currentSupp),~]=quad_min_algo_pool(k,Q(currentSupp,currentSupp),q(currentSupp),c,[],[],beta(currentSupp),Sstruct,iPara,rPara,IotherPara(23),1,0,[]);  % isSoftStop=1, isTF=0 means do not use targetfbest
            else    
               [~,beta(currentSupp),~]=quad_min_algo_pool(k,Q(currentSupp,currentSupp),q(currentSupp),c,lb(currentSupp),ub(currentSupp),beta(currentSupp),Sstruct,iPara,rPara,IotherPara(23),1,0,[]);  % isSoftStop=1, isTF=0 means do not use targetfbest
            end

            df=norm(beta-pbeta,1)/max(norm(beta,1)); % using L_1 norm
%             df=( fx(pbeta,p,Q,q,c)-fx(beta,p,Q,q,c) );
%             fprintf('df=%1.4f \n',df);
%             if isnan(df)
%                beta' 
%             end
            
            if df<epstol
               break; 
            end
            pbeta=beta;
            iter=iter+1;
        end
        if iter==maxiter, stopflag=5;end
        currentVal=fx(beta,p,Q,q,c);
%         funval=y-X*beta;
%         currentVal=funval'*funval;
        if currentVal<bestVal
           bestVal=currentVal;
           betaOut=beta;
        end
        pbeta=beta0+2*rand(p,1)*max(abs(beta0));  % start the next run from a new random point
        iter=0; % reset the counter
    end


    % as suggested in the reference, once the support set of beta gets
    % stabilized, polish the coefficients on the active set
%     activeAlgo=20; % use unbox CG to find the min. in the reduced space
%     [epsMax,xq,fout]=quad_min_algo_pool(k,Q(bestSupp,bestSupp),q(bestSupp),c,[],[],betaOut(bestSupp),iPara,rPara,activeAlgo,1,0,[]);  % isSoftStop=1, isTF=0 means do not use targetfbest
%     betaOut(bestSupp)=xq;


end

% % the above version provides better results
% function  [Mu,betaOut]=firOrdDisNecCond(p,n,A,b,c,k,L,iPara,rPara)
% % Sec3. Algo.1 page 830 from Bertsimas 2016 BSS via a modern optimization lens
% % based on first order discrete necessary conditions.
% 
% % L>=l where || f(x)-f(y) ||<= l || x-y ||, in the ref. p833 l=max. eigenvalue of (X'X)
% % L=max(eig(X'*X)) where Q=X'X for us
% %pbeta is the previous solution if available, else a vector of zeros with first tmax no. of 1
% 
%     maxiter=1000; % as used in their implementation
%     iter=0; % initialization
%     nruns=50; % 
%     epstol=1e-4; % this tolerance is suggested in the reference p833 
%     %1. initial solution b with ||b||_0 <=k
%    if p<n
%       sol=A\(-b);  % threshold least squares coef.
%    else
%       sol=b./(-diag(A));  % threshold marginal regss. coef.  
%    end
%    [~,idxsort]=sort(abs(sol),'descend'); % find the indices of max k elements
%    pbeta=zeros(p,1);pbeta(idxsort(1:k))=sol(idxsort(1:k));
%     
%     beta0=pbeta;
%     bestVal=Inf; % initialization
%     betaOut=beta0; % final solution to output
% %     bestSupp=find(beta0); % current best support of the solution
%     for i=1:nruns
%         %2. run the loop until converges 
%         while iter<=maxiter
%             gradvec=A*pbeta+b;
%             v=pbeta-(1/L)*gradvec;   % v in Hk(v)
%             [~,idx]=sort(abs(v),'descend');
%             
%             beta=zeros(p,1);beta(idx(1:k))=v(idx(1:k));  % new beta
%             currentSupp=idx(1:k);
%             % polish the coefficients
%             [~,beta(currentSupp),~]=quad_min_algo_pool(k,A(currentSupp,currentSupp),b(currentSupp),c,[],[],beta(currentSupp),iPara,rPara,12,1,0,[]);  % active algo. =20 unconst. CG ,isSoftStop=1, isTF=0 means do not use targetfbest
% 
%             df=norm(beta-pbeta,1)/max(norm(beta,1)); % using L_1 norm
% %             df=( fx(pbeta,p,Q,q,c)-fx(beta,p,Q,q,c) );
% %             fprintf('df=%1.4f \n',df);
% %             if isnan(df)
% %                beta' 
% %             end
%             
%             if df<epstol
%                break; 
%             end
%             pbeta=beta;
%             iter=iter+1;
%         end
%         currentVal=fx(beta,p,A,b,c);
%         if currentVal<bestVal
%            bestVal=currentVal;
%            betaOut=beta;
%         end
%         pbeta=beta0+2*rand(p,1)*max(abs(beta0));  % start the next run from a new random point
%         iter=0; % reset the counter
%     end
%     Mu=max(abs(betaOut));
% 
% end