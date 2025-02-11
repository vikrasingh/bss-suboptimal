function [rParaOut,fx,x,fXarray,X]=sffs(p,k,A,d,c,lb,ub,diagInv,iPara,rPara,IotherPara,IstopCondPara,targetfbest)
% last updated 13Dec23, for quadminpackage with c0 and Struct
% version 2 on 19April23
% 24June22, Sequential Forward floating selection, from Pudil 1994Floating seach methods in feature selection

% parameters for quad. min.
todebug=0;
whichalgo=IotherPara(23); % selected algo. to fit the data in the reduced dim.
isTF=0; % whether to use target fbest during quad_min in the reduced space, we are not using target fbest inside quad.min. call
isSoftStop=1; % to use soft stop while computing min. in the reduced space or not
firstDigit=floor( IotherPara(21)/(10^(length(num2str(IotherPara(21)))-1)) );  % first digit of IotherPara(21) will tell use which version to use. 
iPara(2)=IotherPara(21)-firstDigit*( 10^(length(num2str(IotherPara(21)))-1) ); % modify iPara for refinement
stopflag=0; 
timeLimit=abs(IstopCondPara(6)); % hard cputime limit for algo.
itime=cputime; % to save the cputime for the algorithm
iter=0;
psuppk=[];
%X=zeros(1,k);  %set to hold the indices of the selected features 
%fXarray=zeros(1,k); % will store the f value f(X(1:i)) for i=1:k 

%Step 0 : Sequentially select the first 2 features using forward stepwise selection
[~,xhat,fXarray(1),X(1)]=fs1(p,1,A,d,c,lb,ub,diagInv,iPara,rPara,IotherPara,IstopCondPara,targetfbest);
if todebug==1
fprintf('p=%d, tmax=%d ========================= \n',p,k);    
fprintf('X=');printArray(X,'%d');
fprintf('fX=');printArray(fXarray,'%1.5f');
end

if k==1
   x=xhat;fx=fXarray;
   rParaOut.numOfIter=iter;rParaOut.necConMaxVioRefQM=-1;rParaOut.stopflag=stopflag;
   return
end
Xc=setdiff(1:p,X(1)); % complement of the current support set X(1:ki)
[xhat,fXarray(2),X(2)]=mostSigniXc(fXarray(1),X(1),Xc,p,1, A,d,c,lb,ub,diagInv,  iPara,rPara,whichalgo,isSoftStop, isTF,targetfbest);
if todebug==1
fprintf('X=');printArray(X,'%d');
fprintf('fX=');printArray(fXarray,'%1.5f');
end

if k==2
   X=sort(X); 
   x=zeros(p,1);x(X)=xhat;fx=fXarray(2);
   rParaOut.numOfIter=iter;rParaOut.necConMaxVioRefQM=-1;rParaOut.stopflag=stopflag;
   return 
end
ki=2; % current features in the set X

    while ki~=k
        iter=iter+1;
%         if iter==12
%            disp('here'); 
%         end
        if todebug==1
        fprintf('iter=%d; ki=%d ============== \n',iter,ki);
        fprintf('X=');printArray(X(1:ki),'%d');
        fprintf('fX=');printArray(fXarray(1:ki),'%1.9f');
        end
        %Step 1 : Inclusion
        % istarlocal is the variable selected to be in the support in the reduced space
        Xc=setdiff(1:p,X(1:ki)); % complement of the current support set X(1:ki)
        [~,fXkplus1,istar]=mostSigniXc(fXarray(ki),X(1:ki),Xc,p,ki, A,d,c,lb,ub,diagInv,  iPara,rPara,whichalgo,isSoftStop, isTF,targetfbest);
        Xkplus1=union(X(1:ki),istar,'stable');
        if todebug==1
        fprintf('istar=%d\n',istar);
        fprintf('Xkplus1=');printArray(Xkplus1,'%d');
        end
        
        % save the support if X has select k features at some point
        if (ki+1)==k
           psuppk=Xkplus1; 
        end

        %Step 2 : Conditional Exclusion
        [fstar,xr]=leastSigniX(fXkplus1,Xkplus1,p,ki+1, A,d,c,lb,ub,diagInv, iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);
        % here fstar= f(Xkplus1 - xr)
        if istar==xr
           ki=ki+1;X(1:ki)=Xkplus1;fXarray(ki)=fXkplus1; 
           if todebug==1, fprintf('istar=xr \n'); end
           continue; % go to step 1
        else
           Xkdash=setdiff(Xkplus1,xr,'stable');
           if todebug==1, fprintf('Xkdash=');printArray(Xkdash,'%d'); end
           if ki==2
              X(1:2)=Xkdash;fXarray(2)=fstar; 
              continue; % go to step 1
           end
           for id=1:(ki-1)  % 25June  update the fXarray for the previous models that get affected
              if  Xkdash(id)~=Xkplus1(id)
                  nred=id; 
                  sortX=sort(Xkdash(1:id));
                  [Q,diagD]=eig(A(sortX,sortX),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
                  Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
                  Sstruct.Q=Q(:, (nred-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
                  Sstruct.D=diagD( (nred-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
                  Sstruct.diagInv=diagInv(sortX); 
                  [~,~,fXarray(id)]=quad_min_algo_pool(nred,A(sortX,sortX),d(sortX),c,lb(sortX),ub(sortX),zeros(nred,1),Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest); 
              end
           end
           fXarray(ki)=fstar;
           if todebug==1, fprintf('fX=');printArray(fXarray(1:ki),'%1.9f'); end
        end

        % Step 3 : Continue conditinal exclusion
        step3=0;
        while step3==0
            
            [gstar,xs]=leastSigniX(fstar,Xkdash,p,ki, A,d,c,lb,ub,diagInv, iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);
            % gstar= f(Xkdash-xs)
            if gstar>=fXarray(ki-1)  
               X(1:ki)=Xkdash;    %fXarray(ki)=fstar;
               step3=1;% go to step 1
               if todebug==1
                  fprintf('gstar>=fX(ki-1) \n'); 
               end
            else
                if todebug==1, fprintf('Further reduce Xkdash\n'); end
                Xkdash1=setdiff(Xkdash,xs,'stable');
                ki=ki-1;
                if ki==2
                   X(1:2)=Xkdash1;   %fXarray(2)= fstar;
                   step3=1;% go to step 1
                   if todebug==1
                      fprintf('ki=2, continue \n'); 
                   end
                else
                   for id=1:(ki-1)  % 25June  update the fXarray for the previous models that get affected
                      if  Xkdash1(id)~=Xkdash(id)
                          nred=id; 
                          sortX=sort(Xkdash1(1:id));
                          [Q,diagD]=eig(A(sortX,sortX),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
                          Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
                          Sstruct.Q=Q(:, (nred-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
                          Sstruct.D=diagD( (nred-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
                          Sstruct.diagInv=diagInv(sortX); 
                          [~,~,fXarray(id)]=quad_min_algo_pool(nred,A(sortX,sortX),d(sortX),c,lb(sortX),ub(sortX),zeros(nred,1),Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest); 
                      end
                   end
                   Xkdash=Xkdash1;
                   fstar=gstar;
                   fXarray(ki)=fstar;
                   if todebug==1
                      fprintf('ki=%d \n',ki); 
                      fprintf('Xkdash=');printArray(Xkdash,'%d');
                      fprintf('fXarray=');printArray(fXarray(1:ki),'%1.9f');
                   end
                end
            end

            ctime=(cputime-itime)/60; % intermediate check for cputime limit
            if ctime > timeLimit  % if the cputime becomes greater than timeLimit=|IstopCondPara(6)| min, then stop
               stopflag=6;
               rParaOut.numOfIter=iter;rParaOut.necConMaxVioRefQM=-1;rParaOut.stopflag=stopflag;
               break;
            end
            
        end % end while step3=0

        ctime=(cputime-itime)/60; % intermediate check for cputime limit
        if ctime > timeLimit  % if the cputime becomes greater than timeLimit=|IstopCondPara(6)| min, then stop
           stopflag=6;
           rParaOut.numOfIter=iter;rParaOut.necConMaxVioRefQM=-1;rParaOut.stopflag=stopflag;
           break;
        end
        %fprintf('ki=%d, |supp|=%d \n',ki,length(X(X~=0)));
    end % while ki=/k
    
    if stopflag==6
       X=X(X~=0);
       nred=length(X);
       if nred<k && ~isempty(psuppk)  % if selected support is less than k, use previous set with k support if available
          X=psuppk(psuppk~=0);
          nred=length(X);
       end
    else, nred=k;   
    end
    % fit the data
    x=zeros(p,1);
    sortX=sort(X);
    [Q,diagD]=eig(A(sortX,sortX),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
    Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
    Sstruct.Q=Q(:, (nred-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
    Sstruct.D=diagD( (nred-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
    Sstruct.diagInv=diagInv(sortX);
    [~,x(sortX),fx]=quad_min_algo_pool(nred,A(sortX,sortX),d(sortX),c,lb(sortX),ub(sortX),zeros(nred,1),Sstruct,iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest);
    rParaOut.numOfIter=iter;rParaOut.necConMaxVioRefQM=-1;rParaOut.stopflag=stopflag;
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

%% leastSigniX
function [fstar,istar]=leastSigniX(fX,X,p,ki, A,d,c,lb,ub,diagInv,  iPara,rPara,whichalgo,isSoftStop,isTF,targetfbest) 
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
















