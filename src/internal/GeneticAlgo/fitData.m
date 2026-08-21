function [fout,xout]=fitData(n,xbin,A,b,c,lb,ub,x0,diagInv,targetfbest,iPara,rPara,quadMinAlg,isSoftStop)
% A is of order n x n, b and xbin are of n x 1 order
% xbin is a binary vector, we need to fit the data in least squares sense in the reduced dim. determined by the support fo the xbin


    isTF=0; % do not use targetfbest for lb F inside quadratic min.      
    supp=find(xbin);
    newdim=length(supp); % dim of the reduced quadratic problem
    Ared=A(supp,supp); % initialize
    bred=b(supp); % initialize
%         cred=c; % no change
    if ~isempty(lb)
       lbred=lb(supp);ubred=ub(supp); 
    else, lbred=[];ubred=[];
    end
    %x0=zeros(newdim,1);  % initial point for lb f algo.
    
    % 5 Sep23
    [Q,diagD]=eig(Ared,'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
    Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
    Sstruct.Q=Q(:, (newdim-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
    Sstruct.D=diagD( (newdim-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
    Sstruct.diagInv=diagInv(supp);

    [~,xquadmin,fout]=quad_min_algo_pool(newdim,Ared,bred,c,lbred,ubred,x0(supp),Sstruct,iPara,rPara,quadMinAlg,isSoftStop,isTF,targetfbest);  % isSoftStop=1 means soft stop, isTF=0 do not use targetfbest 
    xout=zeros(n,1);xout(supp)=xquadmin; % get solution x in the original space
    %if necConMaxVioInfQM<epsMax, necConMaxVioInfQM=epsMax;end

end