function [rParaOut,supp,xbest,fbest]=geneticAlgBSS(pDim,lb,ub,Anorm,bnorm,cnorm,tmax,normxRelaxedOpt,preXstarNormSpace,diagInv,...
                                                    targetfbest, IotherPara,IstopCondPara,iPara,rPara,textfileName,toDebug)
% 4 March24, genetic algorithm to solve BSS problem

   if ismember(IotherPara(23),12:20)
      isBoxed=0;  % are we using unconstrained or constrained version, if isBoxed=1 constrained, isBoxed=0 unconstrained 
      funhand=@(x1,x2) fitData(x1,x2,Anorm,bnorm,cnorm,[],[],normxRelaxedOpt,diagInv,targetfbest,iPara,rPara,IotherPara(23),0);  % isSoftStop=0
   else
      isBoxed=1;  % are we using unconstrained or constrained version, if isBoxed=1 constrained, isBoxed=0 unconstrained  
      funhand=@(x1,x2) fitData(x1,x2,Anorm,bnorm,cnorm,lb,ub,normxRelaxedOpt,diagInv,targetfbest,iPara,rPara,IotherPara(23),0);  % isSoftStop=0 
   end
   
   
   [iParaGA,rParaGA]=defineGApara(IotherPara);  % read the parameters from a separate file
   firstDigit=floor( IotherPara(21)/(10^(length(num2str(IotherPara(21)))-1)) );  % first digit of IotherPara(21) will tell use which box to use for refinement.
   iParaGA(2)=IotherPara(21)-firstDigit*(10^(length(num2str(IotherPara(21)))-1)); % remaining digits will tell us the no. of iterations to use for refinement
   iParaGA(3)=tmax; iParaGA(4)=0; % toDebug=0 inside the binaryGA1viaCC
   iParaGA(5)=IstopCondPara(5); % max. iter for GA 09/13/2024
   rParaGA(3)=IstopCondPara(6)*60; % max. cputime limit for GA in seconds 10/1/24
   idx1=min(IotherPara(24),pDim-tmax+1); % to pick a new support, starting from idx1Supp 
   [~,initialSuppnRun]=sort(abs(normxRelaxedOpt),'descend');
   cpuStart=cputime;
   initialSupp=initialSuppnRun( idx1:tmax-1+idx1 );
   x0GA=zeros(pDim,1,'logical');
   x0GA(initialSupp)=1; % initial binary vector for genetic algorithm
   [rParaOut,xbest,~] = binaryGA3viaCC(pDim,x0GA,iParaGA,rParaGA,funhand); 
   supp=find(xbest);  % find the indices of the non zero components
   nred=length(supp);   % reduce dim.
   [Q,diagD]=eig(Anorm(supp,supp),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
   Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
   Sstruct.Q=Q(:, (nred-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
   Sstruct.D=diagD( (nred-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
   Sstruct.diagInv=diagInv(supp);
   if isBoxed==1
      [epsMax,xtemp,fbest]=quad_min_algo_pool(nred,Anorm(supp,supp),bnorm(supp),cnorm,lb(supp),ub(supp),normxRelaxedOpt(supp),Sstruct,iPara,rPara,IotherPara(23),1,0,targetfbest);  % isSoftStop=1 soft stop, isTF=0 means do not use targetfbest 
   else
      [epsMax,xtemp,fbest]=quad_min_algo_pool(nred,Anorm(supp,supp),bnorm(supp),cnorm,[],[],normxRelaxedOpt(supp),Sstruct,iPara,rPara,IotherPara(23),1,0,targetfbest);  % isSoftStop=1 soft stop, isTF=0 means do not use targetfbest
   end
   xbest=zeros(pDim,1); 
   xbest(supp)=xtemp;
   rParaOut.necConMaxVioRefQM=epsMax;
   rParaOut.cpusec=cputime-cpuStart;

end