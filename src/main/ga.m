function [rParaOut,supp,xbest,fbest]=ga(pDim,Anorm,bnorm,cnorm,tmax,normxRelaxedOpt,preXstarNormSpace,diagInv,...
                                                    targetfbest, IotherPara,IstopCondPara,iPara,rPara,textfileName,toDebug)
% VS:4 March24, genetic algorithm to solve BSS problem

   
   funhand=@(x1,x2) fitData(x1,x2,Anorm,bnorm,cnorm,[],[],normxRelaxedOpt,diagInv,targetfbest,iPara,rPara,IotherPara(23),0);  % isSoftStop=0
   
   [iParaGA,rParaGA]=defineGApara(IotherPara);  % read the parameters from a separate file
   firstDigit=floor( IotherPara(21)/(10^(length(num2str(IotherPara(21)))-1)) );  % first digit of IotherPara(21) will tell use which box to use for refinement.
   iParaGA(2)=IotherPara(21)-firstDigit*(10^(length(num2str(IotherPara(21)))-1)); % remaining digits will tell us the no. of iterations to use for refinement
   iParaGA(3)=tmax; iParaGA(4)=0; % toDebug=0 inside the binaryGA1viaCC

   nruns=min(IotherPara(24),pDim-tmax+1); % no. of runs with different initial support for GA
   fbestList=zeros(1,nruns);
   prefbest=inf;
   [~,initialSuppnRun]=sort(abs(normxRelaxedOpt),'descend');
   for irun=1:nruns
       cpuStart=cputime;
       initialSupp=initialSuppnRun( irun:tmax-1+irun );
       x0GA=zeros(pDim,1);
       x0GA(initialSupp)=1; % initial binary vector for genetic algorithm
       [rParaOut,xbest,fbest] = binaryGA1viaCC(pDim,x0GA,iParaGA,rParaGA,funhand);
       cpuEnd=cputime-cpuStart;
       fbestList(irun)=fbest;
       if fbest<prefbest
          prexbest=xbest;prerParaOut=rParaOut;precpu=cpuEnd; 
       end
   end
   cpuStart1=cputime;
   xbest=prexbest; rParaOut=prerParaOut;
   supp=find(xbest);  % find the indices of the non zero components
   nred=length(supp);   % reduce dim.
   [Q,diagD]=eig(Anorm(supp,supp),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
   Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
   Sstruct.Q=Q(:, (nred-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
   Sstruct.D=diagD( (nred-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
   Sstruct.diagInv=diagInv(supp);
   [epsMax,xbest(supp),fbest]=quad_min_algo_pool(nred,Anorm(supp,supp),bnorm(supp),cnorm,[],[],normxRelaxedOpt(supp),Sstruct,iPara,rPara,IotherPara(23),1,0,targetfbest);  % isSoftStop=1 soft stop, isTF=0 means do not use targetfbest
  
   rParaOut.necConMaxVioRefQM=epsMax;
   rParaOut.fbestList=fbestList; rParaOut.cpusec=precpu+(cputime-cpuStart1);

end %===================================================