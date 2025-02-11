function [XstarNormSpace,xRelaxedOpt2,outputPara,rParaOut,fx_tilde,x_tilde,Xstar,fbest]=setupForIntvalAlgo(nPts,pDim,yArray,xMatrix,trueb,...
                                                         A,b,c,tmax,targetfbest,IstopCondPara,IotherPara,iPara,rPara,version_flag,xRelaxedOpt,...
                                                               fxRelaxedOpt,textfileName,toDebug)
    
    Mu=Inf; % Mu will store the actual bound used for enlargement to define the initial box
    % For BSS algorithm, check if tmax is 0 or the xRelaxedOpt is already feasible
    if tmax==0
        Xstar=-ones(pDim,1);fbest=-1;normFactor=zeros(pDim,1);xRelaxedOpt2=xRelaxedOpt;xRelaxOptConstVio=0;
        outputPara=(-1)*ones(17,1);volDeleted=-1;cpuAlg=0;XstarNormSpace=Xstar;fx_tilde=-1;x_tilde=Xstar;
        disp('Error: tmax should be > 0.');
        rParaOut.normFactor=normFactor;rParaOut.scaleQP=nan;rParaOut.xRelaxOptConstVio=xRelaxOptConstVio;rParaOut.volDeleted=volDeleted;rParaOut.cpuIntvalAlgo=cpuAlg;rParaOut.tmax=tmax;rParaOut.Mu=Mu;
        return
    end
    numOf0sInxRelaxedOpt=sum(xRelaxedOpt==0); 
    if numOf0sInxRelaxedOpt>=(pDim-tmax)
        Xstar=xRelaxedOpt;fbest=fxRelaxedOpt;fx_tilde=fxRelaxedOpt;x_tilde=xRelaxedOpt;
        outputPara=zeros(17,1);volDeleted=0;cpuAlg=0;XstarNormSpace=xRelaxedOpt;normFactor=zeros(pDim,1);xRelaxedOpt2=xRelaxedOpt;
        xRelaxOptConstVio=0;
        disp('Return: xRelaxedOpt is feasible.');
        rParaOut.normFactor=normFactor;rParaOut.scaleQP=nan;rParaOut.xRelaxOptConstVio=xRelaxOptConstVio;rParaOut.volDeleted=volDeleted;rParaOut.cpuIntvalAlgo=cpuAlg;rParaOut.tmax=tmax;rParaOut.Mu=Mu;
        return
    end
   
    % Initial Box
    xRelaxOptConstVio=inf;
    normFactor=ones(pDim,1); % initialization
    if version_flag==34
        Mu=max(abs(xRelaxedOpt)); % Mu=max( |xRelaxedOpt| )
        Mu=0.01*IotherPara(2)*Mu;   % Enlargement of the box given
        xRelaxedOpt2=xRelaxedOpt+sign(xRelaxedOpt)*Mu; % xRelaxedOpt2 used to enlarge the box
        normFactor=max(abs(xRelaxedOpt2))*ones(pDim,1);  % max( |xRelaxedOpt2| )
               
        % set the box
        boundForBox=(abs(xRelaxedOpt2))./normFactor;
        upbnd=boundForBox;lowbnd=-upbnd;
    end

    % get xRelaxedOpt in norm space
    normxRelaxedOpt=xRelaxedOpt./normFactor;

    % get A,b,c in the normalized space 
    U=diag(normFactor);
    Anorm=U*A*U;      % A_norm,b_norm,c_norm in the normalized space. Will be used to evaluate the inclusion function F and f
    bnorm=(b'*U)';
    cnorm=c;
   
    % scale the data for numerical stability
    scaleQP=max( max(max(abs(Anorm))) , max(abs(bnorm)) );
    scaleQP=scaleQP/10; % we want the values to be less than 10 for better computation results
    Anorm=Anorm/scaleQP;
    bnorm=bnorm/scaleQP;
    cnorm=cnorm/scaleQP;
    targetfbest=targetfbest/scaleQP;  % we need to scale the targetfbest value as well
   
    % find diag inverse of Anorm
    diagA=diag(Anorm);
    diagInv=zeros(pDim,1);
    for idiag=1:1:pDim 
       if diagA(idiag)>rPara(3)
          diagInv(idiag)=1/diagA(idiag);  
       end
    end
    
    
    % main algorithms
    if version_flag==9 % sfs1
       outputPara=nan(17,1);
       timeGivenalgo=cputime;
       [rParaAG9,~,x_tilde,fx_tilde,Xstar,fbest]=AG9(pDim,[],[],Anorm,bnorm,cnorm,tmax,normxRelaxedOpt,diagInv,targetfbest, IotherPara,IstopCondPara,iPara,rPara);
       cpuAlg=(cputime-timeGivenalgo); % in sec
       outputPara(16)=rParaAG9.necConMaxVioRefQM; outputPara(5)=rParaAG9.numOfIter;outputPara(7)=rParaAG9.stopflag;

    elseif version_flag==92  % sfs2
        outputPara=nan(17,1);
        timeGivenalgo=cputime;
        [rParaAG92,~,x_tilde,fx_tilde,Xstar,fbest]=AG92(pDim,Anorm,bnorm,cnorm,tmax,normxRelaxedOpt, diagInv, targetfbest, IotherPara,IstopCondPara,iPara,rPara);
        cpuAlg=(cputime-timeGivenalgo); % in sec
        outputPara(16)=rParaAG92.necConMaxVioRefQM; outputPara(5)=rParaAG92.numOfIter;outputPara(7)=rParaAG92.stopflag; 

    elseif version_flag==93 % sffs
        outputPara=nan(17,1);
        timeGivenalgo=cputime;
        [rParaOut93,fbest,Xstar,~,~]=sffs(pDim,tmax,Anorm,bnorm,cnorm,diagInv,iPara,rPara,IotherPara,IstopCondPara,targetfbest); % the given box is just a dummy input
        %** the Xstar,fbest above are already in the original space as we used original xMatrix,yArray, A,b and c 
        cpuAlg=(cputime-timeGivenalgo); % in sec
        outputPara(7)=rParaOut93.stopflag; outputPara(5)=rParaOut93.numOfIter; 

    elseif version_flag==94 % ga
        outputPara=nan(17,1);
%         timeGivenalgo=cputime;
        [rParaAG94,~,Xstar,fbest]=geneticAlgBSS(pDim,Anorm,bnorm,cnorm,tmax,normxRelaxedOpt,[],diagInv,targetfbest, IotherPara,IstopCondPara,iPara,rPara,textfileName,toDebug);
%         cpuIntvalAlgo=(cputime-timeGivenalgo); % in sec
        cpuAlg=rParaAG94.cpusec;
        outputPara(16)=rParaAG94.necConMaxVioRefQM; outputPara(5)=rParaAG94.numOfIter;outputPara(7)=rParaAG94.stopflag; %printArray(rParaAG94.fbestList,'%1.6f');

    elseif version_flag==35 % fs
        outputPara=nan(17,1);
        timeGivenalgo=cputime;
        [outputPara(7),Xstar,fbest,~]=fs1(pDim,tmax,Anorm,bnorm,cnorm,diagInv,iPara,rPara,IotherPara,IstopCondPara,targetfbest);
        cpuAlg=(cputime-timeGivenalgo); % in sec 

    elseif version_flag==36 % dfo
        stepL=max(eig(Anorm/2)); % step length for projected gradient method
        outputPara=nan(17,1);
        %timeGivenalgo=cputime;
        [rParaOutdfo1,fbest,Xstar]=dfo1wrapper(pDim,nPts,[],[],tmax,Anorm,bnorm,cnorm,diagInv,[],normxRelaxedOpt,stepL,iPara,rPara,IotherPara,IstopCondPara); 
        %cpuIntvalAlgo=(cputime-timeGivenalgo); % in sec
        cpuAlg=rParaOutdfo1.cpusec;
        outputPara(7)=rParaOutdfo1.stopflag;   %printArray(rParaOutdfo1.fbestList,'%1.6f');   

    elseif version_flag==34  % mio
        outputPara=nan(17,1); % initialization
        if min(nPts,pDim)<tmax
            Xstar=-ones(pDim,1);fbest=-1;cpuAlg=0;fx_tilde=-1;x_tilde=-ones(pDim,1);
        else 
           timeGivenalgo=cputime; 
           [outputPara(7),fx_tilde,x_tilde,Xstar,fbest]=mio(pDim,nPts,xMatrix*diag(normFactor),yArray,tmax,Anorm,bnorm,cnorm,lowbnd,upbnd,diagInv,abs(IstopCondPara(6)),iPara,rPara,IotherPara,IstopCondPara,normxRelaxedOpt,targetfbest ); 
           cpuAlg=(cputime-timeGivenalgo); % in sec
        end
        if toDebug==1,fprintf('CPU time for MIO = %1.8f min. \n',cpuAlg/60);end    

    end
 
    % xstar in normed space
    XstarNormSpace=Xstar;
    
    % save violation of the Xstar
    outputPara(8)=sum(logical(Xstar))-tmax;  % get constraint vio by Xstar in the normalized space.
     
    % Convert Xstar and x_tilde back to the original space
    Xstar=Xstar.*normFactor; % or Xstar=U*Xstar , it is the same, no need to update fbest
    if version_flag~=34 %  we are not using scaleQP for AG34
       fbest=scaleQP*fbest; % only for the quadratic function 
       fx_tilde=scaleQP*fx_tilde; % scale the function value of the initial point as well
    end
    x_tilde=x_tilde.*normFactor; 
    
    % set output parameters
    rParaOut.normFactor=normFactor;rParaOut.scaleQP=scaleQP;rParaOut.xRelaxOptConstVio=xRelaxOptConstVio;rParaOut.cpuAlg=cpuAlg;rParaOut.tmax=tmax;rParaOut.Mu=Mu;
          
%===========================================================================================================================================================    
end  % end of the function setupForIntvalAlgo.
%===========================================================================================================================================================
