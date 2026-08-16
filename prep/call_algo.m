function [outputPara,rParaOut,Xstar,fbest]=call_algo(nPts,pDim,yArray,xMatrix,A,b,c,tmax,...
                                                        IstopCondPara,IotherPara,iPara,rPara,version_flag,xRelaxedOpt,textfileName,toDebug)
    
    
    Mu=-1; % dummy initialization

    scaleQP=max( max(max(abs(A))) , max(abs(b)) ); % a scaling factor for quad. min.
    scaleQP=scaleQP/10; % we want the values to be less than 10 for better computation results
    Anorm=A/scaleQP;
    bnorm=b/scaleQP;
    cnorm=c/scaleQP;
        
    diagA=diag(Anorm);
    diagInv=zeros(pDim,1);
    for idiag=1:1:pDim
       if diagA(idiag)>rPara(3)
          diagInv(idiag)=1/diagA(idiag);
       end
    end
    targetfbest=0;
    % main algo. call
    %=======================================================
    if version_flag==9  % sfs1
        outputPara=nan(17,1);timeGivenalgo=cputime;
        [rParaAG9,~,x_tilde,fx_tilde,Xstar,fbest]=sfs1(pDim,Anorm,bnorm,cnorm,tmax,xRelaxedOpt,[],diagInv,targetfbest, IotherPara,IstopCondPara,iPara,rPara,textfileName,toDebug);
        cpuIntvalAlgo=(cputime-timeGivenalgo); % in sec
        outputPara(16)=rParaAG9.necConMaxVioRefQM; outputPara(5)=rParaAG9.numOfIter;outputPara(7)=rParaAG9.stopflag;

    elseif version_flag==92  % sfs2
        outputPara=nan(17,1);timeGivenalgo=cputime;
        [rParaAG92,~,x_tilde,fx_tilde,Xstar,fbest]=sfs2(pDim,Anorm,bnorm,cnorm,tmax,xRelaxedOpt,[], diagInv, targetfbest, IotherPara,IstopCondPara,iPara,rPara,textfileName,toDebug);
        cpuIntvalAlgo=(cputime-timeGivenalgo); % in sec
        outputPara(16)=rParaAG92.necConMaxVioRefQM; outputPara(5)=rParaAG92.numOfIter;outputPara(7)=rParaAG92.stopflag;

    elseif version_flag==93 % sffs
        fx_tilde=nan;
        outputPara=nan(17,1);timeGivenalgo=cputime;
        [rParaOut93,fbest,Xstar,~,~]=sffs(pDim,tmax,Anorm,bnorm,cnorm,diagInv,iPara,rPara,IotherPara,IstopCondPara,targetfbest); % the given box is just a dummy input
        %** the Xstar,fbest above are already in the original space as we used original xMatrix,yArray, A,b and c
        cpuIntvalAlgo=(cputime-timeGivenalgo); % in sec
        outputPara(7)=rParaOut93.stopflag; outputPara(5)=rParaOut93.numOfIter;
        
    elseif version_flag==94 % ga
        outputPara=nan(17,1);fx_tilde=nan;
        [rParaAG94,~,Xstar,fbest]=ga(pDim,Anorm,bnorm,cnorm,tmax,xRelaxedOpt,[],diagInv,targetfbest, IotherPara,IstopCondPara,iPara,rPara,textfileName,toDebug);
        cpuIntvalAlgo=rParaAG94.cpusec;
        outputPara(16)=rParaAG94.necConMaxVioRefQM; outputPara(5)=rParaAG94.numOfIter;outputPara(7)=rParaAG94.stopflag;  %printArray(rParaAG94.fbestList,'%1.6f');
        
    elseif version_flag==34  % my implementation of BSS
        outputPara=nan(17,1); % initialization
        if min(nPts,pDim)<tmax
            Xstar=-ones(pDim,1);fbest=-1;cpuIntvalAlgo=0;fx_tilde=-1;x_tilde=-ones(pDim,1);
        else
            Mu=max(abs(xRelaxedOpt)); % Mu=max( |xRelaxedOpt| )
            Mu=0.01*IotherPara(2)*Mu;   % Enlargement of the box given
            xRelaxedOpt2=xRelaxedOpt+sign(xRelaxedOpt)*Mu; % xRelaxedOpt2 used to enlarge the box
            upbnd=abs(xRelaxedOpt2);lowbnd=-upbnd;
           timeGivenalgo=cputime;
           [outputPara(7),fx_tilde,x_tilde,Xstar,fbest]=mio(pDim,nPts,xMatrix,yArray,tmax,Anorm,bnorm,cnorm,lowbnd,upbnd,diagInv,abs(IstopCondPara(6)),iPara,rPara,IotherPara,IstopCondPara,xRelaxedOpt,targetfbest );
           % x_tilde is the solution from projected gradient descent, used as a starting point for gurobi solver
           cpuIntvalAlgo=(cputime-timeGivenalgo); % in sec
        end
        
    elseif version_flag==35 % fs 
        outputPara=nan(17,1);fx_tilde=nan;
        timeGivenalgo=cputime;
        [outputPara(7),Xstar,fbest,~]=fs(pDim,tmax,Anorm,bnorm,cnorm,diagInv,iPara,rPara,IotherPara,IstopCondPara,targetfbest);
        cpuIntvalAlgo=(cputime-timeGivenalgo); % in sec
        
    elseif version_flag==36 % discrete first order method to solve the BSS
        stepL=max(eig(Anorm/2)); % step length for projected gradient method
        outputPara=nan(17,1);fx_tilde=[];
        [rParaOutdfo1,fbest,Xstar]=dfo(pDim,nPts,tmax,Anorm,bnorm,cnorm,diagInv,xRelaxedOpt,stepL,iPara,rPara,IotherPara);
        cpuIntvalAlgo=rParaOutdfo1.cpusec;
        outputPara(7)=rParaOutdfo1.stopflag;
        
    elseif version_flag==37  % abess   
        outputPara=nan(17,1); fx_tilde=nan; 
        tmax1=min(tmax, min(nPts,round( nPts/(log(pDim)*log(log(nPts))) ))); % in abess smax=min(n, [n/logp*loglogn]) 
        inputFile=sprintf('%s/input.mat',pwd);
        outputFile=sprintf('%s/output.mat',pwd);
        save(inputFile, 'xMatrix', 'yArray', 'tmax1');
        Rscript = '/apps/spack/anvil/apps/r/4.4.1-gcc-11.2.0-kth7vej/bin/Rscript';
        modelUbxR = '/home/x-vsingh7/lls_package/utility/R_files/modelUbx.R';
        cmd = sprintf('%s %s %s %s', Rscript, modelUbxR, inputFile, outputFile);
        %cmd = sprintf('Rscript abessR.R xMatrix yArray tmax lowbnd upbnd');
        [status, result] = system(cmd);
        disp(result)
        if status ~= 0
           error('R failed');
        end
        res = load(outputFile);
        Xstar = res.beta(2:pDim+1,end);
        cpuIntvalAlgo=res.ctime;
        fbest=fx(Xstar,pDim,A,b,c);
        %cpuIntvalAlgo=-1; % measured in R function only
             
    end % main algo. call
    %================================================================================

    if toDebug==1,fprintf('cputime for the interval algo. in min. = %1.8f \n',cpuIntvalAlgo/60);end
    
    %% save violation of the Xstar
    outputPara(8)=length(find(abs(Xstar)>eps))-tmax; % get constraint vio by Xstar in the normalized space.
      
    %% Convert Xstar and x_tilde back to the original space
    if version_flag~=34 % if not using QRD to find fbest, we are not using scaleQP for AG34
      fbest=scaleQP*fbest; % only for the quadratic function
    end
   
    rParaOut.scaleQP=scaleQP;
    rParaOut.cpuIntvalAlgo=cpuIntvalAlgo;
    rParaOut.Mu=Mu;

%===========================================================================================================================================================
end  % end of the function call_algo.
%===========================================================================================================================================================
