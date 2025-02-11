function [bestAlgoFlag,lb,ub,solsOfQuadMin,cpuxRelaxedOpt,xRelaxedOpt,fxRelaxedOpt,A,b,c]=getxRelaxedOpt(pDim,yArray,xMatrix,trueb,iPara,rPara,QuadMinFunPara,IotherPara,folderPath)
%==============================================================================================================================================
% Subroutine to find the parameters A,b,c such that ||y-X*z||^2 = 0.5*z'A z + b'z + c , and 
% to find the relaxed optimal solution using a box or without using a box.
% ===========================================================================================================================================
% Calls the following subsroutines:
% quadMinFun.m
% iterWayToGetRelaxOpt.m
%==========================================================================================================================================
% Called in the following subroutines:
% runLLS.m, runLLSwithParSetLoop.m,runLLSwithParTmLoop,runLLSwithParEgLoop.m, runAnEgLLS.m
%==========================================================================================================================================
% 08/18/23 S is a structure S.Q = Q matrix from spectral decomposition such that A=QDQ'  where Q is an orthonormal matrix and D is diagonal with eig. value in the decreasing order 
%                  S.nNonZeroEig = no. of non zero eigenvalues  
%                  S.D = Diag matrix  of order nNonZeroEig x nNonZeroEig
%                  S.diagInv= 1/diag(A)  nx1 vector


if ~isempty(folderPath)
   fileid=fopen(fullfile(folderPath,sprintf('quadMinFunOutputDim%d.txt',pDim)),'w');
   fprintf(fileid,'%s\n',datetime('now'));
else, fileid=[];   
end

A=2*(xMatrix'*xMatrix);
b=-2*xMatrix'*yArray;
c=yArray'*yArray; 

% for debugging
eigenvalues=eig(A);
% if eigenvalues>=0, disp('PD problem.');
% elseif eigenvalues<=0, disp('ND problem.');
% else, disp('Indefinite problem.');
% end
scaleQP=max( max(max(abs(A))) , max(abs(b)) );
scaleQP=scaleQP/10; % we want the values to be less than 10 for better computation results
Anorm=A/scaleQP;
bnorm=b/scaleQP;
cnorm=c/scaleQP;
[Q,D]=eig(Anorm,'vector');
% vs 08/16/23
diagA=diag(Anorm);
diagInv=zeros(pDim,1);
for idiag=1:1:pDim 
   if diagA(idiag)>rPara(3)
      diagInv(idiag)=1./diagA(idiag);  
   end
end
S.nNonZeroEig=sum(D>eps);
S.Q=Q(:,(pDim-S.nNonZeroEig+1):end );  % Q is in the reduced space
% S.D=diag( D(1:S.nNonZeroEig , 1:S.nNonZeroEig) );  % D vector of eig/diag of D are in decreasing order
S.D=D(pDim-S.nNonZeroEig+1:end); 
S.diagInv=diagInv;

% added on 08/17/23
iParaRelaxX=iPara; rParaRelaxX=rPara; 
iParaRelaxX(2)=iPara(2)*rPara(9); rParaRelaxX(5)=rPara(5)*rPara(8);


if IotherPara(1)==0  % use the whole space, i.e. unconstrained algo. to find xRelaxedOpt
   % for BSS we are using unconstrained minima as xRelaxedOpt
   lb=[];ub=[]; % dummy variables
   select_only_ubox_algo=(QuadMinFunPara(1,:)==1);  % the first row of QuadMinFunpara
   select_only_ubox_algo(1:11)=0; % make the box constrained active algo. inactive
   [bestAlgoFlag,saveAllTheSols,cpuxRelaxedOpt,xRelaxedOpt,fxRelaxedOpt]=quadMinFun(pDim,Anorm,bnorm,cnorm,lb,ub,S,iParaRelaxX,rParaRelaxX,select_only_ubox_algo);
   solsOfQuadMin=saveAllTheSols(1:pDim,sum( select_only_ubox_algo(1:bestAlgoFlag) )  ); % save the best solution to be used in LLS

elseif IotherPara(1)==5  % use sequential/iterative way of determining xRelaxedOpt and fxRelaxedOpt
    
   [bestAlgoFlag,lb,ub,saveAllTheSols,cpuxRelaxedOpt,xRelaxedOpt,fxRelaxedOpt]=iterWayToGetRelaxOpt(pDim,Anorm,bnorm,cnorm,S,iParaRelaxX,rParaRelaxX,QuadMinFunPara,fileid);
   solsOfQuadMin=saveAllTheSols(1:pDim,sum( QuadMinFunPara(1:bestAlgoFlag) )  ); % save the best solution to be used in LLS

else % IotherPara(1)~=5  use the user defined box X to find the xRelaxedOpt and fxRelaxedOpt
    
   if IotherPara(1)==1,lb=-max(abs(trueb))*ones(pDim,1); ub=max(abs(trueb))*ones(pDim,1);            % - max(|b_i|) <= X_i <= max(|b_i|)  box X for quad.Min. same as true b
   elseif IotherPara(1)==2,lb=-2*max(abs(trueb))*ones(pDim,1); ub=2*max(abs(trueb))*ones(pDim,1);    % - 2*max(|b_i|) <= X_i <= 2*max(|b_i|)  box X for quad. Min. 2 times bigger than b 
   elseif IotherPara(1)==3,lb=-3*max(abs(trueb))*ones(pDim,1); ub=3*max(abs(trueb))*ones(pDim,1);    % - 3*max(|b_i|) <= X_i <= 3*max(|b_i|) box X is 3 times bigger than b
   elseif IotherPara(1)==4,lb=-(1/2)*max(abs(trueb))*ones(pDim,1); ub=(1/2)*max(abs(trueb))*ones(pDim,1);    % -(1/2)*max(|b_i|) <= X_i <= (1/2)*max(|b_i|) box X is half of the true b
   end

   [bestAlgoFlag,saveAllTheSols,cpuxRelaxedOpt,xRelaxedOpt,fxRelaxedOpt]=quadMinFun(pDim,Anorm,bnorm,cnorm,lb,ub,S,iParaRelaxX,rParaRelaxX,QuadMinFunPara);
   solsOfQuadMin=saveAllTheSols(1:pDim,sum( QuadMinFunPara(1:bestAlgoFlag) )  ); % save the best solution to be used in LLS   

end
fxRelaxedOpt=scaleQP*fxRelaxedOpt;

if ~isempty(fileid)
    fprintf(fileid,'bestAlgoFlag = %d \n',bestAlgoFlag);
    if pDim<=50
        fprintf(fileid,'Solutions are =\n');printArray(saveAllTheSols(1:pDim,:),'%1.4f',fileid);
        fprintf(fileid,'A=\n');printArray(A,'%1.4f',fileid);
        fprintf(fileid,'b=\n');printArray(b,'%1.4f',fileid);fprintf(fileid,'c=%g \n',c);
    else % if the dim. is large, print only first 100 entries of the first row of A, and only the first 100 entries of vector b
        fprintf(fileid,'Solutions are =\n');printArray(saveAllTheSols(1:50,:),'%1.4f',fileid);fprintf(fileid,'. . . \n');
        fprintf(fileid,'A(1,1:50)=\n');printArray(A(1,1:50),'%1.4f',fileid);fprintf(fileid,'. . . \n');
        fprintf(fileid,'b(1:50)=\n');printArray(b(1:50),'%1.4f',fileid);fprintf(fileid,'. . . \n');
    end
    fprintf(fileid,'scaleQP used to find xRelaxedOpt is %1.8f \n',scaleQP);
    fprintf(fileid,'eig A=\n');
    printArray(eigenvalues,'%1.4f',fileid);
    if ~isempty(xMatrix), fprintf(fileid,'Rank(xMatrix)=%d \n',rank(xMatrix)); end
    fprintf(fileid,'cond. no.= max eig / min eig =%1.8f\n',max(eigenvalues)/min(eigenvalues));
    fprintf(fileid,'1/cond. no.= min eig / max eig =%1.8f\n',min(eigenvalues)/max(eigenvalues));
    fprintf(fileid,'%s\n',datetime('now'));
    
    fclose(fileid);
end







end