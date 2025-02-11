function [rParaOut,bestVal,betaOut]=dfo1wrapper(p,n,X,y,k,Q,q,c,lb,ub,diagInv,XtX,normxRelaxedOpt,L,iPara,rPara,IotherPara,IstopCondPara)
% 23 June24, to use multiple initial points and choose the results corresponding the best solution
   
%   firstDigit=floor( IotherPara(21)/(10^(length(num2str(IotherPara(21)))-1)) );  % first digit of IotherPara(21) will tell use which box to use for refinement.
%   iPara(2)=IotherPara(21)-firstDigit*(10^(length(num2str(IotherPara(21)))-1)); % remaining digits will tell us the no. of iterations to use for refinement
   iPara(2)=IstopCondPara(5); % 10 Oct24
   idx1=min(IotherPara(24),p-k+1); % to pick a new support, starting from idx1Supp 
   [~,initialSupp]=sort(abs(normxRelaxedOpt),'descend');
   cpuStart=cputime;
   pbeta=normxRelaxedOpt;  
   if idx1>1
      pbeta( initialSupp(idx1-1) )=0;  
   end
   [stopflag,bestVal,betaOut]=dfo1(p,n,X,y,k,Q,q,c,lb,ub,diagInv,XtX,pbeta,L,iPara,rPara,IotherPara);
   cpuEnd=cputime;
   rParaOut.stopflag=stopflag; rParaOut.cpusec=cpuEnd-cpuStart;
   %fprintf('DFO cpustart=%1.10f, cpuend=%1.10f \n',cpuStart,cpuEnd);
   
end