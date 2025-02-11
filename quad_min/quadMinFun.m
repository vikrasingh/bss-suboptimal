%{
07 Oct20, function that consists of different versions to solve quadratic minimization problem.
%}
function [bestAlgoFlag,saveAllTheSols,cpuxRelaxedOpt,xRelaxedOpt,fxRelaxedOpt]=quadMinFun(n,A,b,c,lb,ub,S,iPara,rPara,QuadMinFunPara)
% till 08/18/23 [bestAlgoFlag,saveAllTheSols,cpuxRelaxedOpt,xRelaxedOpt,fxRelaxedOpt]=quadMinFun(n,A,b,c,lb,ub,iPara,rPara,QuadMinFunPara)

dummyMatrix=zeros(n+2,20);  % Initialize the matrix 
bestAlgoFlag=0;cpuFlag=inf;
fxRelaxedOpt=inf;   % just for implementation purpose
xRelaxedOpt=[]; % initialization
cpuxRelaxedOpt=0;
isSoftStop=1; % use soft stop for the quadratic min. algo.
isTF=0;targetfbest=[];

if QuadMinFunPara(1)==1  % use minq8
   if isempty(xRelaxedOpt)
      x0=0.5*(lb+ub);
   else,x0=xRelaxedOpt;
   end  
   tstart=tic;
  [~,xminq8,fminq8]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,1,isSoftStop,isTF,targetfbest);  % which_algo=1 
  tend=toc(tstart);
  dummyMatrix(:,1)=[xminq8;fminq8;tend];
  
  if fminq8<fxRelaxedOpt
      xRelaxedOpt=xminq8;fxRelaxedOpt=fminq8;cpuxRelaxedOpt=tend;
      bestAlgoFlag=1;cpuFlag=tend;
  end
end
if QuadMinFunPara(2)==1  % use SD2
   if isempty(xRelaxedOpt)
      x0=(lb+ub)/2;
   else,x0=xRelaxedOpt;
   end
   tstart=tic;
   [~,xPsdSd,fPsdSd]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,2,isSoftStop,isTF,targetfbest);  % which_algo=2 
   tend=toc(tstart);
   dummyMatrix(:,2)=[xPsdSd;fPsdSd;tend];
  if (fPsdSd + eps )<fxRelaxedOpt
     xRelaxedOpt=xPsdSd;fxRelaxedOpt=fPsdSd;cpuxRelaxedOpt=tend; 
     bestAlgoFlag=2;cpuFlag=tend;
  elseif abs(fPsdSd-fxRelaxedOpt)< (10)^(-5)
      if tend<cpuFlag
         bestAlgoFlag=2;cpuFlag=tend;
      end
  end
  
end
if QuadMinFunPara(3)==1 % use LLS
   tstart=tic;
   [~,vLls,fvLls]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,3,isSoftStop,isTF,targetfbest);  % which_algo=3  
   tend=toc(tstart);
   dummyMatrix(:,3)=[vLls;fvLls;tend];
   if (fvLls + eps )<fxRelaxedOpt
       xRelaxedOpt=vLls;fxRelaxedOpt=fvLls;cpuxRelaxedOpt=tend;
       bestAlgoFlag=3;cpuFlag=tend;
   elseif abs(fvLls-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=3;cpuFlag=tend;
       end
   end
   
end
if QuadMinFunPara(4)==1 % use Hbd
   if isempty(xRelaxedOpt)
      x0=(lb+ub)/2;
   else,x0=xRelaxedOpt;
   end
   tstart=tic;
   [~,vHbd,fvHbd]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,4,isSoftStop,isTF,targetfbest);  % which_algo=4 
   tend=toc(tstart);
   dummyMatrix(:,4)=[vHbd;fvHbd;tend];
   if (fvHbd + eps )<fxRelaxedOpt
      xRelaxedOpt=vHbd;fxRelaxedOpt=fvHbd;cpuxRelaxedOpt=tend; 
      bestAlgoFlag=4;cpuFlag=tend;
   elseif abs(fvHbd-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=4;cpuFlag=tend;
       end  
   end
end
if QuadMinFunPara(5)==1 % use SequMiniQuadFun_v11  min function
    tstart=tic;
    [~,xv11,fv11]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,5,isSoftStop,isTF,targetfbest);  % which_algo=5 
    tend=toc(tstart);
   dummyMatrix(:,5)=[xv11;fv11;tend];
   if (fv11 + eps )<fxRelaxedOpt
     xRelaxedOpt=xv11;fxRelaxedOpt=fv11;cpuxRelaxedOpt=tend; 
     bestAlgoFlag=5;cpuFlag=tend;
   elseif abs(fv11-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=5;cpuFlag=tend;
       end
   end
end
if QuadMinFunPara(6)==1 % use convexQuadByParTanCorDesBox  min function 
   if isempty(xRelaxedOpt)
      x0=0.5*(lb+ub);
   else,x0=xRelaxedOpt;
   end 
   tstart=tic;
   [~,xPtCd,fPtCd]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,6,isSoftStop,isTF,targetfbest);  % which_algo=6  
   tend=toc(tstart);
   dummyMatrix(:,6)=[xPtCd;fPtCd;tend];
   if (fPtCd + eps )<fxRelaxedOpt
      xRelaxedOpt=xPtCd;fxRelaxedOpt=fPtCd;cpuxRelaxedOpt=tend; 
      bestAlgoFlag=6;cpuFlag=tend;
   elseif abs(fPtCd-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=6;cpuFlag=tend;
       end
   end
end
if QuadMinFunPara(7)==1 % use matlab quadprog function
  tstart=tic;
  [~,xquadprog,fxquadprog]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,7,isSoftStop,isTF,targetfbest);  % which_algo=7  
  tend=toc(tstart);
  dummyMatrix(:,7)=[xquadprog;fxquadprog;tend];
  fxquadprog=fxquadprog+c;
   if (fxquadprog + eps )<fxRelaxedOpt
     xRelaxedOpt=xquadprog;fxRelaxedOpt=fxquadprog;cpuxRelaxedOpt=tend; 
     bestAlgoFlag=7;cpuFlag=tend;
   elseif abs(fxquadprog-fxRelaxedOpt)< (10)^(-5)
      if tend<cpuFlag
          bestAlgoFlag=7;cpuFlag=tend;
      end
   end
end
if QuadMinFunPara(8)==1 % use SD  
   if isempty(xRelaxedOpt)
      x0=0.5*(lb+ub);
   else,x0=xRelaxedOpt;
   end 
   tstart=tic;
   [~,xSD,fxSD]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,8,isSoftStop,isTF,targetfbest);  % which_algo=8 
   tend=toc(tstart);
   dummyMatrix(:,8)=[xSD;fxSD;tend];
   if (fxSD + eps )<fxRelaxedOpt
      xRelaxedOpt=xSD;fxRelaxedOpt=fxSD;cpuxRelaxedOpt=tend; 
      bestAlgoFlag=8;cpuFlag=tend;
   elseif abs(fxSD-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=8;cpuFlag=tend;
       end
   end
end
if QuadMinFunPara(9)==1 % use Coordinate Descent  
   if isempty(xRelaxedOpt)
      x0=0.5*(lb+ub);
   else,x0=xRelaxedOpt;
   end 
   tstart=tic;
   [~,xCD,fxCD]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,9,isSoftStop,isTF,targetfbest);  % which_algo=9 
   tend=toc(tstart);
   dummyMatrix(:,9)=[xCD;fxCD;tend];
   if (fxCD + eps )<fxRelaxedOpt
      xRelaxedOpt=xCD;fxRelaxedOpt=fxCD;cpuxRelaxedOpt=tend; 
      bestAlgoFlag=9;cpuFlag=tend;
   elseif abs(fxCD-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=9;cpuFlag=tend;
       end
   end
end
if QuadMinFunPara(10)==1 % use Conjugate Gradient  
   if isempty(xRelaxedOpt)
      x0=0.5*(lb+ub);
   else,x0=xRelaxedOpt;
   end 
   tstart=tic;
   [~,xCG,fxCG]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,10,isSoftStop,isTF,targetfbest);  % which_algo=10 
   tend=toc(tstart);
   dummyMatrix(:,10)=[xCG;fxCG;tend];
   if (fxCG + eps )<fxRelaxedOpt
      xRelaxedOpt=xCG;fxRelaxedOpt=fxCG;cpuxRelaxedOpt=tend; 
      bestAlgoFlag=10;cpuFlag=tend;
   elseif abs(fxCG-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=10;cpuFlag=tend;
       end
   end
end
if QuadMinFunPara(11)==1 % use Parallel Tangent  
   if isempty(xRelaxedOpt)
      x0=0.5*(lb+ub);
   else,x0=xRelaxedOpt;
   end 
   tstart=tic;
   [~,xPT,fxPT]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,11,isSoftStop,isTF,targetfbest);  % which_algo=11 
    tend=toc(tstart);
   dummyMatrix(:,11)=[xPT;fxPT;tend];
   if (fxPT + eps )<fxRelaxedOpt
      xRelaxedOpt=xPT;fxRelaxedOpt=fxPT;cpuxRelaxedOpt=tend; 
      bestAlgoFlag=11;cpuFlag=tend;
   elseif abs(fxPT-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=11;cpuFlag=tend;
       end
   end
end
if QuadMinFunPara(12)==1 % use convexQuadByConGraCorDesUnbox for min.   
   if isempty(xRelaxedOpt)
      x0=zeros(n,1);
   else,x0=xRelaxedOpt;
   end 
   tstart=tic;
   [~,xCgCd,fCgCd]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,12,isSoftStop,isTF,targetfbest);  % which_algo=12 
   tend=toc(tstart);
   dummyMatrix(:,12)=[xCgCd;fCgCd;tend];
   if (fCgCd + eps )<fxRelaxedOpt
      xRelaxedOpt=xCgCd;fxRelaxedOpt=fCgCd;cpuxRelaxedOpt=tend; 
      bestAlgoFlag=12;cpuFlag=tend;
   elseif abs(fCgCd-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=12;cpuFlag=tend;
       end
   end
end
if QuadMinFunPara(13)==1 % use unconstrained SD   
   if isempty(xRelaxedOpt)
      x0=zeros(n,1);
   else,x0=xRelaxedOpt;
   end 
   tstart=tic;
   [~,xuSD,fxuSD]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,13,isSoftStop,isTF,targetfbest);  % which_algo=13 
   tend=toc(tstart);
   dummyMatrix(:,13)=[xuSD;fxuSD;tend];
   if (fxuSD + eps )<fxRelaxedOpt
      xRelaxedOpt=xuSD;fxRelaxedOpt=fxuSD;cpuxRelaxedOpt=tend; 
      bestAlgoFlag=13;cpuFlag=tend;
   elseif abs(fxuSD-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=13;cpuFlag=tend;
       end
   end
end
if QuadMinFunPara(14)==1 % use unconstrained PT   
   if isempty(xRelaxedOpt)
      x0=zeros(n,1);
   else,x0=xRelaxedOpt;
   end 
   tstart=tic;
   [~,xuPT,fxuPT]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,14,isSoftStop,isTF,targetfbest);  % which_algo=14 
   tend=toc(tstart);
   dummyMatrix(:,14)=[xuPT;fxuPT;tend];
   if (fxuPT + eps )<fxRelaxedOpt
      xRelaxedOpt=xuPT;fxRelaxedOpt=fxuPT;cpuxRelaxedOpt=tend; 
      bestAlgoFlag=14;cpuFlag=tend;
   elseif abs(fxuPT-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=14;cpuFlag=tend;
       end
   end
end
if QuadMinFunPara(15)==1 % use unconstr. adaptive Method 1   
   if isempty(xRelaxedOpt)
      x0=zeros(n,1);
   else,x0=xRelaxedOpt;
   end 
   tstart=tic;
   [~,xuA1,fxuA1]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,15,isSoftStop,isTF,targetfbest);  % which_algo=15 
   tend=toc(tstart);
   dummyMatrix(:,15)=[xuA1;fxuA1;tend];
   if (fxuA1 + eps )<fxRelaxedOpt
      xRelaxedOpt=xuA1;fxRelaxedOpt=fxuA1;cpuxRelaxedOpt=tend; 
      bestAlgoFlag=15;cpuFlag=tend;
   elseif abs(fxuA1-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=15;cpuFlag=tend;
       end
   end
end
if QuadMinFunPara(16)==1 % use unconstr. adaptive Method 2   
   if isempty(xRelaxedOpt)
      x0=zeros(n,1);
   else,x0=xRelaxedOpt;
   end 
   tstart=tic;
   [~,xuA2,fxuA2]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,16,isSoftStop,isTF,targetfbest);  % which_algo=16 
   tend=toc(tstart);
   dummyMatrix(:,16)=[xuA2;fxuA2;tend];
   if (fxuA2 + eps )<fxRelaxedOpt
      xRelaxedOpt=xuA2;fxRelaxedOpt=fxuA2;cpuxRelaxedOpt=tend; 
      bestAlgoFlag=16;cpuFlag=tend;
   elseif abs(fxuA2-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=16;cpuFlag=tend;
       end
   end
end
if QuadMinFunPara(17)==1 % use unconstr. adaptive Method 3   
   if isempty(xRelaxedOpt)
      x0=zeros(n,1);
   else,x0=xRelaxedOpt;
   end 
   tstart=tic;
   [~,xuA3,fxuA3]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,17,isSoftStop,isTF,targetfbest);  % which_algo=17 
   tend=toc(tstart);
   dummyMatrix(:,17)=[xuA3;fxuA3;tend];
   if (fxuA3 + eps )<fxRelaxedOpt
      xRelaxedOpt=xuA3;fxRelaxedOpt=fxuA3;cpuxRelaxedOpt=tend; 
      bestAlgoFlag=17;cpuFlag=tend;
   elseif abs(fxuA3-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=17;cpuFlag=tend;
       end
   end
end
if QuadMinFunPara(18)==1 % use unconstr. Coordinate Descent   
   if isempty(xRelaxedOpt)
      x0=zeros(n,1);
   else,x0=xRelaxedOpt;
   end 
   tstart=tic;
   [~,xuCD,fxuCD]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,18,isSoftStop,isTF,targetfbest);  % which_algo=18 
   tend=toc(tstart);
   dummyMatrix(:,18)=[xuCD;fxuCD;tend];
   if (fxuCD + eps )<fxRelaxedOpt
      xRelaxedOpt=xuCD;fxRelaxedOpt=fxuCD;cpuxRelaxedOpt=tend; 
      bestAlgoFlag=18;cpuFlag=tend;
   elseif abs(fxuCD-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=18;cpuFlag=tend;
       end
   end
end
if QuadMinFunPara(19)==1 % use unconstr. adaptive Method 4   
   if isempty(xRelaxedOpt)
      x0=zeros(n,1);
   else,x0=xRelaxedOpt;
   end 
   tstart=tic;
   [~,xuA4,fxuA4]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,19,isSoftStop,isTF,targetfbest);  % which_algo=19 
   tend=toc(tstart);
   dummyMatrix(:,19)=[xuA4;fxuA4;tend];
   if (fxuA4 + eps )<fxRelaxedOpt
      xRelaxedOpt=xuA4;fxRelaxedOpt=fxuA4;cpuxRelaxedOpt=tend; 
      bestAlgoFlag=19;cpuFlag=tend;
   elseif abs(fxuA4-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=19;cpuFlag=tend;
       end
   end
end
if QuadMinFunPara(20)==1 % use unconstr. conjugate gradient   
   if isempty(xRelaxedOpt)
      x0=zeros(n,1);
   else,x0=xRelaxedOpt;
   end 
   tstart=tic;
   [~,xuCG,fxuCG]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,20,isSoftStop,isTF,targetfbest);  % which_algo=20 
   tend=toc(tstart);
   dummyMatrix(:,20)=[xuCG;fxuCG;tend];
   if (fxuCG + eps )<fxRelaxedOpt
      xRelaxedOpt=xuCG;fxRelaxedOpt=fxuCG;cpuxRelaxedOpt=tend; 
      bestAlgoFlag=20;cpuFlag=tend;
   elseif abs(fxuCG-fxRelaxedOpt)< (10)^(-5)
       if tend<cpuFlag
          bestAlgoFlag=20;cpuFlag=tend;
       end
   end
end

saveAllTheSols=zeros(n+2,sum(QuadMinFunPara(1:20)));count=0;
for i=1:20
    if dummyMatrix(:,i)~=0
        count=count+1;
        saveAllTheSols(:,count)=dummyMatrix(:,i);
    end
end






end