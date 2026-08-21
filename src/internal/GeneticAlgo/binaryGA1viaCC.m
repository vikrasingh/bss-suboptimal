function [rParaOut,xbest,fbest] = binaryGA1viaCC(dim,x0,iPara,rPara,fun)
%2/3/2024 binaryGA1 via Cardinality Constraint
%feasible binary col vector x0(i)=0 or 1
%iPara, [1]=NP (say 50-100), [2]=max_Iter, tmax=iPara(3)
%rPara, [1]=CR=crossover rate (say 0.80), [2]=MR=mutation rate (say 0.15), 

%not used yet (iout--- [1]=nf (original), [2]=idExit (0 by suffCond(), 1 by nmax(soft), 2 by gg<eps, 3 by abs(alpha)<eps)
%         [3]=max_Iter
%not used yet (rout--- 
% rParaOut is structure with fields, numOfIter, necConMaxVioRefQM, stopflag used
%fun = @jFitnessFunction;  %MS:obj function (dim,x) of binary variables, 

NP=iPara(1); %population size
max_Iter=iPara(2);
tmax=iPara(3);
toDebug=iPara(4); % whether to print more info. in debug mode or not
CR=rPara(1);     %rPara(1); %crossover rate ~ 0.8, 
MR=rPara(2);    %rPara(2); %mutation rate ~ 0.15
stopflag=0;  % =0 converged, =6 max cpu limit, =5 max iter. limit  
rParaOut.numOfIter=0; rParaOut.stopflag=stopflag;

X   = zeros(dim,NP);   %MS:population

for k = 1:NP
  itp=0;  
  for i = 1:dim
    if rand() > 0.5 && itp < tmax
      X(i,k) = 1;  %x(i)=0 or 1
      itp=itp+1;
    end
  end
  while itp<tmax % if X(:,k) does not have exactly tmax 1's
     for i=1:dim
         if X(i,k)==0
            X(i,k)=1;
            itp=itp+1;
            break;
         end
     end
  end 
end
X(:,1) = x0;

fit  = zeros(1,NP);  %MS:row, fitness of population
%fbest = inf;        %MS:fbest w/ corresponding Xgb
for k = 1:NP
  fit(k) = fun(dim,X(:,k));
  %if fit(k) < fbest
  %  fbest = fit(k); 
  %  xbest  = X(:,k);
  %end
end
[~, idx] = sort(fit,'ascend');
fbest=fit(idx(1)); xbest=X(:,idx(1));
fworstG=fit(idx(NP));   %worst of current generation

it = 1;        %MS:iter count
%---Generations start------------------------------------------------
while it <= max_Iter
  %MS:set probability distribution of population based on fitness:  
  Ifit = fworstG-fit; 
  if sum(abs(Ifit))<1e-5
     stopflag=0; rParaOut.numOfIter=it; rParaOut.stopflag=stopflag;  
     return 
  end
  % include early exit if Ifit becomes a zero vector
  prob = Ifit / (1e-30+sum(Ifit)); %%**assume: the better/smaller of fitness, the larger of probability  
  
  %MS:generate offsprings via CO:
  k1  = selectParent(prob); %MS:parent 1
  k2  = selectParent(prob); %MS:parent 2
  X1=zeros(dim,NP);X2=X1;
  X1(:,1)  = X(:,k1);      %zeros(dim,1);  %MS:initialization of 2 base offsprings, may be expanded up to 2*NP offsprings
  X2(:,1)  = X(:,k2);      %zeros(dim,1);
  z   = 1;

  for k = 1:NP
    if rand() < CR
      k1  = selectParent(prob); %MS:parent 1
      k2  = selectParent(prob); %MS:parent 2
      ctr=0;
      while k1==k2 && ctr<20
         k2=selectParent(prob);
         ctr=ctr+1;
      end
      if k1~=k2  % only then do the crossover
          P1  = X(:,k1);
          P2  = X(:,k2);
          ind = randi([1, dim - 1]); %MS: CO site
          X1(:,z) = [P1(1 : ind); P2(ind + 1 : dim)];  %n1+(tmax-n2)
          X2(:,z) = [P2(1 : ind); P1(ind + 1 : dim)];  %n2+(tmax-n1)
          n1=sum(P1(1 : ind));
          n2=sum(P2(1 : ind));
          if n1<n2
             %increase P1(1 : ind) towards P2(1 : ind) by n2-n1
             ii=0;
             for j=1:ind
                 if ii<n2-n1 && P2(j)==1 && P1(j)==0 
                     X1(j,z)=1;
                     X2(j,z)=0;
                     ii=ii+1;
                 end
             end
          elseif n1>n2
             %increase P2(1 : ind) towards P1(1 : ind) by n1-n2
             ii=0;
             for j=1:ind
                 if ii<n1-n2 && P1(j)==1 && P2(j)==0 
                     X2(j,z)=1;
                     X1(j,z)=0;
                     ii=ii+1;
                 end
             end
          end
          z = z + 1;
      end
      
    end

  end

  X1=X1(:,1:z-1);X2=X2(:,1:z-1);
  
  %MS:modify offsprings via MT:
  Xnew = [X1, X2];
  Noff   = size(Xnew,2);  %MS:total # of offsprings
  Fnew = zeros(1,Noff);
  for k = 1:Noff
    for i = 1:2:(dim-1)
      if rand() <= MR && Xnew(i,k)~=Xnew(i+1,k)
        Xnew(i,k) = 1 - Xnew(i,k);    %MS: mutated to the opposite value
        Xnew(i+1,k) = 1 - Xnew(i+1,k); 
      end
    end
  end 
  
  %MS:evaluate all new offsprings:  ok
  for k = 1:Noff
    Fnew(k) = fun(dim,Xnew(:,k));
    if Fnew(k) < fbest
      xbest  = Xnew(:,k); 
      fbest = Fnew(k);
    end
  end 
  
  %MS:use NP best of current population and new offsprings as new generation: ok
  XX  = [X, Xnew]; 
  FF  = [fit, Fnew]; 
  [FF, idx] = sort(FF,'ascend');
  X   = XX(:,idx(1:NP)); 
  fit = FF(1:NP);
  fworstG=fit(NP);
  
  %MS: print 
%   if toDebug>=2
%      fprintf('Iteration %d; fbest= %1.6f;fworst=%1.6f; card=%d \n',it,fbest,fworstG,sum(xbest));
%   end
%   fprintf('Iteration %d; fbest= %f; card=%d \n',it,fbest,sum(xbest));
  
  it = it + 1;
  if it>max_Iter
     stopflag=5; rParaOut.numOfIter=it; rParaOut.stopflag=stopflag;
     break; % stop the while loop and exit
  end
%   sum(Xnew,1)
end

%MS:set output:

end

%************************************************** ok
function Index = selectParent(prob)
%MS: randomly select an index based on the reversed probability "prob"
%%**assume: the better/smaller of fitness, the larger of probability
    Index=1;
    C = cumsum(prob); 
    P = rand();   %MS: r.n. in (0, 1)
    for i = 1:length(C)
        if C(i) >= P
           Index = i;
           break;
        end
    end

end
