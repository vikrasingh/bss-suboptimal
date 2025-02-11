function [out]=data_consistency_check(numOfSets,chooseParaToRun,IotherPara,IstopCondPara,QuadMinFunPara,seedToRunNInstances,seedToRunSetOfEgs,pDim,infoCol,nDcol,nPts,alphaValues,lambdaValues,tmaxValues,snr,...
                              chooseEgsToRun,numOfEgs,giveTruebManually,giveValuesOf_b)
out=1;  

for set=1:numOfSets
    if IotherPara(chooseParaToRun(set,1),15)==51 || IotherPara(chooseParaToRun(set,1),15)==52  % if running portfolio selection algorithm
       return; 
    end
end


%1. chooseEgsToRun should have same length as no. of rows of giveValuesOf_b
if length(chooseEgsToRun)~=numOfEgs 
   error('numOfEgs is not equal to the chooseEgsToRun.'); 
end

if giveTruebManually==1
   [~,numOfcols]=size(giveValuesOf_b(:,chooseEgsToRun)); 
   if numOfcols~=length(chooseEgsToRun)
      error('values defined in chooseEgsToRun para are not equal to the numOfEgs parameter.') 
   end 
   if numOfcols~=numOfEgs,error('numOfEgs is not equal to the chooseParaToRun length'); end
else
   [numOfrows,numOfcols]=size(giveValuesOf_b(:,1:numOfEgs)); 
   if numOfcols~=numOfEgs || numOfrows~=pDim
      error('randomly generated giveValuesOf_b array is not consistent with pDim or number of examples.'); 
   end
end

if length(chooseEgsToRun)>numOfcols
   error('giveValuesOf_b array is not consistent with the values given in chooseEgsToRun. ') 
end


%2. nDcol should be 0 for IotherPara(15)=0,1,2,3 and 5
for set=1:numOfSets
    if IotherPara(chooseParaToRun(set,1),15)~=4 || IotherPara(chooseParaToRun(set,1),15)~=6
       if nDcol~=0
          disp('WARNING: Algo. 1,2,3,4 or 5 has been called with nonzero nDcol.');
          userInput=input('To continue, type 1 and Enter, to stop the computation and redefine the para. type 2 and Enter \n');
          if userInput==1
             continue
          elseif userInput==2
             out=2;return 
          end
       end
    end
end

%3.  To deal with dependent columns in the xMatrix, we can only use DIMENSION REDUCTION method.
for set=1:numOfSets
    if IotherPara(chooseParaToRun(set,1),15)==4 || IotherPara(chooseParaToRun(set,1),15)==6
       if nDcol~=0
          error('Dimension Reduction algorithm IotherPara(15)=4 or 6 has been called with nDcol=0 '); 
          
       end
    end
end

%4. seedForTheEg should be 0 if testing multiple examples.
if ~isequal(length(chooseEgsToRun),length(seedToRunSetOfEgs)) && (seedToRunSetOfEgs~=0)
   error('seedForTheEg should be 0 if want to test multiple examples.') 
end


% if algo. 3 is running with tmax value not 0, it is good to use IotherPara(14)=1 for variable selection
% for tmax not 0 and using normalization approach 2 there is a chance that the initial box X is infeasible for the given value of tmax.
for set=1:numOfSets
    if IotherPara(chooseParaToRun(set,1),15)==3  
       if IotherPara(chooseParaToRun(set,1),14)==2 
          for i=1:length(tmaxValues)
              if tmaxValues(i)> (pDim)
                 error('tmax value should be smaller than pDim value.'); 
              end
          end
       elseif IotherPara(chooseParaToRun(set,1),14)==1
          fprintf('WARNING: Use IotherPara(14)=2 for the inequality constraint algo. 3 \n');
          userInput=input('Type 1 and Enter, to continue. Type 2 and enter, to stop the computation and redefine the para.  \n');
          if userInput==1
             continue
          elseif userInput==2
             out=2;return 
          end
          for i=1:length(tmaxValues)
              if tmaxValues(i)> (pDim) ||  tmaxValues(i)> nPts
                 error('tmax value should be smaller than pDim and nPts value.'); 
              end
          end
       end
    
    end
end


%6. if using semi-interval 2 approach for penalized regss. i.e. IotherPara(15)=2  then the only choices for
%evaluating inclusion function i.e. IotherPara(19)=0,1,2,3  that is no option that find F(V) by deleting and
%adding contribution of some partitioned sub-box is not valid.
for set=1:numOfSets
    if IotherPara(chooseParaToRun(set,1),15)==2 
       if IotherPara(chooseParaToRun(set,1),19)>3
          error('Inclusion function eval. (IotherPara(19) can only be equal to 0,1,2,3) when running semi-interval 2 (IotherPara(15)=2) approach');
       end
    end
end

%7. if using semi-interval 2 approach, then for no. of cuts parameter IotherPara(9) should be equal to 0 or 1m where m cuts (from longest to shortest side of the box) at the midpoint of the box X,
% cut at the feasible point cannot be used as we find feasible point later on.
for set=1:numOfSets
    if IotherPara(chooseParaToRun(set,1),9)==1 && IotherPara(chooseParaToRun(set,1),15)~=77 && IotherPara(chooseParaToRun(set,1),15)~=79
          error('IotherPara(9)=1 is no longer in use.');
    end 
    
    if IotherPara(chooseParaToRun(set,1),15)==2
       firstDigit=floor(IotherPara(chooseParaToRun(set,1),9)/10);numOfCuts=IotherPara(9)-(floor(IotherPara(9)/10))*10; 
       if firstDigit>1 && numOfCuts>1 
          error(' IotherPara(9) can only be equal to 0 or 1m when IotherPara(15)=2');
       end
    end
end



%9. If using algo. 1 or algo. 3 only then keep the dummy array of either  tmaxValues or lambdaValues to have
%entries which are not same (eg. avoid keep array of zeros or ones)
if length(lambdaValues)>1 
   if isequal(lambdaValues,lambdaValues(1)*ones(1,length(lambdaValues)))
      if (length(alphaValues))>1 
         error('Keep the lambdaValues array with all entries as different (like 1:k )');
      end
   end
end
if length(tmaxValues)>1
   if isequal(tmaxValues,tmaxValues(1)*ones(1,length(tmaxValues)))
      if (length(tmaxValues))>1 
         error('Keep the tmaxValues array with all entries as different (like 1:k )'); 
      end
   end
end

%10. IotherPara(18) the plotting control parameter should be same for the given no. of numOfSets for a run
plotPara=IotherPara(chooseParaToRun(1,1),18);
for set=2:numOfSets    
    if IotherPara(chooseParaToRun(set,1),18)~=plotPara
          error('IotherPara(18) should have the same value for all the given numOfSets.');
    end 
end

   
%11. if the IotherPara(9)=-1 or -2 or -3  but versionflag~=7
for set=1:numOfSets
    if IotherPara(chooseParaToRun(set,1),9)==-1 || IotherPara(chooseParaToRun(set,1),9)==-2 || IotherPara(chooseParaToRun(set,1),9)==-3
       if ~ismember(IotherPara(chooseParaToRun(set,1),15),[7 8 81 82 83 84 9 91 92 93 71 72 73 74 75 76 77 78 79 20 21 22 23 24 25 30 31 32 33 34 35 36]) 
          error('IotherPara(9) should be positive for if IotherPara(15) is not 7');
       end
    end 
    if ismember(IotherPara(chooseParaToRun(set,1),5),[5 6 7])  % if the selection option is 5,6,7
       if ~ismember( IotherPara(chooseParaToRun(set,1),15),[7 8 81 82 83 84 9 91 92  93 71 72 73 74 75 76 77 78 79 20 21 22 23 24 25  30 31 32 33 34 35 36] ) 
          error('IotherPara(5) should not be 5,6 or 7, if IotherPara(15) is not 7');
       end
    end
end

%12. check that there is no repetition in the userdefined lambdaValues and tmaxValues array
for set=1:numOfSets
    if IotherPara(chooseParaToRun(set,1),15)==1 || IotherPara(chooseParaToRun(set,1),15)==2
       if (length(lambdaValues))~=length(unique(lambdaValues))
          error('Repeated values in the lambdaValues array.'); 
       end
    else
       if (length(tmaxValues))~=length(unique(tmaxValues))
          error('Repeated values in the tmaxValues array.'); 
       end 
    end
end

%{
% 13. 18Feb22 if v71 has been called with unconstrained refinement or lb f options 
for set=1:numOfSets
   if ismember(IotherPara(chooseParaToRun(set,1),15),[7 8 71 74 75 76]) 
      if ismember(IotherPara(chooseParaToRun(set,1),23),[13:20])
         error('Boxed v7,8 or 71 algo. with unconstr. quad. min. for refinement option. ');
      end
      if ismember( find(QuadMinFunPara(chooseParaToRun(set,3),:)==1,1) ,[13:20])
         error('Boxed v7,8 or 71 algo. with unconstr. quad. min. for lb f option. ');
      end
   end 
end

%14. 18Feb22 if v72 has been called with box constr. refinement or lb f options. 
for set=1:numOfSets
   if ismember(IotherPara(chooseParaToRun(set,1),15),[72 73 77 78])
      if ismember(IotherPara(chooseParaToRun(set,1),23),[1:11]) && IotherPara(chooseParaToRun(set,1),19)==8
         error('V72 /73 /77 /78: Unconstr. algo. with box constr. quad. min. for refinement option. ');
      end
      if ismember( find(QuadMinFunPara(chooseParaToRun(set,3),:)==1,1) ,[1:11] )
         error('V72 /73 /77: Unconstr. algo. with box constr. quad. min. for lb f option. ');
      end
   end 
end
%}

%15. 14Mar22 if delay the F call IotherPara(7)=1 has been used for AG7
for set=1:numOfSets
   if ismember(IotherPara(chooseParaToRun(set,1),7),[1 2]) && IotherPara(chooseParaToRun(set,1),15)==7
      error('AG7 with delay F call IotherPara(7)= 1 or 2 has been called. ');
   end 
end

%16. 16Mar22 if using AG8 stepTm i.e. IotherPara(17) should be a positive integer <pDim
for set=1:numOfSets
   if ~IotherPara(chooseParaToRun(set,1),15)==8  % i.e. AG8 
      if ismember(IotherPara(chooseParaToRun(set,1),17),(1:(pDim-1)) ) 
         error('Error in IotherPara(17) for AG8. ');
      end
   end 
end

%17. 16Mar22 IotherPara(5) 
for set=1:numOfSets
   if ismember(IotherPara(chooseParaToRun(set,1),15),[7 8 81 82 83 84 71 72 73 74 75 76 77 78 79])  %  
      if ismember(IotherPara(chooseParaToRun(set,1),5),[2 3] ) 
         error('Error in IotherPara(5) for the AG used. ');
      end
   end 
end

%18. 25Mar22,if IotherPara(1)=0 then QuadMinFunpara(1,:) should have atleast one unconstrained algo. selected
for set=1:numOfSets
   activeAlgo=(QuadMinFunPara(1,:)==1); 
   if IotherPara(chooseParaToRun(set,1),1)==0  % flag to use the whole space while finding xRelaxedOpt 
      if sum(activeAlgo(12:20))==0
         error('IotherPara(1)=0 option with only box constrained algo. to find xRelaxedOpt has been used. ');
      end
   end 
   if sum(activeAlgo(21:end))>0   % right now only 20 active algo.
      error('QuadMinFunPara(i,:) algo. > 20 is falsely active. '); 
   end
end

%19. 28April22, if using multiple cuts IotherPara(9)= _m , make sure that IotherPara(19)=8 and alpha=0 with AG7 
for set=1:numOfSets 
   if length(num2str(abs( IotherPara(chooseParaToRun(set,1),9) )))>9 && ismember(0,alphaValues) &&  ismember(IotherPara(chooseParaToRun(set,1),15), [7 71 72 73 74 75 76 77 78 79] ) % flag to use the whole space while finding xRelaxedOpt 
      if IotherPara(chooseParaToRun(set,1),19)~=8
         error('Multiple cuts option IotherPara(9) has been called with IotherPara(19)/=8 .');
      end
   end 
end

%20. 13May22, if using QRD appraoch to find lb f, i.e. if IotherPara(19)= 9 or 10, then p<=n
for set=1:numOfSets
   if pDim>nPts && (IotherPara(chooseParaToRun(set,1),19)==9 || IotherPara(chooseParaToRun(set,1),19)==10)
      error('QRD approach, IotherPara(19)=9 or 10 only works when pDim<=nPts') 
   end
end

% %21. 22May22, infoCol==pDim always
% for set=1:numOfSets
%    if infoCol~=pDim 
%       error('infoCol /= pDim.'); 
%    end
% end



end