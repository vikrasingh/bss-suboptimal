%{
26 Jan21, lasso approximation, pick those solutions of lasso which are in the range provided by lbOfError and
ubOfError.
%}
function [ridgeInfo]=getTrainErrorRidge(yArray,xMatrix,tol)

k=0.01:0.01:0.8;
getTime=tic;
[B]=ridge(yArray,xMatrix,k); % lasso call 
cpuRidge=toc(getTime);
n=length(k);
%
numOfZeros=sum(~(B(:,1)));numOfRidgeSol=0;whichSol=zeros(1,n);numOf0sInRidgeSol=10^(12)*ones(1,n); % initialization
for i=1:n
    if sum(~(B(:,i))) == numOfZeros
       numOfZeros=sum(~(B(:,i)));
       numOfRidgeSol=numOfRidgeSol+1;
       whichSol(numOfRidgeSol)=i;
       numOf0sInRidgeSol(numOfRidgeSol)=numOfZeros;
       numOfZeros=numOfZeros+1;
    end
end
id=logical(whichSol);whichSol=whichSol(id);
numOf0sInRidgeSol=numOf0sInRidgeSol(id);

ridgePara=k(whichSol);
ridgeSol=B(:,whichSol);
ErrorMat=xMatrix*ridgeSol-yArray;
trainErrorRidge=10^(12)*ones(1,numOfRidgeSol);

for i=1:numOfRidgeSol
   trainErrorRidge(i)=(norm(ErrorMat(:,i)))^2; 
end


%}
% save all the ouput solutions to a structure variable ridgeInfo
ridgeInfo.sol=ridgeSol;
ridgeInfo.para=ridgePara;
ridgeInfo.te=trainErrorRidge;
ridgeInfo.numOfSol=numOfRidgeSol;
ridgeInfo.cpu=cpuRidge;
ridgeInfo.numOf0s=numOf0sInRidgeSol;













end