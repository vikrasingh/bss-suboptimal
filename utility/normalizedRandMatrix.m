%{
22 May20, to generate xMatrix with desired features
%}
function [X]=normalizedRandMatrix(nPts,pDim,infoCol,mu,sigma,seed)

rng(seed,'twister');
%   X=randn(pDim,nPts);X=X';  % generate a random matrix
 

X=mvnrnd(mu,sigma,nPts);  % updated 16 Feb21

if infoCol>0
   if pDim~=infoCol 
      rightBlock=mvnrnd(zeros(1,pDim-infoCol),eye(pDim-infoCol),nPts); 
      X=[X rightBlock];
   end  
end

% % to standardize the data, check p 217 of An Intro to statistical learning book.
% %
% for j=1:pDim
%     X_jbar=mean(X(:,j));  % take the mean of each predictor/column 
%     std_dev=sqrt( (sum( ( X(:,j)-X_jbar ).^2 ))/nPts );  % standard deviation of each predictor.
%     for i=1:nPts
%         X(i,j)=X(i,j)/std_dev;
%     end
% end

% % 28 July21, standardize matrix X to have mean 0 and unit length(using 2-norm) for each column 
% for i=1:pDim
%     col=X(:,i);
%     stdcol=col - mean(col);  % make mean of each column = 0
%     X(:,i)=stdcol./norm(stdcol);  % make unit length = 1   
% end
X=normalize(X,'center'); % 25Oct22, make mean of columns 0
X=normalize(X,'norm');  % make 2-norm of columns 1

%}


% make some columns with large integer entries
%{
blocked 13th july20

colIndex1=randi([1 nPts],1,floor(nPts/3));  % pick almost a 3rd of columns to change to integer values
X(:,colIndex1)=10*randi([-200 200],pDim,floor(nPts/3)); % assign integer entries to those columns, with integer bw [-200 200], change accordingly


% make some columns of the matrix having 0 entries with a few large integers
colIndex2=setdiff([1:nPts],colIndex1);  % pick the remaining columns
getSize=length(colIndex2);
colIndex2=colIndex2(randperm(getSize, ceil(nPts/4)));   % randomly pick a 4th of columns from the remaining cols ,other than picked earlier in colIndex1
getColumns=zeros(pDim,length(colIndex2));
for i=1:length(colIndex2)
    pickEntries=randi([-500 500],ceil(pDim/2),1);
    pickIndexForThoseEntries=randi([1 pDim],ceil(pDim/2),1);
    col=zeros(pDim,1);col(pickIndexForThoseEntries)=pickEntries;
    getColumns(:,i)=col;
end
X(:,colIndex2)=getColumns;  % assign those cols


% make some entires of the matrix close to zero.
howManyEntries=floor((nPts*pDim)/5);
vec1=randi([1 pDim],1,howManyEntries);vec2=randi([1 nPts],1,howManyEntries);
if pDim==1 || nPts==1
    eps1=zeros(floor(howManyEntries/2),1);
else
    eps1=normrnd(0,0.05,[floor(howManyEntries/2) 1]);  % some really close to zero
end
eps2=normrnd(0,1,[ceil(howManyEntries/2) 1]);  % 
eps=vertcat(eps1,eps2);   % this will have howManyEntries no. of value close to zeros, lets assign these value to matrix X now
                                % where indices will be given by vec1 and vec2
for i=1:howManyEntries
   X(vec1(i),vec2(i))=eps(i); 
end
X=X';
%}


end