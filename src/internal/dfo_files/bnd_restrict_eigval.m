function [Mu,gammak]=bnd_restrict_eigval(pDim,xMatrix,yArray,tmax)
% updated 25 Oct22
% Procedure to find bound Mu where ||beta||_infty <= Mu , using the method given in Sec. 2.3.1 of Bertsimas et.al. 2016, BSS via Modern optimization lens

% reorder the columns of X matrix
absCorrVal=zeros(1,pDim);
for idx=1:pDim
    absCorrVal(idx)=abs(xMatrix(:,idx)'*yArray);
end
[sortedabsCorrVal,sortidx]=sort(absCorrVal,'descend');
xMatrix=xMatrix(:,sortidx);


mu=-inf;
for i=1:pDim
   for j=(i+1):pDim
      corrij=abs( xMatrix(:,i)'*xMatrix(:,j) );
      if corrij > mu
         mu=corrij; 
      end
   end
end
% fprintf('coherence mu=%1.6f \n',mu);
gammak=1-mu*(tmax-1);
% fprintf('gammak=%1.6f \n',gammak);

if gammak>0
   betainf=min( sqrt(  sum( sortedabsCorrVal(1:tmax).^2 )  )/gammak , norm(yArray)/sqrt(gammak) );
   Mu=betainf;
else
   Mu=NaN; % a dummy flag
end

% fprintf('bound for the box Mu= %1.6f \n',Mu);

end