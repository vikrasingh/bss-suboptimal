function  [nDcol,idxnDcol,D_t]=find_nDcol_inData(p,n,xMatrix)
% 20 March21, subroutine to find which and how many columns of xMatrix are dependent and save the dependent columns 
% in the matrix D_t to be used as dependent block later on. 

  nIDcol=0;        % for implementation, to save the no. of independent columns in the xMatrix
  idxnDcol=ones(1,p);  % to store the indices of the dependent columns in the xMatrix
  E=rref(xMatrix); % find the row reduced echelon form of xMatrix
  for i=1:p
      unitVec=zeros(n,1);
      unitVec(nIDcol+1)=1;
      if isequal(E(:,i),unitVec)  % if ith column of echelon form E is not a unit vector
         nIDcol=nIDcol+1;
         idxnDcol(i)=0;
      end 
      if nIDcol==p || nIDcol==n
         break; 
      end
  end
  D_t=E(:,idxnDcol==1) ;  % define the D_t matrix to be dependent columns in the xMatrix 
  nDcol=sum(idxnDcol);


end