% 18 March21, to compare the num of 0s in the solution Xstar with given tolerances

function  [numOfSol0s]=get0sInSol(Xstar,tol)

   numOfSol0s=(-1)*ones(6,1); % Actual zeros of trueb in Xstar.
   tolVec=[tol  2*tol  3*tol  5*tol  10*tol  100*tol];
   for i=1:6      
%       [idx]=find(abs(Xstar) < tolVec(i));
%       numOfSol0s(i)=sum(logical(idx));
      numOfSol0s(i)=sum(abs(Xstar)<tolVec(i));
   end 







end