function [X]=Irandom(r,m,n)
%{
Input: 
r is an array with two components giving the lower and upper bounds from which random values for 
the interval vector/matrix will get picked. eg. r=[-5 5] so that lb and ub of each interval entry will be 
bw -5 and 5.
m x n is the order of the array.

%}
if m==1  % row interval vec with n components.
   lb_X=randi([r(1)  (r(1)+r(2))/2],1,n);ub_X=zeros(1,n);
   for i=1:n
       ub_X(i)=randi([lb_X(i) r(2)],1,1);
   end
   X=BiasHull(lb_X,ub_X);
elseif n==1 % column interval vec with m components.
   lb_X=randi([r(1) (r(1)+r(2))/2],m,1);ub_X=zeros(m,1);
   for i=1:m
       ub_X(i)=randi([lb_X(i) r(2)],1,1);
   end
   X=BiasHull(lb_X,ub_X);
end




end