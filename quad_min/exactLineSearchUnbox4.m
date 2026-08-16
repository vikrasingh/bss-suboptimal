%8/30/22
%********************************************************************%%8/1/23,26
function [isDxSmall,isHhSmall,xnew,Ah,hAh]=exactLineSearchUnbox4(n,A,d,S,x,epsx,maxDx,epsg,epsgg,g,h)
%h: search dir (unidirectional) at x. May not be g-related. May not be descent.
%assume: convex
    xnew=x; 
    isDxSmall=0; isHhSmall=0;
    Ah=zeros(n,1); hAh=0;
    
    hh=-(h'*g);    

    %if hh<0, x is already optimal along h.
    if hh<=epsgg    %no update. When hh=0, h=0. So optimal and stop.
        isHhSmall=1;       
        return; 
    else
        if S.nNonZeroEig<n/2
          Ah=(S.Q)*( (S.D).*( (S.Q)'*h ) );  % A*h;
        else
          Ah=A*h;
        end
        hAh=h'*Ah;   %hAh=h'*A*h;       
    end
    %{
    if hAh<=epsgg   %This is actually impossible now (in unbox case) since hh>epsgg>0. When hAh=0, hh must be 0. So optomal. So block this if-statement
        isHhSmall=1;       
        return; 
    end  
    %}       
    tstar=(hh/hAh)*h;   %dx*
    
    %     
    %---check possible exits:
    aa=max(abs(tstar));
    if aa <=epsx      %8/1/23  now use epsxScaled=epsx*(1+A(1,1))
        isDxSmall=1; 
        return;    
    end
    %

    %if aa>=maxDx    %as protection against dAd=0     8/18/22 blocked it. 
    %end 
      
    xnew=x+tstar;  %fx is not updated yet to save cpu time
               
end