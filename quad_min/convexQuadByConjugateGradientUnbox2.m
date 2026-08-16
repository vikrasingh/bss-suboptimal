%******************************MS:8/30/2022, 7/24/23,26,8/1 , VS:08/18/23,26 
function [isNew,xnew,fxnew,niters,idExit]=convexQuadByConjugateGradientUnbox2(n,A,d,S,x,epsx,maxDx,epsg,nMax,isSoftStop,isTF,targetF)
% till 08/18 [isNew,xnew,fxnew,niters,idExit]=convexQuadByConjugateGradientUnbox2(n,A,d,c,x,epsx,maxDx,epsg,nMax,isSoftStop,isTF,targetF)
% S is a structure S.Q = right-side block of Q matrix from spectral decomposition such that A=QDQ'  where Q is an orthonormal matrix and D is diagonal with eig. value in the increasing order 
%                  S.nNonZeroEig = no. of non zero eigenvalues  
%                  S.D = "2nd part of Diag vector"  of order nNonZeroEig 
%                  S.diagInv= 1/diag(A)  nx1 vector

%dfeps was not imported.
%function [isNew,xnew,fxnew,nfAdd,idExit]=convexQuadByConjugateGradientUnbox2(n,A,d,c,x,epsx,maxDx,epsg,nMax,isSoftStop,isTF,targetF)
%was function [isNew,xnew,fxnew,nfAdd,idExit]=refineByConjugateGradientUnbox2(n,A,d,c,x,epsx,maxDx,epsg,nMax,isSoftStop,isTF,targetF)
%convex Quad (p.s.d) over whole R^n By Conjugate Gradient alg.
%idExit=0  (best)by sufficient cond. i.e. |g_i|<eps for i=1...n
%      =1  the value of x-step size is small enough, So x is not updated.
%      =2  the value of x-step size gets bigger the max. value allowed. So x is not updated. 
%      =3  x is not updated for other reasons. 
%      =4  (worst)hard stop by max. iter. allowed
    nfAdd=0;
    isNew=0;  %overall update of x 
    niter=0;  %#(iter)
    idExit=3;

    niters=0;   %7/26/23
    n2=100;    
   
    epsgg=epsg*epsg;  %MS22h18

    epsxScaled=epsx/(1+A(1,1));  %8/1/23:  
        
    toStop=0;  
    % updated vs 08/21/23
    uplimit=0.5*n;
 while toStop==0   
    %assume x:
    isUpdate=0;   %update of x(j), per sweep   
    % updated VS 08/21/23
    if S.nNonZeroEig < uplimit
       g1=(S.Q)*( (S.D).*( (S.Q)'*x ) );  % reduced dim.
    else
       g1=A*x; 
    end
    g=g1+d;

    nfAdd=1+nfAdd;   %g=A*x+d;
    if (max(abs(g))) <= epsg  
        idExit=0; 
        %toStop=1;
        xnew=x;
        fxnew=0.5*(x'*g1)+d'*x;     %MS:8/25/23 c=0 is assumed
        nfAdd=1+nfAdd; %gnew=g1+d;
        niters=niter;   %7/26/23
        return;       
    end
    %8/25/2022:
    if isTF ~= 0
      if mod(niter,n2)==0                         %7/26/23
        fxnew=0.5*(x'*g1)+d'*x; 
        nfAdd=1+nfAdd;
        if fxnew<=targetF+10^(-15)
            idExit=0; 
            %toStop=1;
            xnew=x; 
            niters=niter;   %7/26/23
            return;       
        end
      end
    end     
            
    if niter==0
      %==== S-D step:
      h=-g;
    else
      %normal CG step:
      beta=(g'*Ah0)/h0Ah0;       %7/26/23  beta=g'*Ah0/h0Ah0;
      h=-g+beta*h0;   
    end
    [isDxSmall,isHhSmall,x,Ah0,h0Ah0]=exactLineSearchUnbox4(n,A,d,S,x,epsxScaled,maxDx,epsg,epsgg,g,h);  %MS:8/26/23 _c0 version
    niter=1+niter;
    %debug:
    %g1=A*x; f1=0.5*(x'*g1)+d'*x+c;
    %fprintf('At iter %d: after SD got f= %d \n',niter,f1);  

    if isHhSmall==1, idExit=0; toStop=1; break; end   
    if isDxSmall==1, idExit=1; toStop=1; break; end   
    
    isUpdate=1;
    isNew=1;
    h0=h;      
    
    if isSoftStop==0 && niter>=nMax
        idExit=5;  
        toStop=1; 
        break;        
    end
    
    if isUpdate==0
        idExit=3;        
        toStop=1; 
        break;
    end                 
            
 end  %while toStop==0

 xnew=x;
 % updated VS 08/21/23,26
if S.nNonZeroEig < uplimit
    Qx=(S.Q)'*x; 
    fxnew=0.5*( (Qx)'*( (S.D).*(Qx) ) ) + d'*x;    %MS8/26/23: assume c=0.
else
   g1=A*x; 
   fxnew=0.5*(x'*g1)+d'*x; 
end

nfAdd=1+nfAdd; %gnew=g+d; 
niters=niter;   %7/26/23

end


%8/30/22 ********************************************************************
function [isDxSmall,isHhSmall,xnew,Ah,hAh]=exactLineSearchUnbox4_notUsed(n,A,d,c,x,epsx,maxDx,epsg,epsgg,g,h)
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
        Ah=A*h;
        hAh=h'*Ah;   %hAh=h'*A*h;       
    end
    if hAh<=epsgg   %This is actually impossible now (in unbox case) since hh>epsgg>0. When hAh=0, hh must be 0. So optomal. So block this if-statement
        isHhSmall=1;       
        return; 
    end  
            
    tstar=(hh/hAh)*h;   %dx*
    
    %---check possible exits:
    aa=max(abs(tstar));
    if aa <=epsx  
        isDxSmall=1; 
        return;    
    end

    %if aa>=maxDx    %as protection against dAd=0     8/18/22 blocked it. 
    %end 
      
    xnew=x+tstar;  %fx is not updated yet to save cpu time
               
end

%MS:3/22/22*************************** new
function [isBreak,isUpdate,isNew,niterNew,idExit,toStop,xnew,gnew]=exactLineSearchUnbox2(n,A,d,c,x,epsx,maxDx,epsg,nMax,isSoftStop,niter,g,h,isSDstep)
    isBreak=0;
    niterNew=niter;
    isUpdate=0; isNew=0;
    idExit=3; toStop=0; 
    xnew=x; gnew=g;
    
    if isSDstep==1 
       hh=h'*h;
    else
       hh=-(h'*g);   %hh=-(h'*gxold);
    end    
    hAh=h'*A*h;           
    tstar=(hh/hAh)*h;   %dx*
    
    %=== Debugging =====================================================================================
    fprintf('Iter=%d\n',niter);fprintf('hAh=%1.4f\n',hAh);fprintf('hh=%1.4f\n',h'*h);fprintf('hh/hAh=%1.4f\n',(h'*h)/hAh );
    disp('tstar=');tstar
    %=======================================================================================================

    %---check possible exits:
    aa=max(abs(tstar));
    if aa <=epsx    %
        idExit=1; 
        toStop=1;
        isBreak=1;  % break;        
    end
    if aa>=maxDx    %as protection against dAd=0      
        idExit=2; 
        toStop=1;
        isBreak=1;  % break;        
    end 
      
    isUpdate=1; isNew=1;
    niterNew=1+niter;
    xnew=x+tstar;  %fx is not updated yet to save cpu time
    %if isSDstep==1
    %   xnew=x+tstar  %fx is not updated yet to save cpu time
    %else
    %   xnew=xold+tstar 
    %end
              
    %---check other possible exits:
    gnew=A*xnew+d;  
    if (max(abs(gnew))) <= epsg           
        idExit=0; 
        toStop=1; 
        isBreak=1;  % break;
    end 
    if isSoftStop==0 && niterNew>=nMax
        idExit=4;  
        toStop=1; 
        isBreak=1;  % break;        
    end
    if isUpdate==0
        idExit=3;        
        toStop=1; 
        isBreak=1;  % break;
    end   
end

%MS:3/22/22***************************
function[isBreak,isUpdate,isNew,niterNew,idExit,toStop,xnew,gnew]=exactLineSearchUnbox(n,A,d,c,x,xold,epsx,maxDx,epsg,nMax,isSoftStop,niter,g,h,gxold,isSDstep)
    isBreak=0;
    niterNew=niter;
    isUpdate=0; isNew=0;
    idExit=3; toStop=0; 
    xnew=x; gnew=g;
    
    if isSDstep==1 
       hh=h'*h;
    else
       hh=-(h'*gxold); 
    end    
    hAh=h'*A*h;           
    tstar=((hh)/hAh)*h;   %dx*
    
    %=== Debugging =====================================================================================
    fprintf('Iter=%d\n',niter);fprintf('hAh=%1.4f\n',hAh);fprintf('hh=%1.4f\n',h'*h);fprintf('hh/hAh=%1.4f\n',(h'*h)/hAh );
    disp('tstar=');tstar
    %=======================================================================================================

    %---check possible exits:
    aa=max(abs(tstar));
    if aa <=epsx    %
        idExit=1; 
        toStop=1;
        isBreak=1;  % break;        
    end
    if aa>=maxDx    %as protection against dAd=0      
        idExit=2; 
        toStop=1;
        isBreak=1;  % break;        
    end 
      
    isUpdate=1; isNew=1;
    niterNew=1+niter;
    if isSDstep==1
       xnew=x+tstar  %fx is not updated yet to save cpu time
    else
       xnew=xold+tstar 
    end
              
    %---check other possible exits:
    gnew=A*xnew+d;  
    if (max(abs(gnew))) <= epsg           
        idExit=0; 
        toStop=1; 
        isBreak=1;  % break;
    end 
    if isSoftStop==0 && niterNew>=nMax
        idExit=4;  
        toStop=1; 
        isBreak=1;  % break;        
    end
    if isUpdate==0
        idExit=3;        
        toStop=1; 
        isBreak=1;  % break;
    end   
end

