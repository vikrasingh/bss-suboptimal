%******************************MS:3/31/2022,4/8,11,8/27,8/6/23,7  , VS:%08/18/23,26
function [isNew,xnew,fxnew,niters,idExit]=convexQuadByParallelTangentBox3(n,A,d,S,low,up,x,epsx,epsg,dfeps,nMax,isSoftStop,isTF,targetF)
% till 08/18/23  [isNew,xnew,fxnew,niters,idExit]=convexQuadByParallelTangentBox3(n,A,d,c,low,up,x,epsx,epsg,dfeps,nMax,isSoftStop,isTF,targetF)
% S is a structure S.Q = right-side block of Q matrix from spectral decomposition such that A=QDQ'  where Q is an orthonormal matrix and D is diagonal with eig. value in the increasing order 
%                  S.nNonZeroEig = no. of non zero eigenvalues  
%                  S.D = "2nd part of Diag vector"  of order nNonZeroEig 
%                  S.diagInv= 1/diag(A)  nx1 vector

%for p.s.d. quad minimization w/ >= 0.
%refine (if possible) the current solution x by Parallel Tangent
%idExit=0  (best)by sufficient cond. i.e. |g_i|<eps for i=1...n, or f<=targetF+10^(-15)
%      =1  the value of x-step size is small enough, So x is not updated.
%      (=2  the value of x-step size gets bigger the max. value allowed. So x is not updated. 
%      =3  x is not updated for other reasons. 
%      =4  (worst)hard stop by max. iter. allowed
% dfeps = tol for acceptable decrease in the fx value

    nfAdd=1;  % 1 due to prefx below
    isNew=0;  %overall update of x 
    niter=0;  %#(iter)
    idExit=3; 
    %prefx=0.5*(x'*(A*x))+d'*x+c; % to check the soft stop

    % updated VS 08/21/23
    uplimit=0.5*n;
    if S.nNonZeroEig < uplimit
       Ax=(S.Q)*( (S.D).*( (S.Q)'*x ) );
    else
       Ax=A*x;    
    end
    prefx=0.5*(x'*Ax)+d'*x ; % to check the soft stop
    g=Ax+d;

    n2=100;    %7/24/23
    niters=0;   %7/26/23
    nitersAdd=0;   %7/31/23
    epsxScaled=epsx/(1+A(1,1));  %8/1/23:  
    
    %checkBds=(x<low | x>up);
    %if sum(checkBds)>0
    %   disp('WARNING: Initial sol x is out of bound in convexQuadByParallelTangentBox().') 
    %end
    %MS8/20/22:
    for i=1:1:n
           if x(i)<low(i)     
             x(i) = low(i); isNew=1;
           elseif x(i)>up(i)  
             x(i) = up(i); isNew=1;
           end
    end %i       
 
    epsgg=epsg*epsg;  %MS22h18
    xold=x;  xold2=x;   %8/1/23   just to define them, their values will be correctly set inside while loopt
               
    toStop=0;   
    
 while toStop==0   
    %assumed: x.
    niter=1+niter;   %7/26/23
    isUpdate=0;   %update of x(j), per sweep   
 
    if niter==1, xold=x;  end   %8/1/23 so that also good under soft stopping reset
    
    %{
    %check g-based exit condition:
    g1=A*x;
    g=g1+d; nfAdd=1+nfAdd;
    if (max(abs(g))) <= epsg  
        idExit=0; toStop=1;
        xnew=x; fxnew=0.5*(x'*g1)+d'*x+c; nfAdd=1+nfAdd; %gnew=g1+d;
        return;       
    end
    %}

    %8/25/2022:
    if isTF ~= 0
      if mod(niter,n2)==0                         %7/24/23
        % updated VS 08/21/23  
        if S.nNonZeroEig < uplimit
           Qx=(S.Q)'*x; 
           fxnew=0.5*( (Qx)'*( (S.D).*(Qx) ) ) + d'*x;  
        else
           fxnew=0.5*(x'*(A*x))+d'*x;
        end

        nfAdd=1+nfAdd;
        if fxnew<=targetF+targetF*10^(-7)
            idExit=0; %toStop=1;
            xnew=x;
            niters=niter+nitersAdd;   %7/26/23
            return;       
        end
      end
    end
           
    %==== S-D step:

    %set a feasible descent search dir: h is the box-truncated version of -g
    %set isBdy:
    [isBdy]=Set_isBdy(n,low,up,epsxScaled,x);     
    if isBdy==0    %(x is in interior of X):
      h=-g;   % is a feasible dir
    else %if isBdy==1 (x is in bdy of X):
      h=-g;   %just initialization
      %project h to box:
      bindingLow=(low-epsxScaled<=x    & x<=low+epsxScaled);   %MS22h21: bindingLow=(low<=x     & x<=low+epsx);   %4/4/22
      bindingUp=(up-epsxScaled<=x & x<=up+epsxScaled);  %MS22h21: bindingUp= (up-epsx<=x & x<=up);         %4/4/22
      for j=1:1:n
           if bindingLow(j)==1 && g(j)>0, h(j)=0; end
           if bindingUp(j)==1 && g(j)<0, h(j)=0; end
      end %j            
    end %if isBdy 
    
    %check (projected g)-based exit condition:
    if (max(abs(h))) <= epsg  
        idExit=0; toStop=1;
        xnew=x; 
        % updated VS 08/21/23
        if S.nNonZeroEig < uplimit
           Qx=(S.Q)'*x; 
           fxnew=0.5*( (Qx)'*( (S.D).*(Qx) ) ) + d'*x; 
        else
           fxnew=0.5*(x'*(A*x))+d'*x; 
        end
        
        nfAdd=1+nfAdd; %gnew=g1+d;
        niters=niter+nitersAdd;
        return;       
    end    
        
    %[isBreak,isUpdate,isNew,niter,idExit,toStop,x,g,isBdy]=lineSearchBoxed(n,A,d,c,low,up,x,epsx,epsg,nMax,isSoftStop,niter,g,h,isNew,isBdy);
    isUpdate1=0;
    %[isDxSmall,isHhSmall,x]=lineSearchBoxed2(n,A,d,c,low,up,x,epsxScaled,epsg,epsgg,g,h);
    [isDxSmall,isHhSmall,x,g]=lineSearchBoxed3(n,A,d,S,low,up,x,epsxScaled,epsg,epsgg,g,h);   %MS8/26/23: c=0 version
    %niter=1+niter;
    if isHhSmall==1, idExit=0; toStop=1; break; end   
    if isDxSmall==1, idExit=1; toStop=1; break; end   %fprintf('idExit=1 at PT step. \r');
    isUpdate1=1;
    %isNew=1;
                     
    if niter==1
        xold2=x;    %gxold2=g;
        continue;   %so far, xold remains unchanged
    end
    
    %{
    [isBdy]=Set_isBdy(n,low,up,epsxScaled,x);     
    if isBdy==1
        xold2=x;    %gxold2=g;
        continue;   %so far, xold remains unchanged
    end
    %}

    %===== The other line search step:  (this is skipped, thus no truncated line search, if current x is on bdy)
    %Now isBdy=0:
    xmid=x;  %at mid of a complete iteration
    %Now assume niter>1:
    h=xmid-xold;                     %must be feasible since isBdy=0
   
    %[isBreak,isUpdate,isNew,niter,idExit,toStop,x,g,isBdy]=lineSearchBoxed(n,A,d,c,low,up,x,epsx,epsg,nMax,isSoftStop,niter,g,h,isNew,isBdy);
    isUpdate2=0;
    %[isDxSmall,isHhSmall,x]=lineSearchBoxed2(n,A,d,c,low,up,x,epsxScaled,epsg,epsgg,g,h);
    [isDxSmall,isHhSmall,x,g]=lineSearchBoxed3(n,A,d,S,low,up,x,epsxScaled,epsg,epsgg,g,h);
    if isHhSmall~=1 && isDxSmall~=1    %fprintf('idExit=1 at PT step. \r');
      isUpdate2=1;
    end
       
    xold=xold2;  
    xold2=x; 
 
    if isUpdate1+isUpdate2 > 0
        isUpdate=1; 
        isNew=1; 
    end
    
    %check stopping conditions:  
    if isSoftStop==0 && niter>=nMax
        idExit=5;  
        toStop=1; 
        break;   
    elseif isSoftStop==1 && niter>=nMax 
       % updated VS 08/21/23 
       if S.nNonZeroEig < uplimit
           Qx=(S.Q)'*x; 
           fx=0.5*( (Qx)'*( (S.D).*(Qx) ) ) + d'*x;
       else
           fx=0.5*(x'*(A*x))+d'*x; 
       end  
       
       nfAdd=1+nfAdd;
       if (prefx-fx)<dfeps*prefx % if improvement of fx is not significant
          idExit=4;  
          toStop=1; 
          break;
       end
       nitersAdd=nitersAdd+niter;   %7/26/23
       prefx=fx;niter=0; % continue with while loop   
    end    
    if isUpdate==0
        idExit=3;        
        toStop=1; 
        break;
    end   
    
 end  %while toStop==0

 xnew=x; 
 % updated VS 08/21/23 
if S.nNonZeroEig < uplimit
   Qx=(S.Q)'*x; 
   fxnew=0.5*( (Qx)'*( (S.D).*(Qx) ) ) + d'*x;
else
   fxnew=0.5*(x'*(A*x))+d'*x; 
end 
 
 nfAdd=1+nfAdd; %gnew=g+d;  
 niters=niter+nitersAdd;   %7/26/23

end

%MS:8/20/22,22*************************** 
function [isBdy]=Set_isBdy_moved(n,low,up,epsx,x)
    %set isBdy:
    isBdy=0; 
    bindingLow=(low-1*epsx<=x    & x<=low+epsx);   %MS22h21: bindingLow=(low<=x    & x<=low+epsx);
    bindingUp=(up-epsx<=x & x<=up+1*epsx);  %MS22h21: bindingUp=(up-epsx<=x & x<=up);
    if sum(bindingLow) > 0, isBdy=1; end
    if sum(bindingUp)  > 0, isBdy=1; end   
end    
    
%MS:3/31/22,4/8,11,8/20*************************** 
function [isDxSmall,isHhSmall,xnew]=lineSearchBoxed2_notUsed(n,A,d,c,low,up,x,epsx,epsg,epsgg,g,h)
%assume: h is a feasible descent search dir at x, epsgg=epsg*epsg
    xnew=x;    
    isDxSmall=0; isHhSmall=0;
    
    hh=-(h'*g); 
    %MS:8/18/22
    if hh<=epsgg    %no update
        isHhSmall=1;       
        return; 
    else
        hAh=h'*A*h;       
    end    
               
    aStar=hh/hAh;    %needed
    tstar=aStar*h;
    xStar=x+tstar;
    
    %---check possible exits:
    aa=max(abs(tstar));
    if aa <=epsx    %
        isDxSmall=1; 

        %for debug:
        %fprintf('g= \r'); g'        
        %fprintf('h= \r'); h'        
        %fprintf('hh= %f hAh= %f aStar=hh/hAh= %f \r', hh,hAh,aStar);
        %fprintf('aa= %f epsx= %f  \r', aa,epsx);
        
        return;    
    end

    %update x:
    if sum(low<=xStar & xStar<=up) == n
        x=xStar;                                        
    else      
      %now xStar is out of box:
      %find alpha in (0, aStar) so that xnew = x + alpha*d is on bdy of box: 
      alpha=aStar;   %initialize    
      for i=1:1:n
        if xStar(i)<low(i)     %h(i) must be nonzero in this case
            atp=(low(i)-x(i))/h(i);  %x(i) + atp*h(i) = low(i)
            if alpha>atp, alpha=atp; end
        elseif xStar(i)>up(i)  %h(i) must be nonzero in this case
            atp=(up(i)-x(i))/h(i);   %x(i) + atp*h(i) = up(i)
            if alpha>atp, alpha=atp; end
        end
      end %i
      dxtp=alpha*h;     %4/8/22
      aa=max(abs(dxtp));   
      if aa <=epsx    
        isDxSmall=1; 
        return;          
      end        
      x = x + dxtp;
      for i=1:1:n
           if x(i)<low(i)     
             x(i) = low(i);
           elseif x(i)>up(i)  
             x(i) = up(i);
           end
      end %i        
    end
    
    xnew=x;           

end

%MS:3/31/22,4/8,11*************************** 
function [isBreak,isUpdate,isNew,niter,idExit,toStop, xnew,gnew,isBdy]=lineSearchBoxed(n,A,d,c,low,up,x,epsx,epsg,nMax,isSoftStop,niter0,g,h,isNew0,isBdy0)
%assume: h is a feasible descent search dir at x
    isBreak=0;            %isBreak=1 to immediately break the while loop in the calling function
    niter=niter0;
    isUpdate=0; isNew=isNew0; 
    isBdy=isBdy0;
    idExit=3; toStop=0;   %toStop=1 to stop (after going through the remaining statements) the while loop in the calling function w/ idExit
    xnew=x; gnew=g;   
    
    hh=-(h'*g);   
    hAh=h'*A*h;           
    aStar=hh/hAh; 
    xStar=x+aStar*h;
    
    %---check possible exits:
    aa=max(abs(aStar*h));
    if aa <=epsx    %
        %idExit=1; 
        %toStop=1;
        %isBreak=1;  %break;  to be determined in the calling function (may take a SD step instead), not here
        return;               
    end

    %update x:
    if sum(low<=xStar & xStar<=up) == n
        x=xStar;                                        
    else      
      %now xStar is out of box:
      %find alpha in (0, aStar) so that xnew = x + alpha*d is on bdy of box: 
      alpha=aStar;   %initialize    
      for i=1:1:n
        if xStar(i)<low(i)     %h(i) must be nonzero in this case
            atp=(low(i)-x(i))/h(i);  %x(i) + atp*h(i) = low(i)
            if alpha>atp, alpha=atp; end
        elseif xStar(i)>up(i)  %h(i) must be nonzero in this case
            atp=(up(i)-x(i))/h(i);   %x(i) + atp*h(i) = up(i)
            if alpha>atp, alpha=atp; end
        end
      end %i
      dxtp=alpha*h;     %4/8/22
      aa=max(abs(dxtp));   
      if aa <=epsx    
        return;          
      end        
      x = x + dxtp;
      if aStar ~= alpha     %MS: 4/11/22    
         for i=1:1:n
           if x(i)<low(i)     
             x(i) = low(i);
           elseif x(i)>up(i)  
             x(i) = up(i);
           end
         end %i      end
      end
    end            
    isUpdate=1+isUpdate; isNew=1;
    niter=1+niter;
    xnew=x;
    %set isBdy:
    isBdy=0; 
    bindingLow=(low<=x    & x<=low+epsx);
    bindingUp=(up-epsx<=x & x<=up);
    if sum(bindingLow) > 0, isBdy=1; end
    if sum(bindingUp)  > 0, isBdy=1; end
            
    %---check other possible exits: 
    gnew=A*xnew+d;  
    if (max(abs(gnew))) <= epsg           
        idExit=0; 
        toStop=1; 
        isBreak=1;  % break;
    end 
    if isSoftStop==0 && niter>=nMax
        idExit=4;  
        toStop=1; 
        isBreak=1;  %break;        
    end
    %if isUpdate==0  %impossible to happen at this point. Let the calling function to decide wether to break
        %idExit=3;        
        %toStop=1; 
        %isBreak=1;  %break;
    %end             

end



