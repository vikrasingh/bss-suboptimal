%MS:3/31/22,4/8,11,8/20,7/20/23,8/7,26 *************************** 
function [isDxSmall,isHhSmall,xnew,gnew]=lineSearchBoxed3(n,A,d,S,low,up,x,epsx,epsg,epsgg,g,h)
%based on function [isDxSmall,isHhSmall,xnew]=lineSearchBoxed2(n,A,d,c,low,up,x,epsx,epsg,epsgg,g,h)
%MS:7/16/23 added as a separate m-file, copied from convexQuadByParallelTangentBox3.m 

%assume: h is a feasible descent search dir at x, epsgg=epsg*epsg
    xnew=x;   gnew=g; 
    isDxSmall=0; isHhSmall=0;
    
    hh=-(h'*g); 
    %MS:8/18/22
    if hh<=epsgg    %no update
        isHhSmall=1;       
        return; 
    else
        %MS:8/26/23:
        if S.nNonZeroEig < n/2
          Ah=(S.Q)*( (S.D).*( (S.Q)'*h ) );  % Ah=A*h;
          hAh=h'*Ah;         
        else
          hAh=h'*A*h;     
        end
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
        if abs(h(i))<epsg, continue; end      %7/20/23
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
    %8/26/23   gnew=A*xnew+d;
    if S.nNonZeroEig < n/2
       Ax=(S.Q)*( (S.D).*( (S.Q)'*x ) );
    else
       Ax=A*x;    
    end
    gnew=Ax+d;

end
