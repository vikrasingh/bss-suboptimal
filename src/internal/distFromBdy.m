function [distBdry]=distFromBdy(pDim,IotherPara,Xstar,xRelaxedOpt2)
% VS: 08/27/22
% if Xstar is outside the box, distBdry is the max distance of Xstar(i) for i=1...pDim where Xstar(i) is outside the box
% if Xstar is inside the box, distBdry is the min. distance of Xstar(i) from its closest face for i=1...pDim
%-------------------------------------------------------------------------------------
% Called by runAnEgLLS.m

distOutBdry=0; % will give the max. distance of all components which are outside the box, sign will be +ive
distInBdry=-inf; % will give the min. distance from the bdry of all the components which are inside the box, sign will be -ive 
    if sign(IotherPara(14))==-1 % all quadrant box
       for i=1:pDim
          if Xstar(i)>abs(xRelaxedOpt2(i)), distOutBdry=max( distOutBdry,Xstar(i)-abs(xRelaxedOpt2(i)) );
          elseif Xstar(i)<-abs(xRelaxedOpt2(i)),  distOutBdry=max(distOutBdry, -abs(xRelaxedOpt2(i))-Xstar(i) );
          elseif Xstar(i)>=0, distInBdry=max(distInBdry, Xstar(i)-abs(xRelaxedOpt2(i)) );
          elseif Xstar(i)<=0, distInBdry=max(distInBdry, -abs(xRelaxedOpt2(i))-Xstar(i) );
          end
       end
      
    else  % if one quadrant box 
       for i=1:pDim
           if sign(xRelaxedOpt2(i))>0
              if Xstar(i)>xRelaxedOpt2(i), distOutBdry=max( distOutBdry,Xstar(i)-xRelaxedOpt2(i) );
              elseif Xstar(i)<0,  distOutBdry=max(distOutBdry, 0-Xstar(i) );
              else %  0<= Xstar(i) <= xRelaxedOpt2(i)
                  distInBdry=max( [distInBdry, -Xstar(i), Xstar(i)-xRelaxedOpt2(i)] );
              end
           else   % if signFlag(i)<0,
               if Xstar(i)<xRelaxedOpt2(i), distOutBdry=max( distOutBdry,xRelaxedOpt2(i)-Xstar(i) );
               elseif Xstar(i)>0, distOutBdry=max( distOutBdry,Xstar(i)-0 );
               else % - |xRelaxedOpt2(i)| <= Xstar(i) <=0
                   distInBdry=max( [distInBdry, Xstar(i) , xRelaxedOpt2(i)-Xstar(i)] ); 
               end
           end
       end
    end
   
    if distOutBdry>0,distBdry=distOutBdry;
    else, distBdry=distInBdry; 
    end


end