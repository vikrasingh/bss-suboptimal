%MS:8/20/22,22*************************** 
function [isBdy]=Set_isBdy(n,low,up,epsx,x)
    %set isBdy:
    isBdy=0; 
    bindingLow=(low-1*epsx<=x    & x<=low+epsx);   %MS22h21: bindingLow=(low<=x    & x<=low+epsx);
    bindingUp=(up-epsx<=x & x<=up+1*epsx);  %MS22h21: bindingUp=(up-epsx<=x & x<=up);
    if sum(bindingLow) > 0, isBdy=1; end
    if sum(bindingUp)  > 0, isBdy=1; end   
end   