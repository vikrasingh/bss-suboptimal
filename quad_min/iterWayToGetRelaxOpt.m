function [bestAlgoFlag,lb,ub,saveAllTheSols,cpuxRelaxedOpt,xRelaxedOpt,fxRelaxedOpt]= iterWayToGetRelaxOpt(n,A,b,c,S,choseniPara,rPara,QuadMinFunPara,fileid)    
    
    key=0;    % will help to check if there is a previously saved solution or not.
    flag=0;   % will help us to identify bw first and second iter. solutions
    maxIter=50;  % max no. of times we are going to use quadMinFun for different initial boxes
    fxRelaxedOpt=10^(12);  % for implementation purpose
    enlargeBoxSlightly=0; % initilize
    cpuxRelaxedOpt=0;
    for i=1:maxIter
        lb=-i*10*ones(n,1);ub=i*10*ones(n,1);
        if flag==0    % this loop will only get exceuted in the first iteration
            
           [bestAlgoFlag,saveAllTheSols,cputimeSec,xRelaxedOpt,fxRelaxedOpt]=quadMinFun(n,A,b,c,lb,ub,S,choseniPara,rPara,QuadMinFunPara);
           cpuxRelaxedOpt=cpuxRelaxedOpt+cputimeSec;
           % save the solution to compare with the solution from the next call of quadMinFun 
           savefxRelaxOpt=fxRelaxedOpt;
           flag=1;
           if ~isempty(fileid)
               fprintf(fileid,'No. of times we use quadMinFun for different box X= %d \n',i);
               fprintf(fileid,'lb for solution box X = \n');printArray(-i*10*ones(n,1),'%g',fileid);fprintf(fileid,'ub for solution box X = \n');printArray(i*10*ones(n,1),'%g',fileid);
               fprintf(fileid,'final sol. xRelaxedopt =\n');printArray(xRelaxedOpt,'%g',fileid);
               fprintf(fileid,'final sol. fxRelaxedOpt=%g\n',fxRelaxedOpt);
           end
           
        else
            
           if key==1  % if its not the first iteration compare the result with the previous saved solution
              if abs(fxRelaxedOpt-savefxRelaxOpt)<10^(-4)  % stop. criteria 
                 if ~isempty(fileid)
                     fprintf(fileid,'No. of times we use quadMinFun for different box X= %d \n',i);
                     fprintf(fileid,'lb for solution box X = \n');printArray(-i*10*ones(n,1),'%g',fileid);fprintf(fileid,'ub for solution box X = \n');printArray(i*10*ones(n,1),'%g',fileid);
                     fprintf(fileid,'final sol. xRelaxedopt =\n');printArray(xRelaxedOpt,'%g',fileid);
                     fprintf(fileid,'final sol. fxRelaxedOpt=%g\n',fxRelaxedOpt);
                 end
                 break
              else
                 savefxRelaxOpt=fxRelaxedOpt;  % if stop criteria does not holds, make save current fxRelaxedOpt
              end
           end
        
           [bestAlgoFlag,saveAllTheSols,cputimeSec,xRelaxedOpt,fxRelaxedOpt]=quadMinFun(n,A,b,c,lb,ub,S,choseniPara,rPara,QuadMinFunPara);
           cpuxRelaxedOpt=cpuxRelaxedOpt+cputimeSec;
%            fprintf('Iter. for quadMinFun sequential approach = %d \n',i);xRelaxedOpt,fxRelaxedOpt
           key=1;
           
        end
        
    end
%}
    
%{    
    for i=1:maxIter
        
        if enlargeBoxSlightly==0  
            lb=-i*10*ones(n,1);ub=i*10*ones(n,1);  % enlarge the box twice
        else
            lb=-((4*i)/3)*10*ones(n,1);ub=((4*i)/3)*10*ones(n,1); % enlarge the box slightly
        end
        
       
        [bestAlgoFlag,saveAllTheSols,cpuxRelaxedOpt,xRelaxedOpt,fxRelaxedOpt]=quadMinFun(n,A,b,c,lb,ub,choseniPara,rPara,QuadMinFunPara);
        fprintf('Iter. for quadMinFun sequential approach = %d \n',i);xRelaxedOpt,fxRelaxedOpt
      
        % optimality criteria's check          
        nonoptimalflag=0;  % to make sure that xRelaxedOpt is not on the boundary of the box X
        for j=1:n
           if ub(j)==xRelaxedOpt(j) || lb(j)==xRelaxedOpt(j)
              nonoptimalflag=1;
              break;
           end
        end
        
        stopsign2=1; % initilization, for gradient condition check
        gradf=abs(A*xRelaxedOpt+b); % to check if the gradient of the f is 0 at xRelaxedOpt
        for k=1:n
           if gradf(k)>10^(-4) 
              stopsign2=0;break
           end 
        end
                                 
        if stopsign2==1 && nonoptimalflag==0 % if the grad f(xRelaxedOpt) =0  and xRelaxedOpt is not on the bdry.
           fprintf('No. of times we use quadMinFun for different box X= %d \n',i-1);
           disp('lb for solution box X = '),-(i)*10*ones(n,1),disp('ub for solution box X = '),(i)*10*ones(n,1)
           disp('final sol. xRelaxedopt ='),xRelaxedOpt,disp('final sol. fxRelaxedOpt ='),fxRelaxedOpt
           break
        elseif stopsign2==1 && nonoptimalflag==1 % if grad f(xRelaxedOpt)=0 but xRelaxedOpt is on the bdry.
           disp('xRelaxedOpt is on the boundary, enlarge the box');
        elseif stopsign2==0  % grad f(xRelaxedOpt) is not 0
           disp('Enlarge the box slightly.'); enlargeBoxSlightly=1;                 
        end
        
       
    end    
%}    

    
    
end 