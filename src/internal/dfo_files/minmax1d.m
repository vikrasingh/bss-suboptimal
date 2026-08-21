function [Mu]=minmax1d(p,Q,b,c)
% Procedure to find bound Mu where ||beta||_infty <= Mu , using the method given in Sec. 2.3.2 of Bertsimas et.al. 2016, BSS via Modern optimization lens

    x0=zeros(p,1);
    ub=c;
    Mu=- (10^9);
    for i=1:p
        fun=@(x)objfun(x,i);
    
        constr=@(x)quadconstr(x,Q,b,c-ub);
        % options=optimoptions(@fmincon,'Algorithm','interior-point','SpecifyObjectiveGradient',true);
        options=optimoptions('fmincon','Algorithm','interior-point','MaxIterations',15000,'MaxFunctionEvaluations',100000,'Display','off');
        [~,fval]=fmincon(fun,x0,[],[],[],[],[],[],constr,options);
        gfun=@(x)-objfun(x,i);
        [~,gval]=fmincon(gfun,x0,[],[],[],[],[],[],constr,options);
        Mi=max( abs(fval),abs(gval) );
        if Mi>Mu
            Mu=Mi;
        end
    
    end


end

function [ineq,eq]=quadconstr(x,Q,b,c)
  
   ineq=0.5*(x'*Q*x)+b'*x+c;
   eq=[];
end

function f=objfun(y,i)
   f=y(i);
end

%***************************************************************************************************