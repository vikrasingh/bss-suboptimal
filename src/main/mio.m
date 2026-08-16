function [stopflag,fx_tilde,x_tilde,xstar,fstar]=mio(p,n,X,y,k,Q,q,c,Ml,Mu,diagInv,timelimit,iPara,rPara,IotherPara,IstopCondPara,xRelaxedOpt,targetfbest)
% VS:11/10/22, MATLAB implementation of the bss in R package, following the code from Ryan Tibshirani's github page
% https://github.com/ryantibs/best-subset/blob/master/bestsubset/R/bs.R
% which further uses the setup from Bertsimas et.al.(2015) BSS via modern optimization lens. reference


% definition of each stopflag status 
% https://www.gurobi.com/documentation/10.0/refman/optimization_status_codes.html#sec:StatusCodes

%  4 Sep 2023, adding the structure below to be passed to quadratic minimization package
%  S is a structure S.Q = Q matrix from spectral decomposition such that A=QDQ'  where Q is an orthonormal matrix and D is diagonal with eig. value in the decreasing order 
%                  S.nNonZeroEig = no. of non zero eigenvalues  
%                  S.D = diag vector  of order nNonZeroEig x 1 of the D matrix 
%                  S.diagInv= 1/diag(A)  nx1 vector


% Input: Mu=xRelaxedopt like the way we calculate in setupforIntervalAlgo.m filr
   % X'X matrix to be used for step length in projected gradient method
   itime=cputime; % record the initial cputime 
   XtX=X'*X;
   if IotherPara(22)==0 % just use truncated xRelaxed as initial point
      x_tilde=zeros(p,1); 
      [~,insupp]=maxk(abs(xRelaxedOpt),k);
      x_tilde(insupp)=xRelaxedOpt(insupp);
      fx_tilde=fx(x_tilde,p,Q,q,c);
   else
      tempIstopCondPara=IstopCondPara;tempIstopCondPara(6)=0.1*timelimit; % give ag9 not more than 10% of the total cputime limit for finding x0
      [~,~,~,~,x_tilde,fx_tilde]=AG9(p,Ml,Mu,Q,q,c,k,xRelaxedOpt,[],diagInv,targetfbest,IotherPara,tempIstopCondPara,iPara,rPara,[],0);  % ixstar=[]; textfilename=[]; toDebug=0; 
   end
   
   I=eye(p); % identity matrix of order pxp
   if p<=n  % Form 1, using eq. 2.5 on p822
      rvec=[zeros(p,1);ones(p,1)];
      model.A=sparse([cat(2,I,-Mu.*I); cat(2,-I,Ml.*I); rvec' ]); % for non unifrom box, make Mu a px1 vector |xRelaxedopt| and change the operation to ,
      model.sense=repmat('<',2*p+1,1);
      model.rhs=[zeros(2*p,1) ; k];
      model.lb=[Ml.*ones(p,1);zeros(p,1) ];
      model.ub=[Mu.*ones(p,1);ones(p,1) ];
      model.obj=[-2*(X'*y); zeros(p,1)];
      model.Q=sparse( blkdiag(XtX,zeros(p,p)) );
      model.vtype=[repmat('C',p,1); repmat('B',p,1)]; % variable type cont. or binary
      model.modelsense='min'; % no need to provide, as default is min.
      model.start=[x_tilde; logical(x_tilde) ]; % beta=x_tilde
   else % form 2, using eq 2.6 on p823
      rvec=[zeros(p,1);ones(p,1);zeros(n,1)];
      model.A=sparse( [cat(2,I,-Mu.*I,zeros(p,n)); cat(2,-I,Ml.*I,zeros(p,n)); rvec'; cat(2,X,zeros(n,p),-1*eye(n)) ] );
      model.sense=[repmat('<',2*p+1,1); repmat('=',n,1)];
      model.rhs=[zeros(2*p,1) ; k; zeros(n,1)];
      sortedMatrix=sort(abs(X),2,'descend'); % sort the rows abs(X) matrix in descending order  
      Muzeta=max(sum(sortedMatrix(:,1:k),2))*max( max(abs(Mu)),max(abs(Ml)) );   % using theorem 2.1 part d on page825
      model.lb=[Ml.*ones(p,1); zeros(p,1) ; -Muzeta*ones(n,1)]; % -Mu^zeta<= psi_j <= Mu^zeta for j=1...n
      model.ub=[Mu.*ones(p,1); ones(p,1); Muzeta*ones(n,1) ];
      model.obj=[-2*(X'*y); zeros(p+n,1)];
      model.Q=sparse( blkdiag(zeros(2*p,2*p), eye(n)) );
      model.vtype=[repmat('C',p,1); repmat('B',p,1); repmat('C',n,1)]; % variable type cont. or binary
      model.start=[x_tilde;logical(x_tilde); (X*x_tilde)]; % beta=x_tilde
   end

   params.TimeLimit=timelimit*60 -(cputime-itime);  % TimeLimit parameter for gurobi in seconds
   params.Threads=1; 
   params.BestObjStop=targetfbest-(y'*y);   % if objective fn. value of a feasible point is <= targetfbest return
   params.outputflag=0; % to suppress the intermediate from gurobi
   if isunix  % if running this file on HPC cluster, need to add the following commands
      addpath('/apps/anvil/external/apps/gurobi/9.5.1/matlab');
      setenv('GRB_LICENSE_FILE', '/apps/anvil/external/apps/gurobi/9.5.1/license/gurobi.lic');
     % addpath('/apps/anvil/external/apps/gurobi/9.5.1/matlab/'); % will connect to the folder where gurobi is installed on cluster
     % run('/apps/anvil/external/apps/gurobi/9.5.1/matlab/gurobi_setup.m'); % will run gurobi_setup.m file to connect MATLAB with gurobi
   end
   gurout=gurobi(model,params);
   if isfield(gurout,'x')
      xstar=gurout.x(1:p);fstar=fx(xstar,p,Q,q,c);
      [stopflag]=stopStatusToFlag(gurout.status); 
      if stopflag==2 % which is the optimal flag for gurobi
         stopflag=0; % because we are using 0 as the optimal flag, else output the gurobi output status flag 
      end
   else
      xstar=x_tilde;fstar=fx_tilde;stopflag=-2; % dummy stopflag    
   end 

end % bssMATLAB ===================================================================================================

function [stopCode]=stopStatusToFlag(str)
% definition of each stopflag status 
% https://www.gurobi.com/documentation/10.0/refman/optimization_status_codes.html#sec:StatusCodes
stopCode=-1;   
   if strcmp(str,'LOADED'), stopCode=1;
   elseif strcmp(str,'OPTIMAL'), stopCode=2;
   elseif strcmp(str,'INFEASIBLE'), stopCode=3;
   elseif strcmp(str,'INF_OR_UNBD'), stopCode=4;
   elseif strcmp(str,'UNBOUNDED'), stopCode=5;
   elseif strcmp(str,'CUTOFF'), stopCode=6;
   elseif strcmp(str,'ITERATION_LIMIT'), stopCode=7;
   elseif strcmp(str,'NODE_LIMIT'), stopCode=8;
   elseif strcmp(str,'TIME_LIMIT'), stopCode=9;
   elseif strcmp(str,'SOLUTION_LIMIT'), stopCode=10;
   elseif strcmp(str,'INTERRUPTED'), stopCode=11;
   elseif strcmp(str,'NUMERIC'), stopCode=12;
   elseif strcmp(str,'SUBOPTIMAL'), stopCode=13;
   elseif strcmp(str,'INPROGRESS'), stopCode=14; 
   elseif strcmp(str,'USER_OBJ_LIMIT'), stopCode=15;
   elseif strcmp(str,'WORK_LIMIT'), stopCode=16;    
   end

end %==================================================

