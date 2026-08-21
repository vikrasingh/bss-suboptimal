function [stopCode]=gurobiStopStatusToFlag(str)
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

end