function pairExSumm(fbestDataAllInpFil,cpuDataAllInpFil,InputFiles,varargin)
% pairExSumm( , , )  the first 3 arg. are necessary
% pairExSumm( , , , dfbest) after that we can give dfbest if not using default value
% pairExSumm( , , ,[],dcpu)  will set user defined dcpu with default dfbest
% pairExSumm( , , , dfbest, dcpu) will change both dfbest and dcpu from the default values
% idx: variable that will list the no. of egs in each of the input files, in the same order as given in the main input file 
% this file is for 56 pair input files setup only, need some modification for other types
%{
 Dec13, 22
    if setuptype==1 % empty as of now
    elseif setuptype==2
%         p56PPair_Box7a92c3d3c_28Ex20t2kBd2TmCcRh0p8.m
%         p56PPair_Box7a92c3d3c_28Ex20t2kBd2TmCcRh0p8F0B1.m
%         p56PPair_Box7a92c3d3c_28Ex20t2kBd2TmCcRh0p8F1B1.m
%         p56PPair_Ubx7b92c3d3c_28Ex20t2kBd2TmCcRh0p8.m
%         p56PPair_Ubx7b92c3d3c_28Ex20t2kBd2TmCcRh0p8F0U1.m
%         p56PPair_Ubx7b92c3d3c_28Ex20t2kBd2TmCcRh0p8F1U1.m

        %    small                  medium                                 large
        idx=[8*ones(1,8),   4*ones(1,8), 4*ones(1,8),    2*ones(1,8), 2*ones(1,8), 2*ones(1,8), 2*ones(1,8) ];  % total 192 egs
       
    elseif setuptype==3
%         p56PPair_Box7a92c3d3c_Two28Ex20t2kBd2TmCcRh0p8.m
%         p56PPair_Ubx7b92c3d3c_Two28Ex20t2kBd2TmCcRh0p8.m
        %         small     medium     large
        idx=repmat([8       4   4     2  2 2 2] , 1 , 8);  
        
    elseif setuptype==4
%         p14PSet_Mix5a_Ran12Ex20t2kOd2TmCcRh0p9Sn1.m
        idx=repmat([8       4   4     2  2 2 2] , 1 , 2); 
       
    elseif setuptype==5
%         p14PSet_Mix5b_Ran12Ex20t2kOd2TmCcRh0p9Sn1.m
        idx=[8 8 4 4 4 4 2 2 2 2 2 2 2 2];
    elseif setuptype==6
%         p28PSet_Mix5b_Ran12Ex20t2kOdUd2TmCcRh0p9Sn1.m
       idx=[8*ones(1,4)    4*ones(1,8)   2*ones(1,16)];
    end
%}    

   

    nvarargs=length(varargin);
    if nvarargs==0, dfbest=2;dcpu=1e-5; % default
    elseif nvarargs==1   % we are definitely providing dfbest if varargin has one input variable
       dfbest=varargin{1};dcpu=1e-5;
    elseif nvarargs==2
       if isempty(varargin{1})
          dfbest=2;dcpu=varargin{2};
       else, dfbest=varargin{1};dcpu=varargin{2};   
       end
    end
    % coverted all the main input files in this format now, Dec13,22
    %         small     medium     large
    idx=repmat([8       4   4     2  2 2 2] , 1 , 8); 
    %=========== Here on fixed for all kind of main input files =========================================================================================================
    ninpfil=length(idx);
    fbestBlock=zeros(ninpfil,6); cpuBlock=zeros(ninpfil,6);  % 6 is the no. of algo we used in the algo. input file  
    idx_i_1=1;
    for i=1:ninpfil  % 192 is the no. of egs for each tmax, each pDim, each snr,...
       fbestBlock(i,:)=sum(  fbestDataAllInpFil(idx_i_1 : sum(idx(1:i)),:)  ,1 )/idx(i);
       cpuBlock(i,:)=sum(  cpuDataAllInpFil(idx_i_1:sum(idx(1:i)),:)  ,1 )/idx(i);
       idx_i_1=idx_i_1+idx(i);
    end
    fbestBlock
    cpuBlock
%     % for the 5 set algo. file,   7-34_9-34_7-23_34_9-33
%     diffbestBlock=fbestBlock;
%     diffbestBlock(:,1)=fbestBlock(:,1)-fbestBlock(:,4) ;   % first column 7-34 data
%     diffbestBlock(:,2)=fbestBlock(:,2)-fbestBlock(:,4) ;   % second column 9-34 data
%     diffbestBlock(:,3)=fbestBlock(:,1)-fbestBlock(:,3) ;   % third column 7-23 data
%     diffbestBlock(:,4)=fbestBlock(:,4) ;                      % fourth column  34 data
%     diffbestBlock(:,5)=fbestBlock(:,2)-fbestBlock(:,5) ;   % fifth column 9-33 data
%     sdiffbestBlock=sign(diffbestBlock); 
%     idxsmallvalue=(abs(diffbestBlock)<dfbest);  % find the indices of the values smaller than given tol. dfbest 
%     sdiffbestBlock(idxsmallvalue)=0;  % and set them equal to 0
% 
%     difcpuBlock=cpuBlock;
%     difcpuBlock(:,1)=cpuBlock(:,1)-cpuBlock(:,4) ;   % first column 7-34 data
%     difcpuBlock(:,2)=cpuBlock(:,2)-cpuBlock(:,4) ;   % second column 9-34 data
%     difcpuBlock(:,3)=cpuBlock(:,1)-cpuBlock(:,3) ;   % third column 7-23 data
%     difcpuBlock(:,4)=cpuBlock(:,4) ;                      % fourth column  34 data
%     difcpuBlock(:,5)=cpuBlock(:,2)-cpuBlock(:,5) ;   % fifth column 9-33 data
%     sdifcpuBlock=sign(difcpuBlock);
%     idxsmallvalue=(abs(difcpuBlock)<dcpu);  % find the indices of the values smaller than the given tol. dcpu 
%     sdifcpuBlock(idxsmallvalue)=0;  % and set them equal to 0
%  
%     
%     lefttable=cell2table(InputFiles,'VariableNames',{'alg.file','ex.file'});
%   
%     % Sheet 1  average values for the fbest for each input file pairs
%     middletableSheet1=array2table(fbestBlock,'VariableNames',{'Ag7 fxAvg','Ag9 fxAvg','Ag23 fxAvg','Ag34 fxAvg','Ag35 fxAvg'});
%     righttableSheet1=array2table(cpuBlock,'VariableNames',{'Ag7 cpuAvg','Ag9 cpuAvg','Ag23 cpuAvg','Ag34 cpuAvg','Ag35 cpuAvg'});
%     writetable([lefttable,middletableSheet1,righttableSheet1],'summary5alg.xlsx','Sheet','Avg_val');
%     
%     % Sheet 2 difference in the average values fbest and cpu
%     middletableSheet2=array2table(diffbestBlock,'VariableNames',{'Ag7-Ag34 fxAvg','Ag9-Ag34 fxAvg','Ag7-Ag23 fxAvg','Ag34 fxAvg','Ag9-Ag35 fxAvg'});
%     righttableSheet2=array2table(difcpuBlock,'VariableNames',{'Ag7-Ag34 cpuAvg','Ag9-Ag34 cpuAvg','Ag7-Ag23 cpuAvg','Ag34 cpuAvg','Ag9-Ag35 cpuAvg'});
%     writetable([lefttable,middletableSheet2,righttableSheet2],'summary5alg.xlsx','Sheet','7-34_9-34_7-23_34_9-35');
%     
%     % Sheet 3 sign of the difference in the average values fbest and cpu
%     middletableSheet3=array2table(sdiffbestBlock,'VariableNames',{'Ag7-Ag34 fxAvg','Ag9-Ag34 fxAvg','Ag7-Ag23 fxAvg','Ag34 fxAvg','Ag9-Ag35 fxAvg'});
%     righttableSheet3=array2table(sdifcpuBlock,'VariableNames',{'Ag7-Ag34 cpuAvg','Ag9-Ag34 cpuAvg','Ag7-Ag23 cpuAvg','Ag34 cpuAvg','Ag9-Ag35 cpuAvg'});
%     writetable([lefttable,middletableSheet3,righttableSheet3],'summary5alg.xlsx','Sheet','sign7-34_9-34_7-23_34_9-35');
% 
%     %% Extract groups of examples together
% 
%     % save OD and UD egs. together
%     iod=[1:28]; % indices of over-determined egs. avg.
%     iud=[29:56]; % indices of under-determined egs. avg.
% 
%     writetable([lefttable(iod,:),middletableSheet1(iod,:),righttableSheet1(iod,:)],'summary5alg.xlsx','Sheet','OD egs','Range','A1');
%     writetable([lefttable(iod,:),middletableSheet2(iod,:),righttableSheet2(iod,:)],'summary5alg.xlsx','Sheet','OD egs','Range','N1');
%     writetable([lefttable(iod,:),middletableSheet3(iod,:),righttableSheet3(iod,:)],'summary5alg.xlsx','Sheet','OD egs','Range','AA1');
% 
%     writetable([lefttable(iud,:),middletableSheet1(iud,:),righttableSheet1(iud,:)],'summary5alg.xlsx','Sheet','UD egs','Range','A1');
%     writetable([lefttable(iud,:),middletableSheet2(iud,:),righttableSheet2(iud,:)],'summary5alg.xlsx','Sheet','UD egs','Range','N1');
%     writetable([lefttable(iud,:),middletableSheet3(iud,:),righttableSheet3(iud,:)],'summary5alg.xlsx','Sheet','UD egs','Range','AA1');
% 
%     % save egs. with same dim. group together
%     ismall=[1:7:56]; % indices of small dim. egs. avg.
%     imed=sort([2:7:56, 3:7:56]); % indices of medium dim. egs. avg.
%     ilarge=sort([4:7:56, 5:7:56, 6:7:56, 7:7:56]); % indices of large dim. egs. avg.
% 
%     % sheet for small dim. egs. avg.
%     writetable([lefttable(ismall,:),middletableSheet1(ismall,:),righttableSheet1(ismall,:)],'summary5alg.xlsx','Sheet','small egs','Range','A1');
%     writetable([lefttable(ismall,:),middletableSheet2(ismall,:),righttableSheet2(ismall,:)],'summary5alg.xlsx','Sheet','small egs','Range','N1');
%     writetable([lefttable(ismall,:),middletableSheet3(ismall,:),righttableSheet3(ismall,:)],'summary5alg.xlsx','Sheet','small egs','Range','AA1');
% 
%     % sheet for medium dim. egs. avg.
%     writetable([lefttable(imed,:),middletableSheet1(imed,:),righttableSheet1(imed,:)],'summary5alg.xlsx','Sheet','med egs','Range','A1');
%     writetable([lefttable(imed,:),middletableSheet2(imed,:),righttableSheet2(imed,:)],'summary5alg.xlsx','Sheet','med egs','Range','N1');
%     writetable([lefttable(imed,:),middletableSheet3(imed,:),righttableSheet3(imed,:)],'summary5alg.xlsx','Sheet','med egs','Range','AA1');
% 
%     % sheet for larget dim. egs. avg.
%     writetable([lefttable(ilarge,:),middletableSheet1(ilarge,:),righttableSheet1(ilarge,:)],'summary5alg.xlsx','Sheet','large egs','Range','A1');
%     writetable([lefttable(ilarge,:),middletableSheet2(ilarge,:),righttableSheet2(ilarge,:)],'summary5alg.xlsx','Sheet','large egs','Range','N1');
%     writetable([lefttable(ilarge,:),middletableSheet3(ilarge,:),righttableSheet3(ilarge,:)],'summary5alg.xlsx','Sheet','large egs','Range','AA1');


end


