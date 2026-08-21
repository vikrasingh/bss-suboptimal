
function  [out]=gen_model_para(p_dim,n_pts,k_values,num_instances,numOfSets,numOfEgs,seed_for_num_instances,seed_for_set_of_exs,...
                            seed_for_ex,chooseEgsToRun,snr,mu,sigma, giveValuesOf_b,chooseParaToRun,IotherPara,IstopCondPara, ...
                            QuadMinFunPara,rPara,iPara,inp_fil_ctr,p_ctr,toUseRealDataFile,toDebug,...
                            tod_date,level8dirpath)
% last updated VS 11/16/23
%  for a given no. of nInstances, seedToRunNInstances will generate 1 x nInstances  dim array 
%  seedForEachInstance is the seed to generate every instance of the run
%  for every instance seedToRunSetOfEgs 1 x numOfEgs dim array  which will be used to generate xMatrix,yArray and epsilon for each example
%                                    seedToRunNInstances
%                seedToRunNInstances(1)               seedToRunNInstances(2)                               
%                 seedToRunSetOfEgs                     seedToRunSetOfEgs
%           seedForTheEg=seedToRunSetOfEgs(1)       seedForTheEg=seedToRunSetOfEgs(1)
%           seedForTheEg=seedToRunSetOfEgs(2)       seedForTheEg=seedToRunSetOfEgs(2)
%---------------------------------------------------------------------------------------------------------------------
% output: out is a structure with the following fields-
% out.topDirString: save the topmost folder name 
% out.bestfbestoverSet: save the min(fbest) for all sets for all tmax values to be provided in the example
% input file. Note: bestfbestoverSet is the min fbest values for each tmax in increasing order.
%---------------------------------------------------------------------------------------------------------------------
    
    % initialization of targetcputime, we will change each value to max. hard cpulimit in the paraid loop below
    level6timefileid=fopen(sprintf('timestampL6pDim%d.txt',p_dim),'a+');intimelevel6=datetime('now'); % timestamp Level 6 for the starting time
    fprintf(level6timefileid,'Level 6 \n');fprintf(level6timefileid,'start time=%s \n',intimelevel6); incputimelevel6=cputime;fclose(level6timefileid); 
    
    fprintf('pDim=%d\n',p_dim);          
    if seed_for_num_instances==0
        seed_for_num_instances=randi([1 200],1,num_instances);
    end

    %% Loop to go over all the instances
    %==========================================================================================================================================================================================
    for num_inst_ctr=1:num_instances
        
        seed_for_each_inst=seed_for_num_instances(num_inst_ctr);
        rng(seed_for_each_inst,'twister');
          
        %% Define the seed to generate random examples
        if seed_for_set_of_exs==0
            seed_for_set_of_exs=randi([1 500],numOfEgs,1);
        end
        % topmost directory folder
        top_dir_str=sprintf('np%dInst%dDim%d_',p_ctr,num_inst_ctr,p_dim); % this format is same as used in main input file, imp for ranking 
        out.topDirString=top_dir_str;
        level4dirpath=pwd;
        mkdir(top_dir_str);
        cd(top_dir_str);   % change current directory to the new one

        level5timefileid=fopen('timestampL5.txt','a+');intimelevel5=datetime('now'); % timestamp Level 5 for the starting time
        fprintf(level5timefileid,'Level 5 \n');fprintf(level5timefileid,'start time=%s \n',intimelevel5); incputimelevel5=cputime;fprintf(level5timefileid,'Set;wallclocktime(min);cputime(min) \n');fclose(level5timefileid); 

        mkdir('SummaryEachSetAllEgs'); % create a new directory
        mkdir('SummaryEachEgAllSets');
        
        totalClockTime=tic;  % save the time for whole code.
        totalCPUtime=cputime; % save the cputime for the whole code.
        bestfbestoverSet=inf(1,length(k_values)); % initialization, to save best fbest over all sets 
        saveTableAllEgs1Set=array2table(zeros(0));saveTable1Eg1Set=array2table(zeros(0));saveTableAllEgsAllSets=array2table(zeros(0));% dummy variables, just to initialize
        
        %% SETS LOOP
        [alphabets]=letters; % to save table in excel files
        for set_ctr=1:numOfSets  % outermost loop to go over the no. of sets given as input
                         
            numOfParaVal=length(k_values);
            
            setFolderLevel4=sprintf('Set%d',set_ctr);mkdir(fullfile(pwd,setFolderLevel4)); % make the outermost folder name to be Set_i
            outerFolderNameArray=cell(1,numOfEgs);  % initialize the caption array

            level4timefileid=fopen(fullfile(pwd,setFolderLevel4,'timestampL4.txt'),'a+');intimelevel4=datetime('now'); % timestamp Level 4 for the starting time
            fprintf(level4timefileid,'Level 4 \n');fprintf(level4timefileid,'start time=%s \n',intimelevel4); incputimelevel4=cputime;fclose(level4timefileid); 
            
            %% EXAMPLE LOOP
            %=================================================================================================================================================
            for ex_ctr=1:numOfEgs   % 2nd loop for the given no. of examples
            %=====================================================================================================================================================
                
                % set the seed for each example
                if isequal(seed_for_ex,0)  
                   seed_for_each_ex=seed_for_set_of_exs(ex_ctr);
                else   % i.e. re-running one particular example  
                   seed_for_each_ex=seed_for_ex(ex_ctr); 
                end
                rng(seed_for_each_ex,'twister'); % give each example data a seed to reproduce the same example later on

                
                % generate model parameters, y_obs, X_des and epsilon
                if isempty(toUseRealDataFile) % if the toUseRealDataFile parameter is empty, generate random xMatrix 
                    X_des=mvnrnd(mu,sigma,n_pts);
                    X_des=normalize(X_des,'center'); % make mean of columns 0
                    X_des=normalize(X_des,'norm');  % make 2-norm of columns 1
                    
                else  % if running real data set  
                    if isequal(toUseRealDataFile{1},'OzoneData')
                       X_des=readmatrix( fullfile(level8dirpath,'OzoneData.xlsx') , 'Range','A3:AR332'  );                          
                       X_des=normalize(X_des,'center'); % make mean of columns 0
                       X_des=normalize(X_des,'norm');  % make 2-norm of columns 1
                    elseif isequal(toUseRealDataFile{1},'LeukemiaData')   
                       X_des=readmatrix( fullfile(level8dirpath,'LeukemiaData.xlsx') , 'Range','A3:EGI74'  );
                       X_des=normalize(X_des,'center'); % make mean of columns 0
                       X_des=normalize(X_des,'norm');  % make 2-norm of columns 1   
                    else, error('X_des is not setup for the realdata ex');   
                    end

                end
               
                trueb=giveValuesOf_b(:,chooseEgsToRun(ex_ctr)); % will use the trueb value defined in the chooseEgsToRun array
                Xtrueb=X_des*trueb; 
                sigmaNoise=sqrt( (Xtrueb'*Xtrueb)/snr ); % standard deviation in the noise vector for the linear model
                epsilon=randn(n_pts,1)*sigmaNoise; % noise vector
                    
                if isempty(toUseRealDataFile) % if the toUseEgDataFile parameter is empty               
                   y_obs=X_des*trueb+epsilon;
                else % if testing real data sets
                    if isequal(toUseRealDataFile{1},'OzoneData')
                         y_obs=readmatrix( fullfile(level8dirpath,'OzoneData.xlsx'), 'Range','A333:A662'  );     
                    elseif isequal(toUseRealDataFile{1},'LeukemiaData')
                         y_obs=readmatrix( fullfile(level8dirpath,'LeukemiaData.xlsx') , 'Range','A75:A146'  );       
                    end
                    y_obs=normalize(y_obs,'center'); % make mean of columns 0
                    y_obs=normalize(y_obs,'norm');  % make 2-norm of columns 1
                end
                     
                % Define the 2nd level outer folder name containing ex id
                exFolderLevel3=sprintf('Eg%dp%d_',chooseEgsToRun(ex_ctr),p_dim);
                level5dirpath=fullfile(pwd,setFolderLevel4,char(erase(exFolderLevel3,'_')));
                mkdir(level5dirpath)
                level3timefileid=fopen(fullfile(level5dirpath,'timestampL3.txt'),'a+');intimelevel3=datetime('now');fprintf(level3timefileid,'Level 3 \n');fprintf(level3timefileid,'start time=%s \n',intimelevel3);incputimelevel3=cputime;fclose(level3timefileid); % timestamp Level 3 for the starting time
                fprintf('Ctr=%d;Starting InpFilPair=%d;pDim=%d;Instance=%d;Set=%d;Example=%d;datetime=%s; ... \n',set_ctr,inp_fil_ctr,p_dim,num_inst_ctr,set_ctr,ex_ctr,intimelevel3);
               
                % Calculate x_ols using the quad. min. 
                [~,~,~,~,cputime_ols,x_ols,~,A,b,c]=getxRelaxedOpt(p_dim,y_obs,X_des,trueb,iPara(:,chooseParaToRun(set_ctr,4)),rPara,...
                          QuadMinFunPara(1,:), IotherPara(chooseParaToRun(set_ctr,1),:) ,IotherPara(chooseParaToRun(set_ctr,1),:) , level5dirpath );  % level5dirpath=[]
               
                summary_table=cell(1,numOfParaVal);  % Initialization
                %% k loop
                %==============================================================================================================================================================================
                for paraCtr=1:numOfParaVal    % 4th loop for all k values
                %==============================================================================================================================================================================    

                    level1_k_dir=sprintf('Eg%d_%d_k%g',p_dim,chooseEgsToRun(ex_ctr),k_values(paraCtr)); 

                    level1_k_dir_path=fullfile(level5dirpath,level1_k_dir);
                    mkdir(level1_k_dir_path)   % to create a separate folder for the example inside the outer folder,  in the current directory

                    % save the timestamp 
                    level1timefileid=fopen(fullfile(level1_k_dir_path,'timestampL1.txt'),'a+');
                    intimelevel1=datetime('now');fprintf(level1timefileid,'Level 1 \n');fprintf(level1timefileid,'start time=%s \n',intimelevel1);incputimelevel1=cputime;fclose(level1timefileid); % timestamp Level 1 for the starting time


                    [rParaOut,numOfRows,~,T,~,f_star]=run_instance(n_pts,p_dim,y_obs,X_des,trueb,A,b,c,k_values,cputime_ols,x_ols,set_ctr,paraCtr,...
                                               IotherPara(chooseParaToRun(set_ctr,1),:),IstopCondPara(chooseParaToRun(set_ctr,2),:),...
                                               iPara(:,chooseParaToRun(set_ctr,4)),rPara,IotherPara(chooseParaToRun(set_ctr,1),15),...
                                               level1_k_dir,level1_k_dir_path,tod_date,toDebug);

                       
                    if f_star<bestfbestoverSet(paraCtr), bestfbestoverSet(paraCtr)=f_star; % update the best fbest over all sets for a particular tmax value
                    end
                    
                    % save the table in the excel file
                    save_table_in_excel(T,level1_k_dir,level1_k_dir_path,set_ctr);
                    
                    summary_table{paraCtr}=T;
                    
                    level1timefileid=fopen(fullfile(level1_k_dir_path,'timestampL1.txt'),'a+');
                    outtimelevel1=datetime('now');fprintf(level1timefileid,'End time=%s \n',outtimelevel1);fprintf(level1timefileid,'Time taken for the run=%1.8f min \n',minutes(outtimelevel1-intimelevel1)); % timestamp level0 at the end
                    outcputimelevel1=(cputime-incputimelevel1)/60;fprintf(level1timefileid,'cputime for the run=%1.6f min \n',outcputimelevel1);fprintf(level1timefileid,'fbest for tmax=%d is %1.8f \n',k_values(paraCtr),f_star);
                    fprintf(level1timefileid,'scaleQP=%1.8f \n',rParaOut.scaleQP);
                    
                    fclose(level1timefileid);

                end  % for TMAX LOOP paraid=1:numOfParaVal
                    
                %% write the output data to the 2nd level excel file.
                summary_table=cat(2,summary_table{:}); % horizontal concatenation of all the tables with different tmax values
                [num_entries,row_indices]=extract_row_indices(numOfParaVal);
                excel_file = sprintf('%sSet%d.xlsx',exFolderLevel3,set_ctr);
                excel_file_path1=fullfile(level5dirpath,excel_file);
                for idx5=1:numOfParaVal % to write table for each k value in the 2nd level excel file.
                    writetable(summary_table( (1+num_entries*(idx5-1)):(num_entries+num_entries*(idx5-1)),: ),excel_file_path1,'WriteRowNames',true,'Range',sprintf('A%d', (idx5-1)*(num_entries+3)+(num_entries+4) ));
                end
                 
                %% write the Summary-same para. differ Egs.  and Summary-same Egs differ para. set  excel files
                row_id_to_sum=[1:4]; % index of table rows to select the parameters which we want to sum over all tmax/lambda values
                para_sum_col=zeros(numOfRows,1);
                para_sum_col(row_id_to_sum,1)=sum( summary_table{row_id_to_sum,3*(1:numOfParaVal)},2,'omitnan' ); % get the last column, which saves the cputime it takes for the code to run all lambda/tmax values.

                level1_summary_table=[summary_table(:,[1 2 row_indices+2])   array2table(para_sum_col,'VariableNames',{sprintf('ParaSum (%d)',set_ctr)} )];
                excel_file_table2 = sprintf('egsWithSet%d.xlsx',set_ctr);

                excel_file_path2=fullfile(pwd,'SummaryEachSetAllEgs',excel_file_table2); rangeCell=append(alphabets{ 1+(5+numOfParaVal)*(ex_ctr-1) },'1');
                writetable(level1_summary_table,excel_file_path2,'WriteRowNames',true,'Range',rangeCell);

                % For excel files which show summary of one example with different parameter sets.
                excel_file_table3 = sprintf('%s.xlsx',exFolderLevel3);
                excel_file_path3=fullfile(pwd,'SummaryEachEgAllSets',excel_file_table3);rangeCell=append(alphabets{ 1+(5+numOfParaVal)*(set_ctr-1) },'1');                 
                writetable(level1_summary_table,excel_file_path3,'WriteRowNames',true,'Range',rangeCell);
                    
                %% saveTable1Eg1Set data to create plots in the end
                table_temp=level1_summary_table(:,(1:numOfParaVal+2)); 
                saveTable1Eg1Set=table2array(table_temp); outerFolderNameArray{ex_ctr}={exFolderLevel3};
                if ex_ctr==1
                   rowNamesForInstanceTable=table_temp.Properties.RowNames; % extract rownames for avg and med data for nInstances               
                end

                if ex_ctr==1
                   saveTableAllEgs1Set=saveTable1Eg1Set;
                else
                   saveTableAllEgs1Set=[saveTableAllEgs1Set;saveTable1Eg1Set];
                end
            
                level3timefileid=fopen(fullfile(level5dirpath,'timestampL3.txt'),'a+');    
                outtimelevel3=datetime('now');fprintf(level3timefileid,'End time=%s \n',outtimelevel3);fprintf(level3timefileid,'Time taken for the run=%1.8f min \n',minutes(outtimelevel3-intimelevel3) ); % timestamp level 3 at the end
                outcputimelevel3=(cputime-incputimelevel3)/60;fprintf(level3timefileid,'cputime for the run=%1.6f min \n',outcputimelevel3);fclose(level3timefileid);
                fprintf('Ctr=%d;Ending pDim=%d;Instance=%d;Set=%d;Example=%d;datetime=%s;...;Clocktimetaken=%1.6f min.;cputimetaken=%1.6f min. \n',set_ctr,p_dim,num_inst_ctr,set_ctr,ex_ctr,outtimelevel3,minutes(outtimelevel3-intimelevel3),outcputimelevel3);
            %=======================================================================================================================================================================    
            end  % end EXAMPLE LOOP, egid=1:howManyEgs
            %======================================================================================================================================================================

            if set_ctr==1   % to save data of same egs with different para(lambda) and different sets
               saveTableAllEgsAllSets=saveTableAllEgs1Set;
            else
               saveTableAllEgsAllSets=[saveTableAllEgsAllSets  saveTableAllEgs1Set];
            end
 
            level4timefileid=fopen(fullfile(pwd,setFolderLevel4,'timestampL4.txt'),'a+');    
            outtimelevel4=datetime('now');fprintf(level4timefileid,'End time=%s \n',outtimelevel4);fprintf(level4timefileid,'Time taken for the run=%1.8f min \n',minutes(outtimelevel4-intimelevel4)); % timestamp level 4 at the end
            outcputimelevel4=(cputime-incputimelevel4)/60;fprintf(level4timefileid,'cputime for the run=%1.6f min \n',outcputimelevel4);fclose(level4timefileid);
             % save the level4 timestamp file info. in level 5 timestamp file, to compare the results more effectively
            level5timefileid=fopen('timestampL5.txt','a+');fprintf(level5timefileid,'%d;%1.8f;%1.8f \n',set_ctr,minutes(outtimelevel4-intimelevel4),outcputimelevel4);fclose(level5timefileid);
            
        %===========================================================================================================================================================================    
        end % end SET LOOP, set=1:numOfSets
        %============================================================================================================================================================================
        
        
        numOfParaInEachSet=zeros(1,numOfSets);  % vector which will save the length(tmaxValues or lambdaValues) in each set
        truebindx=zeros(1,numOfSets);truebindx(1)=1; % to delete the trueb columns in the saveTableAllEgsAllSets below
        for i=1:numOfSets
           numOfParaInEachSet(i)=length(k_values);
           if i>1
              truebindx(i)=truebindx(i-1)+numOfParaInEachSet(i-1)+1;
           end
        end
        
        saveTableAllEgsAllSetsRemovingTruebCols=saveTableAllEgsAllSets;saveTableAllEgsAllSetsRemovingTruebCols(:,[truebindx  truebindx+1])=[]; % delete the trueb and xRelaxedOpt cols
        
        fileid=fopen('InputParaForThisRun.txt','a+');  % open the file again
        totalWallClockTime=toc(totalClockTime)/60;fprintf(fileid,'total wall clock time = %1.6f min\n',totalWallClockTime);  %fprintf('total wall clock time = %1.6f min\n',totalWallClockTime);
        totalcputime=(cputime-totalCPUtime)/60;fprintf(fileid,'total cputime = %1.6f min\n',totalcputime);       %fprintf('total cputime = %1.6f min\n',totalcputime);
        fclose(fileid);

        cpuDataForPlot=zeros(numOfEgs*numOfParaVal,numOfSets);
        fbestDataForPlot=zeros(numOfEgs*numOfParaVal,numOfSets);
        niterDataForPlot=zeros(numOfEgs*numOfParaVal,numOfSets);
        stopflagData=zeros(numOfEgs*numOfParaVal,numOfSets);
%                                 cpu set1 |  cpu set2|  . . . 
%         cpuDataForPlot=Eg 1 Tm1    -            -
%                        Eg 1 Tm2    -            -
%                        ...
%                        Eg n Tm1    -            -
%                        Eg n Tm2    -            -
        for seti=1:numOfSets
            for egj=1:numOfEgs
                for tmaxk=1:numOfParaVal
                    cpuDataForPlot( (egj-1)*numOfParaVal +tmaxk, seti )=saveTableAllEgsAllSets( (egj-1)*numOfRows + 3, 2+tmaxk+(seti-1)*(numOfParaVal+2) );   %3 because cputime is the 3rd row in the ouput table
                    fbestDataForPlot( (egj-1)*numOfParaVal +tmaxk, seti )=saveTableAllEgsAllSets( (egj-1)*numOfRows + 1, 2+tmaxk+(seti-1)*(numOfParaVal+2) );   % 1 because fbest is the 1st row in the ouput table
                    niterDataForPlot( (egj-1)*numOfParaVal +tmaxk, seti )=saveTableAllEgsAllSets( (egj-1)*numOfRows + 4, 2+tmaxk+(seti-1)*(numOfParaVal+2) );   % 4 because no. of iter. is the 4th row in the ouput table
                    stopflagData( (egj-1)*numOfParaVal +tmaxk, seti )=saveTableAllEgsAllSets( (egj-1)*numOfRows + 2, 2+tmaxk+(seti-1)*(numOfParaVal+2) );   % 2 because stopflag is the 2nd row in the ouput table
                end
            end    
        end

        out.cpuDataForPlot=cpuDataForPlot;out.fbestDataForPlot=fbestDataForPlot;out.niterDataForPlot=niterDataForPlot;out.stopflagData=stopflagData;
        
        level5timefileid=fopen('timestampL5.txt','a+');
        outtimelevel5=datetime('now');fprintf(level5timefileid,'End time=%s \n',outtimelevel5);fprintf(level5timefileid,'Time taken for the run=%1.8f min \n',minutes(outtimelevel5-intimelevel5)); % timestamp level 5 at the end
        outcputimelevel5=(cputime-incputimelevel5)/60;fprintf(level5timefileid,'cputime for the run=%1.6f min \n',outcputimelevel5);fprintf(level5timefileid,append('min(fbest) over all sets for tmax =',repmat('%d;',1,length(k_values)),'\n' ),flip(k_values) );fprintf(level5timefileid, append(repmat('%1.8f ',1,length(k_values) ),'\n') ,flip(bestfbestoverSet));out.bestfbestoverSet=flip(bestfbestoverSet);
        fclose(level5timefileid);
        cd(level4dirpath);   % change the directory back
        
        if num_inst_ctr==1
            saveEveryInstPlotData=struct('data',saveTableAllEgsAllSetsRemovingTruebCols);  % to save the final plotting data for all instances, use a structure saveEveryInstPlotData 
        else
            saveEveryInstPlotData(num_inst_ctr).data=saveTableAllEgsAllSetsRemovingTruebCols;
        end
        fprintf('\n');
    %=============================================================================================================================================================================================    
    end  % end INSTANCES LOOP,  nInstCtr=1:nInstances
    %===========================================================================================================================================================================================
    
    finalSolsFornInstances(p_ctr,p_dim,num_instances,numOfEgs,k_values,numOfSets,IotherPara,chooseParaToRun,saveEveryInstPlotData,outerFolderNameArray,numOfRows,rowNamesForInstanceTable,numOfParaInEachSet,tod_date)
    
    level6timefileid=fopen(sprintf('timestampL6pDim%d.txt',p_dim),'a+');
    outtimelevel6=datetime('now');fprintf(level6timefileid,'End time=%s \n',outtimelevel6);fprintf(level6timefileid,'Time taken for the run=%1.8f min \n',minutes(outtimelevel6-intimelevel6)); % timestamp level 6 at the end
    outcputimelevel6=(cputime-incputimelevel6)/60;fprintf(level6timefileid,'cputime for the run=%1.6f min \n',outcputimelevel6);fclose(level6timefileid);

    fclose('all');
    fprintf('total wall clock time for pDim %d run =%1.6f min. \n',p_dim,toc(totalClockTime)/60);
    fprintf('total cputime for pDim %d run=%1.6f min \n',p_dim,(cputime-totalCPUtime)/60);
%     fprintf('saving ranking excel files ... \n');
    
    
end % gen_model_para function================================================

function  save_table_in_excel(T,innerFolderName,innerFolderPath,set)
% 2July22, save the deepest level table, the #0 in tol. table and the stopFlag reference table

%  8/22/22   str2=append(innerFolderName,sprintf('Set%d-',set),toDate);
    str2=append(innerFolderName,sprintf('Set%d-',set));
    excelFile = sprintf('%s.xlsx',str2);excelFilePath=fullfile(innerFolderPath,excelFile);

    writetable(T,excelFilePath,'WriteRowNames',true,'Range','A1');
 

end %=================================================================

function [num_entries,row_indices]=extract_row_indices(numOfParaVal)
 % the subroutine to get the topmost table of the second level excel file, by manipulating the data                                                         

     num_entries=3;
     row_indices=zeros(1,numOfParaVal); % rowid1, rowid2 indices to be used to extract the rows
     for idx1=1:numOfParaVal   % extract the first 3 rows for different k values 
         row_indices(idx1)=1+num_entries*(idx1-1);
     end
 
end %==========================
