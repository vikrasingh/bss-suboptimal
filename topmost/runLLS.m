
function  [out]=runLLS(pDim,nPts,infoCol,nDcol,alphaValues,lambdaValues,tmaxValues,nInstances,numOfSets,numOfEgs,seedToRunNInstances,seedToRunSetOfEgs,...
                            seedForEg,chooseEgsToRun,snr,mu,sigma, giveValuesOf_b,seedForTruebAllEgs,check_nDcol,targetfbestValues,chooseParaToRun,IotherPara,IstopCondPara, ...
                            delCondPara,QuadMinFunPara,rPara,iPara,InpFilCtr,pDimCtr,saveIntermOutput,inputFileNameLocation,toUseRealDataFile,wantToRe_runAnEg,whichRe_run,toDebug,...
                            toDate,level7dirpath,level8dirpath)
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
    targetcpuValues=ones(1,length(targetfbestValues))*abs( IstopCondPara(chooseParaToRun(1,2),6) ); % hard cpulimit for the first set will be the targetcpuValues for all tmax values all pDim, one example
    numOfBoxDelAllSets=zeros(8,1);  % dim should be same as numOfBoxDel variable in setupForIntvlAlgo.m line 31
    level6timefileid=fopen(sprintf('timestampL6pDim%d.txt',pDim),'a+');intimelevel6=datetime('now'); % timestamp Level 6 for the starting time
    fprintf(level6timefileid,'Level 6 \n');fprintf(level6timefileid,'start time=%s \n',intimelevel6); incputimelevel6=cputime;fclose(level6timefileid); 
    
    fprintf('pDim=%d\n',pDim);          
    if seedToRunNInstances==0
        seedToRunNInstances=randi([1 200],1,nInstances);
    end

    %% Loop to go over all the instances
    %==========================================================================================================================================================================================
    for nInstCtr=1:nInstances
        
        seedForEachInstance=seedToRunNInstances(nInstCtr);
        rng(seedForEachInstance,'twister');
          
        %% Define the seed to generate random examples
        if seedToRunSetOfEgs==0
            seedToRunSetOfEgs=randi([1 500],numOfEgs,1);
        end
        % topmost directory folder
        topDirString=sprintf('np%dInst%dDim%d_',pDimCtr,nInstCtr,pDim); % this format is same as used in main input file, imp for ranking 
        out.topDirString=topDirString;
        level4dirpath=pwd;
        mkdir(topDirString);
        cd(topDirString);   % change current directory to the new one

        level5timefileid=fopen('timestampL5.txt','a+');intimelevel5=datetime('now'); % timestamp Level 5 for the starting time
        fprintf(level5timefileid,'Level 5 \n');fprintf(level5timefileid,'start time=%s \n',intimelevel5); incputimelevel5=cputime;fprintf(level5timefileid,'Set;wallclocktime(min);cputime(min) \n');fclose(level5timefileid); 
        
        totalClockTime=tic;  % save the time for whole code.
        totalCPUtime=cputime; % save the cputime for the whole code.
        bestfbestoverSet=inf(1,length(tmaxValues)); % initialization, to save best fbest over all sets 
        maxdistBdryoverSet=-inf; % Initialization, to save max(distBdry) over all sets
        maxVioRefQMoverSet=-inf; % initialize, max vio of the necess. cond. for QM during refinement over all sets
        maxVioInfQMoverSet=-inf; % initialize, max vio of the necess. cond. for QM during Inf F over all sets
        saveTableAllEgsAllSets=cell(1,numOfSets); %initialization
        saveTableAllEgs1Set=cell(1,numOfEgs);

        %% SETS LOOP
        for setCtr=1:numOfSets  % outermost loop to go over the no. of sets given as input
                         

            numOfParaVal=length(tmaxValues);
            setFolderLevel4=sprintf('Set%d',setCtr);mkdir(fullfile(pwd,setFolderLevel4)); % make the outermost folder name to be Set_i
            
            level4timefileid=fopen(fullfile(pwd,setFolderLevel4,'timestampL4.txt'),'a+');intimelevel4=datetime('now'); % timestamp Level 4 for the starting time
            fprintf(level4timefileid,'Level 4 \n');fprintf(level4timefileid,'start time=%s \n',intimelevel4); incputimelevel4=cputime;fclose(level4timefileid); 
            
            %% EXAMPLE LOOP
            %=================================================================================================================================================
            for egCtr=1:numOfEgs   % 2nd loop for the given no. of examples
            %=====================================================================================================================================================
                
                % set the seed for each example
                if isequal(seedForEg,0)  
                   seedForThisEg=seedToRunSetOfEgs(egCtr);
                else   % i.e. re-running one particular example  
                   seedForThisEg=seedForEg(egCtr); 
                end
                rng(seedForThisEg,'twister'); % give each example data a seed to reproduce the same example later on
    
                   
                %% Define the input data nPts, true b, yArray, xMatrix and epsilon
                % xMatrix
                xMatrix=normalizedRandMatrix(nPts,pDim,infoCol,mu,sigma,seedForThisEg); % define the xMatrix using multivariate normal random dist.
                
                % D_t the block of dependent column in xMatrix
                if nDcol==0 && check_nDcol==1  % if user provided nDcol=0 and check for the dependent block of the xMatrix is on
                    [nDcol,idxnDcol,D_t]=find_nDcol_inData(pDim,nPts,xMatrix); % when user does not know nDcols
                elseif nDcol>0 % if user gave nDcol>0
                    D_t=randi([-5  5],pDim-nDcol,nDcol);  % random integer matrix, where each row will give the coefficients of linear combination for the first pDim-nDcol to create each col of the depBlockMatrix
                    idxnDcol=[zeros(1,pDim-nDcol) ones(1,nDcol)];  % when user gave nDcols, make the last nDcols dependent
                    depBlockxMatrix=xMatrix(:,(idxnDcol==0))*D_t;
                    xMatrix=[xMatrix(:,(idxnDcol==0)) depBlockxMatrix];   % now xMatrix is of order nPts by pDim, with each column for the last nDcol depends on the first pDim-nDcol
                else % nDcol=0 and check_nDcol=0
                    idxnDcol=zeros(1,pDim);D_t=[]; % dummy variable
                end
                
                trueb=giveValuesOf_b(:,chooseEgsToRun(egCtr)); % will use the trueb value defined in the chooseEgsToRun array
                Xtrueb=xMatrix*trueb; 
                sigmaNoise=sqrt( (Xtrueb'*Xtrueb)/snr ); % standard deviation in the noise vector for the linear model
                epsilon=randn(nPts,1)*sigmaNoise; % noise vector
                yArray=xMatrix*trueb+epsilon;
                       
                    
                %% Define the 2nd level outer folder name containing egid
                exFolderLevel3=sprintf('Eg%dp%d_',chooseEgsToRun(egCtr),pDim);
                level5dirpath=fullfile(pwd,setFolderLevel4,char(erase(exFolderLevel3,'_')));
                mkdir(level5dirpath)
                level3timefileid=fopen(fullfile(level5dirpath,'timestampL3.txt'),'a+');intimelevel3=datetime('now');fprintf(level3timefileid,'Level 3 \n');fprintf(level3timefileid,'start time=%s \n',intimelevel3);incputimelevel3=cputime;fclose(level3timefileid); % timestamp Level 3 for the starting time
                fprintf('Ctr=%d;Starting InpFilPair=%d;pDim=%d;Instance=%d;Set=%d;Example=%d;datetime=%s; ... \n',setCtr,InpFilCtr,pDim,nInstCtr,setCtr,egCtr,intimelevel3);
    
                %% Calculate xRelaxedOpt using the quad. min. approach with different versions to choose from
                [~,lb,ub,solsOfQuadMin,cpuxRelaxedOpt,xRelaxedOpt,fxRelaxedOpt,A,b,c]=getxRelaxedOpt(pDim,yArray,xMatrix,trueb,iPara(:,chooseParaToRun(setCtr,4)),rPara,...
                          QuadMinFunPara(1,:), IotherPara(chooseParaToRun(setCtr,1),:) , level5dirpath );  % level5dirpath=[]
                      
                % assign values of eta, lambda and tmax
                alpha=alphaValues(1);
                maxVioRefQMoverTm=-inf; % initialize, to save max vio. QM during refinement for all tmax values
                maxVioInfQMoverTm=-inf; % initialize, to save max vio. QM during Inf F for all tmax values
                preXstarNormSpace=[];  % 16 July21, for penalized regss. with a particular value of alpha, solution for lambda1 can be used as a warm start for lambda2 run where lambda1>=lambda2
                
                %% add the alpha extension to the 2nd level outer folder name
                setNum1=sprintf('_Set%d_AG%d_',setCtr,IotherPara(chooseParaToRun(setCtr,1),15));
                alphaFolderLevel2=append(exFolderLevel3,sprintf('AH%g',alpha),setNum1);
                level2FolderPath=fullfile(pwd,setFolderLevel4,char(erase(exFolderLevel3,'_')),alphaFolderLevel2);
                mkdir(level2FolderPath);
                level2timefileid=fopen( fullfile(level2FolderPath,'timestampL2.txt'),'a+');intimelevel2=datetime('now');incputimelevel2=cputime;fprintf(level2timefileid,'Level 2 \n');fprintf(level2timefileid,'start time=%s \n',intimelevel2);fclose(level2timefileid); % timestamp Level 2 for the starting time
    
                % Initialization
                numOfCols=2+numOfParaVal;
                numOfRows=5+pDim;
                Tdata=zeros(numOfRows,numOfCols);
                header=cell(1,numOfCols);
                % Col. 1: data for true b used in the linear model y=Xb+e
                header{1}=sprintf('trueb(%d)',setCtr); 
                Tdata(:,1)=[trueb'*( (0.5*A)*trueb+b) + c;missing;missing;missing;missing;trueb;];
                % Col. 2: data for OLS solution
                header{2}=sprintf('xRelax(%d)',setCtr);
                Tdata(:,2)=[fxRelaxedOpt ;missing;cpuxRelaxedOpt/60;missing;missing; xRelaxedOpt];
                rowNamesTail=cell(1,pDim);
                for j=1:pDim, rowNamesTail{j}=sprintf('x(%d)',j);end
                rowNames=['fbest','stopFlag','cputime min','niter','totalFcalls',rowNamesTail];


                %% TMAX VALUES LOOP
                %==============================================================================================================================================================================
                for paraCtr=1:numOfParaVal    % 4th loop for all tmax/lambda values
                %==============================================================================================================================================================================    
                    targetfbest=targetfbestValues(paraCtr);  % targetfbest for a particular tmax value 
                    ldTmFolderLevel1=sprintf('Eg0%d0%d_TM%g',pDim,chooseEgsToRun(egCtr),tmaxValues(paraCtr));   
                    level1FolderPath=fullfile(level2FolderPath,ldTmFolderLevel1);
                    mkdir(level1FolderPath)   % to create a separate folder for the example inside the outer folder,  in the current directory
    
                    % save the timestamp 
                    level1timefileid=fopen(fullfile(level1FolderPath,'timestampL1.txt'),'a+');
                    intimelevel1=datetime('now');fprintf(level1timefileid,'Level 1 \n');fprintf(level1timefileid,'start time=%s \n',intimelevel1);incputimelevel1=cputime;fclose(level1timefileid); % timestamp Level 1 for the starting time
                            
                    [header{2+paraCtr},rParaOut,numOfBoxDel,Tdata(:,2+paraCtr),XstarNormSpace,outputPara,fbest]=runAnEgLLS(nPts,pDim,yArray,xMatrix,trueb,A,b,c,tmaxValues(paraCtr),targetfbest,xRelaxedOpt,fxRelaxedOpt,preXstarNormSpace,InpFilCtr,pDimCtr,egCtr,setCtr,paraCtr,...
                                                   IotherPara(chooseParaToRun(setCtr,1),:),IstopCondPara(chooseParaToRun(setCtr,2),:),...
                                    iPara(:,chooseParaToRun(setCtr,4)),rPara,IotherPara(chooseParaToRun(setCtr,1),15),delCondPara(chooseParaToRun(setCtr,5),:),...
                                    ldTmFolderLevel1,level1FolderPath,toDate,saveIntermOutput,toDebug,level8dirpath);
                    
                    preXstarNormSpace=XstarNormSpace;  % to be used as a warm start for the next value of lambda, assuming it is smaller than the lamdba(paraid) used in the above run
                    if fbest<bestfbestoverSet(paraCtr), bestfbestoverSet(paraCtr)=fbest; % update the best fbest over all sets for a particular tmax value
                    end
                    
                    if maxVioRefQMoverTm<outputPara(16), maxVioRefQMoverTm=outputPara(16); end
                    if maxVioInfQMoverTm<outputPara(17),maxVioInfQMoverTm=outputPara(17); end
                    numOfBoxDelAllSets=numOfBoxDelAllSets+numOfBoxDel;
                    
                    level1timefileid=fopen(fullfile(level1FolderPath,'timestampL1.txt'),'a+');
                    outtimelevel1=datetime('now');fprintf(level1timefileid,'End time=%s \n',outtimelevel1);fprintf(level1timefileid,'Time taken for the run=%1.8f min \n',minutes(outtimelevel1-intimelevel1)); % timestamp level0 at the end
                    outcputimelevel1=(cputime-incputimelevel1)/60;fprintf(level1timefileid,'cputime for the run=%1.6f min \n',outcputimelevel1);fprintf(level1timefileid,'fbest for tmax=%d is %1.8f \n',tmaxValues(paraCtr),fbest); 
                    fprintf(level1timefileid,'maxVioRefQM=%1.8f and maxVioInfQM=%1.8f for tmax=%d \n',outputPara(16), outputPara(17),tmaxValues(paraCtr) );
                    fprintf(level1timefileid,'normFactor=');printArray(rParaOut.normFactor,'%1.8f',level1timefileid);fprintf(level1timefileid,'scaleQP=%1.8f \n',rParaOut.scaleQP);
                    fprintf(level1timefileid,'numOfBoxDel=');printArray(numOfBoxDel','%d',level1timefileid);
                    fclose(level1timefileid);
                
                %=============================================================================================================================================================================    
                end  % for TMAX LOOP paraid=1:numOfParaVal
                %=============================================================================================================================================================================
                
                % save the table in the excel file
                saveTableInExcelFile(Tdata,header,rowNames,level5dirpath,setCtr);
                
                        
                %% saveTable1Eg1Set data to create plots in the end                           
         
                level2timefileid=fopen( fullfile(level2FolderPath,'timestampL2.txt'),'a+');
                outtimelevel2=datetime('now');fprintf(level2timefileid,'End time=%s \n',outtimelevel2);fprintf(level2timefileid,'Time taken for the run=%1.8f min \n',minutes(outtimelevel2-intimelevel2)); % timestamp level 2 at the end
                outcputimelevel2=(cputime-incputimelevel2)/60;fprintf(level2timefileid,'cputime for the run=%1.6f min \n',outcputimelevel2);
                fprintf(level2timefileid,'maxVioRefQM=%1.8f and maxVioInfQM=%1.8f for all tmax values. \n',maxVioRefQMoverTm,maxVioInfQMoverTm);
                fclose(level2timefileid);
                
                % Also, save maxdistBdryoverTm in the level 5 timestamp to make it easy to check the violations for all sets
                level5timefileid=fopen('timestampL5.txt','a+');
                fprintf(level5timefileid,'maxVioRefQM=%1.8f and maxVioInfQM=%1.8f for all tmax values Set=%d \n',maxVioRefQMoverTm,maxVioInfQMoverTm,setCtr);fclose(level5timefileid);
                
                saveTableAllEgs1Set{egCtr}=Tdata;
                % if egCtr==1
                %    saveTableAllEgs1Set=saveTable1Eg1Set;
                % else
                %    saveTableAllEgs1Set=[saveTableAllEgs1Set;saveTable1Eg1Set];
                % end
                
                level3timefileid=fopen(fullfile(level5dirpath,'timestampL3.txt'),'a+');    
                outtimelevel3=datetime('now');fprintf(level3timefileid,'End time=%s \n',outtimelevel3);fprintf(level3timefileid,'Time taken for the run=%1.8f min \n',minutes(outtimelevel3-intimelevel3) ); % timestamp level 3 at the end
                outcputimelevel3=(cputime-incputimelevel3)/60;fprintf(level3timefileid,'cputime for the run=%1.6f min \n',outcputimelevel3);fclose(level3timefileid);
                fprintf('Ctr=%d;Ending pDim=%d;Instance=%d;Set=%d;Example=%d;datetime=%s;...;Clocktimetaken=%1.6f min.;cputimetaken=%1.6f min. \n',setCtr,pDim,nInstCtr,setCtr,egCtr,outtimelevel3,minutes(outtimelevel3-intimelevel3),outcputimelevel3);
            %=======================================================================================================================================================================    
            end  % end EXAMPLE LOOP, egid=1:howManyEgs
            %====================================================================================================================================================================== 
            saveTableAllEgsAllSets{setCtr}=cat(2,saveTableAllEgs1Set{:});

            % if setCtr==1   % to save data of same egs with different para(lambda) and different sets
            %    saveTableAllEgsAllSets=saveTableAllEgs1Set;
            % else
            %    saveTableAllEgsAllSets=[saveTableAllEgsAllSets  saveTableAllEgs1Set];
            % end
 
            if maxVioRefQMoverSet<maxVioRefQMoverTm, maxVioRefQMoverSet=maxVioRefQMoverTm; end
            if maxVioInfQMoverSet<maxVioInfQMoverTm, maxVioInfQMoverSet=maxVioInfQMoverTm; end
            
            level4timefileid=fopen(fullfile(pwd,setFolderLevel4,'timestampL4.txt'),'a+');    
            outtimelevel4=datetime('now');fprintf(level4timefileid,'End time=%s \n',outtimelevel4);fprintf(level4timefileid,'Time taken for the run=%1.8f min \n',minutes(outtimelevel4-intimelevel4)); % timestamp level 4 at the end
            outcputimelevel4=(cputime-incputimelevel4)/60;fprintf(level4timefileid,'cputime for the run=%1.6f min \n',outcputimelevel4);fclose(level4timefileid);
             % save the level4 timestamp file info. in level 5 timestamp file, to compare the results more effectively
            level5timefileid=fopen('timestampL5.txt','a+');fprintf(level5timefileid,'%d;%1.8f;%1.8f \n',setCtr,minutes(outtimelevel4-intimelevel4),outcputimelevel4);fclose(level5timefileid);
            
        %===========================================================================================================================================================================    
        end % end SET LOOP, set=1:numOfSets
        %============================================================================================================================================================================
        saveTableAllEgsAllSets=cat(2,saveTableAllEgsAllSets{:});
        
        fileid=fopen('InputParaForThisRun.txt','a+');  % open the file again
        totalWallClockTime=toc(totalClockTime)/60;fprintf(fileid,'total wall clock time = %1.6f min\n',totalWallClockTime);  %fprintf('total wall clock time = %1.6f min\n',totalWallClockTime);
        totalcputime=(cputime-totalCPUtime)/60;fprintf(fileid,'total cputime = %1.6f min\n',totalcputime);       %fprintf('total cputime = %1.6f min\n',totalcputime);
        fclose(fileid);

        %saveListOfSols(numOfEgs,alphaValues,lambdaValues,tmaxValues,numOfSets,IotherPara,chooseParaToRun,saveTableAllEgsAllSets,outerFolderNameArray,numOfRows,rowNamesForInstanceTable,numOfParaInEachSet)
        
        % The file below work only for BSS algorithms when we are running only constrained regss. algo. with
        % tmax values only, not penalized regss. with lambda values
%         save('data.mat','numOfEgs','alphaValues','lambdaValues','tmaxValues','numOfSets','IotherPara','chooseParaToRun','saveTableAllEgsAllSets','outerFolderNameArray','numOfRows','rowNamesForInstanceTable','numOfParaInEachSet')
        
        cpuDataForPlot=zeros(numOfEgs*numOfParaVal,numOfSets);
        fbestDataForPlot=zeros(numOfEgs*numOfParaVal,numOfSets);
        FcallsDataForPlot=zeros(numOfEgs*numOfParaVal,numOfSets);
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
                    fbestDataForPlot( (egj-1)*numOfParaVal +tmaxk, seti )=saveTableAllEgsAllSets( (egj-1)*numOfRows + 1, 2+tmaxk+(seti-1)*(numOfParaVal+2) );   % 1 because fbest is the 1st row in the ouput table
                    stopflagData( (egj-1)*numOfParaVal +tmaxk, seti )=saveTableAllEgsAllSets( (egj-1)*numOfRows + 2, 2+tmaxk+(seti-1)*(numOfParaVal+2) );       % 2 because stopflag is the 2nd row in the ouput table
                    cpuDataForPlot( (egj-1)*numOfParaVal +tmaxk, seti )=saveTableAllEgsAllSets( (egj-1)*numOfRows + 3, 2+tmaxk+(seti-1)*(numOfParaVal+2) );     % 3 because cputime is the 3rd row in the ouput table
                    FcallsDataForPlot( (egj-1)*numOfParaVal +tmaxk, seti )=saveTableAllEgsAllSets( (egj-1)*numOfRows + 5, 2+tmaxk+(seti-1)*(numOfParaVal+2) );  % 5 because totalF is the 5th row in the ouput table
                end
            end    
        end

        out.cpuDataForPlot=cpuDataForPlot;out.fbestDataForPlot=fbestDataForPlot;out.targetfbestValues=targetfbestValues;out.targetcpuValues=targetcpuValues;out.FcallsDataForPlot=FcallsDataForPlot;out.stopflagData=stopflagData;
        
        level5timefileid=fopen('timestampL5.txt','a+');
        outtimelevel5=datetime('now');fprintf(level5timefileid,'End time=%s \n',outtimelevel5);fprintf(level5timefileid,'Time taken for the run=%1.8f min \n',minutes(outtimelevel5-intimelevel5)); % timestamp level 5 at the end
        outcputimelevel5=(cputime-incputimelevel5)/60;fprintf(level5timefileid,'cputime for the run=%1.6f min \n',outcputimelevel5);fprintf(level5timefileid,append('min(fbest) over all sets for tmax =',repmat('%d;',1,length(tmaxValues)),'\n' ),flip(tmaxValues) );fprintf(level5timefileid, append(repmat('%1.8f ',1,length(tmaxValues) ),'\n') ,flip(bestfbestoverSet));out.bestfbestoverSet=flip(bestfbestoverSet);
        fprintf(level5timefileid,'max(distBdry) over all sets with boxed algorithms is %1.8f \n',maxdistBdryoverSet);fprintf(level5timefileid,'maxVioRefQM=%1.8f and maxVioInfQM=%1.8f for all sets.', maxVioRefQMoverSet,maxVioInfQMoverSet);out.maxVioRefQMoverSet=maxVioRefQMoverSet; out.maxVioInfQMoverSet=maxVioInfQMoverSet;
        fclose(level5timefileid);
        cd(level4dirpath);   % change the directory back
        fprintf('\n');
    %=============================================================================================================================================================================================    
    end  % end INSTANCES LOOP,  nInstCtr=1:nInstances
    %===========================================================================================================================================================================================
    
    level6timefileid=fopen(sprintf('timestampL6pDim%d.txt',pDim),'a+');
    outtimelevel6=datetime('now');fprintf(level6timefileid,'End time=%s \n',outtimelevel6);fprintf(level6timefileid,'Time taken for the run=%1.8f min \n',minutes(outtimelevel6-intimelevel6)); % timestamp level 6 at the end
    outcputimelevel6=(cputime-incputimelevel6)/60;fprintf(level6timefileid,'cputime for the run=%1.6f min \n',outcputimelevel6);fprintf(level6timefileid,'numOfBoxDelAllSets=');printArray(numOfBoxDelAllSets','%d',level6timefileid);fclose(level6timefileid);

    fclose('all');
    fprintf('total wall clock time for pDim %d run =%1.6f min. \n',pDim,toc(totalClockTime)/60);
    fprintf('total cputime for pDim %d run=%1.6f min \n',pDim,(cputime-totalCPUtime)/60);
%     fprintf('saving ranking excel files ... \n');
    
    
end % run LLS function

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

function  saveTableInExcelFile(T,header,rowNames,folderPath,set)
% 2July22, save the deepest level table, the #0 in tol. table and the stopFlag reference table

%  8/22/22   str2=append(innerFolderName,sprintf('Set%d-',set),toDate);
    str2=append(folderPath,sprintf('Set%d',set));
    excelFilePath = sprintf('%s.xlsx',str2); 
    Tab=array2table(T,"RowNames",rowNames,"VariableNames",header);
    format longE
    % write the deepest level table
    writetable(Tab,excelFilePath,'WriteRowNames',true,'Range','A1');
    format short

end %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

function struct2vars(s)
    %STRUCT2VARS Extract values from struct fields to workspace variables
    names = fieldnames(s);
    for i = 1:numel(names)
        assignin('caller', names{i}, s.(names{i}));
    end

end %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
