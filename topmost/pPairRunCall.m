
function  pPairRunCall(InputFiles,paraStruct)
% 12/3/22   

% Implemented par Inp file pair 16July22
% Generic main input file to run input file pair in parllel loop

% delete(gcp('nocreate')); % if there is a pool already running,kill it first
if isempty(gcp('nocreate'))  % if no parpool running, execute the following commands
   parpool('local'); % will start the parallel pool for all the available workers
%parpool('local',numOfWorkers)  will start the parpool on local machine with the given numOfWorkers no. of workers
end
mainfile=paraStruct.mainfile;
toDebug=paraStruct.toDebug;

diary('diary.txt');  
formatOut='yymmmdd';toDate=datestr(now,formatOut); % todate string to be used in the captionformatOut='ddmmmyy';toDate=datestr(now,formatOut); % todate string to be used in the caption
%% Define the Level 1 output directory
% topdirpath=pwd;
% if isfolder(append('output_',toDate))==false
%     mkdir(append('output_',toDate));cd(append('output_',toDate));  % create the output_date folder and then change the current folder to that
% else % if folder aleardy exists, do not overwrite, make a new folder with output(num) data format  
%     num=1;
%     while isfolder(append(sprintf('output(%d)_',num),toDate))==true
%        num=num+1;
%     end
%     mkdir(append(sprintf('output(%d)_',num),toDate));cd(append(sprintf('output(%d)_',num),toDate)); % create the output_date folder and then change the current folder to that
% end
% 
% topfileid=fopen('timestamp.txt','w+');intimetopdir=datetime('now'); % timestamp Level 0 for the starting time
% fprintf(topfileid,'start time=%s \n',intimetopdir);  

startclock=tic;
datetime1=datetime('now');
cputime1=cputime;

nInpFil=size(InputFiles,1); % no. of input files to run 
% to save fbest,cpu,Fcalls data for all the input files combined
fbestDataAllInpFil=cell(1,nInpFil);cpuDataAllInpFil=cell(1,nInpFil);FcallsAllInpFil=cell(1,nInpFil);stopflagAllInpFil=cell(1,nInpFil);
targetfbestAllInpFil=cell(1,nInpFil);targetcpuAllInpFil=cell(1,nInpFil);
% Create a directory to save matlab data files for the option to resume an algorithm
mkdir('saveListData');  % Directory name should be the same as used in each algorithm file which has the option to resume from the previous run

%$$$$$$$$$$$$$$$ End of user defined input $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


eval( append('[s.numOfSets,s.chooseParaToRun,s.IotherPara,s.IstopCondPara,s.delCondPara,s.QuadMinFunPara,s.rPara,s.iPara,s.intermFbestIter,s.intermMaxIter,s.intermCPU,s.algoParafile]=',...
             InputFiles{1,1},';') )
eval( append('[s.npDim,s.pDim_values,s.nPts_values,s.nEgs_eachpDim,s.tmaxVal_eachpDim,s.nDcol_values,s.infoCol_values,s.targetfbest_eachpDim,s.giveCorrParaManually_option,',...
              's.corrflag_option,s.rho_values,s.snr_values,s.n_non_zeros_trueb_values,s.trueb_type_option,s.chooseEgsToRun_eachpDim,s.egParafile,s.toUseRealDataFile]=',...
              InputFiles{1,2},';') )
inputFileNameLocation.mainfile=mainfile;
inputFileNameLocation.egParafile=s.egParafile;
inputFileNameLocation.algoParafile=s.algoParafile;
numOfSetsAllInpFil=zeros(1,nInpFil);numOfSetsAllInpFil(1)=s.numOfSets;
for preInpFilCtr=2:nInpFil
   eval( append('[t.numOfSets,t.chooseParaToRun,t.IotherPara,t.IstopCondPara,t.delCondPara,t.QuadMinFunPara,t.rPara,t.iPara,t.intermFbestIter,t.intermMaxIter,t.intermCPU,t.algoParafile]=',...
             InputFiles{preInpFilCtr,1},';') )
   eval( append('[t.npDim,t.pDim_values,t.nPts_values,t.nEgs_eachpDim,t.tmaxVal_eachpDim,t.nDcol_values,t.infoCol_values,t.targetfbest_eachpDim,t.giveCorrParaManually_option,',...
              't.corrflag_option,t.rho_values,t.snr_values,t.n_non_zeros_trueb_values,t.trueb_type_option,t.chooseEgsToRun_eachpDim,t.egParafile,t.toUseRealDataFile]=',...
              InputFiles{preInpFilCtr,2},';') )
   inputFileNameLocation(preInpFilCtr).mainfile=mainfile;
   inputFileNameLocation(preInpFilCtr).egParafile=t.egParafile;
   inputFileNameLocation(preInpFilCtr).algoParafile=t.algoParafile;
   numOfSetsAllInpFil(preInpFilCtr)=t.numOfSets;
   s(preInpFilCtr)=t;
   clear t
end


                
% 1 :==============================================================================================================================================================
parfor InpFilCtr=1:nInpFil % go over each pair of example and algorithm input file
%=================================================================================================================================================================    
        level8fileid=fopen(sprintf('timestampL8Ctr%d.txt',InpFilCtr),'a+');intimelevel8=datetime('now');incputimelevel8=cputime; % timestamp for the starting time
        fprintf(level8fileid,sprintf('Level 8, Input file pair %d \n',InpFilCtr));fprintf(level8fileid,'start time=%s \n',intimelevel8);fclose(level8fileid); 
       
        wantToRe_runAnEg=0;   % MAKE IT 1, IF WANT TO RE-RUN an example with your current folder as the one in which binary files of the initial run are saved.  
        whichRe_run='0'; % if running the examples for the first time, no need to change. Otherwise replace 0 with the integer representing which re-run it is
        maxVioRefQMoverpDim=-inf; % initialize max vio. of the necess. cond for quad. min sol. during refinement for all pDim values used
        maxVioInfQMoverpDim=-inf; % initialize max vio. of the necess. cond for quad. min sol. during Inf F for all pDim values used
        
        %% Define example and algorithm input file
        numOfSets=s(InpFilCtr).numOfSets;
        chooseParaToRun=s(InpFilCtr).chooseParaToRun;
        IotherPara=s(InpFilCtr).IotherPara;
        IstopCondPara=s(InpFilCtr).IstopCondPara;
        delCondPara=s(InpFilCtr).delCondPara;
        QuadMinFunPara=s(InpFilCtr).QuadMinFunPara;
        rPara=s(InpFilCtr).rPara;
        iPara=s(InpFilCtr).iPara;
        intermFbestIter=s(InpFilCtr).intermFbestIter;
        intermMaxIter=s(InpFilCtr).intermMaxIter;
        intermCPU=s(InpFilCtr).intermCPU;

        npDim=s(InpFilCtr).npDim;
        pDim_values=s(InpFilCtr).pDim_values;
        nPts_values=s(InpFilCtr).nPts_values;
        nEgs_eachpDim=s(InpFilCtr).nEgs_eachpDim;
        tmaxVal_eachpDim=s(InpFilCtr).tmaxVal_eachpDim;
        nDcol_values=s(InpFilCtr).nDcol_values;
        infoCol_values=s(InpFilCtr).infoCol_values;
        targetfbest_eachpDim=s(InpFilCtr).targetfbest_eachpDim; 
        giveCorrParaManually_option=s(InpFilCtr).giveCorrParaManually_option;
        corrflag_option=s(InpFilCtr).corrflag_option;
        rho_values=s(InpFilCtr).rho_values;
        snr_values=s(InpFilCtr).snr_values;
        n_non_zeros_trueb_values=s(InpFilCtr).n_non_zeros_trueb_values;
        trueb_type_option=s(InpFilCtr).trueb_type_option;
        chooseEgsToRun_eachpDim=s(InpFilCtr).chooseEgsToRun_eachpDim;
        toUseRealDataFile=s(InpFilCtr).toUseRealDataFile ;
        
        %% Define the output directory
        level8dirpath=pwd;
%  2Aug23       outDirName=append( InputFiles{InpFilCtr,1},'_',InputFiles{InpFilCtr,2},'_',toDate);
        outDirName=append( InputFiles{InpFilCtr,2},'_',toDate );
        if isfolder( outDirName )==false
            mkdir( outDirName );cd( outDirName );  % create the output_date folder and then change the current folder to that
        else % if folder aleardy exists, do not overwrite, make a new folder with output(num) data format  
            num=0;
            while isfolder( outDirName )==true
               num=num+1;
               outDirName=append( InputFiles{InpFilCtr,1},'_',InputFiles{InpFilCtr,2},'_',toDate,sprintf('(%d)',num) ); 
            end
            mkdir( outDirName );cd( outDirName ); % create the output_date folder and then change the current folder to that
        end

  
        cpuDataForPlot=cell(1,npDim); % variable to gather data for perprof for all tmax, all pDim together
        fbestDataForPlot=cell(1,npDim); % variable to gather data for perprof for all tmax, all pDim together
        targetfbestAllpDim=cell(1,npDim); % variable to gather target fbest as a row for all tmax all pDim
        targetcpuAllpDim=cell(1,npDim); % variable to gather target CPU as a row for all tmax all pDim
        FcallsDataForPlot=cell(1,npDim); % variable to gather Fcalls as a row for all tmax all pDim
        stopflagData=cell(1,npDim); % variable to save stopflags for all tmax all pDim together
        %% Loop to go over different pDim values
        %========================================================================================================================================================================
        for pDimCtr=1:npDim   % outermost loop to run multiple pDim runs for the parameters defined above
        %========================================================================================================================================================================    
                pDim=pDim_values(pDimCtr);    % give the dimension of the predictors
                nPts=nPts_values(pDimCtr);    %10*(pDim+nDcol)   % define the nPts which is the no. of observations
                level7dirpath=pwd;
                level7fileid=fopen(sprintf('timestampL7pDim%d.txt',pDim),'a+');
                intimelevel7=datetime('now'); % timestamp for the starting time
                incputimelevel7=cputime; % save cputime for level8 in timestampL8
                fprintf(level7fileid,'Level 7 \n');fprintf(level7fileid,'start time=%s \n',intimelevel7);fclose(level7fileid);
               
                %Input Part 1. simulation type=feasibilitySamplings under alg 3: -------------------------------------------
                % alpha =-1 and 0, works for IotherPara(15)=1,3 only   
                alphaValues=[0];  % give the alpha values we want to test, always keep the values 2 and 1 in the array while testing.
                %{
                            alpha > 0 use regular penalty p(x)= |x_1|^alpha + ... + |x_pDim|^alpha   where alpha>0
                            alpha = 0 runs exact L_0 penalty  p(x)=0 if x is 0 and p(x)=1 if x is not 0
                            alpha = -1 runs perturbed L_0 penalty, need to define the eps0 for this penalty in IotherPara(22) below
                                            = 0   if |x|<=eps0
                                       p(x) = ( |x|/eps0 - 1 )  if eps0<|x|<2*eps0
                                            = 1   if |x|>=2*eps0 
                %}
                lambdaValues=[10];  % lambda ,the penalty coefficient, in the penalized problem.
                tmaxValues=tmaxVal_eachpDim(pDimCtr,:);  % define tmaxValues in the decreasing order. Length of lambdaValues and tmaxValues should be the same.
                tmaxValues=flip(tmaxValues(tmaxValues~=0)); % remove the 0 flags in the end of the row
                targetfbestValues=targetfbest_eachpDim(pDimCtr,:);
                targetfbestValues=flip(targetfbestValues(tmaxValues~=0));  % discard those targetfbest values for which tmax is 0
    
                % NOTE: the value of tmax that we provide here is already in the reduced space as we take this exact value to be RHS in the eq/ineq constraint.
                % -----when running PENALIZED regression:
                % if alpha>0 i.e. for L_alpha Ridge penalty
                %      tmax= 0  to set tmax= |xRelaxOptNormSpace_1|^alpha + . . . + |xRelaxOptNormSpace_p|^alpha  and use above defined lambda 
                % elseif alpha==0 and -1
                %      tmax= user defined value, should be any positive no.  and use above defined lambda 
                % 
                % ------when running QUADRATIC MIN. wrt EQUALITY or INEQUALITY constraint
                %      tmax= any user defined scalar value,will be the RHS parameter of the eq/ineq constraint
                
                
                
                %Input Part 2. Define the lagorithm-parameters below---------------------------------------------------
%                 toDebug=0;  % =1 print the box flag for AG7, =0 do not print.                                
                                                
                % give data for examples below---------------
%                 toUseRealDataFile={} ; % '0', then generate data synthetically, using the para. given below,
                % if non zero , then the given no. will represent the unique id of the example file name to read from text file.
                nInstances=1;    %$1 number of instances to run 
                seedToRunNInstances=1;   % change 0 to other seed if you want to re-run a particular instance, otherwise no need to alter
                seedToRunSetOfEgs=[0];  % change 0 to other seed if you want to re-run a PARTICULAR set of examples, otherwise no need to alter.
                seedForEg=0 ; % give seed/seeds of an example or a subset of examples we want to re-run from a previous run , else no need to alter.
                
                seedForTrueb=[0]; % seed to be used when we are generating trueb randomly, otherwise will not be used
                 
                giveCorrParaManually=giveCorrParaManually_option(pDimCtr);      % 0 means generate covariance matrix randomly, 1 means generate manually
                nDcol=nDcol_values(pDimCtr);     % give the no. of dependent columns we want to add to xMatrix
                infoCol=infoCol_values(pDimCtr);  % how many informative/correlated columns we want to add to the xMatrix. Only makes the first infoCol no. of cols of xMatrix correlated
        
                % corrflag and rho only get used when generating covariance matrix sigma MANUALLY 
                corrflag=corrflag_option(pDimCtr);    % 1 means constant correlation, 2 means exponential correlation, 3 means take the user defined matrix in the corrMatrix_MeanVecPool 
                rho=rho_values(pDimCtr);    % when corrflag=1, sigma(i,j)=rho when i/=j and sigma(i,i)=1
                          % when corrflag=2, sigma(i,j)=rho^(|i-j|) and sigma(i,i)=1
                          
                 
                % If infoCol is 0 then degOfCorr parameter will not be used. if giveCorrMatManually=1, this parameter will not be used.
                % D_t defines a matrix of order pDim by nDcol with independent columns with integer entries in the range [-2 2], check line 286
                snr=snr_values(pDimCtr);   % signal to noise ratio, will be used to calculate the sd for the noise epsilon
                
                check_nDcol=0;  % if nDcol is non-zero as defined above then check_nDcol is taken to be 1 internally
                %  however, if nDcol=0 and check_nDcol=1 then find the Dependent col columns of xMatrix and use dimension reduction
                numOfEgs=nEgs_eachpDim(pDimCtr);  % define how many examples we want to run, should be same as the length of chooseEgsToRun array below
                n_non_zeros_trueb=n_non_zeros_trueb_values(pDimCtr);  % define the number of 0s in randomly generated trueb
                
                
                trueb_type=trueb_type_option(pDimCtr);           % different options for defining trueb
                
                if trueb_type==-1 % generate the trueb values randomly
                    chooseEgsToRun=(1:numOfEgs);    % If randomly generating trueb keep chooseEgsToRun=[1 2 ... numOfEgs] as a dummy variable
                    [seedForTruebAllEgs,giveValuesOf_b]=randomTrueb(pDim,numOfEgs,pDim-n_non_zeros_trueb,seedForTrueb,seedToRunNInstances);
                    % seedForTheEg will generate another seed of arrays to be used to create random trueb for given no. of examples
                    % to regenerate trueb for an example, provide 
        
                elseif trueb_type==0  % read trueb values from the truebPool.m
                    seedForTruebAllEgs=zeros(numOfEgs,1);  % dummy variable
                    % define the values of true parameter b for all examples, each column represent a value of b for each example.
                    chooseEgsToRun=chooseEgsToRun_eachpDim(pDimCtr,:);    % give which trueb values we want to test from the array below. If randomly generating trueb keep chooseEgsToRun=[1 2 ... numOfEgs]
                    chooseEgsToRun=chooseEgsToRun(chooseEgsToRun~=0);    % remove the dummy 0 entries in the end of the row
                    [truebPoolFilePath,giveValuesOf_b]=truebPool(pDim);
                    % if using trueb values from the pool, save a copy of truebPool
                    getFileName=strsplit( truebPoolFilePath,filesep  );backupFile=fullfile( pwd,append(getFileName{length(getFileName)},'_copy.m') );copyfile(strcat( truebPoolFilePath,'.m' ),backupFile);
        
                elseif ismember(trueb_type,(1:5))  % define trueb using different options
                    chooseEgsToRun=(1:numOfEgs);   % dummy variable
                    seedForTruebAllEgs=zeros(numOfEgs,1);  % dummy variable
                    giveValuesOf_b=zeros(pDim,numOfEgs);
                    for i=1:numOfEgs
                        giveValuesOf_b(:,i)=trueb_type_pool(pDim,numOfEgs,n_non_zeros_trueb,trueb_type);
                    end
        
                end
                
                if giveCorrParaManually==1 && infoCol>0
                    [corrMatrix_MeanVecPoolFilePath,mu,sigma]=corrMatrix_MeanVecPool(infoCol,rho,corrflag);  % defined at the end of this subroutine, set the covariance matrix in that subroutine accordingly
                    % if using corr. matrix and mean vector from corrMatrix_MeanVecPool then save a copy of the file 
                    getFileName=strsplit( corrMatrix_MeanVecPoolFilePath,filesep  );backupFile=fullfile( pwd,append(getFileName{length(getFileName)},'_copy.m') );copyfile(strcat( corrMatrix_MeanVecPoolFilePath,'.m' ),backupFile);
        
                elseif giveCorrParaManually==0 && infoCol>0
                    [~,sigma]=randomCorrMatrix_MeanVec(infoCol); % get mu and sigma randomly, seed=0 has been used to generate mu and sigma 
                     mu=zeros(1,infoCol);
                else % if infoCol=0 
                    sigma=eye(pDim);mu=zeros(1,pDim);
                end
                                                                
                [~]=data_consistency_check(numOfSets,chooseParaToRun,IotherPara,IstopCondPara,QuadMinFunPara,seedToRunNInstances,seedToRunSetOfEgs,pDim,infoCol,nDcol,nPts,alphaValues,lambdaValues,tmaxValues,snr,...
                chooseEgsToRun,numOfEgs,trueb_type,giveValuesOf_b);
%                 if out==2,return;end
                
                saveIntermOutput=struct('FbestIter',intermFbestIter,'MaxIter',intermMaxIter,'CPU',intermCPU);                
                [out]=runLLS(pDim,nPts,infoCol,nDcol,alphaValues,lambdaValues,tmaxValues,nInstances,numOfSets,numOfEgs,seedToRunNInstances,seedToRunSetOfEgs,...
                            seedForEg,chooseEgsToRun,snr,mu,sigma, giveValuesOf_b,seedForTruebAllEgs,check_nDcol,targetfbestValues,chooseParaToRun,IotherPara,IstopCondPara, ...
                            delCondPara,QuadMinFunPara,rPara,iPara,InpFilCtr,pDimCtr,saveIntermOutput,inputFileNameLocation(InpFilCtr),toUseRealDataFile,wantToRe_runAnEg,whichRe_run,toDebug,...
                            toDate,level7dirpath,level8dirpath);

                % to get avg. data using partial and full ranking for each pDim value
                for nInstCtr=1:nInstances
                    pathToRankFiles=fullfile(level7dirpath,out.topDirString,'ListOfSolutions');
                    if nInstCtr==1
                       [leftTablePR,rightTablePR,rightSumDataTablePR]=rankOutputDataPR(pathToRankFiles); 
                       [~,rightTableFR]=rankOutputDataFR(pathToRankFiles);
                       avgDataPR=table2array(rightTablePR);  % 8 and 9 because only displaying 2+6 para.captions as given in rankOutputData11.m line 18,
                       avgDataFR=table2array(rightTableFR); 
                       avgSumDataTablePR=table2array(rightSumDataTablePR);
                    else
                       [~,rightTablePR,rightSumDataTablePR]=rankOutputDataPR(pathToRankFiles);  % data from partial ranking 
                       [~,rightTableFR]=rankOutputDataFR(pathToRankFiles); % data for full ranking
                       avgDataPR = avgDataPR + table2array(rightTablePR); 
                       avgDataFR = avgDataFR + table2array(rightTableFR);
                       avgSumDataTablePR=avgSumDataTablePR+table2array(rightSumDataTablePR);
                    end
                end
%    8/22/22 no need to make a new vars.             avgRanks1pDimTablePR=( avgDataPR./nInstances );
%                 avgRanks1pDimTableFR=( avgDataFR./nInstances );
%                 if pDimCtr>1
%                    avgRanksAllpDimPR=avgRanks1pDimTablePR+avgRanksAllpDimPR;
%                    avgRanksAllpDimFR=avgRanks1pDimTableFR+avgRanksAllpDimFR;
%                 else
%                    avgRanksAllpDimPR=avgRanks1pDimTablePR;
%                    avgRanksAllpDimFR=avgRanks1pDimTableFR;
%                 end
                if pDimCtr>1
                   avgRanksAllpDimPR=avgDataPR./nInstances+avgRanksAllpDimPR;
                   avgRanksAllpDimFR=avgDataFR./nInstances+avgRanksAllpDimFR;
                   avgAllpDimSumDataTablePR=( avgSumDataTablePR./nInstances ) + avgAllpDimSumDataTablePR;
                else
                   avgRanksAllpDimPR=avgDataPR./nInstances;
                   avgRanksAllpDimFR=avgDataFR./nInstances;
                   avgAllpDimSumDataTablePR=( avgSumDataTablePR./nInstances );
                end
                cpuDataForPlot{pDimCtr}=out.cpuDataForPlot;fbestDataForPlot{pDimCtr}=out.fbestDataForPlot;FcallsDataForPlot{pDimCtr}=out.FcallsDataForPlot;stopflagData{pDimCtr}=out.stopflagData;
                targetfbestAllpDim{pDimCtr}=out.targetfbestValues; targetcpuAllpDim{pDimCtr}=out.targetcpuValues;

                if maxVioRefQMoverpDim<out.maxVioRefQMoverSet,maxVioRefQMoverpDim=out.maxVioRefQMoverSet; end
                if maxVioInfQMoverpDim<out.maxVioInfQMoverSet,maxVioInfQMoverpDim=out.maxVioInfQMoverSet; end
        
                level7fileid=fopen(sprintf('timestampL7pDim%d.txt',pDim),'a+');
                outtimelevel7=datetime('now');outcputimelevel7=(cputime-incputimelevel7)/60;fprintf(level7fileid,'start time=%s \n',outtimelevel7);fprintf(level7fileid,'Wall clock time taken for the run = %s \n',outtimelevel7-intimelevel7); % timestamp level0 at the end
                fprintf(level7fileid,'cputime for matlab worker, to run pDim=%d from input file pair %d is %1.8f min\n',pDimCtr,InpFilCtr,outcputimelevel7); fprintf(level7fileid,'maxVioRefQM=%1.8f and maxVioInfQM=%1.8f for pDim=%d \n', out.maxVioRefQMoverSet,out.maxVioInfQMoverSet,pDim_values(pDimCtr));
                fclose(level7fileid);
                
                %save min(fbest) over all sets for given tmax values, to be used in example input files to provide targetfbest
                targetfbestfileid=fopen('targetfbest.txt','a+');
                if pDimCtr==1,fprintf(targetfbestfileid,'targetfbest_eachpDim=[');end
                fprintf(targetfbestfileid, append(repmat('%1.8f ',1,length(tmaxValues) ),'\n') ,out.bestfbestoverSet);
                if pDimCtr==npDim,fprintf(targetfbestfileid,']; \n');end
                fclose(targetfbestfileid);

        %======================================================================================================================================================================
        end % for pDimCtr=1:npDim
        %======================================================================================================================================================================
    
        % save partial and full ranking avg. data files
        warning('off','MATLAB:xlswrite:AddSheet');
        writetable( [leftTablePR , array2table(avgRanksAllpDimPR./npDim,'VariableNames',rightTablePR.Properties.VariableNames) ],'avgRankAllpDimPR.xlsx','Range','A1','AutoFitWidth',false ); % final avg. Partial ranking for each pDim 
        warning('off','MATLAB:xlswrite:AddSheet');
        writetable( [leftTablePR , array2table(avgRanksAllpDimFR./npDim,'VariableNames',rightTableFR.Properties.VariableNames) ],'avgRankAllpDimFR.xlsx','Range','A1','AutoFitWidth',false ); % final avg. Full ranking for each pDim
        warning('off','MATLAB:xlswrite:AddSheet');
        writetable( [leftTablePR , array2table(avgAllpDimSumDataTablePR./npDim,'VariableNames',rightSumDataTablePR.Properties.VariableNames) ],'avgRankAllpDimPR.xlsx','Sheet','rankingData','Range','A1','AutoFitWidth',false ); % sum of the data for all tmax ,all examples, all pDim values

        cd(level8dirpath);  % change to the level 1 folder

        cpuDataForPlot=cat(1,cpuDataForPlot{:});fbestDataForPlot=cat(1,fbestDataForPlot{:});FcallsDataForPlot=cat(1,FcallsDataForPlot{:});stopflagData=cat(1,stopflagData{:});  % concatenate the data from all the cell arrays row wise
        targetfbestAllpDim=cat(2,targetfbestAllpDim{:});targetcpuAllpDim=cat(2,targetcpuAllpDim{:});
        %==================================================================================================================
%         nproblems=length(targetcpuAllpDim); % or length(targetfbestAllpDim) or nrows of cpuDataForPlot or nrows of fbestDataForPlot
%         perprof(cpuDataForPlot,ones(nproblems,numOfSets),'\tau','Prob.','cpu w/o fbest constraint'); % create performance profile graph and save in the current folder
%         saveas(gcf,sprintf('perProCPUinpFil%dwoFbest',InpFilCtr));
%         perprof(fbestDataForPlot,ones(nproblems,numOfSets ),'\tau','Prob.','fx w/o cpu constraint'); % create performance profile graph and save in the current folder
%         saveas(gcf,sprintf('perProFxInpFil%dwoCPU',InpFilCtr));
        %{
        There are 3 blocked places, 1st Blocked on 27June23, to reduce the execution time while debugging
        perprof( cpuDataForPlot,fbestDataForPlot,targetfbestAllpDim,'\tau','Prob.','cpu',[] ); % create performance profile graph with fbest constraint and save in the current folder
        saveas(gcf,sprintf('perProCPUinpFil%d',InpFilCtr));
        perprof( fbestDataForPlot,fbestDataForPlot,targetfbestAllpDim,'\tau','Prob.','fx',[] ); % create performance profile graph with cputime constraint and save in the current folder
        saveas(gcf,sprintf('perProFxInpFil%d',InpFilCtr));
        lineplots(FcallsDataForPlot,'Examples','Fcalls','F calls for different algo.');saveas(gcf,sprintf('lineplotFcallsInpFil%d',InpFilCtr)); % make line plots for Fcalls
        %}
        parforsave(sprintf('dataPerProfInpFil%d.mat',InpFilCtr),cpuDataForPlot,fbestDataForPlot,targetcpuAllpDim,targetfbestAllpDim,FcallsDataForPlot,stopflagData);
        %===================================================================================================================


        level8fileid=fopen(sprintf('timestampL8Ctr%d.txt',InpFilCtr),'a+');
        outtimelevel8=datetime('now');outcputimelevel8=(cputime-incputimelevel8)/60;fprintf(level8fileid,'End time=%s \n',outtimelevel8);fprintf(level8fileid,'Wall clock time taken for the run = %s \n',outtimelevel8-intimelevel8); % timestamp level8 at the end
        fprintf(level8fileid,'cputime for matlab worker, to run input file pair %d is %1.8f min.\n',InpFilCtr,outcputimelevel8); fprintf(level8fileid,'maxVioRefQM=%1.8f and maxVioInfQM=%1.8f for Input file pair=%d over all pDimvalues. \n', maxVioRefQMoverpDim,maxVioInfQMoverpDim,InpFilCtr);
        fclose(level8fileid);

        % update fbest,cpu and Fcalls for all input files
        fbestDataAllInpFil{InpFilCtr}=fbestDataForPlot;cpuDataAllInpFil{InpFilCtr}=cpuDataForPlot;FcallsAllInpFil{InpFilCtr}=FcallsDataForPlot;stopflagAllInpFil{InpFilCtr}=stopflagData;
        targetfbestAllInpFil{InpFilCtr}=targetfbestAllpDim; targetcpuAllInpFil{InpFilCtr}=targetcpuAllpDim;
%end 1:===================================================================================================================================================
end  %  input file pair loop InpFilCtr=1:#input file pairs
%==================================================================================================================================================

% topfileid=fopen('timestamp.txt','a+');
% outtimetopdir=datetime('now');fprintf(topfileid,'start time=%s \n',outtimetopdir);fprintf(topfileid,'Time taken for the run=%s \n',outtimetopdir-intimetopdir); % timestamp level0 at the end
% fclose(topfileid);
% 
% cd(topdirpath); % change directory back to the current
maxNumOfSets=max(numOfSetsAllInpFil);
if isequal(numOfSetsAllInpFil,numOfSetsAllInpFil(1)*ones(1,nInpFil))  % if all the input files have the same no. of sets
    fbestDataAllInpFil=cat(1,fbestDataAllInpFil{:});cpuDataAllInpFil=cat(1,cpuDataAllInpFil{:});FcallsAllInpFil=cat(1,FcallsAllInpFil{:});stopflagAllInpFil=cat(1,stopflagAllInpFil{:}); 
    targetfbestAllInpFil=cat(2,targetfbestAllInpFil{:});targetcpuAllInpFil=cat(2,targetcpuAllInpFil{:});
    %{
    There are 3 blocked places, 1st Blocked on 27June23, to reduce the execution time while debugging
    perprof( cpuDataAllInpFil,fbestDataAllInpFil,targetfbestAllInpFil,'\tau','Prob.','cpu',[] ); % create performance profile graph with fbest constraint and save in the current folder
    saveas(gcf,'perProCPUAllInpFiles');
    perprof( fbestDataAllInpFil,fbestDataAllInpFil,targetfbestAllInpFil,'\tau','Prob.','fx',[] ); % create performance profile graph with cputime constraint and save in the current folder
    saveas(gcf,'perProFxAllInpFiles');
    %}
    perprof2( cpuDataAllInpFil,'\tau','Prob.','cpu all input files',cellfun(@(u) append('Set', num2str(u) ), num2cell(1:maxNumOfSets) ,'UniformOutput',false) ); % create performance profile graph with fbest constraint and save in the current folder
    saveas(gcf,'perProCPUAllInpFiles');
    perprof2( fbestDataAllInpFil,'\tau','Prob.','fx all input files',cellfun(@(u) append('Set', num2str(u) ), num2cell(1:maxNumOfSets) ,'UniformOutput',false) ); % create performance profile graph with cputime constraint and save in the current folder
    saveas(gcf,'perProFxAllInpFiles');
    boxplotfun(fbestDataAllInpFil); saveas(gcf,'boxPlotFxAllInpFiles');
else % if no. of sets are different for all input files

    collectfbest=cell(1,nInpFil);collectcpu=cell(1,nInpFil);collectFcalls=cell(1,nInpFil);collectstopflag=cell(1,nInpFil);
    for iInpFil=1:nInpFil
       [nproblems,~]=size(fbestDataAllInpFil{iInpFil});
       rightfill=nan(nproblems,maxNumOfSets-numOfSetsAllInpFil(iInpFil));
       collectfbest{iInpFil}=cat(2,fbestDataAllInpFil{iInpFil},rightfill);
       collectcpu{iInpFil}=cat(2,cpuDataAllInpFil{iInpFil},rightfill);
       collectFcalls{iInpFil}=cat(2,FcallsAllInpFil{iInpFil},rightfill);
       collectstopflag{iInpFil}=cat(2,stopflagAllInpFil{iInpFil},rightfill);
    end
    fbestDataAllInpFil=cat(1,collectfbest{:});cpuDataAllInpFil=cat(1,collectcpu{:});FcallsAllInpFil=cat(1,collectFcalls{:});stopflagAllInpFil=cat(1,collectstopflag{:}); 
    targetfbestAllInpFil=cat(2,targetfbestAllInpFil{:});targetcpuAllInpFil=cat(2,targetcpuAllInpFil{:});
end
save('dataForAllInpFiles.mat','fbestDataAllInpFil','cpuDataAllInpFil','FcallsAllInpFil','targetfbestAllInpFil','targetcpuAllInpFil','stopflagAllInpFil');

totalTimeAllInputFilePairs=(toc(startclock))/60;
fprintf('Wall clock time for all input file pairs = %1.6f min \n',totalTimeAllInputFilePairs);
clockTextFileName=sprintf('clocktime_%s.txt',toDate);
if isfile( clockTextFileName )==true
    num=1;
    clockTextFileName=sprintf('clocktime_%s(%d).txt',toDate,num);
    while isfolder( clockTextFileName )==true
       num=num+1;
    end
end


totaltimefileid=fopen( clockTextFileName,'w+' );  % to write the total time for the run in a text file
fprintf(totaltimefileid,'total wallclocktime to run all the input file pairs = %1.8f min. \n',totalTimeAllInputFilePairs);
datetime2=datetime('now');
fprintf(totaltimefileid,'start time=%s \n', datetime1);
fprintf(totaltimefileid,'end time=%s \n', datetime2);
fprintf(totaltimefileid,'Time increment = %s \n', datetime2-datetime1);
totalTime=(cputime- cputime1)/60;
fprintf(totaltimefileid,'total cputime = %1.8f min. \n',totalTime);
fclose(totaltimefileid);
diary off   

% pairExSumm(fbestDataAllInpFil,cpuDataAllInpFil,InputFiles); % create the avg. fbest and cpu data for all input files combined

end % end of runParInpFilPairLoop

function parforsave(filename,cpuDataForPlot,fbestDataForPlot,targetcpuAllpDim,targetfbestAllpDim,FcallsDataForPlot,stopflagData)
% a remedy for the fact that we cannot call save function inside a parfor loop    
   save(filename,'cpuDataForPlot','fbestDataForPlot','targetcpuAllpDim','targetfbestAllpDim','FcallsDataForPlot','stopflagData');
end



