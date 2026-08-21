function finalSolsFornInstances(pDimRun,pDim,nInstances,numOfEgs,tmaxValues,numOfSets,IotherPara,chooseParaToRun,S,outerFolderNameArray,numOfdataRows,rowNames,numOfParaInEachSet,toDate)
% S is the structure which saves the data for all the instances
    % directory name
    folderName=sprintf('np%dInstData_dim%d_%degs_',pDimRun,pDim,numOfEgs);  %folderName=append(sprintf('np%dInstData_dim%d_%degs_',pDimRun,pDim,numOfEgs),toDate);
    mkdir(folderName);folderPath=fullfile(pwd,folderName);
    
    save(fullfile(folderPath,'inputData'),'nInstances','numOfEgs','tmaxValues','numOfSets','IotherPara','chooseParaToRun','S','outerFolderNameArray','numOfdataRows','rowNames','numOfParaInEachSet');
    % get numOfdataCols for each instance 
 
    numOfdataCols=sum(numOfParaInEachSet);

    % get avgData and medianData of each instance
    rowsForEachInstance=numOfEgs*numOfdataRows;
    stackedDataOfAllInstances=zeros(nInstances*rowsForEachInstance,numOfdataCols);   % initialization
    for idx=1:nInstances
        stackedDataOfAllInstances([1+rowsForEachInstance*(idx-1) : rowsForEachInstance*idx] , : )=S(idx).data;   % stack each instance data in one matrix
    end
    avgDataForAllInstances=zeros(rowsForEachInstance,numOfdataCols);
    medDataForAllInstances=zeros(rowsForEachInstance,numOfdataCols);
    dummyVec=zeros(1,nInstances);
    for i=1:rowsForEachInstance
        for j=1:numOfdataCols
           for k=1:nInstances
               dummyVec(k)=stackedDataOfAllInstances( i+rowsForEachInstance*(k-1) , j);
           end
           avgDataForAllInstances(i,j)=mean(dummyVec);
           medDataForAllInstances(i,j)=median(dummyVec);
        end       
    end
    
    
    % get row and column names for stacked data for each instance and write data in a excel file
    
    excelfilepath1=fullfile(folderPath,'avgDataForNInstances.xlsx');
    excelfilepath2=fullfile(folderPath,'medianDataForNInstances.xlsx');
    
    
    headerNamesForInstanceTable=cell(1,sum(numOfParaInEachSet));
    for egctr=1:numOfEgs
       sumPreNumOfParaVal=0; 
       egName=outerFolderNameArray{egctr}; 
       for i=1:numOfSets
           if i>1,sumPreNumOfParaVal=sum(numOfParaInEachSet(1:(i-1)));end
           for j=1:numOfParaInEachSet(i)
              headerNamesForInstanceTable{ j+sumPreNumOfParaVal }=append( erase(egName{1},'_'),sprintf(' Set%d TM%g',i,tmaxValues(j)) );
           end
       end
       avgDataTable=array2table(avgDataForAllInstances( (1+numOfdataRows*(egctr-1):numOfdataRows*egctr),: ),'VariableNames',headerNamesForInstanceTable,'RowNames',rowNames);
       medDataTable=array2table(medDataForAllInstances( (1+numOfdataRows*(egctr-1):numOfdataRows*egctr),: ),'VariableNames',headerNamesForInstanceTable,'RowNames',rowNames);
       writetable(avgDataTable,excelfilepath1,'WriteRowNames',true,'Range',sprintf('A%d',1+(numOfdataRows+2)*(egctr-1) ) );
       writetable(medDataTable,excelfilepath2,'WriteRowNames',true,'Range',sprintf('A%d',1+(numOfdataRows+2)*(egctr-1) ) );
    end
    
    
    
    
    
    
    
%     plots of the avg and median data  (not ready yet)
%     % save plots for average data of each instance     
%     for egCtr=1:numOfEgs
%         numRowsForEachEg=numOfdataRows*length(alphaValues);  % because for each eg. we are saving 11 rows
%         rowCtr=[1+numRowsForEachEg*(egCtr-1):numRowsForEachEg + numRowsForEachEg*(egCtr-1)];
%         %     fileName=outerFolderNameArray{egCtr};tableDataEveryIter=saveTableAllEgsAllSets(rowCtr,:);
%         %     save(append(pwd,'\saveDataToReRun_'),'numOfParaVal','alphaValues','lambdaValues','tmaxValues','numOfSets','IotherPara','chooseParaToRun','tableDataEveryIter','fileName');
%         
%         compareResultsWithPlots11(alphaValues,lambdaValues,tmaxValues,numOfSets,IotherPara,chooseParaToRun,avgDataForAllInstances(rowCtr,:),outerFolderNameArray{egCtr},numOfdataRows);
%     end
%     
%     % save plots for median data of each instance     
%     for egCtr=1:numOfEgs
%         numRowsForEachEg=numOfdataRows*length(alphaValues);  % because for each eg. we are saving 11 rows
%         rowCtr=[1+numRowsForEachEg*(egCtr-1):numRowsForEachEg + numRowsForEachEg*(egCtr-1)];
%         %     fileName=outerFolderNameArray{egCtr};tableDataEveryIter=saveTableAllEgsAllSets(rowCtr,:);
%         %     save(append(pwd,'\saveDataToReRun_'),'numOfParaVal','alphaValues','lambdaValues','tmaxValues','numOfSets','IotherPara','chooseParaToRun','tableDataEveryIter','fileName');
%         
%         compareResultsWithPlots11(alphaValues,lambdaValues,tmaxValues,numOfSets,IotherPara,chooseParaToRun,medDataForAllInstances(rowCtr,:),outerFolderNameArray{egCtr},numOfdataRows);
%     end
%     
    
    
    
end