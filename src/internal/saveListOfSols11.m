function saveListOfSols11(pDimCtr,numOfEgs,alphaValues,lambdaValues,tmaxValues,numOfSets,IotherPara,QuadMinFunPara,chooseParaToRun,D,outerFolderNameArray,numofdatarows,rowNames,numOfParaInEachSet)
    %  26 May22   function to make special excel files which contains list of solutions for 1 eg and 1 tmax and all sets together
    % After running this file, run rankOutputDataPR(pwd) and rankOutputDataFR(pwd) in the command window, to rank all the sets with all the tmax values in one Egsummary.xlsx file. 

    % structure of D: say 2egs, 3 set with AG3,5,1 resp. 2 tmaxvalues, 3 lambdavalues,2 alphavalues
    %                  Set 1                                       Set 2                                     Set 3
    %       |     trueb  xRelaxOpt  Tm1  Tm2  |     trueb  xRelaxOpt  Ld1   Ld2    Ld3   |     trueb  xRelaxOpt  Ld1    Ld2   Ld3
    %  Eg1  | Ah1                             | Ah1                                      | Ah1
    %       | Ah2                             | Ah2                                      | Ah2
    %        --------------------------------------------------------------------------------------------------------------------------
    %       |     trueb  xRelaxOpt  Tm1  Tm2  |     trueb  xRelaxOpt  Ld1   Ld2    Ld3   |     trueb  xRelaxOpt  Ld1    Ld2   Ld3
    %  Eg2  | Ah1                             | Ah1                                      | Ah1
    %       | Ah2                             | Ah2                                      | Ah2
    %
    
    
    
    
    % directory name
    folderName='ListOfSolutions';mkdir('ListOfSolutions');folderPath=fullfile(pwd,folderName);
    
    % save the input of the subroutine
%     save(append(folderPath,'\inputData'),'numOfEgs','alphaValues','lambdaValues','tmaxValues','numOfSets','IotherPara','chooseParaToRun','D','outerFolderNameArray','numofdatarows','rowNames','numOfParaInEachSet');
    
    %===================================================================================================================
    %% save data for every example in the excel file
    %===================================================================================================================
    for egctr=1:numOfEgs
        
        egName=outerFolderNameArray{egctr};egName=erase(egName{1},'_');
        headerNamesFor1Eg1Tm=cell( 1, 1+1+numOfSets );
        headerNamesFor1Eg1Tm{2}='xRelaxOpt';headerNamesFor1Eg1Tm{1}='trueb';

        for tmaxctr=1:length(tmaxValues)

            extractAlgoCtr=cell(1,numOfSets);
            numofcols=1+1+numOfSets*length(alphaValues);
            row2ExcelFile=cell(1,1+ numofcols );
            % loop to go every set for each tmax value together in an excel file SolListTmCtr%dFR.xlsx and SolListTmCtr%dPR.xlsx
            for setctr=1:numOfSets  
        
                if setctr==1
                   % D_1Eg1TmAllSets contains output including trueb columns, for the first set 
                   D_1Eg1TmAllSets=D(( 1+numofdatarows*length(alphaValues)*(egctr-1) : numofdatarows*length(alphaValues)*egctr ) ,   [1 2  2+tmaxctr]  ); 
                else
                   D_1Eg1TmAllSets=[D_1Eg1TmAllSets  D(( 1+numofdatarows*length(alphaValues)*(egctr-1) : numofdatarows*length(alphaValues)*egctr ), ( 2 + (2+length(tmaxValues))*(setctr-1)+tmaxctr ))]; 
                end
                
                headerNamesFor1Eg1Tm{2+setctr}=sprintf('Tm=%g Set%d AG%d',tmaxValues(tmaxctr),setctr,IotherPara(chooseParaToRun(setctr,1),15));
                extractAlgoCtr{setctr}=sprintf('%d',IotherPara(chooseParaToRun(setctr,1),15));
                excelFileNamePR=append( egName,sprintf('SolListTmCtr%dPR.xlsx',tmaxctr ) ); 
                excelFileNameFR=append( egName,sprintf('SolListTmCtr%dFR.xlsx',tmaxctr ) );
                idxoftruebvalues=zeros(1,length(alphaValues)); % in the matlab table
%                 row2ExcelFile(1+setctr)={sprintf('Tm%g Set%d AG%d',tmaxValues(tmaxctr),setctr,IotherPara(chooseParaToRun(setctr,1),15))};

            end  % end setctr=1:numOfSets
             

            row1ExcelFile=cell(1,1+ numofcols );row1ExcelFile(:)={egName};  % egid for  1st col (row names) + 2nd col (trueb) + remaining cols(interval sols) + last 5 cols

            excelfilepathPR=fullfile(folderPath,excelFileNamePR );
            Table=array2table(D_1Eg1TmAllSets,'VariableNames',headerNamesFor1Eg1Tm,'RowNames',rowNames);
            [alphabetCells]=letters;  % will give us the excel file column names, atmost 260 cols can be handled
            
            %******************* we want to give a threshold of 1% to a number of parameters as default
            % use rows (1)fbest,(3)general err., (4)pred. err.,(11)stop flag,(12)cputime in min.,(13)#iter, (25)#incl.fun.call, (26)#feval,
            % (30) lastfbestUpdateIter,(35)cpu refinement,(36) cpu lb F
            thresCol=nan(numofdatarows,1);thresCol([1 3 4 12 13 25 26 30 35 36])=1;thresCol(11)=0; % default value for some para. in threshold %
            
            modiTable=[array2table(thresCol,'VariableNames',{'Rank.Thres.%'} )  Table];
            writetable(modiTable,excelfilepathPR,'WriteRowNames',true,'Range','A4' );
            writecell(row1ExcelFile,excelfilepathPR,'Range','B2');  % write 1st row at the TOP containing egname in the excel file
            writecell(row1ExcelFile,excelfilepathPR,'Range',sprintf('B%d',numofdatarows+6) ); % write 1st row at the BOTTOM of the table as well
%                 writecell(row2ExcelFile,excelfilepath1,'Range','B3');  % write 2nd row at the TOP of the table
%                 writecell(row2ExcelFile,excelfilepath1,'Range',sprintf('B%d',numofdatarows+7)); % write 2nd row at the BOTTOM of the table as well
            writecell(Table.Properties.VariableNames,excelfilepathPR,'Range','C3');  % write 2nd row at the TOP of the table
            writecell(Table.Properties.VariableNames,excelfilepathPR,'Range',sprintf('C%d',numofdatarows+7)); % write 2nd row at the BOTTOM of the table as well
            
            % write row to write the flag of 0 and 1 to tell which column to include for ranking
            % default is keep 0 for trueb and xRelaxedOpt and 1 for all the other columns.
            writematrix([0 0 ones(1,numofcols-2)],excelfilepathPR,'Range','C1');
             
            % make a copy of the excel file with FR 
            copyDataCellArray=readcell(excelfilepathPR);  % read the data in the excel file again
            idxCharCells=cellfun('isclass',copyDataCellArray,'char');
            tempDataCellArray=copyDataCellArray;tempDataCellArray(idxCharCells)={0}; % change all the cells with char. to numeric 0 temperorary
            
            idxMissingValues=cellfun(@ismissing,tempDataCellArray); % logical index of the missing values
            copyDataCellArray(idxMissingValues)={NaN}; % change missing entries to NaN and then write in the excel file
            excelfilepathFR=fullfile(folderPath,excelFileNameFR );
            writecell(copyDataCellArray,excelfilepathFR, 'Range','A1' ); % make a copy of the excel file and rename as FR instead of PR


            % add a special sign $$ to the "first row to the left of the last column" and "first col. at the
            % bottom of the last row" in the table
            writecell({'$$'},excelfilepathPR,'Range',sprintf('A%d',5+numofdatarows));
            writecell({'$$'},excelfilepathPR,'Range',append(alphabetCells{3+numofcols},'1'));
            writecell({'$$'},excelfilepathFR,'Range',sprintf('A%d',5+numofdatarows));
            writecell({'$$'},excelfilepathFR,'Range',append(alphabetCells{3+numofcols},'1'));
            
            
            % the following saved matlab data file will be used to rank the cols after adding ranking threshold to the excel file  
            save( fullfile(folderPath,sprintf('DataEg%dTm%d',egctr,tmaxctr) ),'excelFileNamePR','excelFileNameFR','Table','headerNamesFor1Eg1Tm','rowNames','idxoftruebvalues',...
                           'tmaxValues','numOfParaInEachSet','QuadMinFunPara','IotherPara','chooseParaToRun','excelfilepathPR','pDimCtr','numOfEgs','numOfSets','numofdatarows','numofcols',...
                           'alphabetCells','row1ExcelFile','row2ExcelFile','folderPath','egName'  );

            
%                 addColorsToExcelSheet(excelfilepathPR,IotherPara(chooseParaToRun(setctr,1),15),alphaValues,numofdatarows,numOfSets,setctr,alphabetCells,idxoftruebvalues);
                

             
        end % end tmaxctr=1:numOfParaInEachSet
       
    %===================================================================================================================    
    end % egctr=1:numOfEgs
    %===================================================================================================================
    
    
    
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%--------------------------------------------------------------------------------------------

function addColorsToExcelSheet(excelfilepath1,version_flag,alphaValues,numOfdataRows,numOfSets,setctr,alphabetCells,idxoftruebvalues)

    
    Ex=actxserver('Excel.application');   % Connect to Excel
    Wb=Ex.Workbooks.Open(excelfilepath1,0,false);  % Get Workbook object
    
% block of code below was for coloring the header of last 5 columns, in the summaryDataFor1Eg1set_Ah   
%     % get range of the cells
%     if version_flag==3
%        cellrange5bestParaHeader=sprintf('%s4:%s4', alphabetCells{3+(numOfParaInEachSet(setctr))*length(alphaValues)+1},alphabetCells{3+(numOfParaInEachSet(setctr))*length(alphaValues)+5} );
%     else
%        cellrange5bestParaHeader=sprintf( '%s4:%s4', alphabetCells{2+(1+numOfParaInEachSet(setctr))*length(alphaValues)+1},alphabetCells{2+(1+numOfParaInEachSet(setctr))*length(alphaValues)+5} ); 
%     end
%     Wb.Worksheets.Item(1).Range( cellrange5bestParaHeader ).Interior.ColorIndex = 6; % Set the color of the 5 Best parameter cells of Sheet 1 to yellow
    
    % color the trueb columns separately
%     if version_flag==3 
       Wb.Worksheets.Item(1).Range( sprintf('C1:C%d',numOfdataRows+4) ).Interior.ColorIndex = 34;  % trueb columns color light cyan 
%     else
%        Wb.Worksheets.Item(1).Range( sprintf('%s1:%s%d', alphabetCells{2+idxoftruebvalues(1)}, alphabetCells{2+idxoftruebvalues(1)},numOfdataRows+4 ) ).Interior.ColorIndex = 8;  % trueb columns color cyan  
%        for i=2:length(idxoftruebvalues)
%            Wb.Worksheets.Item(1).Range( sprintf('%s1:%s%d', alphabetCells{2+idxoftruebvalues(i)} ,  alphabetCells{2+idxoftruebvalues(i)} , numOfdataRows+4 ) ).Interior.ColorIndex = 8;  % trueb columns color cyan 
%        end
%     end
    
  
        
    cellrangealphavalueoutput=sprintf('E1:%s%d', alphabetCells{ 5+numOfSets-1}, numOfdataRows+4 );    
    Wb.Worksheets.Item(1).Range( cellrangealphavalueoutput ).Interior.ColorIndex = 35;  % light green 
       
   
    
    Wb.Save();  % Save Workbook
    Wb.Close(); % Close Workbook
    Ex.Quit(); % Quit Excel
    %       Ex=actxGetRunningServer('Excel.Application');Ex.Quit   % to get the currently running excelfile  and close it       
    %  system('taskkill /F /IM EXCEL.EXE');  % to kill the open excel file
   
    

end


%%%%%%%%%%%%%%%%%%---------------------------------------------------------------------------------------------------------------------------





