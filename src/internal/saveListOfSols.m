function saveListOfSols(numOfEgs,alphaValues,lambdaValues,tmaxValues,numOfSets,IotherPara,chooseParaToRun,D,outerFolderNameArray,numofdatarows,rowNames,numOfParaInEachSet)
    %  13 July21   function to make special excel files which contains list of solution for 1 eg and all tmax/lambda values altogether
    
    % structure of D: say 2egs, 3 set with AG3,5,1 resp. 2 tmaxvalues, 3 lambdavalues,2 alphavalues
    %                  Set 1                        Set 2                         Set 3
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
    save(fullfile(folderPath,'inputData'),'numOfEgs','alphaValues','lambdaValues','tmaxValues','numOfSets','IotherPara','chooseParaToRun','D','outerFolderNameArray','numofdatarows','rowNames','numOfParaInEachSet');
    numOfParaInEachSetWithtrueb=numOfParaInEachSet+1+1;
    
    %===================================================================================================================
    %% save data for every example in the excel file
    %===================================================================================================================
    for egctr=1:numOfEgs
        
        egName=outerFolderNameArray{egctr};egName=erase(egName{1},'_');
        extractAlgoCtr=cell(1,numOfSets);
        
        % loop to go every set
        for setctr=1:numOfSets          
            
            % D_1Eg1Set_Ah contains output including trueb columns
            D_1Eg1Set_Ah=D(( 1+numofdatarows*length(alphaValues)*(egctr-1) : numofdatarows*length(alphaValues)*egctr ), ( 1+sum( numOfParaInEachSetWithtrueb(1:(setctr-1)) ) : sum(numOfParaInEachSetWithtrueb(1:setctr)) )); 
            
            extractAlgoCtr{setctr}=sprintf('%d',IotherPara(chooseParaToRun(setctr,1),15));
            idxoftruebvalues=zeros(1,length(alphaValues)); % in the matlab table
            if abs( IotherPara( chooseParaToRun(setctr,1),15 ) )==3 || ismember( IotherPara(chooseParaToRun(setctr,1),15),[7 8 9 71 72 73 74 20 21 22 23 24 25 30 31 32] )
               excelFileName=append( egName,sprintf('SolListSet%dAh%dTm%dAG%d.xlsx',setctr,length(alphaValues),length(tmaxValues),IotherPara(chooseParaToRun(setctr,1),15) ) ); 
            else
               excelFileName=append( egName,sprintf('SolListSet%dAh%dLd%dAG%d.xlsx',setctr,length(alphaValues),length(lambdaValues),IotherPara(chooseParaToRun(setctr,1),15) ) );
            end
            numofcols=(1+1+numOfParaInEachSet(setctr))*length(alphaValues);
            headerNamesFor1Eg1Set_Ah=cell( 1, numofcols );
            dummyTable1Eg1Set_Ah=zeros( numofdatarows, numofcols );
            dummyheader=cell(1,1+1+numOfParaInEachSet(setctr));
            for alphactr=1:length(alphaValues)
                dummyTable1Eg1Set_Ah(:, ( 1+(numOfParaInEachSet(setctr)+1+1 )*(alphactr-1) :( numOfParaInEachSet(setctr)+1+1 )*alphactr ) )=D_1Eg1Set_Ah( (1+numofdatarows*(alphactr-1):numofdatarows*alphactr),(1:(1+1+numOfParaInEachSet(setctr))) );
                dummyheader{1}=sprintf('trueb Ah%g Set%d',alphaValues(alphactr),setctr);
                dummyheader{2}=sprintf('xRelaxOpt Ah%g Set%d',alphaValues(alphactr),setctr);
                for paractr=1:numOfParaInEachSet(setctr)
                    if abs(IotherPara(chooseParaToRun(setctr,1),15))==3 || ismember( IotherPara(chooseParaToRun(setctr,1),15),[7 8 9 71 72 73 74 20 21 22 23 24 25 30 31 32] )
                       dummyheader{1+1+paractr}=sprintf('Ah%g Set%d Tm%g',alphaValues(alphactr),setctr,tmaxValues(paractr) );
                    else
                       dummyheader{1+1+paractr}=sprintf('Ah%g Set%d Ld%g',alphaValues(alphactr),setctr,lambdaValues(paractr) );
                    end
                end
                headerNamesFor1Eg1Set_Ah( 1+(1+1+numOfParaInEachSet(setctr))*(alphactr-1):(1+1+numOfParaInEachSet(setctr))*alphactr )=dummyheader;
            end
                
            row1ExcelFile=cell(1,1+ numofcols );row1ExcelFile(:)={egName};  % egid for  1st col (row names) + 2nd col (trueb) + remaining cols(interval sols) + last 5 cols
            row2ExcelFile=cell(1,1+ numofcols );
            row2ExcelFile(:)={sprintf('AG%d',IotherPara(chooseParaToRun(setctr,1),15))};row2ExcelFile{1}=[];
                
                
            for excelcolctr=1:length(alphaValues)
                idxoftruebvalues(excelcolctr)=1+(1+1+numOfParaInEachSet(setctr))*(excelcolctr-1);
                row2ExcelFile{ 2+(1+1+numOfParaInEachSet(setctr))*(excelcolctr-1) }=sprintf('Set%d trueb',setctr);
            end



            excelfilepath1=fullfile(folderPath,excelFileName );
            Table=array2table(dummyTable1Eg1Set_Ah,'VariableNames',headerNamesFor1Eg1Set_Ah,'RowNames',rowNames);
            [alphabetCells]=letters;  % will give us the excel file column names, atmost 260 cols can be handled
            
            % we want to give a threshold of 1% to a number of parameters as default
            % use rows (1)fbest,(3)general err., (4)pred. err.,(12)cputime in min.,(13)#iter, (25)#incl.fun.call, (26)#feval,
            % (28) lastfbestUpdateIter,(33)cpu refinement,(34) cpu lb F
            thresCol=nan(numofdatarows,1);thresCol([1 3 4 12 13 25 26 28 33 34])=1; % default value for some para. in threshold %
            
            modiTable=[array2table(thresCol,'VariableNames',{'Rank.Thres.%'} )  Table];
            writetable(modiTable,excelfilepath1,'WriteRowNames',true,'Range','A4' );
            writecell(row1ExcelFile,excelfilepath1,'Range','B2');  % write 1st row at the TOP containing egname in the excel file
            writecell(row1ExcelFile,excelfilepath1,'Range',sprintf('B%d',numofdatarows+6) ); % write 1st row at the BOTTOM of the table as well
            writecell(row2ExcelFile,excelfilepath1,'Range','B3');  % write 2nd row at the TOP of the table
            writecell(row2ExcelFile,excelfilepath1,'Range',sprintf('B%d',numofdatarows+7)); % write 2nd row at the BOTTOM of the table as well
            
            % write row to write the flag of 0 and 1 to tell which column to include for ranking
            % default is keep 0 for trueb and xRelaxedOpt and 1 for all the other columns.
            writematrix([0 0 ones(1,numofcols-2)],excelfilepath1,'Range','C1');
             
            
            % add a special sign $$ to the "first row to the left of the last column" and "first col. at the
            % bottom of the last row" in the table
            writecell({'$$'},excelfilepath1,'Range',sprintf('A%d',5+numofdatarows));
            writecell({'$$'},excelfilepath1,'Range',append(alphabetCells{3+numofcols},'1'));
            
            % the following saved matlab data file will be used to rank the cols after adding ranking threshold to the excel file  
            save( fullfile(folderPath,sprintf('rankingdataeg%dset%d',egctr,setctr) ),'excelFileName','Table','headerNamesFor1Eg1Set_Ah','rowNames','idxoftruebvalues',...
                           'numOfParaInEachSet','IotherPara','chooseParaToRun','excelfilepath1','numOfEgs','numOfSets','numofdatarows','numofcols',...
                           'alphabetCells','row1ExcelFile','row2ExcelFile','folderPath','egName'  );

            
            addColorsToExcelSheet(excelfilepath1,IotherPara(chooseParaToRun(setctr,1),15),alphaValues,numofdatarows,numOfParaInEachSet,setctr,alphabetCells,idxoftruebvalues);
            

            
        end  % end setctr=1:numOfSets
        
       
    %===================================================================================================================    
    end % egctr=1:numOfEgs
    %===================================================================================================================
    
    
    
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%--------------------------------------------------------------------------------------------

function addColorsToExcelSheet(excelfilepath1,version_flag,alphaValues,numOfdataRows,numOfParaInEachSet,setctr,alphabetCells,idxoftruebvalues)

    
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
    if version_flag==3 
       Wb.Worksheets.Item(1).Range( sprintf('C1:C%d',numOfdataRows+4) ).Interior.ColorIndex = 34;  % trueb columns color light cyan 
    else
       Wb.Worksheets.Item(1).Range( sprintf('%s1:%s%d', alphabetCells{2+idxoftruebvalues(1)}, alphabetCells{2+idxoftruebvalues(1)},numOfdataRows+4 ) ).Interior.ColorIndex = 8;  % trueb columns color cyan  
       for i=2:length(idxoftruebvalues)
           Wb.Worksheets.Item(1).Range( sprintf('%s1:%s%d', alphabetCells{2+idxoftruebvalues(i)} ,  alphabetCells{2+idxoftruebvalues(i)} , numOfdataRows+4 ) ).Interior.ColorIndex = 8;  % trueb columns color cyan 
       end
    end
    
    % color each alpha value output differently
    for alphactr=1:length(alphaValues)
        if version_flag==3
           cellrangealphavalueoutput=sprintf('%s1:%s%d', alphabetCells{ 5+( numOfParaInEachSet(setctr) )*(alphactr-1)},alphabetCells{ 5+( numOfParaInEachSet(setctr) )-1+( numOfParaInEachSet(setctr) )*(alphactr-1)}, numOfdataRows+4 );
        else
           cellrangealphavalueoutput=sprintf('%s1:%s%d', alphabetCells{ 5+( 1+1+numOfParaInEachSet(setctr) )*(alphactr-1)},alphabetCells{ 5+(numOfParaInEachSet(setctr))-1+( 1+1+numOfParaInEachSet(setctr) )*(alphactr-1)}, numOfdataRows+4 );
        end
        if rem(alphactr,2)==1  % alternate colors for each alpha value output
           Wb.Worksheets.Item(1).Range( cellrangealphavalueoutput ).Interior.ColorIndex = 35;  % light green 
        else
           Wb.Worksheets.Item(1).Range( cellrangealphavalueoutput ).Interior.ColorIndex = 37;  % light blue 
        end
    end
    
    Wb.Save();  % Save Workbook
    Wb.Close(); % Close Workbook
    Ex.Quit(); % Quit Excel
    %       Ex=actxGetRunningServer('Excel.Application');Ex.Quit   % to get the currently running excelfile  and close it       
    %  system('taskkill /F /IM EXCEL.EXE');  % to kill the open excel file
   
    

end


%%%%%%%%%%%%%%%%%%---------------------------------------------------------------------------------------------------------------------------





