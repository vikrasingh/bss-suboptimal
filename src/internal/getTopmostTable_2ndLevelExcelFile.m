% 8th June21, the subroutine to get the topmost table of the second level excel file, by manipulating the data
function [topMostTable,summaryTableData,topTableDataLevel2ExcelFile,countEntries,rowid1,...
                                rowNamesForTopMostTable]=getTopmostTable_2ndLevelExcelFile(summaryTable1,numOfParaVal,set,...
                                                                          tmaxValues)

     summaryTableData=(table2array(summaryTable1))'; % to form a special table on the top of the 2nd level excel files
     
     countEntries=3;
     
     topTableDataLevel2ExcelFile=zeros(countEntries,3*numOfParaVal); % get the data for the special table at the top of 2nd level excel  file
     rowid1=zeros(1,numOfParaVal);rowid2=zeros(1,numOfParaVal); % row1d1,rowid2 indices to be used to extract the rows
     for idx1=1:numOfParaVal   % extract the first 3 columns fbest,t.e.,g.e. for different lambda/tmax values and put them sidewise for all given alpha values.
         rowid1(idx1)=1+countEntries*(idx1-1);rowid2(idx1)=countEntries+countEntries*(idx1-1);
         topTableDataLevel2ExcelFile(:,(1+3*(idx1-1):3+3*(idx1-1)))=summaryTableData((rowid1(idx1):rowid2(idx1)),(1:3));       
     end
     rowNamesForTopMostTable=cell(countEntries,1); % row names for the special table
     rowNamesForTopMostTable{1}='trueb';rowNamesForTopMostTable{2}=sprintf('x_ols(%d)',set);rowNamesForTopMostTable{3}='C.Regss.'; 
     
    topMostTableHeader=cell(1,3);
    for idx4=1:numOfParaVal
        topMostTableHeader{1+3*(idx4-1)}=sprintf('Exact L_0 TM%g fbest',tmaxValues(idx4));topMostTableHeader{2+3*(idx4-1)}=sprintf('Exact L_0 TM%g t.e.',tmaxValues(idx4));topMostTableHeader{3+3*(idx4-1)}=sprintf('Exact L_0 TM%g g.e.',tmaxValues(idx4));    
    end
    
     topMostTable=array2table(topTableDataLevel2ExcelFile,'VariableNames',topMostTableHeader,'RowNames',rowNamesForTopMostTable);

     









end