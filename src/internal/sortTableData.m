%{
 01 Jan20, to sort the table data in an increasing order of number of zeros in the solution for Table Output.
For solutions with same no. of zeros sort such as the lower training error column comes first among cols with
same no. of zeros.
%}
function [colOrder]=sortTableData(summaryTableData,saveThirdRow)

  
  [~,numOfCols]=size(summaryTableData);
  colOrder=(1:numOfCols);
  [numOf0sWithoutRepetitions]=unique(saveThirdRow);
 
  for i=1:length(numOf0sWithoutRepetitions)
      [index]=find(saveThirdRow==numOf0sWithoutRepetitions(i));
      firstRow=summaryTableData(1,index);
      [~,sortingIndex]=sort(firstRow);
      colOrder(index)=index(sortingIndex);
  end





end