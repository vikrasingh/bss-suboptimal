%{
      2 Feb21
      Ranking in columns B,C and J is straightforward, for columns D:J there are two values i,j assigned as rank
      >> i represents rank wrt no. of zeros(more zeros lower the rank) in the above tol column, for solutions with same no. of zeros tie
      breaker is train error(less train. error lower the rank).
      >> j represents rank wrt no. of zeros(more zeros lower the rank) in the above tol column, for solutions with same no. of zeros tie
      breaker is general error(less general error. lower the rank).
    
      Ranking criteria for D:J columns:  1st value = rank using # of 0s, tie breaker t.e.
                                         2nd value = rank using # of 0s, tie breaker g.e.
 %}


function   saveRankInfo_2ndLevelExcelFile(dataToRank,numOfLassoSol,numOfRidgeSol,saveAlphaValues,newRowNames,excelFilePath1)


    nrows=3+length(saveAlphaValues)+numOfLassoSol+numOfRidgeSol;lowerPartData_left=zeros(nrows,3);lowerPartData_middle=cell(nrows,6);lowerPartData_right=zeros(nrows,1);indexVec=(1:nrows);
    [c1,idx]=sort(dataToRank(:,1));c1(idx)=indexVec;lowerPartData_left(:,1)=c1;
    [c2,idx]=sort(dataToRank(:,2));c2(idx)=indexVec;lowerPartData_left(:,2)=c2;
    [c3,idx]=sort(dataToRank(:,3));c3(idx)=indexVec;lowerPartData_left(:,3)=c3;
    [c10,idx]=sort(dataToRank(:,10));c10(idx)=indexVec;lowerPartData_right(:,1)=c10;   
    for i=1:6
        vec=cell(nrows,1);
        saveCount=0;
        for j=(max(dataToRank(:,3+i))+1):-1:1
            [sameNumOf0sIndices]=find(dataToRank(:,3+i)==(j-1));
            [index_te]=c2(sameNumOf0sIndices);[vec2,rankIndex_te]=sort(index_te);vec2(rankIndex_te)=(1:length(sameNumOf0sIndices));vec2=saveCount+vec2;
            [index_ge]=c3(sameNumOf0sIndices);[vec3,rankIndex_ge]=sort(index_ge);vec3(rankIndex_ge)=(1:length(sameNumOf0sIndices));vec3=saveCount+vec3;
            
            
            for k=1:length(sameNumOf0sIndices)
               vec{sameNumOf0sIndices(k)}=sprintf('%d,%d',vec2(k),vec3(k));
            end
            saveCount=length(sameNumOf0sIndices)+saveCount;
        end
        lowerPartData_middle(:,i)=vec;
    end
    lowerPartTable_left=array2table(lowerPartData_left,'VariableNames',{'rank using fbest','rank using t.e.','rank using g.e.'},'RowNames',newRowNames);
    lowerPartTable_middle=cell2table(lowerPartData_middle,'VariableNames',{'rank using tol','rank using 2*tol','rank using 3*tol','rank using 5*tol','rank using 10*tol','rank using 100*tol'},'RowNames',newRowNames);
    lowerPartTable_right=array2table(lowerPartData_right,'VariableNames',{'rank using cpu'},'RowNames',newRowNames);
    lowerPartTable=[lowerPartTable_left lowerPartTable_middle lowerPartTable_right];
    writetable(lowerPartTable,excelFilePath1,'WriteRowNames',true,'Range',sprintf('A%d',nrows+3));
    footnote={'Ranking Criteria for E:J cols.','1st value=rank using','#0s,tie breaker t.e.','2nd value=rank using','#0s,tie breaker g.e.'};
    writecell(footnote,excelFilePath1,'Range',sprintf('A%d',2*nrows+5));






end