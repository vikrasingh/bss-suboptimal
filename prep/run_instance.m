function  [rParaOut,numOfRows,numOfCols,T,outputPara,f_star]=run_instance(nPts,pDim,yArray,xMatrix,trueb,A,b,c,k_values,cputime_ols,...
                           x_ols,set_ctr,paraCtr,...
                           IotherPara,IstopCondPara,iPara,rPara,version_flag,...
                           level1_k_dir,level1_k_dir_path,tod_date,...
                           toDebug)
     
        tmax=k_values(paraCtr);
              
        % Define no. of rows and columns of the final output table according to the input data
        numOfRows=4+pDim;
        numOfCols=3;
        dataMat=ones(numOfRows,numOfCols);header=cell(1,numOfCols);
        header{1}=sprintf('k=%g trueb(%d)',tmax,set_ctr);
       
        %% Loop to call different versions for above defined example data
        
        ctr1=2;colsForQuadMin=0;
        colIndex=ctr1;  % 20 because we have 20 algo. so far to find xRelaxedOpt
        
        header{2+colsForQuadMin}=sprintf('k=%g xols(%d)',tmax,set_ctr);
        
        
        % BSS using exact or heuristic methods
        if version_flag==9
           if toDebug>=1,str=append('sfs1',level1_k_dir,tod_date);diaryFilePath=fullfile(level1_k_dir_path,sprintf('%s.txt',str));fileid=fopen(diaryFilePath,'a+');fprintf(fileid,'Date time = %s \n',datetime('now'));fprintf(fileid,'setCtr = %d \n',set_ctr);end
           string=sprintf('k=%g sfs1(%d)',tmax,set_ctr);
        elseif version_flag==92
           if toDebug>=1,str=append('sfs2',level1_k_dir,tod_date);diaryFilePath=fullfile(level1_k_dir_path,sprintf('%s.txt',str));fileid=fopen(diaryFilePath,'a+');fprintf(fileid,'Date time = %s \n',datetime('now'));fprintf(fileid,'setCtr = %d \n',set_ctr);end
           string=sprintf('k=%g sfs2(%d)',tmax,set_ctr); 
        elseif version_flag==93
           if toDebug>=1,str=append('sffs',level1_k_dir,tod_date);diaryFilePath=fullfile(level1_k_dir_path,sprintf('%s.txt',str));fileid=fopen(diaryFilePath,'a+');fprintf(fileid,'Date time = %s \n',datetime('now'));fprintf(fileid,'setCtr = %d \n',set_ctr);end
           string=sprintf('k=%g sffs(%d)',tmax,set_ctr);
        elseif version_flag==94
           if toDebug>=1,str=append('ga',level1_k_dir,tod_date);diaryFilePath=fullfile(level1_k_dir_path,sprintf('%s.txt',str));fileid=fopen(diaryFilePath,'a+');fprintf(fileid,'Date time = %s \n',datetime('now'));fprintf(fileid,'setCtr = %d \n',set_ctr);end
           string=sprintf('k=%g ga(%d)',tmax,set_ctr);   
        elseif version_flag==34 
           if toDebug>=1,str=append('mio',level1_k_dir,tod_date);diaryFilePath=fullfile(level1_k_dir_path,sprintf('%s.txt',str));fileid=fopen(diaryFilePath,'a+');fprintf(fileid,'Date time = %s \n',datetime('now'));fprintf(fileid,'setCtr = %d \n',set_ctr);end
           string=sprintf('k=%g mio(%d)',tmax,set_ctr);   
        elseif version_flag==35   
           if toDebug>=1,str=append('sfs',level1_k_dir,tod_date);diaryFilePath=fullfile(level1_k_dir_path,sprintf('%s.txt',str));fileid=fopen(diaryFilePath,'a+');fprintf(fileid,'Date time = %s \n',datetime('now'));fprintf(fileid,'setCtr = %d \n',set_ctr);end
           string=sprintf('TM=%g sfs(%d)',tmax,set_ctr);   
        elseif version_flag==36   
           if toDebug>=1,str=append('dfo',level1_k_dir,tod_date);diaryFilePath=fullfile(level1_k_dir_path,sprintf('%s.txt',str));fileid=fopen(diaryFilePath,'a+');fprintf(fileid,'Date time = %s \n',datetime('now'));fprintf(fileid,'setCtr = %d \n',set_ctr);end
           string=sprintf('TM=%g dfo(%d)',tmax,set_ctr); 
        elseif version_flag==37
           if toDebug>=1,str=append('abess',level1_k_dir,tod_date);diaryFilePath=fullfile(level1_k_dir_path,sprintf('%s.txt',str));fileid=fopen(diaryFilePath,'a+');fprintf(fileid,'Date time = %s \n',datetime('now'));fprintf(fileid,'setCtr = %d \n',set_ctr);end
           string=sprintf('k=%g abess(%d)',tmax,set_ctr);  
        end
            
        textfileIntermOut=fopen( append(level1_k_dir_path,'\summary_',sprintf('k%g',tmax),'.txt' ),'a');  % save another text file for intermediate results for stepTm values
        fprintf(textfileIntermOut,'%s \n',datetime('now'));
        if ismember(version_flag,[9 91 92 93]),fprintf(textfileIntermOut,'Iter ;fbest ;cpuSec ;fbest reduction due to refine. within that iter. ; \n');end  % for the first iteration only

        [outputPara,rParaOut,x_star,f_star]=call_algo(nPts,pDim,yArray,xMatrix,A,b,c,tmax,IstopCondPara,IotherPara,iPara,rPara,version_flag,x_ols,textfileIntermOut,toDebug);
        fprintf(textfileIntermOut,'%s \n',datetime('now'));
        fclose(textfileIntermOut); 
        colIndex=colIndex+1;
        % Calculate the distance bw the boundary of the box and Xstar
        dataCol=[f_star;outputPara(7);rParaOut.cpu_elapsed/60;outputPara(5);x_star];  % for table format
        dataMat(:,colIndex)=dataCol;  % for table format
        header{colIndex}=string;      % header will give us the top row of the table.
        if toDebug>=1,fprintf(fileid,'Date time =%s \n',datetime('now'));fclose(fileid);end

        % to display xRelaxedOpt column in the output table
        dataMat(:,2)=[fx(x_ols,pDim,A,b,c) ;missing;cputime_ols/60;missing;x_ols];
        
        % To display true b data column to the output table
        fbestTrueb=fx(trueb,pDim,A,b,c);
        truebData=[fbestTrueb;missing;missing;missing;trueb];
        
        dataMat(:,1)=truebData;
        
        %% To prepare the output table and write it into an excel file
        % block below is to create the Table format---
        rownames=cell(1,pDim+4);
        rownames(1:4)={'rss','stop flag','CPU time min','num iter'};
        for j=1:pDim
            rownames{j+4}=sprintf('x(%d)',j);
        end
        format short
        T=array2table(dataMat,'VariableNames',header,'RowNames',rownames);

        
        
end
