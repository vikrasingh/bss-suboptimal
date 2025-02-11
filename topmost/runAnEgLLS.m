function  [header,rParaOut,numOfBoxDel,dataMat,XstarNormSpace,outputPara,...
        fbest]=runAnEgLLS(nPts,pDim,yArray,xMatrix,trueb,A,b,c,tmax,targetfbest,...
                           xRelaxedOpt,fxRelaxedOpt,preXstarNormSpace,InpFilCtr,pDimCtr,egCtr,setCtr,paraCtr,...
                           IotherPara,IstopCondPara,choseniPara,rPara,version_flag,delCondPara,...
                           level1FolderName,level1FolderPath,toDate,saveIntermOutput,toDebug,level8dirpath)
   
             
        %% Loop to call different versions for above defined example data
        
        colIndex=2;  % Col 1. trueb , Col 2. xRelaxedOpt
        
        str1=append('IntermSol.',sprintf('setCtr%d-',setCtr),level1FolderName);saveFilePathIntermOut=fullfile(level1FolderPath,sprintf('%s.txt',str1));
        str2=append('lbFafbest_',sprintf('setCtr%d-',setCtr),level1FolderName);saveFilePathlbFandfbest=fullfile(level1FolderPath,sprintf('%s.txt',str2));
       
       
        %% BSS using exact or heuristic methods
        if ismember(version_flag,[23 34])  % to test hybrid/relaxed interval algo. with equality or inequality constraint
            if version_flag==23 % BB
               header=sprintf('TM=%g BB(%d)',tmax,setCtr); 
               if toDebug>=1,str=append('BB',level1FolderName,toDate);diaryFilePath=fullfile(level1FolderPath,sprintf('%s.txt',str));fileid=fopen(diaryFilePath,'a+');fprintf(fileid,'Date time = %s \n',datetime('now'));fprintf(fileid,'setCtr = %d \n',setCtr);end
            elseif version_flag==34 % bssMATLAB
                header=sprintf('TM=%g MIO(%d)',tmax,setCtr);
               if toDebug>=1,str=append('MIO',level1FolderName,toDate);diaryFilePath=fullfile(level1FolderPath,sprintf('%s.txt',str));fileid=fopen(diaryFilePath,'a+');fprintf(fileid,'Date time = %s \n',datetime('now'));fprintf(fileid,'setCtr = %d \n',setCtr);end   
            end
            
            textfileIntermOut=fopen( append(level1FolderPath,'\summary_',sprintf('TM%g',tmax),'.txt' ),'a');  % save another text file for intermediate results for stepTm values
            fprintf(textfileIntermOut,'%s \n',datetime('now'));
            [XstarNormSpace,xRelaxedOpt2,outputPara,numOfBoxDel,rParaOut,fx_tilde,x_tilde,Xstar,fbest]=setupForIntvalAlgo(nPts,pDim,yArray,xMatrix,trueb,...
                                                   A,b,c,tmax,targetfbest,IstopCondPara,IotherPara,choseniPara,rPara,delCondPara,saveIntermOutput,version_flag,xRelaxedOpt,...
                                                        fxRelaxedOpt,preXstarNormSpace,InpFilCtr,pDimCtr,egCtr,setCtr,paraCtr,textfileIntermOut,toDebug,level8dirpath);
            fprintf(textfileIntermOut,'%s \n',datetime('now'));
            fclose(textfileIntermOut); 
            dataMat=[fbest;outputPara(7);rParaOut.cpuIntvalAlgo/60;outputPara(5);outputPara(15);Xstar];  % for table format
            if toDebug>=1,fprintf(fileid,'Date time =%s \n',datetime('now'));fclose(fileid);end
            
        end %
        
        
        %% IBB+
        if version_flag==71 
            if toDebug>=1,str=append('IBB+ ',level1FolderName,toDate);diaryFilePath=fullfile(level1FolderPath,sprintf('%s.txt',str));fileid=fopen(diaryFilePath,'a+');fprintf(fileid,'Date time =%s \n',datetime('now'));fprintf(fileid,'setCtr = %d \n',setCtr);end
            
            header=sprintf('TM=%g IBB+(%d)',tmax,setCtr);
            if toDebug==2 % save both intermediate updates and boxDeletionFlag output
               textfileBoxFlag=fopen( fullfile(level1FolderPath,sprintf('boxDeleteFlag Tm%d setCtr%d.txt',tmax,setCtr) ) ,'w');
               textfileIntermOut=fopen(saveFilePathIntermOut,'w');fprintf(textfileIntermOut,'%s \n',datetime('now'));
               fprintf(textfileIntermOut,header);fprintf(textfileIntermOut,'\n'); 
            elseif toDebug==1 % only save the intermediate results
               textfileBoxFlag=[];
               textfileIntermOut=fopen(saveFilePathIntermOut,'w');fprintf(textfileIntermOut,'%s \n',datetime('now'));
               fprintf(textfileIntermOut,header);fprintf(textfileIntermOut,'\n'); 
            else % neither intermediate results nor boxDeletionFlag is saved
               textfileIntermOut=[];textfileBoxFlag=[];
            end
            textfilelbFandfbest=fopen(saveFilePathlbFandfbest,'w');

            textfileName=struct('interm_out',textfileIntermOut,'box_flag',textfileBoxFlag,'lbFandfbest',textfilelbFandfbest);
               [XstarNormSpace,xRelaxedOpt2,outputPara,numOfBoxDel,rParaOut,fx_tilde,x_tilde,Xstar,fbest]=setupForIntvalAlgo(nPts,pDim,yArray,xMatrix,trueb,...
                                                   A,b,c,tmax,targetfbest,IstopCondPara,IotherPara,choseniPara,rPara,delCondPara,saveIntermOutput,version_flag,xRelaxedOpt,...
                                                        fxRelaxedOpt,preXstarNormSpace,InpFilCtr,pDimCtr,egCtr,setCtr,paraCtr,textfileName,toDebug,level8dirpath); 
            if toDebug==2                                     
                fprintf(textfileIntermOut,'%s \n',datetime('now'));
                fclose(textfileIntermOut);fclose(textfileBoxFlag);
            elseif toDebug==1, fclose(textfileIntermOut);    
            end  
            fclose(textfileName.lbFandfbest);                                                                              
            dataMat=[fbest;outputPara(7); rParaOut.cpuIntvalAlgo/60 ;outputPara(5);outputPara(15);Xstar];  % for table format
            if toDebug==1,fprintf(fileid,'Date time =%s \n',datetime('now'));fclose(fileid);end
  
        end % end flag which version to run using interval algo
        
        
        
     

%         % to increase the precision of the entries in T
%         new_T = varfun(@(x) num2str(x, ['%' sprintf('.%df',8)]), T);
%         % preserve the variable names and the row names in the original table
%         new_T.Properties.VariableNames = T.Properties.VariableNames;
%         new_T.Properties.RowNames = T.Properties.RowNames;
%         T=new_T;

        
        
        
        
end
