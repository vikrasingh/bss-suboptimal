function readRawOzoneData
% to read data from raw.xlsx
% Col 1 is yArray
% Col 2 to 9 are x1,...,x8 
% Col 10 is just the day of the year, skip it

pDim=44;
nPts=330;
xMatrix=zeros(nPts,pDim);  % initialization
rawData=readmatrix( 'raw.xlsx'  );
yArray=rawData(:,1);
xMatrix(:,1:8)=rawData(:,2:9);
idx=8; % of the column
i=1;

while idx<pDim
    for j=1:i
       idx=idx+1;
       fprintf('x%d%d, ',j,i);
       newCol=xMatrix(:,j).*xMatrix(:,i);
       xMatrix(:,idx)=newCol;
    end
    i=i+1;
end


% write the data in an excel file "OzoneData"
writematrix(pDim,'OzoneData.xlsx','Range','A1');
writematrix(nPts,'OzoneData.xlsx','Range','A2');
writematrix(xMatrix,'OzoneData.xlsx','Range','A3:AR332');
writematrix(yArray,'OzoneData.xlsx','Range','A333:A662');


end