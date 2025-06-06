function [npDim,pDim_values,nPts_values,nEgs_eachpDim,tmaxVal_eachpDim,nDcol_values,infoCol_values,targetfbest_eachpDim,...
          giveCorrParaManually_option,corrflag_option,rho_values,snr_values,n_non_zeros_trueb_values,trueb_type_option,...
          chooseEgsToRun_eachpDim,egParafile,toUseRealDataFile]...
         =One1Ex2000Ud2TmCcRh0p8Sn0p05
% Ud, 
% snr=0.05, rho=0.8, Constant Correlation, trueb type =1

npDim=1   % number of different pDim sets to run
pDim_values=[2000]      % values of pDim to run
nPts_values=[300]      % nPts corresponding to each pDim

% set tmax para.
tmaxVal_eachpDim=[5  10;
                  5  10] % each row corresponds to the tmax values of given pDim run     

targetfbest_eachpDim=[0  0]

% set correlation para.
nDcol_values=[0 ]    % no need to change
infoCol_values=pDim_values    % how many correlated cols for each pDim run
giveCorrParaManually_option=[1  ]    % manually or randomly assigned for each pDim run
corrflag_option=[1  ]    % constant or exponential corr. for each pDim run
rho_values=[0.8  ]     % value of rho para for each pDim run
snr_values=0.05*[1  ]     % value of snr para for each pDim run

% trueb para.
n_non_zeros_trueb_values=[10 ]     % no. of non-zeros in trueb we want for each pDim run
trueb_type_option=[1 ]     % 0 or 1 for each pDim run

chooseEgsToRun_eachpDim=[1  2];  % if giveTruebManually_option=0 , will not get used
                            % give which trueb values (i.e. which examples) we want to test from the array below.


nEgs_eachpDim=[1] % nEgs for each run of pDim
toUseRealDataFile={}; % empty if running synthetic examples
for ip=1:npDim   
   if trueb_type_option(ip)==0  
      Egs1=chooseEgsToRun_eachpDim(ip,:);
      nEgs1=length(Egs1(Egs1~=0));
      nEgs_eachpDim(ip)=nEgs1;       % nEgs for each run of pDim  
   end
end



egParafile=mfilename('fullpath'); % (No need to change) to save a copy of this file in the output folder.

end