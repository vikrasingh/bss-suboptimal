function [npDim,pDim_values,nPts_values,nEgs_eachpDim,tmaxVal_eachpDim,nDcol_values,infoCol_values,targetfbest_eachpDim,...
               giveCorrParaManually_option,corrflag_option,rho_values,snr_values,n_non_zeros_trueb_values,trueb_type_option,...
               chooseEgsToRun_eachpDim,egParafile,toUseRealDataFile]=OzoneData1ExOd43Tm
% VS 3/7/26
% Special example input file to run real Ozone data set

npDim=1    % number of different pDim sets to run
pDim_values=[44  ]      % values of pDim to run: 
nPts_values=[330 ]      % nPts corresponding to each pDim:   -------------------- 10x, (1/5)x.

% set 3 tmax para.
tmaxVal_eachpDim=[1:43]   %---------------------------20%, 80%, other 1-2 small about 5 (up to 10%)
                 
targetfbest_eachpDim=zeros(1,43)
toUseRealDataFile={'OzoneData'} 

%% all the parameters below are just dummy variable for the run       
% set correlation para.
nDcol_values=zeros(1,npDim)    % dummy output for real data
infoCol_values=pDim_values    % dummy output for real data
giveCorrParaManually_option=ones(1,npDim)    % dummy output for real data, manually or randomly assigned for each pDim run
corrflag_option=1*ones(1,npDim)    % dummy output for real data, constant or exponential corr. for each pDim run
rho_values=1*ones(1,npDim)     % dummy output for real data, value of rho para for each pDim run
snr_values=1*ones(1,npDim)     % dummy output for real data, value of snr para for each pDim run

% trueb para.
n_non_zeros_trueb_values=5*ones(1,npDim)     %dummy output for real data, no. of 0s in trueb we want for each pDim run  ----20% of pDim & upto 20
trueb_type_option=ones(1,npDim)     %dummy output for real data, -1 to 5 for each pDim run

%If trueb_type_option=0, then also use this:
chooseEgsToRun_eachpDim=repmat([1 2],npDim,1);  %dummy output for real data, if giveTruebManually_option=0 , will not get used
                           % give which trueb values (i.e. which examples) we want to test from the array below.

nEgs_eachpDim=ones(1,npDim) %dummy output for real data, nEgs for each run of pDim

for ip=1:npDim   
   if trueb_type_option(ip)==0 
      Egs1=chooseEgsToRun_eachpDim(ip,:);
      nEgs1=length(Egs1(Egs1~=0));
      nEgs_eachpDim(ip)=nEgs1;       % nEgs for each run of pDim  
   end
end


egParafile=mfilename('fullpath'); % (No need to change) to save a copy of this file in the output folder.

end
