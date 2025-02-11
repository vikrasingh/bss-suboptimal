function [giveValuesOf_b]=trueb_type_pool(pDim,numOfEgs,n_non_zeros_trueb,trueb_type)
% 18 April22    
% trueb_type= -1  define randomly the way we are doing it now
%           = 0 read from the trueb_pool
%           = 1 this is beta type 1 in the Hastie et.al. (2020 IMS) BSS,FS and LASSO comparison. page 586
%           = 2 this is beta type 2 in the Hastie et.al. (2020 IMS) BSS,FS and LASSO comparison. page 586 
%           = 3 this is beta type 3 in the Hastie et.al. (2020 IMS) BSS,FS and LASSO comparison. page 586 
%           = 4 this is beta type 4 in the Hastie et.al. (2020 IMS) BSS,FS and LASSO comparison. page 586 
%           = 5 this is beta type 5 in the Hastie et.al. (2020 IMS) BSS,FS and LASSO comparison. page 586

    giveValuesOf_b=zeros(pDim,1); % initialization
    if trueb_type==1  % beta_type 1    
       idx=linspace(1,pDim,n_non_zeros_trueb); % find n_non_zeros_Trueb equally spaced points in 1:pDim
       idx=ceil(idx); % round all the non-integer values to the nearest bigger integer
       giveValuesOf_b(idx)=1; % make all the indices in the idx 1, leave rest to be 0   
    
    elseif trueb_type==2  % beta_type 2
       giveValuesOf_b(1:n_non_zeros_trueb)=1;  % change the first n_non_zeros_Trueb entries to be 1
    
    elseif trueb_type==3  % beta_type 3
       equally_spaced_values=flip( linspace(0.5,10,n_non_zeros_trueb) ); % num0strueb entries equally spaced from 10 to 0.5 
       giveValuesOf_b(1:n_non_zeros_trueb)=equally_spaced_values;  % change the first n_non_zeros_Trueb entries
    
    elseif trueb_type==4  % beta_type 4
       equally_spaced_values=linspace(-10,10,n_non_zeros_trueb); % n_non_zeros_Trueb entries equally spaced from -10 to 10 
       giveValuesOf_b(1:n_non_zeros_trueb)=equally_spaced_values;  % change the first n_non_zeros_Trueb entries
    
    elseif trueb_type==5  % beta_type 5
       giveValuesOf_b(1:n_non_zeros_trueb)=1;  % change the first n_non_zeros_Trueb entries to be 1
       for i=(n_non_zeros_trueb+1):pDim
           giveValuesOf_b(i)=0.5^( i - (n_non_zeros_trueb)); % change the remaining components to be decaying exp. to 0
       end
    
    end





end