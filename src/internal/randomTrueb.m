function [seedForEgs,giveValuesOf_b]=randomTrueb(pDim,numOfEgs,num0sTrueb,seedForTrueb,seedToRunNInstances)
% to generate synthetic data for trueb
% if seed=0 and same numOfEgs, will generate same trueb values
% if seed/=0 then re-generate trueb for one particular example
% Restriction : if seed/=0 then numOfEgs=1


    if seedForTrueb==0  
       rng(seedToRunNInstances,'twister');  % use seedForNinstances as a control parameter for the randomness of the trueb
       seedForEgs=randi([0,10000],1,numOfEgs);
    else % if running one particular example
       seedForEgs=seedForTrueb; 
    end
    
    giveValuesOf_b=zeros(pDim,numOfEgs);
    for i=1:numOfEgs
        giveValuesOf_b(:,i)=getRand_b(pDim,num0sTrueb,seedForEgs(i));
    end
end

function [X]=getRand_b(pDim,num0sTrueb,seed)
    
    rng(seed,'twister');
    X=ones(pDim,1);  % generate a vector of ones
    
    % make some entries as large integer values
    entries1=randperm(pDim,num0sTrueb);  % pick random #num0sTrueb entries to be set to 0
    X(entries1)=0;
    
    
    % make some entries of the array having 0 entries
    entries2=setdiff(1:pDim,entries1);  % pick the remaining entries, and assign them random integer values bw -5 and 5
    getSize=length(entries2);
    X( entries2 )=randi([1 5],getSize,1); % assign integer values to those entries, with integer bw [-5 5], change accordingly
  
       
    if X==0  % protection
        [X]=getRand_b(pDim,num0sTrueb,seed+1);   % just to make sure that true parameter b is not a zero vector, change the seed by 1
        if X==0  % second protection
            [X]=getRand_b(pDim,num0sTrueb,seed+2);     % again to add some more security, change the seed by 2
        end
    end
    
end


% function [X]=getRand_b(pDim,num0sTrueb,seed)
%     
%     rng(seed,'twister');
%     X=randn(pDim,1);  % generate a random vector from standard normal dist.
%     
%     % make some entries as large integer values
%     entries1=randperm(pDim,num0sTrueb);  % pick random #num0sTrueb entries to be set to 0
%     X(entries1)=0;
%     
%     
%     % make some entries of the array having 0 entries
%     entries2=setdiff(1:pDim,entries1);  % pick the remaining entries
%     getSize=length(entries2);
%     % entries2(randperm(getSize, floor((1/2)*getSize)));   % randomly pick half of the entries from the remaining entries2 ,other than picked earlier in entries1 to assign integer values
%     X( entries2(randperm(getSize, floor((1/2)*getSize))) )=randi([-5 5],length( entries2(randperm(getSize, floor((1/2)*getSize))) ),1); % assign integer values to those entries, with integer bw [-20 20], change accordingly
%   
%        
%     if X==0  % protection
%         [X]=getRand_b(pDim,num0sTrueb,seed+1);   % just to make sure that true parameter b is not a zero vector, change the seed by 1
%         if X==0  % second protection
%             [X]=getRand_b(pDim,num0sTrueb,seed+2);     % again to add some more security, change the seed by 2
%         end
%     end
%     
% end