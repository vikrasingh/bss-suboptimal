function [mu,sigma]=randomCorrMatrix_MeanVec(infoCol)
    % Define the data correlation parameters mu and sigma randomly, mu is the 1 by pDim mean array used to generate multivariate normal random vectors for the predictor matrix xMatrix.
    rng(0,'twister'); % seed is 0 to generate mu and sigma
    mu=normrnd(0,1,[1,infoCol]);
    % sigma is the symmetric pDim by pDim positive semi-definite correlation matrix to generate multivariate normal random vectors for the predictor matrix xMatrix.
    sigma=normrnd(0,1,[infoCol,infoCol]);  
    for i=1:infoCol
        sigma(:,i)=sigma(:,i)./sigma(i,i);
    end
    sigma=(sigma+sigma')/2;sigma=(sigma*sigma)/2;
   
%     df=infoCol+1;sigma=wishrnd(sigma,df)/df; % df ?? generating sigma from Wishart distribution
     
end