function  [epsMax,x,fval]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,S,iPara,rPara,which_algo,isSoftStop,isTF,targetfbest)
% till 08/18/23 [epsMax,x,fval]=quad_min_algo_pool(n,A,b,c,lb,ub,x0,iPara,rPara,which_algo,isSoftStop,isTF,targetfbest)
%  S is a structure S.Q = Q matrix from spectral decomposition such that A=QDQ'  where Q is an orthonormal matrix and D is diagonal with eig. value in the decreasing order 
%                  S.nNonZeroEig = no. of non zero eigenvalues  
%                  S.D = diag vector  of order nNonZeroEig x 1 of the D matrix 
%                  S.diagInv= 1/diag(A)  nx1 vector

% 8 April22, 7/3/23, MS:7/20/23(rPara(2)=epsx, rPara(3)=epsg,optUsingMinQ), 24
% A common algo. pool for all the algorithms where we are using quadratic min. options either for refinement or for Inf(F(Y))
% Input:=======================================================================================================================
% n,A,b,c are the quadratic algo. parameters
% iPara,rPara are user defined integer and real parameters resp.
% which_algo flag to choose the algo.
% isSoftStop = 1 means use soft stopping criteria for the algo. (for Inf(F(Y)) purpose)
%            = 0 means use hard stopping criteria (for refinement purpose)
% 28Aug22, isTF= 1 use targetfbest provided
%              = 0 do not use the target fbest
%=============================================================================================================================
% Called by the following subroutines
% getAFeasiblePt_L0andPertL0.m  
% lbF_envelope.m\lbUsingSDWorkingBox.m
% AG9.m, AG91.m, AG92.m, AG93.m
% BB_AG71.m, BB_AG72.m, BB_AG73.m, BB_AG74.m
% IBBplus_min.m, BBPPplus_min.m, FBBplus_min.m
%==============================================================================================================================

    epsMax=-inf; % initialization suppressed necessaryCondForBoxQuadMinFull on 23May23
    if which_algo==1  % use minq8
        [x,fval,~,~,~]=optUsingMinQ8(n,A,b,c,lb,ub,x0,rPara);     %7/20/23  [x,fval,g,~,~]=optUsingMinQ8(n,A,b,c,lb,ub,rPara);
%         [~,epsMax]=necessaryCondForBoxQuadMinFull(n,lb,ub,rPara(2),rPara(3),x,g);
    elseif which_algo==2  % use SD2
        [x,fval,~,~]=MiniAllQuadFunOverBox_bySD2(n,A,b,c,lb,ub,x0,iPara,rPara);
%         [~,epsMax]=necessaryCondForBoxQuadMinFull(n,lb,ub,rPara(2),rPara(3),x,A*x+b);
    elseif which_algo==3 % use LLS
        [x,fval,~,~]=MiniBsdQuadFun_byLLS_Wrap(n,A,b,c,lb,ub,iPara,rPara);
%         [~,epsMax]=necessaryCondForBoxQuadMinFull(n,lb,ub,rPara(2),rPara(3),x,A*x+b);
    elseif which_algo==4 % use Hbd
        [x,fval,~,~]=MiniQuadFunOverBox_byHbd_Wrap(n,A,b,c,lb,ub,x0,iPara,rPara);
%         [~,epsMax]=necessaryCondForBoxQuadMinFull(n,lb,ub,rPara(2),rPara(3),x,A*x+b);
    elseif which_algo==5 % use SequMiniQuadFun_v11  min function
        [x,fval,~,~,~,~]=SequMiniQuadFun_v11(n,A,b,c,lb,ub,iPara,rPara,0);
%         [~,epsMax]=necessaryCondForBoxQuadMinFull(n,lb,ub,rPara(2),rPara(3),x,A*x+b);
    elseif which_algo==6 % use parallel tangent + coord. descent  min function
        [~,x,fval,~,~]=convexQuadByParTanCorDesBox(n,A,b,S,lb,ub,x0,rPara(2),rPara(3),rPara(5),iPara(2),rPara(8),rPara(9),isSoftStop,isTF,targetfbest);
        fval=fval+c;
%{      
   08/13/2023  id 12 is used for convexQuadByParTanCorDesBox.m now 
        [x,fval,~,~,~,~]=SequMiniQuadFun_v12(n,A,b,c,lb,ub,iPara,rPara,0);
%         [~,epsMax]=necessaryCondForBoxQuadMinFull(n,lb,ub,rPara(2),rPara(3),x,A*x+b);
%}
    elseif which_algo==7 % use matlab quadprog function
        options=optimset('Display','off');
        [x,fval] = quadprog(A,b,[],[],[],[],lb,ub,x0,options);
        fval=fval+c;
%         [~,epsMax]=necessaryCondForBoxQuadMinFull(n,lb,ub,rPara(2),rPara(3),x,A*x+b);
    elseif which_algo==8  % use SD
        [~,x,fval,~,~]=convexQuadBySteepestDescentBox3(n,A,b,S,lb,ub,x0,rPara(2),rPara(3),rPara(5),iPara(2),isSoftStop,isTF,targetfbest);
        fval=fval+c;
%      used till 25June23   [~,x,fval,~,~]=convexQuadBySteepestDescentBox2(n,A,b,c,lb,ub,x0,rPara(3),rPara(3),iPara(2),isSoftStop,isTF,targetfbest); 
%         [~,epsMax]=necessaryCondForBoxQuadMinFull(n,lb,ub,rPara(2),rPara(3),x,A*x+b);
%    used till 28Aug22     [~,x,fval,~,~]=convexQuadBySteepestDescentBox(n,A,b,c,lb,ub,x0,rPara(3),rPara(3),iPara(2),isSoftStop);
%    earlier version     [x,fval,~,~]=MiniAllQuadFunOverBox_bySD(n,A,b,c,lb,ub,x0,iPara,rPara);
    elseif which_algo==9 % use CD coordinate descent
      [~,x,fval,~,~]=convexQuadByCoordDescentBox3d(n,A,b,S,lb,ub,x0,rPara(2),rPara(3),rPara(5),iPara(2),isSoftStop,isTF,targetfbest); 
      fval=fval+c;
%     used till 25June23  [~,x,fval,~,~]=convexQuadByCoordDescentBox2(n,A,b,c,lb,ub,x0,rPara(3),rPara(3),iPara(2),isSoftStop,isTF,targetfbest); 
%       [~,epsMax]=necessaryCondForBoxQuadMinFull(n,lb,ub,rPara(2),rPara(3),x,A*x+b);
%    used till 28Aug22     [~,x,fval,~,~]=convexQuadByCoordDescentBox(n,A,b,c,lb,ub,x0,rPara(3),rPara(3),iPara(2),isSoftStop);
%    earlier version     [x,fval,~]=MinimizeQuadraticFun_byCD(n,A,b,c,lb,ub,x0,iPara,rPara);
    elseif which_algo==10  % use CG conjugate gradient
        [x,fval,~]=MinimizeQuadraticFun_byCG(n,A,b,c,lb,ub,x0,iPara,rPara);    %MS7/24/23: not done
%         [~,epsMax]=necessaryCondForBoxQuadMinFull(n,lb,ub,rPara(2),rPara(3),x,A*x+b);
    elseif which_algo==11  % use PT parallel tangent
        [~,x,fval,~,~]=convexQuadByParallelTangentBox3(n,A,b,S,lb,ub,x0,rPara(2),rPara(3),rPara(5),iPara(2),isSoftStop,isTF,targetfbest);
        fval=fval+c;
%   used till 25June23      [~,x,fval,~,~]=convexQuadByParallelTangentBox2(n,A,b,c,lb,ub,x0,rPara(3),rPara(3),iPara(2),isSoftStop,isTF,targetfbest);
%         [~,epsMax]=necessaryCondForBoxQuadMinFull(n,lb,ub,rPara(2),rPara(3),x,A*x+b);
%   used till 28Aug22      [~,x,fval,~,~]=convexQuadByParallelTangentBox(n,A,b,c,lb,ub,x0,rPara(3),rPara(3),iPara(2),isSoftStop);
%    earlier version     [x,fval,~]=MinimizeQuadraticFun_byPT(n,A,b,c,lb,ub,x0,iPara,rPara);
    elseif which_algo==12 % conjugate gradient + Coord. descent solution
        [~,x,fval,~,~]=convexQuadByConGraCorDesUnbox(n,A,b,S,x0,rPara(2),rPara(3),rPara(5),iPara(2),rPara(7),rPara(8),rPara(9),isSoftStop,isTF,targetfbest);
        fval=fval+c;
%{
  08/13/2023  id 12 is used for convexQuadByConGraCorDesUnbox.m now 
        warning('off','MATLAB:nearlySingularMatrix');
        x=A\(-b);  % use the backslash operator which uses QR decomposition with pivoting to solve the
                   % first order optimality condition Ax=-b
        fval=fx(x,n,A,b,c); 
%         [~,epsMax]=necessaryCondForUnBoxQuadMinFull(n,rPara(3),x,A*x+b);
%}
    elseif which_algo==13 % unconstrained SD
        [~,x,fval,~,~]=convexQuadBySteepestDescentUnbox3(n,A,b,S,x0,rPara(2),rPara(7),rPara(3),rPara(5),iPara(2),isSoftStop,isTF,targetfbest);
        fval=fval+c;
%  used till 23July23        [~,x,fval,~,~]=refineBySteepestDescentUnbox4(n,A,b,c,x0,rPara(2),rPara(7),rPara(3),rPara(5),iPara(2),isSoftStop,isTF,targetfbest);
%  used till 25June23        [~,x,fval,~,~]=refineBySteepestDescentUnbox3(n,A,b,c,x0,rPara(3),rPara(7),rPara(3),iPara(2),isSoftStop,isTF,targetfbest);
%  used till 10June23        [~,x,fval,~,~]=refineBySteepestDescentUnbox2(n,A,b,c,x0,rPara(3),rPara(7),rPara(3),iPara(2),isSoftStop,isTF,targetfbest);
%         [~,epsMax]=necessaryCondForUnBoxQuadMinFull(n,rPara(3),x,A*x+b);
%   used till 28Aug22        [~,x,fval,~,~]=refineBySteepestDescentUnbox(n,A,b,c,x0,rPara(3),rPara(7),rPara(3),iPara(2),isSoftStop);
%    earlier version         [x,fval,~,~]=MiniAllQuadFunUnconst_bySD(n,A,b,c,x0,iPara,rPara,isSoftStop);   % isSoftStop=1 means use soft stopping criteria
      
    elseif which_algo==14 % unconstrained PT
          [~,x,fval,~,~]=convexQuadByParallelTangentUnbox3(n,A,b,S,x0,rPara(2),rPara(7),rPara(3),rPara(5),iPara(2),isSoftStop,isTF,targetfbest);
          fval=fval+c;
%  used till 23July23          [~,x,fval,~,~]=refineByParallelTangentUnbox4(n,A,b,c,x0,rPara(2),rPara(7),rPara(3),rPara(5),iPara(2),isSoftStop,isTF,targetfbest);
%  used till 25June23    [~,x,fval,~,~]=refineByParallelTangentUnbox3(n,A,b,c,x0,rPara(3),rPara(7),rPara(3),iPara(2),isSoftStop,isTF,targetfbest);
%  used till 10June23      [~,x,fval,~,~]=refineByParallelTangentUnbox2(n,A,b,c,x0,rPara(3),rPara(7),rPara(3),iPara(2),isSoftStop,isTF,targetfbest);
%         [~,epsMax]=necessaryCondForUnBoxQuadMinFull(n,rPara(3),x,A*x+b);
%   used till 28Aug22      [~,x,fval,~,~]=refineByParallelTangentUnbox(n,A,b,c,x0,rPara(3),rPara(7),rPara(3),iPara(2),isSoftStop);
%   earlier version     [x,fval,~,~]=MiniAllQuadFunUnconst_byPT(n,A,b,c,x0,iPara,rPara,isSoftStop);   % isSoftStop=1 means use soft stopping criteria 
        
    elseif which_algo==15 % unconstr. Adaptive Method 1
        [x,fval,~,~]=MiniQuadFunUnconstAdapt1(n,A,b,c,x0,iPara,rPara,isSoftStop); 
%         [~,epsMax]=necessaryCondForUnBoxQuadMinFull(n,rPara(3),x,A*x+b);
    elseif which_algo==16 % unconstr. Adaptive Method 2
        [x,fval,~,~]=MiniQuadFunUnconstAdapt2(n,A,b,c,x0,iPara,rPara,isSoftStop);
%         [~,epsMax]=necessaryCondForUnBoxQuadMinFull(n,rPara(3),x,A*x+b);
    elseif which_algo==17 % unconstr. Adaptive Method 3    
        [x,fval,~,~]=MiniQuadFunUnconstAdapt3(n,A,b,c,x0,iPara,rPara,isSoftStop);
%         [~,epsMax]=necessaryCondForUnBoxQuadMinFull(n,rPara(3),x,A*x+b);

    elseif which_algo==18 % unconstr. CD
        [~,x,fval,~,~]=convexQuadByCoordDescentUnbox3(n,A,b,S,x0,rPara(2),rPara(7),rPara(3),iPara(2),isSoftStop,isTF,targetfbest);
        fval=fval+c;
%  used till 23July23        [~,x,fval,~,~]=refineByCoordDescentUnbox2(n,A,b,c,x0,rPara(2),rPara(7),rPara(3),iPara(2),isSoftStop,isTF,targetfbest);
%         [~,epsMax]=necessaryCondForUnBoxQuadMinFull(n,rPara(3),x,A*x+b);
%  used till 28Aug22       [~,x,fval,~,~]=refineByCoordDescentUnbox(n,A,b,c,x0,rPara(3),rPara(7),rPara(3),iPara(2),isSoftStop);
%   earlier version      [~,x,fval,~,~]=MiniAllQuadFunUnconst_byCD(n,A,b,c,x0,rPara);
    elseif which_algo==19
        [x,fval,~,~]=MiniQuadFunUnconstAdapt4(n,A,b,c,x0,iPara,rPara,isSoftStop);
%         [~,epsMax]=necessaryCondForUnBoxQuadMinFull(n,rPara(3),x,A*x+b);
    elseif which_algo==20
        [~,x,fval,~,~]=convexQuadByConjugateGradientUnbox2(n,A,b,S,x0,rPara(2),rPara(7),rPara(3),iPara(2),isSoftStop,isTF,targetfbest);
        fval=fval+c;
%  used till 23July23        [~,x,fval,~,~]=refineByConjugateGradientUnbox2(n,A,b,c,x0,rPara(2),rPara(7),rPara(3),iPara(2),isSoftStop,isTF,targetfbest);
%         [~,epsMax]=necessaryCondForUnBoxQuadMinFull(n,rPara(3),x,A*x+b);
    end
       
end