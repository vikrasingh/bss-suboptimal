function  [numOfSets,chooseParaToRun,IotherPara,IstopCondPara,delCondPara,QuadMinFunPara,rPara,iPara,intermFbestIter,intermMaxIter,intermCPU,algoParafile]...
    =Box7a2c3dV9
%  VS:05March24
% define the parameters to be used in all the stopping criteria of the Interval algo below---------------------------------
        %               1    2         3          4          5           6           7-10 ONLY FOR SEMI-INTERVAL 2
        %                 (fixed)   (fixed)         SEMI-ITVL 1                      7         8        9     10        11
        %               f*   dF      dX        It_Fbest   max.Iter      min.         X0     Iter_Lt     dt   Iter-Lx   cpu_Fbest
        IstopCondPara=[ 0    0.1     0.002      100000    2000000       1             0.3       20       0.2    30      0;    % 0  < pDim < 10
                        0    0.1     0.002      400000    8000000       5             0.3       30       0.2    40      0;    % 10 < pDim < 50  since 10/15/23
                        0    0.1     0.002      2000000   20000000      10           0.3       40       0.2    50       0;    % 50 < pDim < 200 since 10/15/23
                        0    0.1     0.002      500       5000          0.5          0.3       20       0.2    25       0;  
                        
                        0    0.1     0.002      100       2000000       15           0.3       20       0.2    30         2;   % 0  < pDim < 10  since 02/19/24
                        0    0.1     0.002      500       8000000       30           0.3       30       0.2    40         5;    % 10 < pDim < 50  since 02/19/24
                        0    0.1     0.002      1000      20000000      60           0.3       40       0.2    50         10;    % 50 < pDim < 200 since 02/19/24
                        0    0.1     0.002      2000      40000000      180          0.3       40       0.2    50         15;    % 200 < pDim   since 02/19/24

                        0    0.1     0.002      500       1000000       10           0.3       30       0.2    40         5];
         
        %{
        IstopCondPara(1)=1 means check if fbest-targetfbest < tol=1e-16), if yes stop, line 309, IntvalAlgo and IstopCondPara(4)=0 means do not check this condition
                            also a similar threshold has been used in checking the deletion condition if fbest < lb F(V) - 0.001  for a working box V, discard it
        IstopCondPara(2) defines the tol s.t. width(F(Z))< tol for every box Z in the list L. (NOTE: this stopping will be checked only when IotherPara(19)=0 or 1)
        IstopCondPara(3) defines the tol which ensures that the box selected from the list L for new iteration has width greater than this tol
        IstopCondPara(4) (ONLY FOR SEMI-INTERVAL 1)defines the no. of iterations for which if fbest does not improve then stop the interval algo.
        IstopCondPara(5) (ONLY FOR SEMI-INTERVAL 1)defines the max no. of iterations
        IstopCondPara(6) defines the cpu time limit in min for each run. 
                       
        IstopCondPara(7) (ONLY FOR SEMI-INTERVAL 2) defines the stopping tol for max. width of a box while initial partition of box X.
        IstopCondPara(8) (ONLY FOR SEMI-INTERVAL 2) defines the max. no. of iterations for the list L_T partition.
        IstopCondPara(9) (ONLY FOR SEMI-INTERVAL 2) defines the tolerance such that if the max. width of the box in the list L_T < tol then stop the algorithm.
        IstopCondPara(10) (ONLY FOR SEMI-INTERVAL 2) defines the max iterations while partitioning X box and adding subboxes in the list L_X
        IstopCondPara(11) (only for L0 constraint algo.) defined the cputime in min. such that if fbest does not improve for the given limit stop
        %}
       %            DC1   DC2   DC3    DC4    DC5    DC6   DC7   DC8 
        delCondPara=[1      1     1       0     1      1     0    1; 
                     1      1     1       0     1      1     1    1 ;
                     1      1     1       0     1      1     0    0 ;
                     1      1     0       0     1      1     0    0 ;
                     1      1     0       0     1      1     0    0 ;
                     0      1     1       0     1      1     0    0 ]
        %{         
                 0 means we are not checking that deletion condition, 1 means checking
        DC1 --- fbest-based. For alg7, we may choose not to use F(X) as we discussed on 8/23/21. In that case, we will simply bypass DC1.
        DC2 --- total infeasibility (supp>Tm)
        DC3 --- there is no acceptable direction to partition
        DC4 --- based on a saved deletion pattern (DP), proposed on 8/23/21
        DC5 --- impossible to get an acceptable feasible sample in the modified midpoint-based feasibility sampling procedure proposed on 8/24/21 although DC2 is not satisfied.
        DC6 --- there is no need for further partition, proposed on 8/24/21
    
        %}

numOfSets=3-1;
        %                IoP ISP QMP iP delP(1=safe,2/all&fair)
        chooseParaToRun=[
                         1   9  1+20   9  4; 
                         2   9  1+20   9  4; 
                         3   9  1+20   9  4; 
                         ]                                     
% Define the parameters in IotherPara below which are different options to be used in the Interval Algo--------------------            
%  1   2     3    4     5    6    *7    8     9     10      11      12     13   14     *15     16      17    18  19    20     21     22    23      24   
% X0QM %ELX0 '-"  X*   sel save F_delay f-s  cut  x*-zero X-Fdelay  dh    ucIvl scaleX alg     eta  dTmAG8  plot  F  dX-frac #SD1n2n FS-L0 refiAlg eps0
IotherPara= [
   5   100   2   NaN    4    0    1    NaN   -81  0.001   NaN     NaN    NaN    -3     71   NaN    NaN     0   8     NaN    420     7      11      NaN ;  % IBB+ 
   5   100   0   NaN    0    0    0    NaN    0  0.001    NaN     NaN    NaN    -3     23   NaN    NaN     0   8     NaN    420     7      11      NaN ;  % BB
   5   100   0   NaN    0    0    0    NaN    0  0.001    NaN     NaN    NaN    -3     34   NaN    NaN     0   8     NaN    420     7      11      NaN ;  % MIO 
   ]   

%{
        IotherPara(1) define the box XquadMinBox for Quad. Min. Fun. to get the relaxed optimal point
         = 0, use the whole space to find xRelaxedOpt 
         = 1,box just containing the true para b,  m=max(|b|) then a=-m and b=m where a<= XquadMinBox <=b
         = 2,box bigger than option 1, m=2*max(|b|) then a=-m and b=m where a<= XquadMinBox <=b
         = 3,box bigger than option 2, m=3*max(|b|) then a=-m and b=m where a<= XquadMinBox <=b
         = 4,box smaller than option 1, m=(1/2)*max(|b|) then a=-m and b=m where a<= XquadMinBox <=b
         = 5,WITHOUT USING true b , use iterative way of determining the solution using quadMinFun without using true b
IotherPara(2)= delta>=0 , Used only when IotherPara(14)=3 or -3
            where delta defines the % enlargement in the box X defined by xRelaxedOpt, if 0 then no enlargement the box defined by xRelaxedOpt will be used. (3 Aug21)      
                       normxRelaxedOpt(i)=xRelaxedOpt(i)/m    for i=1...pDim  and m=Max(abs(xRelaxedOpt))
                       box X is defined as follows:
                       if sign(xRelaxedOpt(i))>=0 , X(i)=[0 , normxRelaxedOpt(i)+0.01*delta*normxRelaxedOpt(i)]
                       elseif sign(xRelaxedOpt(i))<0 , X(i)=[normxRelaxedOpt(i)-0.01*delta*normxRelaxedOpt(i),0]
             
                       when delta=0, we get the box defined by xRelaxedOpt No enlargement
                       when delta=100, we get DOUBLE the size of the box defined by xRelaxedOpt
                       when delta=200, we get TRIPLE the size of the box defined by xRelaxedOpt
             
IotherPara(3)=(only for AG71) 1 use the acc. technique to reduce inclusion fun. eval. based on truncated xRelaxedOpt, for box
                               Y if truncate xRelaxedOpt to get xhat in Y, if f(xhat)<fbest , skip the lb f(Y) eval.
         = 0 means do not use it    
IotherPara(4)= 1 means select the final box Y_tilde from the final list on the basis of feasiblePt with the lowest function value.
      and IotherPara(1)=2 means select the final box Y_tilde from the final list on the basis of the lowest lb F(box) value.
IotherPara(5)= 1 means use the age of the box criteria, that is select the first box in the list to process,
         = 2 means select the maximum width box present in the list L to process.
         = 3 means select the box whose farPoint is at the max. distance from the origin.
         = 4 means pick that box which has the least lb F(X) for any box X in the list
         = 5 (only for AG7) pick that box from the list which has max. entries away from 0. i.e. box with
               max 2 in the box flag 112200
         = 6 (only for AG7) pick that box from the list which has max. no. of sub-boxes containing 0 but are not degenerate  
         = 7 (only for AG7) pick that box from the list which has max. no. of degenerate sub-boxes
         = 8 (only for AG7) best first search
         = 9 pick the last box in the list

IotherPara(6)= 1 to save the info. in the binary file to re-run the algo. later on.
        and 0 do not save the info. for re-run
IotherPara(7)= options for the inclusion function F
         = 0, standard way, evaluate at every iter.
         = 1, delay the call of F.
         = 2, use recursive approach to find F.(HAS NOT BEEN SET YET)
         = 3, do not call inclusion function F at all.
         = 4,(implemented in 72 so far)(only for AG7,71,72) use a prediction mechanism, find lb f using hard stop with no. of iter coming from second digit of IotherPara(21), if for some working box
               fbest <lb f, then verify this by finding actual lb f using soft stop, and if still fbest < lb f, discard that box.

IotherPara(8)= which SAMPLING PROCEDURE to use. Below are the options to choose from:
             = 0 (for AG77) skip sampling for subsequence of flag0 and flag2 boxes
             = -1m (for AG77) use sampling only for the last m child boxes for subsequence of flag0 and flag2 boxes 
           IotherPara(8)==9  uses the inequality constraint for sampling.
           IotherPara(8)==10  uses the equality constraint as given in the paper :  Ying M, Sun M, J-GOP(2015) Some feasibility sampling procedures in interval methods for C-GOP
           IotherPara(8)==11  (06 Nov20 RECURSIVE) sampling using 50-50 weights for bisection over the diagonal joining close and far point of the box.
           IotherPara(8)==12  (RECURSIVE) sampling using 40-60 weights for bisection over the diagonal joining close and far point of the box.
           IotherPara(8)==13  (RECURSIVE) sampling using 30-70 weights for bisection over the diagonal joining close and far point of the box.
           IotherPara(8)==14  (RECURSIVE) sampling using ADAPTIVE weights for bisection over the diagonal joining close and far point of the box.
           IotherPara(8)==15  (14 Nov20) uses face of the box X. Find an adjacent vertex to the close point which is above the constraint and solve for the parameter t to find the feasible point.
           IotherPara(8)==16  (20 Nov20) will use the homogenous ray from origin to the far point of the box. Solve for parameter t to find the intersection point of ray and the equality constraint.
           IotherPara(8)==17  (3 May21) will use homo. ray from origin to the far point of the box, find the feasible by reducing the dim of the problem.
           IotherPara(8)==18  (02 Dec20) solve for paramter t for the case of alpha=1 and 2 to find the feasible point and for rest of alpha values call sampling procedure 16
           IotherPara(8)==19  (04 Dec20) using partial homogenous and solve for paramter t to find the intersection point.

IotherPara(9)= defines how many cuts we want and on which point midpoint/feasiblePt.
         = 0, 1 cut at the midpoint the box X.
         = 1,(NO LONGER IN USE 21 Mar21)(DO NOT USE THIS IF RUNNING BRIDGE REGSS. i.e if IotherPara(15)=1 below) 2n cuts for the box X at the feasible point and add the sub-boxes to the list L w/o calling F.
         = 2, 1 cut at the feasible point along the longest side of the box X.
         = 1m, m cuts (from longest to shortest side of the box) at the midpoint of the box X
         = 2m, m cuts (from longest to shortest side of the box) at the feasible point of the box,
          (e.g. IotherPara(9)= 23 means 3 cuts at the feasible point of the box)
         = -1m, (only for AG7) use m cuts at '0' for the component i with max( width(X_i)) for i in I={j:X_j contains 0}
         = -2m, (only for AG7) use m cuts at '0' for the first possible component i
         = -3m, (only for AG7) use m cuts at '0' for the last possible component i
         = -4m, (only for AG7 and AG8) use m cuts at '0' for that component i which gives the most reduction in the quad. term b'(X'X)b  from the reference value  xols'(X'X)xols where
                b is the  truncated OLS solution after changing one variable to  0 at a time.
         = -5m, (only for AG7 and AG8) use m cuts at '0' for that component i which gives the most reduction in the quad term  b'(X'X)b  from the reference value xols'(X'X)xols where 
                b is the OLS solution in the reduced pDim-1 dim.
         = -6m, (only for AG7 and AG8) use m cuts at '0' for that component i which gives the least reduction in the quad. term b'(X'X)b  from the reference value  xols'(X'X)xols where
                b is the  truncated OLS solution after changing one variable to  0 at a time.
         = -7m, (NOT READY YET)(only for AG7,AG8,AG71,AG72) use m cuts at '0' for that component i given by the criteria in forward stepwise for a variable to enter in the support as given in 
                                           Hastie et.al. 2020 BSS, FS and Lasso comparison.
         (if = -1,-2,... then no. of cuts = 1 )

             
IotherPara(10)= define the tolerance less than which we will consider the entry in the final Xstar to be 0. Used to count no. of zeros in Xstar for output file.
IotherPara(11)= tol (if IotherPara(7)=1, DELAY F call)define the tolerance for the width of the box,once every box in the list become smaller than this tolerance then start calling inclusion function F.
          
IotherPara(12)= define the tolerance to be treated as 0 for h(x)=0 , while finding the feasible point.
IotherPara(13) (Suppressed right now)option to control the initial feasible point for AG7
          = 1 means use feasiblity sampling as in IotherPara(22)
          = 2 means use AG3
             
IotherPara(14)=  4, (1 Quadrant) Normalize using beta from projected gradient descent method as described in "Bertsimas et.al.(2016)BSS Modern optim.lens"
          =   3, (1 Quadrant) Normalize using xRelaxedOpt, will use IotherPara(2) as defined above
          = 2, (1 Quadrant) Normalize using MinMax1d approach sec. 2.3.2 in "Bertsimas et.al.(2016)BSS Modern optim.lens"
          = 1, (1 Quadrant) Do Not Normalize the data. 
          =  -4, (All Quadrant) Normalize using beta from projected gradient descent method as described in "Bertsimas et.al.(2016)BSS Modern optim.lens"
          = -3, (All Quadrant) Normalize using xRelaxedOpt, will use IotherPara(2) as defined above
          = -2, (All Quadrant) Normalize using MinMax1d approach sec. 2.3.2 in "Bertsimas et.al.(2016)BSS Modern optim.lens"
          = -1, (All Quadrant) Do Not Normalize the data.

IotherPara(15)= 0 , run constrained regression problem using equality constraint
          = 1 , run Penalized regression problem using semi-interval 1 approach   8 Feb21
          = 2 , run using the second approach, invloving 2 different lists   9 Feb21
          = 3 , 25 Feb21. Solve quadratic problem subject to inequality constraint. Sampling option given in
                IotherPara(8) will be used. For this option set the RHS parameter for the ineq
                constraint as tmax, see IotherPara(18) below
                Also, it only works for IotherPara(14)=2;
         = 4 , run DIMENSION REDUCTION method using nDcol given 8 March21
         = 5 , penalized regression with using unconstrained interval algo. 13 March21
         = 6 , Phase 1: use penalized regss on indcols, Phase 2: use penalty for both indcols and depcols 29March21
         = 7 ,(BEST SUB.SELECT. using machine eps) ONLY WORKS FOR Ah=0. Uses L0 + sd/ols  depending on the IotherPara(21) (11 Aug21)
         = 8 , (BSS with tmax solution path using AG7) 21 Sep21
         = 9 , (NON INTERVAL ALGO.) Quad. Min. Fun algo. for BSS  06 Oct21 
         = 71, (NON INTERVAL ALGO.) AG7 16 Oct21
         = 20, IBB from Somol et al (2004)Fast B&B for optimal feature selection  10 Nov21
         = 21, BBPP
         = 22, FBB
         = 23, IBB+  which is IBB with the 'minimum solution tree' approach and same for options below
         = 24, BBPP+
         = 25, FBB+
         = 30, bestsubset package in R from Hastie et.al.2020 BSS,FS and LASSO comparison
         = 31, forward stepwise from bestsubset package
         = 32, relaxed lasso 

         NOTE : negative of some options above will run the Itvl algo. which uses machine eps to take care
         of the discontinuity at '0' in L0 and pert. L0 penalty.
         
IotherPara(16) (NO LONGER IN USE) = eta , the parameter used in the computation of t_max for the rhs parameter bound T=[t_min , t_max]
      
IotherPara(17)= (ONLY FOR AG8) define the stepTm value to be used (21 Sep21)
         
IotherPara(18) For plotting (NOTE : use the same option 0,1,2 or 3 for all the sets you are running,i.e. keep IotherPara(18) same for all rows)
          = 0 , do not save any plots for the results.
          = 1 , level 1, save the alpha solution path plot and the lambda/tmax solution path plot,showing shrinkage of the coeff.
          = 2 , level 1, save the plots in option 1 above and along with that save bar plots comparing fbest, t.e. for different sets with differ alpha values for each eg.
          = 3 , level 1 and 2, save all the plots as in option 1 and 2 above and along with those save
                 plots for each examples, comparing interval para. and fbest., t.e at level 2
                      
IotherPara(19) = [ Options to find lb F(X) ]
             5, evaluate F(X) as in option 4 but by evaluating only one bound (NOTE:FOR IotherPara(15)=1 due to implementation reasons options 5 and 2 has been used together)
             6, evaluate F(X) using the EXACT lower bound of (0.5*X^2+JX) by adding and subtracting the components of X_j0 where j0 is the component of X box that we are bisecting
             4, evaluate F(X) by adding and subtracting the components of X_j0 where j0 is the component of X box that we are bisecting.
             2, evaluate only the lower bound of natural interval extension F(X) added on 31 Mar21
             0, evaluate F(X) using the natural interval extension of f
             1, evaluate F(X) using the mean value form F(X)= f(c)+(X-c)' gradF(X),  where c is any point in X, gradF is the inclusion fn. of f'
             3, evaluate only the lower bound of mean value form F(X)
             8, evaluate lb F(X) using quadratic min. algo. given by the first option selected in
                   QuadMinFunPara below
             7, (NOT READY YET)evaluate F(X) using 2 cuts and finding exact lb of quad. term in 2D  (12 May21)
             9, (14 Jan22) evaulate lb F(X) using QRD approach 
             10,(only implemented for AG72 so far) (28 Feb22) evaluate lb F(X) and refinement of the feasible pt. using the QRD approach
Conclusion: cpu5 < cpu6 <cpu 4 < cpu 2 < cpu 0      and  cpu 3 < cpu 1  and cpu 0 ~ cpu 1
        
IotherPara(20) = n ,(ONLY FOR SEMI-INTERVAL 2) where n is defined as the factor times IstopCondPara(7) which gives the new Width tolerance when we have to partition X box again
            i.e.  after putting the initial boxes in the list L_X with a width tol defined as
            IstopCondPara(7), we want to decrease the width tolerance for the second call of adding more
            boxes in the list L_X. So, new width will be  n*IstopCondPara(7)
             
IotherPara(21) (to select the box for refinement and max. no. of iter for refinement) 
          NOTE: only for L0 constraint regss. such as AG3 and AG7, this option will not control the max. iter. for S-D call at the end of AG3. (that is controlled by iPara(2) below)
           = 0 , no refinement
           = 1n ,where integer n > 0 ,use the WORKING box to improve the feasible pt. without changing the support.(26 July21)                                                 
           = 2n ,where integer n > 0 , use the ORIGINAL box to improve the feasible pt. without changing the support. (2 Aug21)      
           = 3n , where integer n>0 , use the WHOLE BDD. box. 
           = 4n , where integer n>0 , use the WHOLE UNBDD. box.
                
IotherPara(22) NOTE: only for L0 constraint regss. such as AG3 and AG7
           = 1 , use original feasibility sampling for exact L0 and pert. L0
           = 2 , use improved(correlation in A based) feasibility sampling for exact L0 and pert. L0 (26 July21)
           = 3 , use Modified feasibility sampling 1  (23 Aug21)
           = 4 , use Modified feasibility sampling 2  (24 Aug21)
           = 5 , use Modified feasibility sampling 3  (8 Sep21)
           = 6 , use Modified feasibility sampling (12 Sep21)

IotherPara(23) will be used when IotherPara(21)/=-1, select the quadratic min fun options for the refinement for AG7, AG8 and AG9
           = 1 means use minq8 algo. 
           = 2 means use SD2 algo.
           = 3 means use LLS algo. 
           = 4 means use Hbd algo.
           = 5 means use SequMiniQuadFun_v11 
           = 6 means use SequMiniQuadFun_v12 
           = 7 means use quadprog 
           = 8 means use SD 
           = 9 means use Coordinate descent   
           = 10 means use Conjugate gradient 
           = 11 means use Parallel Tangent
           = 12 means ols solution
           = 13 means use SD over whole space 12/14/2021
           = 14 means use PT over whole space 12/14/2021
           = 15 means use unconstr. adaptive Method 1  2/21/22
           = 16 means use unconstr. adaptive Method 2
           = 17 means use unconstr. adaptive Method 3
           = 18 means use CD over whole space 03/16/2022
           = 19 means use unconstr. adaptive Method 4,    04/03/22
             
IotherPara(24) = tolerance eps0 for the perturbed penalty IotherPara(21)==3 only, will not get used for other penalty functions.
%}
                
                
                           
        
        % the following 3 parameters will be used only if IstopCondPara(4 or 5 or 6) is negative as defined above.
        intermFbestIter=[4  10 50 100]    % define the intermediate values corresponding to IstopCondpara(4), strictly less than IstopCondpara(4) above
        intermMaxIter=[100 200 500]      % define the intermediate values corresponding to IstopCondpara(5), strictly less than IstopCondpara(5) above
        intermCPU=[0.2 1]      % define the intermediate values corresponding to IstopCondpara(6), strictly less than IstopCondpara(6) above       
                       
                    
        % define which quadratic min. approach we want to use below---------------------------------------------------
        %               1   2   3    4    5    6    7   8   9  10  11   12   13  14   15   16   17  18  19   20  ||  21   22  23
        QuadMinFunPara=[0   0   0    0    0    0    0   0   0   0   0    0    0   0    0    0    0   0   0    1       0    0   0;
            % 1st row algo. for xRelaxedOpt=========================================================================================
                        1   0   0    0    0    0    0   0   0   0   0    0    0   0    0    0    0   0   0    0       0    0   0;
                        0   1   0    0    0    0    0   0   0   0   0    0    0   0    0    0    0   0   0    0       0    0   0;
                        0   0   1    0    0    0    0   0   0   0   0    0    0   0    0    0    0   0   0    0       0    0   0;
                        0   0   0    1    0    0    0   0   0   0   0    0    0   0    0    0    0   0   0    0       0    0   0;
                        0   0   0    0    1    0    0   0   0   0   0    0    0   0    0    0    0   0   0    0       0    0   0;
                        0   0   0    0    0    1    0   0   0   0   0    0    0   0    0    0    0   0   0    0       0    0   0;
                        0   0   0    0    0    0    1   0   0   0   0    0    0   0    0    0    0   0   0    0       0    0   0;
                        0   0   0    0    0    0    0   1   0   0   0    0    0   0    0    0    0   0   0    0       0    0   0;
                        0   0   0    0    0    0    0   0   1   0   0    0    0   0    0    0    0   0   0    0       0    0   0;
                        0   0   0    0    0    0    0   0   0   1   0    0    0   0    0    0    0   0   0    0       0    0   0;
                        0   0   0    0    0    0    0   0   0   0   1    0    0   0    0    0    0   0   0    0       0    0   0;
                        0   0   0    0    0    0    0   0   0   0   0    1    0   0    0    0    0   0   0    0       0    0   0;
                        0   0   0    0    0    0    0   0   0   0   0    0    1   0    0    0    0   0   0    0       0    0   0;
                        0   0   0    0    0    0    0   0   0   0   0    0    0   1    0    0    0   0   0    0       0    0   0;
                        0   0   0    0    0    0    0   0   0   0   0    0    0   0    1    0    0   0   0    0       0    0   0;  %since 9/29/22
                        0   0   0    0    0    0    0   0   0   0   0    0    0   0    0    1    0   0   0    0       0    0   0;
                        0   0   0    0    0    0    0   0   0   0   0    0    0   0    0    0    1   0   0    0       0    0   0;
                        0   0   0    0    0    0    0   0   0   0   0    0    0   0    0    0    0   1   0    0       0    0   0;
                        0   0   0    0    0    0    0   0   0   0   0    0    0   0    0    0    0   0   1    0       0    0   0;
                        0   0   0    0    0    0    0   0   0   0   0    0    0   0    0    0    0   0   0    1       0    0   0;
                        0   0   0    0    0    0    0   0   0   0   0    0    0   1    0    0    0   0   0    1       0    0   0;  %1st custom choice
                       ]
              
        % NOTE:the first active algo. in a row is the one used to get lb F(X), if IotherPara(19)=8 above
                    
        %{
        QuadMinFunPara(1) =1 means use minq8 algo. to get xRelaxedOpt
        QuadMinFunPara(2) =1 means use SD2 algo. to get xRelaxedOpt
        QuadMinFunPara(3) =1 means use LLS algo. to get xRelaxedOpt
        QuadMinFunPara(4) =1 means use Hbd algo. to get xRelaxedOpt
        QuadMinFunPara(5) =1 means use SequMiniQuadFun_v11 to get xRelaxedOpt
        QuadMinFunPara(6) =1 means use SequMiniQuadFun_v12 to get xRelaxedOpt
        QuadMinFunPara(7) =1 means use quadprog to get xRelaxedOpt
        QuadMinFunPara(8) =1 means use SD to get xRelaxedOpt
        QuadMinFunPara(9) =1 means use Coordinate descent to get xRelaxedOpt                
        QuadMinFunPara(10) =1 means use Conjugate gradient to get xRelaxedOpt
        QuadMinFunPara(11) =1 means use Parallel Tangent to get xRelaxedOpt 
        QuadMinFunPara(12) = 1 means use OLS i.e. the first order optimality condition with MATLAB 'mldivide' which uses QR decomposition to solve Ax=b system                 
        QuadMinFunPara(13) = 1 means use unconstrained SD  12/14/2021
        QuadMinFunPara(14) = 1 means use unconstrained PT 12/14/2021
        QuadMinFunPara(15) = 1 means use unconstr. adaptive Method 1  2/21/22
        QuadMinFunPara(16) = 1 means use unconstr. adaptive Method 2
        QuadMinFunPara(17) = 1 means use unconstr. adaptive Method 3
        QuadMinFunPara(18) = 1 means use unconstr. CD
        QuadMinFunPara(19) = 1 means use unconstr. adaptive Method 4
        QuadMinFunPara(20) = 1 means use unconstr. CG
        *QuadMinFunPara(21) =1 means include LASSO results in the output
        *QuadMinFunPara(22) =1 means include Ridge results in the output
        *QuadMinFunPara(23) =1 means include OLS results in the output  
        QuadMinFunPara(24) = 1 means use unconstr. minq8 algorithm 07/19/23 
        *controls whether to include these algo. for comparison in the final output or not
        NOTE: if adding more quad. min. algo. i.e. extending 1:14 then change , runLLS line 207 indexOfIntvalSol  and runAnEgLLS line 63 colIndex, BB_AG71 line 1130,             
        %}
                        
                        %
                        % Define input parameters for Quadratic minimization function below-------------------------
                        rPara=[1e-14      %rPara(1)=10^(-14)=eps
                               1e-9       %rPara(2)=10^(-9)=epsPivot used in matrix factorization
                               1e-9       %rPara(3)=10^(-6)=epsSD   10^(-3) is not good
                               0.5*1e-3   %rPara(4)=0.5*10^(-3), MS9/9/19 relative error of f-value
                               1e-9       % used in MiniQuadFunUnconstAdapt1.m and MiniQuadFunUnconstAdapt3 , if -df < rPara(5)*|fx0| , not acceptable progress
                               0.5        % used in MiniQuadFunUnconstAdapt2.m  , if  -gdx*rPara(6) > -df , not acceptable progress  
                               1e18       % maxDx, max value allowed for a step size in uSD,uPT,uCD
                               1e5        %rparaFac0=1e+5, since 8/16/23
                               1e-1]      % iparaFac0=1e-1    
                            
                        %        1    2     3     4      5    6       7     8     9
                        iPara=[  0,   0,    0,    0,     0,   0,      0,    0     0        %iPara(1)=option used in LDLtOfSymMatrix21.m  for prearrangement
                               1000,1000, 1000, 1000, 1000, 1000,   1000, 1000  500      %iPara(2)=max # of iterations for S-D alg
                                 0,   0,    0,    0,     0,   0,      0,    0     0        %iPara(3)= MS:9/8/19 toDebug
                                21,  21,   21,   21,    21,  21,     21,   21    21        %iPara(4)= MS:9/9/19   #-1 is max length of printed out vector/matrix
                                 0,   0,    0,    0,     0,   0,      1,    2     3        %iPara(5)=0 means use full recursion, iPara(5)=1 means use 1-recursion + steepest descent to solve instead of full recursion
                                 0,   0,    0,    0,     0,   0,      0,    0     0        %iPara(6)=1 means use SD2 in the framework of v41
                                 0,   0,    0,    0,     0,   0,      0,    0     0        %iPara(7)= MS:9/15/19 as option to pre-permutation of variables based on domain size
                                 1,   1,    1,    1,     1,   1,      1,    1     1        %iPara(8)= MS:10/20/19 as option to check suff cond of optimality for earlier exit (1=y, 0=n)
                                 0,   0,    0,    0 ,    0,   0,      0,    0     0        %iPara(9)= MS:10/31/19 as option to use 0=LDU (original/default), 1=QDQ' ([Q,D]=eig(A); Q=o.n. e.v. for symmetric case, enough!)
                                 0,   0,    0,    0,     0,   0,      0,    0     0        %iPara(10)=MS:10/31/19 as option to use different choices of Xstart[]/Ystart[]: 0=Box(original/default), 1=space-Box(hybrid)
                                 5,   5,    5,    5,     5,   5,      5,    5     5        %*  iPara(11)=VS:11/7/19, select bdy face 0=original, 1,2,3 for 3 versions of min dist from a line, 4 for min dist from a face approach 11/20
                                 0,   0,    0,    0,     0,   0,      0,    0     0        %iPara(12)=tie breaker flag to select only one bdy face.
                                 0,   1,    2,    3,     0,   0,      0,    0     0];      %iPara(13)=as a flag controlling a pre-search procedure! It applies to: Sd,AllSd,PsdSd,
                        %                                        *     *     * (reserved last 3 cols)
                     
                        % end of input parameters required for recursive min function--------------------------------
                        %}
                        % end of defining input parameters-----

       algoParafile=mfilename('fullpath'); % (No need to change) to save a copy of this file in the output folder.

end