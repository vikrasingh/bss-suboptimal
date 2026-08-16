function   [truebPoolFilePath,giveValuesOf_b]=truebPool(pDim)
% define trueb manually for each given pDim
    
    truebPoolFilePath=mfilename('fullpath'); % (No need to change) to save a copy of this file in the output folder.
    if pDim==2
        giveValuesOf_b=[-3  2
                       -0.5  4]';
    
    elseif pDim==3, giveValuesOf_b=[ -2   3   5;   %1a
            0  -2   3;   %2a
            0 -0.2  3;   %3a
            0 -0.02 3;   %4a
            0   0  -3;   %5a
            0   0  -0.3; %6a
            0   0  -0.03;  %7a
            -5   2   4;   %1b
            0  -5   4;   %2b
            0 -0.5   4;  %3b
            0 -0.05  4;  %4b
            0   0   4;   %5b
            0   0  0.4;  %6b
            0   0  0.04]';  %7b  [pDim+nDcol]   Quick 3-dim
    
    elseif pDim==5,   giveValuesOf_b=[3     -4    -2    4   -3;
            0.2  0.4   0.6  0.8    1;
            5    4.8   4.5  4.2    4;
            3     3     3     3    3;
            1     0.8   0.7   1.2  1.5;
            3      1     0     0    0;
            0      0     2     0    2]';
    elseif pDim==6,  giveValuesOf_b=[0 0 0 -2 3 5]';
    
    elseif pDim==12
        giveValuesOf_b=[ -2      3      5    -2      3     5  -2     3     5    -2     3     5;   %1a
            0     -2      3     0     -2     3   0    -2     3     0    -2     3;   %2a
            0      0     -3     0      0    -3   0     0    -3     0     0    -3;   %3a
            0   -0.2      3     0   -0.2     3   0  -0.2     3     0  -0.2     3;   %4a
            0      0   -0.3     0      0  -0.3   0     0   -0.3    0     0  -0.3;  %5a
            0  -0.02      3     0  -0.02     3   0  -0.02    3     0 -0.02     3;   %6a
            0      0  -0.03     0      0  -0.03  0     0  -0.03    0     0  -0.03;    %7a
            -5      2      4    -5      2     4  -5     2     4    -5     2     4;   %1b
            0     -5      4     0     -5     4   0    -5     4     0    -5     4;   %2b
            0      0      4     0      0     4   0     0     4     0     0     4;   %3b
            0   -0.5      4     0   -0.5     4   0   -0.5    4     0  -0.5     4;  %4b
            0      0    0.4     0      0   0.4   0     0    0.4    0     0    0.4;  %5b
            0  -0.05      4     0  -0.05     4   0  -0.05    4     0 -0.05     4;  %6b
            0      0   0.04     0      0  0.04   0     0   0.04    0     0    0.04;  %7b  [pDim+nDcol]  Quick 12-dim
            5      5     5      5      5     0    0     0     0     0     0      0;  % 15
            3      0     3      0      3     0    3     0     0     0     0      0;  % 16
            0      0     3      0      3     0    3     0     3     0     0      0;
            0      0     0      0      3     0    3     0     3     0     3      0;
            0      0     0      0      0     0    0     0     3     3     3      3]';
    
    elseif pDim==20
        giveValuesOf_b=[1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
                        1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0]';
    elseif pDim==30
        giveValuesOf_b= [3  2  1     0.2  0    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0.2  1  2  3;  %1
            2  0  0     0    2    0  0  0  2  0  0  0  2  0  0  0  2  0  0  0  2  0  0  0  2  0  0  0  2  0;    %2
            2  0  1.75  0    1.5  0  1.25  0  1  0  0.75  0  0.5  0  0.25  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %3
            5  0  3     0    2    0  4  0  6  0  0  3  0  0  2  0  0  1  0  0  3  0  0  5  0  0  0  0  0  0]';  %4
    
    end


end