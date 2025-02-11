
function  p28pairBox7a2c3dOne24Ex20t2kUd2Tm_V9

if isunix
   addpath(genpath('/home/vsingh12/LLSpackage14Apr24/')); 
   addpath(genpath('/home/vsingh12/MinimizeQuadFunOverBox_c0_2023-08-21/'));
end
% Each row of InputFiles below represents a pair of algo. and example input file to run a test
% Algo. input file first and example input file second, order of input files in a row matters
paraStruct.mainfile=mfilename('fullpath');
paraStruct.toDebug=0;  % =0 no intermediate output, =1 some intermediate output, =2 the most extensive debugging
InputFiles={
 %{  
  %1 OD, snr 0.05
  'Box7a2c3dV9',  'One4Ex20t80Od2TmCcRh0p8Sn0p05'
  'Box7a2c3dV9',  'One2Ex200a300Od2TmCcRh0p8Sn0p05'
  'Box7a2c3dV9',  'One2Ex400a500Od2TmCcRh0p8Sn0p05'
  'Box7a2c3dV9',  'One1Ex800Od2TmCcRh0p8Sn0p05'
  'Box7a2c3dV9',  'One1Ex1000Od2TmCcRh0p8Sn0p05'
  'Box7a2c3dV9',  'One1Ex1500Od2TmCcRh0p8Sn0p05'
  'Box7a2c3dV9',  'One1Ex2000Od2TmCcRh0p8Sn0p05'
 %
  %2 OD,  snr 0.5
  'Box7a2c3dV9',  'One4Ex20t80Od2TmCcRh0p8Sn0p5'
  'Box7a2c3dV9',  'One2Ex200a300Od2TmCcRh0p8Sn0p5' 
  'Box7a2c3dV9',  'One2Ex400a500Od2TmCcRh0p8Sn0p5'
  'Box7a2c3dV9',  'One1Ex800Od2TmCcRh0p8Sn0p5'
  'Box7a2c3dV9',  'One1Ex1000Od2TmCcRh0p8Sn0p5'
  'Box7a2c3dV9',  'One1Ex1500Od2TmCcRh0p8Sn0p5'
  'Box7a2c3dV9',  'One1Ex2000Od2TmCcRh0p8Sn0p5'
  %
  %3 OD, snr 1
  'Box7a2c3dV9',  'One4Ex20t80Od2TmCcRh0p8Sn1'
  'Box7a2c3dV9',  'One2Ex200a300Od2TmCcRh0p8Sn1'
  'Box7a2c3dV9',  'One2Ex400a500Od2TmCcRh0p8Sn1'
  'Box7a2c3dV9',  'One1Ex800Od2TmCcRh0p8Sn1'
  'Box7a2c3dV9',  'One1Ex1000Od2TmCcRh0p8Sn1'
  'Box7a2c3dV9',  'One1Ex1500Od2TmCcRh0p8Sn1'
  'Box7a2c3dV9',  'One1Ex2000Od2TmCcRh0p8Sn1'

 %4 OD, snr 5 
  'Box7a2c3dV9',  'One4Ex20t80Od2TmCcRh0p8Sn5'
  'Box7a2c3dV9',  'One2Ex200a300Od2TmCcRh0p8Sn5'
  'Box7a2c3dV9',  'One2Ex400a500Od2TmCcRh0p8Sn5'
  'Box7a2c3dV9',  'One1Ex800Od2TmCcRh0p8Sn5'
  'Box7a2c3dV9',  'One1Ex1000Od2TmCcRh0p8Sn5'
  'Box7a2c3dV9',  'One1Ex1500Od2TmCcRh0p8Sn5'
  'Box7a2c3dV9',  'One1Ex2000Od2TmCcRh0p8Sn5'
 %}
 %5 UD, snr 0.05
  'Box7a2c3dV9',  'One4Ex20t80Ud2TmCcRh0p8Sn0p05'
  'Box7a2c3dV9',  'One2Ex200a300Ud2TmCcRh0p8Sn0p05'
  'Box7a2c3dV9',  'One2Ex400a500Ud2TmCcRh0p8Sn0p05'
  'Box7a2c3dV9',  'One1Ex800Ud2TmCcRh0p8Sn0p05'
  'Box7a2c3dV9',  'One1Ex1000Ud2TmCcRh0p8Sn0p05'
  'Box7a2c3dV9',  'One1Ex1500Ud2TmCcRh0p8Sn0p05'
  'Box7a2c3dV9',  'One1Ex2000Ud2TmCcRh0p8Sn0p05'
 %
 %6 UD, snr 0.5
  'Box7a2c3dV9',  'One4Ex20t80Ud2TmCcRh0p8Sn0p5'
  'Box7a2c3dV9',  'One2Ex200a300Ud2TmCcRh0p8Sn0p5'
  'Box7a2c3dV9',  'One2Ex400a500Ud2TmCcRh0p8Sn0p5'
  'Box7a2c3dV9',  'One1Ex800Ud2TmCcRh0p8Sn0p5'
  'Box7a2c3dV9',  'One1Ex1000Ud2TmCcRh0p8Sn0p5'
  'Box7a2c3dV9',  'One1Ex1500Ud2TmCcRh0p8Sn0p5'
  'Box7a2c3dV9',  'One1Ex2000Ud2TmCcRh0p8Sn0p5'
 %
 %7 UD, snr 1 
  'Box7a2c3dV9',  'One4Ex20t80Ud2TmCcRh0p8Sn1'
  'Box7a2c3dV9',  'One2Ex200a300Ud2TmCcRh0p8Sn1'
  'Box7a2c3dV9',  'One2Ex400a500Ud2TmCcRh0p8Sn1'
  'Box7a2c3dV9',  'One1Ex800Ud2TmCcRh0p8Sn1'
  'Box7a2c3dV9',  'One1Ex1000Ud2TmCcRh0p8Sn1'
  'Box7a2c3dV9',  'One1Ex1500Ud2TmCcRh0p8Sn1'
  'Box7a2c3dV9',  'One1Ex2000Ud2TmCcRh0p8Sn1'

 %8 UD, snr 5 
  'Box7a2c3dV9',  'One4Ex20t80Ud2TmCcRh0p8Sn5'
  'Box7a2c3dV9',  'One2Ex200a300Ud2TmCcRh0p8Sn5'
  'Box7a2c3dV9',  'One2Ex400a500Ud2TmCcRh0p8Sn5'
  'Box7a2c3dV9',  'One1Ex800Ud2TmCcRh0p8Sn5'
  'Box7a2c3dV9',  'One1Ex1000Ud2TmCcRh0p8Sn5'
  'Box7a2c3dV9',  'One1Ex1500Ud2TmCcRh0p8Sn5'
  'Box7a2c3dV9',  'One1Ex2000Ud2TmCcRh0p8Sn5'
 %}
           };   
% ppDimRunCall(InputFiles,paraStruct);  % call parallel pair run for LLS 
pPairRunCall(InputFiles,paraStruct);
% serialRunCall(InputFiles,paraStruct);
% pSetRunCall(InputFiles,paraStruct);

end % end 




