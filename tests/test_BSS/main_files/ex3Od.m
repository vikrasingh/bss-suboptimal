
function  [fbest_all_inp_files,cpu_all_inp_files,niter_all_inp_files,stopflag_all_inp_files]=ex3Od
% 03/20/24, give the input file pairs to run below
if isunix % if running on cluster
   addpath(genpath('/home/x-vsingh7/bss_suboptimal/'));  
   addpath(genpath('/home/x-vsingh7/MinimizeQuadFunOverBox_c0_2023-08-21/'));
end
% Each row of InputFiles below represents a pair of algo. and example input file to run a test
% Algo. input file first and example input file second, order of input files in a row matters
para_struct.mainfile=mfilename('fullpath');
para_struct.to_debug=0;  % =0 no intermediate output, =1 some intermediate output, =2 the most extensive debugging
para_struct.id='Ex3Od';
InputFiles={
 %   
  %1 OD, snr 0.05  , 32 egs for p 20 to 5k for Sn 0p05
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran4Ex20t80Od2TmCcRh0p8Sn0p05'    
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran2Ex200a300Od2TmCcRh0p8Sn0p05'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran2Ex400a500Od2TmCcRh0p8Sn0p05'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex800Od2TmCcRh0p8Sn0p05v2'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex1000Od2TmCcRh0p8Sn0p05v2'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex1500Od2TmCcRh0p8Sn0p05v2'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex2000Od2TmCcRh0p8Sn0p05v2'
  
  %2 OD,  snr 0.5
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran4Ex20t80Od2TmCcRh0p8Sn0p5'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran2Ex200a300Od2TmCcRh0p8Sn0p5' 
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran2Ex400a500Od2TmCcRh0p8Sn0p5'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex800Od2TmCcRh0p8Sn0p5v2'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex1000Od2TmCcRh0p8Sn0p5v2'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex1500Od2TmCcRh0p8Sn0p5v2'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex2000Od2TmCcRh0p8Sn0p5v2'
  
  %3 OD, snr 1
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran4Ex20t80Od2TmCcRh0p8Sn1'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran2Ex200a300Od2TmCcRh0p8Sn1'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran2Ex400a500Od2TmCcRh0p8Sn1'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex800Od2TmCcRh0p8Sn1v2'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex1000Od2TmCcRh0p8Sn1v2'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex1500Od2TmCcRh0p8Sn1v2'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex2000Od2TmCcRh0p8Sn1v2'
  
 %4 OD, snr 5 
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran4Ex20t80Od2TmCcRh0p8Sn5'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran2Ex200a300Od2TmCcRh0p8Sn5'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran2Ex400a500Od2TmCcRh0p8Sn5'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex800Od2TmCcRh0p8Sn5v2'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex1000Od2TmCcRh0p8Sn5v2'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex1500Od2TmCcRh0p8Sn5v2'
  'Mix99b3e9c9d3f3g3dV9Od',  'Ran1Ex2000Od2TmCcRh0p8Sn5v2'
  
           };   
% ppDimRunCall(InputFiles,paraStruct) % call parallel pair run for LLS 
% pPairRunCall(InputFiles,para_struct);
% pSetRunCall(InputFiles,paraStruct);
[fbest_all_inp_files,cpu_all_inp_files,niter_all_inp_files,stopflag_all_inp_files]=get_ex_data(InputFiles,para_struct);
end % end 




