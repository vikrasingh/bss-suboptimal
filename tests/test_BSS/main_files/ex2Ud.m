
function  [fbest_all_inp_files,cpu_all_inp_files,niter_all_inp_files,stopflag_all_inp_files]=ex2Ud
% 03/20/24, give the input file pairs to run below

% Each row of InputFiles below represents a pair of algo. and example input file to run a test
% Algo. input file first and example input file second, order of input files in a row matters
para_struct.mainfile=mfilename('fullpath');
para_struct.toDebug=0;  % =0 no intermediate output, =1 some intermediate output, =2 the most extensive debugging
para_struct.id='Ex2Ud';
InputFiles={
 %5 UD, snr 0.05
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two4Ex20t80Ud2TmCcRh0p8Sn0p05'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two2Ex200a300Ud2TmCcRh0p8Sn0p05'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two2Ex400a500Ud2TmCcRh0p8Sn0p05'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex800Ud2TmCcRh0p8Sn0p05'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex1000Ud2TmCcRh0p8Sn0p05'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex1500Ud2TmCcRh0p8Sn0p05'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex2000Ud2TmCcRh0p8Sn0p05'
  
 %6 UD, snr 0.5
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two4Ex20t80Ud2TmCcRh0p8Sn0p5'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two2Ex200a300Ud2TmCcRh0p8Sn0p5'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two2Ex400a500Ud2TmCcRh0p8Sn0p5'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex800Ud2TmCcRh0p8Sn0p5'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex1000Ud2TmCcRh0p8Sn0p5'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex1500Ud2TmCcRh0p8Sn0p5'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex2000Ud2TmCcRh0p8Sn0p5'
 
 %7 UD, snr 1 
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two4Ex20t80Ud2TmCcRh0p8Sn1'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two2Ex200a300Ud2TmCcRh0p8Sn1'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two2Ex400a500Ud2TmCcRh0p8Sn1'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex800Ud2TmCcRh0p8Sn1'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex1000Ud2TmCcRh0p8Sn1'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex1500Ud2TmCcRh0p8Sn1'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex2000Ud2TmCcRh0p8Sn1'
  
 %8 UD, snr 5 
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two4Ex20t80Ud2TmCcRh0p8Sn5'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two2Ex200a300Ud2TmCcRh0p8Sn5'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two2Ex400a500Ud2TmCcRh0p8Sn5'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex800Ud2TmCcRh0p8Sn5'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex1000Ud2TmCcRh0p8Sn5'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex1500Ud2TmCcRh0p8Sn5'
  'Mix99b3e9c9d3f3g3dV9Ud',  'Two1Ex2000Ud2TmCcRh0p8Sn5'
  
           };   

[fbest_all_inp_files,cpu_all_inp_files,niter_all_inp_files,stopflag_all_inp_files]=get_ex_data(InputFiles,para_struct);

end % end 




