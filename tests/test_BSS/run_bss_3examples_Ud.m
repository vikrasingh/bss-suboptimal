function run_bss_3examples_Ud
   if isunix
      addpath(genpath('/home/x-vsingh7/bss_suboptimal_29July26/')); 
      addpath(genpath('/home/x-vsingh7/MinimizeQuadFunOverBox_c0_2023-08-21/'));
   end
   fbest=cell(1,3); % collect fx for the three examples
   cpu_data=cell(1,3);
   num_iter=cell(1,3);
   stop_flag=cell(1,3);
   
   current_dir=pwd;
   % ex1Ud
   new_dir='ex1Ud_output';
   mkdir(new_dir);cd(new_dir)
   [fbest{1},cpu_data{1},num_iter{1},stop_flag{1}]=ex1Ud;
   cd(current_dir);
   
   % ex2Ud
   new_dir='ex2Ud_output';
   mkdir(new_dir);cd(new_dir)
   [fbest{2},cpu_data{2},num_iter{2},stop_flag{2}]=ex2Ud;
   cd(current_dir);
   
   % ex3Ud
   new_dir='ex3Ud_output';
   mkdir(new_dir);cd(new_dir)
   [fbest{3},cpu_data{3},num_iter{3},stop_flag{3}]=ex3Ud;
   cd(current_dir);
   
   % concatenate the data for all three examples
   fbest=cat(1,fbest{:});cpu_data=cat(1,cpu_data{:});num_iter=cat(1,num_iter{:});stop_flag=cat(1,stop_flag{:});
   
   % create box plots for rel gap % for all the three examples in each dim.
   % regimes
   upYLimVec={400,50,20};
   boxPlotsAllDimRegimes(fbest, upYLimVec, 'Ud')
   
   % create performance profiles of CPU times
   perProf(cpu_data,'\tau','P(r_{p,s} \leq \tau:s\in S )','OD examples',{'SFS1','SFS2','FS','SFFS','GA','DFO','ABESS','MIO'},[],'Od8Alg')
   
   perProf(cpu_data(:,[1 3 4 6 7]),'\tau','P(r_{p,s} \leq \tau:s\in S )','OD examples',{'SFS1','FS','SFFS','DFO','ABESS'}, [2 5 8], 'Od5Alg')


end