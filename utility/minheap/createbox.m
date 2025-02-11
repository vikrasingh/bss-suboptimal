classdef createbox < handle
   % A data class to store information about a multi-dimensional box
   % VS 28Dec23, same as boxdata.m but with N1FLAG  as a node attribute as well 
   properties
      BOX;    % array of unsigned integers uint8
      LBFBOX; % floating point value to save lb f
      LBFPT;  % floating point array to save the point x, such that f(x)=lb f(X)
      N2FLAG; % no. of 2 flags in the box
      N0FLAG; % no. of 0 flags in the box
      N1FLAG; % no. of 1 flags in the box
      FLAG;   % only used for ag76,ag78 to keep track of whether the box is marked for deletion for a particular tmax value or not
      sumFLAG; % only used for ag76,ag78 to store sum of the entries in FLAG
   end

   methods
      function data = createbox(BOX,LBFBOX,LBFPT,N2FLAG,N0FLAG,N1FLAG,FLAG,sumFLAG)
         % Construct a data object
         if nargin > 0
            data.BOX = BOX;
            data.LBFBOX=LBFBOX;
            data.LBFPT=LBFPT;
            data.N2FLAG=N2FLAG;
            data.N0FLAG=N0FLAG;
            data.N1FLAG=N1FLAG;
            data.FLAG=FLAG;
            data.sumFLAG=sumFLAG;
         end
      end
   end
end