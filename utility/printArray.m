function printArray(A,varargin)
% printArray(A) default is %f format, printArray(A,'%d'),printArray(A,'%d',fileid)
% 9 July20, to print a vector or matrix using fprintf,
% printArray(A) print the array in [%f] form
% printArray(A,'%d') print the array in [%d] form
% printArray(A,'%d',filename) print the array in [%d] form in the filename text file
%   filename is the id after using fopen to open some text file  

[m,n]=size(A);
switch nargin
    case 1
        if m==0 && n==0
            fprintf('\n');
        elseif m==1 && n==1
            fprintf('[%f]\n',A);
        elseif m==1 && n>1 % A is a row vector
            string='[%f ';
            for i=2:(n-1)
                string=[string   '%f '];
            end
            string=[string   '%f]\n'];
            fprintf(string,A);
        elseif n==1 && m>1 % A is a column vector
            string='[%f ';
            for i=2:(m-1)
                string=[string   '%f '];
            end
            string=[string   '%f]\n'];
            fprintf(string,A);
        elseif m>1 && n>1  % A is an mxn matrix
            string='%f ';
            for i=2:(n-1)
                string=[string   '%f '];
            end
            string=[string  '%f;\n' ];
            fprintf('[');
            fprintf(string,A');
            fprintf('] \n');
        end
       return
    case 2
        if m==0 && n==0
            fprintf('\n');
        elseif m==1 && n==1
            fprintf(append('[',varargin{1},']\n'),A);
        elseif m==1 && n>1 % A is a row vector
            string=append('[',varargin{1},' ');
            for i=2:(n-1)
                string=[string   append(varargin{1},' ')];
            end
            string=[string   append(varargin{1},']\n') ];
            fprintf(string,A);
        elseif n==1 && m>1 % A is a column vector
            string=append('[',varargin{1},' ');
            for i=2:(m-1)
                string=[string   append(varargin{1},' ')];
            end
            string=[string   append(varargin{1},']\n') ];
            fprintf(string,A);
        elseif m>1 && n>1  % A is an mxn matrix
            string=append(varargin{1},' ');
            for i=2:(n-1)
                string=[string   append(varargin{1},' ')];
            end
            string=[string   append(varargin{1},';\n') ];
            fprintf('[')
            fprintf(string,A');
            fprintf('] \n');
        end
    case 3
        if m==0 && n==0
            fprintf(varargin{2},'\n');
        elseif m==1 && n==1
            fprintf(varargin{2},append('[',varargin{1},']\n'),A);
        elseif m==1 && n>1 % A is a row vector
            string=append('[',varargin{1},' ');
            for i=2:(n-1)
                string=[string   append(varargin{1},' ')];
            end
            string=[string   append(varargin{1},']\n') ];
            fprintf(varargin{2},string,A);
        elseif n==1 && m>1 % A is a column vector
            string=append('[',varargin{1},' ');
            for i=2:(m-1)
                string=[string   append(varargin{1},' ')];
            end
            string=[string   append(varargin{1},']\n') ];
            fprintf(varargin{2},string,A);
        elseif m>1 && n>1  % A is an mxn matrix
            string=append(varargin{1},' ');
            for i=2:(n-1)
                string=[string   append(varargin{1},' ')];
            end
            string=[string   append(varargin{1},';\n') ];
            fprintf(varargin{2},'[');
            fprintf(varargin{2},string,A');
            fprintf(varargin{2},'] \n');
        end
end

end
