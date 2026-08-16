function inprint(x,s,fmt)
% x is the vector/matrix that we want to print without breaking it into multiple rows in command window
% s is the string of the variable name
% fmt is data format like %f, %i , %1.3f ... that we want to use to print entries given in the array x

[~,n]=size(x);  % order of x
fprintf(s);fprintf('=');
fprintf(1, [repmat(append('  ',fmt), 1, n) '\n'], x');

    
    
end