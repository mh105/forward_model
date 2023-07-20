function [ free_column, free_column1, free_column2, free_column3 ] = fix2freeindexing( fixed_column )
%Function used to convert the column indices of sources in a fixed
%orientation leadfield matrix to column indices of the same sources in a
%free orientation leadfield matrix.

free_column1 = (fixed_column-1)*3+1;
free_column2 = free_column1+1;
free_column3 = free_column1+2;

free_column = sort([free_column1, free_column2, free_column3], 'ascend');

end

