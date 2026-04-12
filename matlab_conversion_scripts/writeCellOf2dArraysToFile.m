function writeCellOf2dArraysToFile(A, t, filename)
% writeArrayToFile Append array A of doubles to filename as [a, b, c] on a new line.
% The doubles array's formatting correspond to a f64 futhark array.
% Usage:
%   writeArrayToFile(A, 'myfile.txt') % appends to specified file

fid = fopen(filename, 'a');
if fid == -1
    error('Cannot open file: %s', filename);
end

% Ensure start on a new line
% fprintf(fid, '\n');

% Handle empty array
if isempty(A)
    fprintf(fid, '[]\n');
    fclose(fid);
    return
end

[m,n] = size(A{1})
fprintf(fid, '[');
for i = 1:length(A)
    fprintf(fid, '[');
    for j = 1:m
        fprintf(fid, ['[%g' t], A{i}(j,1));
        for k = 2:n
            fprintf(fid, [', %g' t], A{i}(j,k));
        end
        if j == m
            fprintf(fid, ']');
        else 
            fprintf(fid, '], ');
        end
    end
    if i == length(A)
        fprintf(fid, ']');
    else
        fprintf(fid, '], ');
    end
end
fprintf(fid, ']\n');
fclose(fid);
end