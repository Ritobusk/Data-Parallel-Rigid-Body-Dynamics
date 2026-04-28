function write2dArrayToFile(A, t, filename)
% writeArrayToFile Append array A of type t to filename as [a, b, c] on a new line.
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

[m,n] = size(A);

fprintf(fid, '[');

for j = 1:m
    % Print opening bracket for this row
    if t == "f64"
        fprintf(fid, '[%.4f', A(j,1));
    else
        fprintf(fid, '[%g', A(j,1));
    end
    fprintf(fid, t);

    % Print remaining elements of the row
    for k = 2:n
        if t == "f64"
            fprintf(fid, ', %.4f', A(j,k));
        else
            fprintf(fid, ', %g', A(j,k));
        end
        fprintf(fid, t);
    end

    % Close this row's bracket and separate rows with comma+space except last
    if j == m
        fprintf(fid, ']');
    else
        fprintf(fid, '], ');
    end
end
fprintf(fid, ']\n');

fclose(fid);
end