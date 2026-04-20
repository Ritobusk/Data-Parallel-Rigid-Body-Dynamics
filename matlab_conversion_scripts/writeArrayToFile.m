function writeArrayToFile(A, t, filename)
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

% Write as [a, b, c]
if t == "f64"
    fprintf(fid, ['[%.4f' t], A(1));
else
    fprintf(fid, ['[%g' t], A(1));
end

for k = 2:numel(A)
    if t == "f64"
        fprintf(fid, [', %.4f' t], A(k));
    else
        fprintf(fid, [', %g' t], A(k));
    end
    
end
fprintf(fid, ']\n');

fclose(fid);
end