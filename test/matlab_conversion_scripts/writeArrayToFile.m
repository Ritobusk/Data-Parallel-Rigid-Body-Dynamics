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
fprintf(fid, ['[%g' t], A(1));
for k = 2:numel(A)
    fprintf(fid, [', %g' t], A(k));
end
fprintf(fid, ']\n');

fclose(fid);
end