function CellOf1dArraysToFile(A, t, filename)

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

m = length(A{1});

fprintf(fid, '[');
for i = 1:length(A)
    fprintf(fid, ['[%g' t], A{i}(1));
    for k = 2:m
        fprintf(fid, [', %g' t], A{i}(k));
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