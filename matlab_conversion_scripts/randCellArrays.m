
function C = randCellArrays(N, M)
% randCellArrays Return 1xN cell of double arrays with 2 decimal places.
%   C = randCellArrays(N, M) where
%     N - number of cells (scalar)
%     M - scalar length for each cell OR a 1xN vector of lengths
%
% Example:
%   C = randCellArrays(4,3)     % 1x4 cell, each cell is 1x3 double
%   C = randCellArrays(3,[1 2 4]) % cells lengths 1,2,4

if nargin < 2
    error('Both N and M must be provided');
end

if isscalar(M)
    lengths = repmat(M, 1, N);
else
    lengths = M(:).';      % ensure row
    if numel(lengths) ~= N
        error('If M is a vector, it must have length N.');
    end
end

C = cell(1,N);
for k = 1:N
    % generate and round to 2 decimals
    C{k} = round(rand(1, lengths(k)), 2)';
end
end