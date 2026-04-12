function writeAutoTreeToFile(N, bf, skew, taper, filename)

%outfile = filename(1:end-2) + 'out'
outfile = extractBefore(filename, strlength(filename)-2+1) + "out";
fidin = fopen(filename, 'a');
fprintf(fidin, '%gi64 \n', N);
fprintf(fidin, '%gf64 \n', bf);
fprintf(fidin, '%gf64 \n', skew);
fprintf(fidin, '%gf64 \n', taper);
fclose(fidin);

model = autoTree(N, bf, 0, 1);

writeArrayToFile(model.parent, 'i64', outfile);
writeCellOf2dArraysToFile(model.Xtree, 'f64', outfile);
writeCellOf2dArraysToFile(model.I, 'f64', outfile);
end