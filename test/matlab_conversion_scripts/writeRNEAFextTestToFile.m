
function writeRNEAFextTestToFile(N, bf, filename)

outfile = extractBefore(filename, strlength(filename)-2+1) + "out";
fid = fopen(filename, 'a');
fprintf(fid, '%gf64 \n', bf);
fclose(fid);
q = round(rand(N,1), 2);
qd = round(rand(N,1), 2);
qdd = round(rand(N,1), 2);

fext = randCellArrays(N, 6)

model = autoTree(N, bf, 0, 1)

tau = ID(model, q, qd, qdd, fext)

writeArrayToFile(q, 'f64', filename)
writeArrayToFile(qd, 'f64', filename)
writeArrayToFile(qdd, 'f64', filename)
CellOf1dArraysToFile(fext, 'f64', filename)
writeArrayToFile(tau, 'f64', outfile)
end