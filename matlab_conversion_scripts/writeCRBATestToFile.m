
function writeCRBATestToFile(N, bf, filename)

outfile = extractBefore(filename, strlength(filename)-2+1) + "out";
fid = fopen(filename, 'a');
fprintf(fid, '%gf64 \n', bf);
fclose(fid);
q = round(rand(N,1), 2);
qd = round(rand(N,1), 2);

model = autoTree(N, bf, 0, 1);

[H, C] = HandC(model, q, qd);

writeArrayToFile(q, 'f64', filename)
writeArrayToFile(qd, 'f64', filename)
writeArrayToFile(C, 'f64', outfile)
write2dArrayToFile(H, 'f64', outfile)
end