function RNEA_bench(Ns, bf, filename_ending)
% You do not write the joint_types to file. That should be done in futhark.

skew = 0; % Initialize skew variable
taper = 1; % Initialize taper variable
outfile = "bf" + num2str(bf) + "_results.out";

fid_out = fopen(outfile, 'a');

for nb = Ns
    disp(nb)
    model = autoTree(nb, bf, skew, taper);
    input_file = "N" + num2str(nb) + "_bf" + num2str(bf) + "_" + filename_ending
    fid = fopen(input_file, 'a');
    fprintf(fid, '%di64 \n', nb);
    fprintf(fid, '%gf64 \n', bf);
    fprintf(fid, '%gf64 \n', skew);
    fprintf(fid, '%gf64 \n', taper);
    fclose(fid);
    %writeArrayToFile(model.parent, 'i64', input_file);
    %writeCellOf2dArraysToFile(model.I, 'f64', input_file);
    %writeCellOf2dArraysToFile(model.Xtree, 'f64', input_file);

    q = round(rand(nb,1), 2);
    qd = round(rand(nb,1), 2);
    qdd = round(rand(nb,1), 2);

    writeArrayToFile(q, 'f64', input_file);
    writeArrayToFile(qd, 'f64', input_file);
    writeArrayToFile(qdd, 'f64', input_file);

    
    M = 10;
    t = zeros(1,M);
    for k = 1:M
        t(k) = timeit(@() ID(model, q, qd, qdd));  
    end
    if fid_out == -1
        error('Cannot open file: %s', filename);
    end
    fprintf(fid_out, 'Nb: %d, bf: %.6f Mean: %.6f s, Std: %.6f s\n', nb, bf, mean(t), std(t));
    

end
fclose(fid_out);
end