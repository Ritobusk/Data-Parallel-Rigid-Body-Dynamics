function writeAllRNEATestToFile(N, bf, filename)

outfile = extractBefore(filename, strlength(filename)-2+1) + "out";
fid = fopen(filename, 'a');
fprintf(fid, '%gf64 \n', bf);
fclose(fid);
q = round(rand(N,1), 2);
qd = round(rand(N,1), 2);
qdd = round(rand(N,1), 2);


writeArrayToFile(q, 'f64', filename)
writeArrayToFile(qd, 'f64', filename)
writeArrayToFile(qdd, 'f64', filename)

model = autoTree(N, bf, 0, 1)

a_grav = get_gravity(model);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd(i) + crm(v{i})*vJ;
  end
  f{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
end

writeCellOf2dArraysToFile(Xup, 'f64', outfile);

for i = model.NB:-1:1
  tau(i,1) = S{i}' * f{i};
  if model.parent(i) ~= 0
    f{model.parent(i)} = f{model.parent(i)} + Xup{i}'*f{i};
  end
end



writeArrayToFile(tau, 'f64', outfile);
end
