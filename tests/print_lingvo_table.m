function print_lingvo_table()
  mynorm = 2;
  olddigits = mp.Digits(34);
  epsilons = [mp('0.01'), mp('0.1'), mp('1')];
  filename = 'lingvo_matrices_table.dat';
  fid = fopen(filename, 'w');

  fprintf(fid, '\\begin{tabularx}{\\textwidth}{cccccXccccc}\n\\toprule\n');
  fprintf(fid, '& $n$ & $\\lambda_{\\min}$ & $\\lambda_{\\max}$ & $r$\n');
  fprintf(fid, '    & & & $\\mytol_1$ & $t_1$  &  $\\mytol_2$  &  $t_2$ \\\\\n');
  fprintf(fid, '\\midrule\n');

  for i = 1:3
    filename = sprintf('data/test_mat_%d.mat', i);
    load(filename);
    B = data;
    [n, ~] = size(B);

    % The unnecessary computation of the eigensystem is to obtain the values that
    % are used in test_Shampoo.m, where the whole eigensystem is needed. Using
    %   eigB = eig(B);
    % gives slightly different values for t_1 and t_2.
    [~, D] = eig(B);
    eigB = diag(D);

    % General information about the matrices.
    fprintf(fid, '$B_{%d}$ & %s\n    & %s & %s & %s & & $\\wt B_{%d}$',...
            i, latex_int_fill(n, 4),...
            latex_sci_not(min(eigB), 1, 1),...
            latex_sci_not(max(eigB), 1, 1),...
            latex_int_fill(rank(B), 3), i);
    % Compute the order of approximations.
    for tol = [0.1, n^(3/2)*eps('single')/2]
      fprintf(fid, '\n    & %s & %s',...
              latex_sci_not(tol, 1, 1),...
              latex_int_fill(sum(eigB >= tol), 3));
    end
    fprintf(fid, '\\\\\n');
  end
  fprintf(fid, '\\bottomrule\n\\end{tabularx}');
  fclose(fid);

  mp.Digits(olddigits);
end
