function test_lingvo()
  scale = true; % Use scaled iterations.

  Bcell = cell(3,1);

  truncat_t = zeros(3,2); % stores the truncation parameters t
  load('data/test_mat_1.mat');
  Bcell{1} = data;

  load('data/test_mat_2.mat');
  Bcell{2} = data;

  load('data/test_mat_3.mat');
  Bcell{3} = data;

  u32 = eps('single')/2;
  tols = zeros(3,2); % the tolerances
  for i=1:3
    tols(i,1) = 0.1;
    tols(i,2) = (size(Bcell{i},1))^(3/2)*u32;
  end

  alphas = single([1e-6 1e-3 1e0]);
  tol_iter = 2^3 * eps('single');

  maxit = 50;

  nmatrices = 3;
  nranks = 2;
  nalphas = length(alphas);

  % format: matrix/rank/alpha
  res_unstruct_schur = zeros(nmatrices, nranks, nalphas);
  res_direct_schur = zeros(nmatrices, nranks, nalphas);
  res_direct_dbp = zeros(nmatrices, nranks, nalphas);
  res_struct_dbp = zeros(nmatrices, nranks, nalphas);

  err_unstruct_schur = zeros(nmatrices, nranks, nalphas);
  err_direct_schur = zeros(nmatrices, nranks, nalphas);
  err_direct_dbp = zeros(nmatrices, nranks, nalphas);
  err_struct_dbp = zeros(nmatrices, nranks, nalphas);

  time_unstruct_schur = zeros(nmatrices, nranks, nalphas);
  time_direct_schur = zeros(nmatrices, nranks, nalphas);
  time_direct_dbp = zeros(nmatrices, nranks, nalphas);
  time_struct_dbp = zeros(nmatrices, nranks, nalphas);

  alphau = zeros(nmatrices, nranks, nalphas); % stability factor ALPHA
  condu = zeros(nmatrices, nranks, nalphas); % condition number

  for i = 1:nmatrices
    B = Bcell{i};
    [V, D] = eig(B);
    V = flip(V,2);
    D = diag(flip(diag(D)));
    diag_D = diag(D); % the vector of descendingly ordered eigenvalues.
    for j = 1:nranks
      r = length(diag_D(diag_D>=tols(i,j)));
      truncat_t(i,j) = r;
      U = V(:,1:r)*sqrt(D(1:r,1:r));
      % rank(U*U')
      % norm(B-U*U') / norm(B)
      for k = 1:nalphas
        alpha = alphas(k);
        n = size(B,1);
        % A = alpha * I + U*U'
        A = U*U' + alpha*eye(n,class(alpha));

        % Compute \alpha_2(X), cond, and reference solution in double precision.
        [X_exact, alphatwoX, condx] = sqrtm(double(A));
        % alphau(i,j,k) = alphatwoX * eps();
        % condu(i,j,k) = condx * eps();

        % Unstructured Schur.
        [X_unstruct, time] = smart_timer(@sqrtm, A);
        res_unstruct_schur(i,j,k) = compute_residual(A, X_unstruct);
        time_unstruct_schur(i,j,k) = time;
        err_unstruct_schur(i,j,k) = compute_rel_err(X_exact, X_unstruct);

        % Direct formula and Schur.
        [X_dir, time] = smart_timer(@sqrtm_dir_schur, alpha, U, U);
        res_direct_schur(i,j,k) = compute_residual(A, X_dir);
        time_direct_schur(i,j,k) = time;
        err_direct_schur(i,j,k) = compute_rel_err(X_exact, X_dir);

        % Direct formula and Denman-Beavers iteration in product form.
        [X_db, time] = smart_timer(@sqrtm_dir_db_prod, alpha, U, U);
        res_direct_dbp(i,j,k) = compute_residual(A, X_db);
        time_direct_dbp(i,j,k) = time;
        err_direct_dbp(i,j,k) = compute_rel_err(X_exact, X_db);

        % Structured Denman-Beavers iteration in product form.
        [X_dbp, time] = smart_timer(@sqrtm_s_db_prod, alpha, U, U, scale, tol_iter, maxit);
        res_struct_dbp(i,j,k) = compute_residual(A, X_dbp);
        time_struct_dbp(i,j,k) = time;
        err_struct_dbp(i,j,k) = compute_rel_err(X_exact, X_dbp);
      end
    end
  end

  %% Save results to file.
  methodnames = {'\unstsqrtm', '\dirsqrtm', '\dirdbp', '\stdbp'};
  nmethods = length(methodnames);

  res_global = zeros(nmatrices, nranks, nalphas, nmethods);
  res_global(:,:,:,1) = res_unstruct_schur;
  res_global(:,:,:,2) = res_direct_schur;
  res_global(:,:,:,3) = res_direct_dbp;
  res_global(:,:,:,4) = res_struct_dbp;

  err_global = zeros(nmatrices, nranks, nalphas, nmethods);
  err_global(:,:,:,1) = err_unstruct_schur;
  err_global(:,:,:,2) = err_direct_schur;
  err_global(:,:,:,3) = err_direct_dbp;
  err_global(:,:,:,4) = err_struct_dbp;

  time_global = zeros(nmatrices, nranks, nalphas, nmethods);
  time_global(:,:,:,1) = time_unstruct_schur;
  time_global(:,:,:,2) = time_direct_schur;
  time_global(:,:,:,3) = time_direct_dbp;
  time_global(:,:,:,4) = time_struct_dbp;

  if scale
    suffix = '_scaled';
  else
    suffix = '';
  end
  filename = sprintf('lingvo_table%s.dat', suffix);
  fid = fopen(filename, 'w');
  fprintf(fid, ['\\begin{tabularx}{\\textwidth}{ccX cc cc cc}\n',...
                '\\toprule\n',...
                '& \\multirow{2}{*}{$t$} &',...
                '\\multicolumn{1}{c}{\\multirow{2}{*}{Method}} &',...
                '\\multicolumn{2}{c}{$\\alpha=10^{-6}$} & ',...
                '\\multicolumn{2}{c}{$\\alpha=10^{-3}$} & ',...
                '\\multicolumn{2}{c}{$\\alpha=1$} \\\\\n',...
                '\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}\n',...
                '& & & ',...
                'res & time & ',...
                'res & time & ',...
                'res & time \\\\\n',...
                '\\midrule\n']);
  for i = 1:nmatrices
    for j = 1:nranks
      for l = 1:nmethods % One line per method
        if l == 1 && j == 1              % Print matrix.
          fprintf(fid, '$B_{%d}$ & $%d$ & ',...
                  i, truncat_t(i, j));
        elseif l == 1                    % Print t.
          fprintf(fid, '& $%d$ & ',...
                  truncat_t(i, j));
        else                             % Print method.
          fprintf(fid, '             & & ');
        end
        fprintf(fid,'%s & ', methodnames{l});
        for k = 1:nalphas % Two columns per value of alpha | res | time |
          fprintf(fid, '%s & %s',...
                  latex_sci_not(res_global(i,j,k,l),0,2),...
                  latex_sci_not(time_global(i,j,k,l),0,2));
          if k ~= nalphas
            fprintf(fid, ' & ');
          end
        end
        if i*j*l ~= nmatrices*nranks*nmethods
          if j*l == nranks*nmethods
            fprintf(fid, '\\\\\n\\midrule\n');
          elseif l == nmethods
            fprintf(fid, '\\\\[5pt]\n');
          else
            fprintf(fid, '\\\\\n');
          end
        else
          fprintf(fid, '\\\\\n\\bottomrule\n\\end{tabularx}');
        end
      end
    end
  end
  fclose(fid);

  %% Plot
  lw = 2;

  clf
  subplot(1,3,1)
  semilogy(res_unstruct_schur(:), 'LineWidth', lw)
  hold on
  semilogy(res_direct_schur(:), 'LineWidth', lw)
  semilogy(res_direct_dbp(:), 'LineWidth', lw)
  semilogy(res_struct_dbp(:), 'LineWidth', lw)
  title('Residual')
  legend('u\_schur', 'dir\_schur', 'dir\_db\_prod', 's\_db\_prod',...
         'interpreter', 'latex')

  subplot(1,3,2)
  semilogy(err_unstruct_schur(:), 'LineWidth', lw)
  hold on
  semilogy(err_direct_schur(:), 'LineWidth', lw)
  semilogy(err_direct_dbp(:), 'LineWidth', lw)
  semilogy(err_struct_dbp(:), 'LineWidth', lw)
  title('Error')
  legend('u\_schur', 'dir\_schur', 'dir\_db\_prod', 's\_db\_prod',...
         'interpreter', 'latex')

  subplot(1,3,3)
  semilogy(time_unstruct_schur(:), 'LineWidth', lw)
  hold on
  semilogy(time_direct_schur(:), 'LineWidth', lw)
  semilogy(time_direct_dbp(:), 'LineWidth', lw)
  semilogy(time_struct_dbp(:), 'LineWidth', lw)
  title('Time')
  legend('u\_schur', 'dir\_schur', 'dir\_db\_prod', 's\_db\_prod',...
         'interpreter', 'latex')
end