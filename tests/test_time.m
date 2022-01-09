function test_time(n, alpha, ratios, scale, tol, rankdef, nonsym, nonnegmat)
% Compare structured methods and diagonalization in terms of execution time.
% The quantity that varies, in this script, is the ration n/k.
  rng(0);

  fastmode = false;       % Skip the expensive (unstructed) methods.
  nsizes = length(ratios);
  sizes = zeros(1, nsizes);

  % Check whether relative residual is below threshold.
  debugmode = false;
  threshold = 1e-10;

  time_unstruct_schur = zeros(1, nsizes);
  time_direct_schur = zeros(1, nsizes);
  time_direct_dbp = zeros(1, nsizes);
  time_struct_dbp = zeros(1, nsizes);

  for i = 1:nsizes
    sizes(i) = round(n*ratios(i));

    k = sizes(i);

    if ~nonsym % the symmetric case
        if nonnegmat
            U = rand(n,k) / n;
            V = U;
        else
            U = randn(n,k) / n;
            V = U;
        end
        if rankdef && k > 1
            U(:,end) = U(:,1:end-1)*ones(k-1,1)/(k-1);
        end
    else % nonsymmetric case: U and V are different from the same distribution
        if nonnegmat
            U = rand(n,k) / n;
            V = rand(n,k) / n;
        else
            U = randn(n,k) / n;
            V = randn(n,k) / n;
        end
        if rankdef && k > 1
            U(:,end) = U(:,1:end-1)*ones(k-1,1)/(k-1);
        end
    end
    A = alpha*eye(n) + U*V';

    if ~fastmode
      % Unstructured Schur algorithm.
      [X, time] = smart_timer(@sqrtm, A);
      time_unstruct_schur(i) = time;
      check_residual(A, X, threshold, debugmode);

    end

    % Direct formula with Schur method.
    [Y, time] = smart_timer(@sqrtm_dir_schur, alpha, U, V);
    time_direct_schur(i) = time;
    check_residual(A, Y, threshold, debugmode);

    % Direct formula with Denman-Beavers iteration in product form.
    [Y, time] = smart_timer(@sqrtm_dir_db_prod, alpha, U, V);
    time_direct_dbp(i) = time;
    check_residual(A, Y, threshold, debugmode);

    % Structured Denman-Beavers iteration in product form.
    [Z, time] = smart_timer(@sqrtm_s_db_prod, alpha, U, V, scale, tol);
    time_struct_dbp(i) = time;
    check_residual(A, Z, threshold, debugmode)

  end

  % Save results to file.
  if nonsym
      suffix = '_nonsym';
  else
      suffix = '';
  end
  if nonnegmat
      suffix = sprintf('%s_nonneg', suffix);
  end
  if scale
      suffix = sprintf('%s_scaled', suffix);
  end
  if rankdef
      suffix = sprintf('%s_rd', suffix);
  end

  filename = sprintf('time_%.3f_%05d%s.dat', alpha, n, suffix);
  fid = fopen(filename, 'w');
  for i = 1:nsizes
    fprintf(fid,... % '%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n',...
            '%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n',...
            ratios(i), time_unstruct_schur(i),...
            time_direct_schur(i), time_direct_dbp(i), time_struct_dbp(i));
  end

  % Plot results.
  clf
  semilogy(ratios, time_unstruct_schur,'v-','MarkerSize',8,'LineWidth',1.5);
  hold on;
  semilogy(ratios, time_direct_schur,'s-','MarkerSize',8,'LineWidth',1.5);
  semilogy(ratios, time_direct_dbp,'o-','MarkerSize',8,'LineWidth',1.5);
  semilogy(ratios, time_struct_dbp,'d-','MarkerSize',8,'LineWidth',1.5);

  lgd = legend('unstruct\_schur', 'dir\_schur',...
        'dir\_dbp', 'struct\_dbp', 'interpreter', 'latex');
  lgd.FontSize = 12;
  set(gca,'linewidth',1.5)
  set(gca,'fontsize',12)
  xlabel('$k/n$','interpreter','latex','FontWeight','normal','fontsize',18)
  ylabel('time')
  box on
  title(sprintf('$n=%d$', n),...
        'interpreter', 'latex',...
        'FontWeight','normal',...
        'fontsize',18)

  % set(gcf, 'Color', 'w'); str = ('../figs/newfig.pdf'); export_fig(str)

  %% LOCAL FUNCTIONS

  function check_residual(A, X, threshold, debugmode)
    if debugmode
      n = size(A, 1);
      norm(X^2-A, 1) / norm(A, 1)
      assert(norm(X^2-A, 1) / norm(A, 1) < threshold);
    end
  end

end