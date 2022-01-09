function test_accuracy(n, epsilon, ratios, scale, tol, rankdef, nonsym, nonnegmat)
% Compare structured methods and diagonalization in terms 1-norm residual.
% The quantity that varies, in this script, is the ration n/k.
  rng(0);

  % If FALSE, the function doesn't compute the reference solution.
  % If the reference solution is not available, the forward error is set to NaN.
  compute_fwd_err = false;

  % If FALSE, residual and forward error are computed in the working precision.
  % Otherwise, higher precision (using the Advanpix Multiprecision Toolbox is
  % used).
  use_high_prec = true;

  % Norm used in the computation of residual and forward error.
  mynorm = 2;

  nsizes = length(ratios);
  sizes = zeros(1, nsizes);
  alphau = zeros(1, nsizes);

  % currently four methods are being tested.
  residual_unstruct_schur = zeros(1, nsizes);
  residual_direct_schur = zeros(1, nsizes);
  residual_direct_dbp = zeros(1, nsizes);
  residual_struct_dbp = zeros(1, nsizes);

  condu = zeros(1, nsizes);
  rel_err_unstruct_schur = zeros(1, nsizes);
  rel_err_direct_schur = zeros(1, nsizes);
  rel_err_direct_dbp = zeros(1, nsizes);
  rel_err_struct_dbp = zeros(1, nsizes);

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
    A = epsilon*eye(n) + U*V';

    % Compute reference solution.
    if compute_fwd_err
      mpdigits = mp.Digits(34);
      Xref = double(sqrtm(mp(A)));
      mp.Digits(mpdigits);
    else
      Xref = nan;
    end

    % Condition numbers.
    if compute_fwd_err % then also compute conditioning.
      [~, alpha, condx] = sqrtm(A);
      condu(i) = condx * eps();
    else % compute only alpha.
      [~, alpha, ~] = sqrtm(A);
      condu(i) = nan;
    end
    alphau(i) = alpha * eps();

    % Unstructured Schur algorithm.
    X = sqrtm(A); % A is now passed for the unstructured sqrtm
    residual_unstruct_schur(i) = compute_residual(A, X, mynorm, use_high_prec);
    rel_err_unstruct_schur(i) = compute_rel_err_if_needed(Xref, X, mynorm, use_high_prec);

    % Direct formula and Schur method.
    X = sqrtm_dir_schur(epsilon, U, V);
    residual_direct_schur(i) = compute_residual(A, X, mynorm, use_high_prec);
    rel_err_direct_schur(i) = compute_rel_err_if_needed(Xref, X, mynorm, use_high_prec);

    % Direct formula and product Denman-Beavers iteration.
    X = sqrtm_dir_db_prod(epsilon, U, V);
    residual_direct_dbp(i) = compute_residual(A, X, mynorm, use_high_prec);
    rel_err_direct_dbp(i) = compute_rel_err_if_needed(Xref, X, mynorm, use_high_prec);

    % Structured Denman-Beavers iteration in product form.
    X = sqrtm_s_db_prod(epsilon, U, V, scale, tol);
    residual_struct_dbp(i) = compute_residual(A, X, mynorm, use_high_prec);
    rel_err_struct_dbp(i) = compute_rel_err_if_needed(Xref, X, mynorm, use_high_prec);

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

  filename = sprintf('accuracy_%.3f_%05d%s.dat', epsilon, n, suffix);
  fid = fopen(filename, 'w');
  for i = 1:nsizes
    % ['%.5e\t',...
    %  '%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t',...
    %  '%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t\\\\\n'],...
    fprintf(fid,...
            ['%.5e\t',...
             '%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t',...
             '%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t\\\\\n'],...
            ratios(i),...
            alphau(i), residual_unstruct_schur(i), residual_direct_schur(i),...
            residual_direct_dbp(i), residual_struct_dbp(i),...
            condu(i), rel_err_unstruct_schur(i), rel_err_direct_schur(i),...
            rel_err_direct_dbp(i), rel_err_struct_dbp(i));
  end

  %% Plot results.
  clf

  % Residual.
  if compute_fwd_err
    subplot(1,2,1);
  end
  semilogy(ratios, residual_unstruct_schur,'v-','MarkerSize',8,'LineWidth',1.5);
  hold on;
  semilogy(ratios, residual_direct_schur,'s-','MarkerSize',8,'LineWidth',1.5);
  semilogy(ratios, residual_direct_dbp,'o-','MarkerSize',8,'LineWidth',1.5);
  semilogy(ratios, residual_struct_dbp,'d-','MarkerSize',8,'LineWidth',1.5);

  axis([0, 1, 1e-16, 1e-12]);

  lgd = legend('unstruct\_schur', 'dir\_schur',...
      'dir\_dbp', 'struct\_dbp', 'interpreter', 'latex');
  lgd.FontSize = 12;
  set(gca,'linewidth',1.5)
  set(gca,'fontsize',12)
  xlabel('$k/n$','interpreter','latex','FontWeight','normal','fontsize',18)
  ylabel('Relative residual')
  box on
  title(sprintf('$n=%d$', n), 'interpreter', 'latex', 'FontWeight','normal',...
        'fontsize',18)

  % Relative error.
  if compute_fwd_err
    subplot(1,2,2);

    semilogy(ratios, rel_err_unstruct_schur,'v-','MarkerSize',8,'LineWidth',1.5);
    hold on;
    semilogy(ratios, rel_err_direct_schur,'s-','MarkerSize',8,'LineWidth',1.5);
    semilogy(ratios, rel_err_direct_dbp,'o-','MarkerSize',8,'LineWidth',1.5);
    semilogy(ratios, rel_err_struct_dbp,'d-','MarkerSize',8,'LineWidth',1.5);
    axis([0, 1, 1e-16, 1e-12]);

    lgd = legend('unstruct\_schur', 'dir\_schur',...
        'dir\_dbp', 'struct\_dbp', 'interpreter', 'latex');
    lgd.FontSize = 12;
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',12)
    xlabel('$k/n$','interpreter','latex','FontWeight','normal','fontsize',18)
    ylabel('Relative forward error')
    box on
    title(sprintf('$n=%d$', n),...
          'interpreter', 'latex',...
          'FontWeight','normal',...
          'fontsize',18)
  end

  %% LOCAL FUNCTIONS
  function err = compute_rel_err_if_needed(Xref, X, mynorm, use_high_prec)
    if compute_fwd_err
      err = compute_rel_err(Xref, X, mynorm, use_high_prec);
    else
      err = nan;
    end
  end

end
