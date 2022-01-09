function test_direct_formulae(n, k, alpha, ratios, conds)

  mynorm = 2;
  use_high_prec = true;

  %% Varying rank.
  nsizes = length(ratios);
  sizes = zeros(1, nsizes);
  residual_direct_naive = zeros(1, nsizes);
  residual_direct_new = zeros(1, nsizes);
  condVU = zeros(1, nsizes);

  for i = 1:nsizes
    sizes(i) = round(n*ratios(i));
    kloc = sizes(i);
    U = randn(n,kloc) / n;

    V = U;
    A = alpha*eye(n) + U*V';
    [X_direct, condVU(i)] = sqrtm_direct_naive(alpha, U, V);

    X_direct_new = sqrtm_direct_new(alpha, U, V);

    residual_direct_naive(i) = compute_residual(A, X_direct, mynorm, use_high_prec);
    residual_direct_new(i) = compute_residual(A, X_direct_new, mynorm, use_high_prec);
  end

  % Draw plot.
  subplot(2, 1, 1)
  semilogy(ratios, residual_direct_naive,...
           'v-', 'LineWidth', 2, 'MarkerSize', 10)
  hold on
  semilogy(ratios, residual_direct_new,...
           'v-', 'LineWidth', 2, 'MarkerSize', 10)
  semilogy(ratios, condVU*eps('double'),...
           'k--', 'LineWidth', 2, 'MarkerSize', 10)
  hold off

  lgd = legend('direct\_naive', 'direct\_new', 'condu', 'Location', 'northwest');
  lgd.FontSize = 12;
  set(gca, 'linewidth', 1.5)
  set(gca, 'fontsize', 12)
  xlabel('$k/n$', 'interpreter', 'latex', 'FontWeight', 'normal', 'fontsize', 18)
  ylabel('relative residual', 'fontsize', 18)

  % Save data.
  filename = sprintf('formulae_comparison_size.dat');
  fid = fopen(filename, 'w');
  for i = 1:nsizes
    fprintf(fid, '%.5e\t%.5e\t %.5e\t%.5e\t\\\\\n',...
            ratios(i),...
            residual_direct_naive(i), residual_direct_new(i),...
            condVU(i)*eps());
  end

  %% Varying condition number.
  nconds = length(conds);
  residual_direct_naive = zeros(1, nconds);
  residual_direct_new = zeros(1, nconds);
  condVU = zeros(1, nconds);

  for i = 1:nconds
    S = logspace(-log10(conds(i)),0,k);
    U = orth(randn(n,k));
    V = U .* S;

    A = alpha*eye(n) + U*V';
    [X_direct, condVU(i)] = sqrtm_direct_naive(alpha, U, V);

    X_direct_new = sqrtm_direct_new(alpha, U, V);

    residual_direct_naive(i) = compute_residual(A, X_direct);
    residual_direct_new(i) = compute_residual(A, X_direct_new);
  end

  % Draw plot.
  subplot(2,1,2)
  loglog(condVU, residual_direct_naive,...
         'v-', 'LineWidth', 2, 'MarkerSize', 10)
  hold on
  loglog(condVU, residual_direct_new,...
         'v-', 'LineWidth', 2, 'MarkerSize', 10)
  loglog(condVU, condVU*eps('double'),...
         'k--', 'LineWidth', 2, 'MarkerSize', 10)
  hold off

  axis([1 1e16 1e-16 1])
  lgd = legend('direct\_naive', 'direct\_new', 'condu',...
               'Location', 'northwest');
  lgd.FontSize = 12;
  set(gca, 'linewidth', 1.5)
  set(gca, 'fontsize', 12)
  xlabel('cond', 'interpreter', 'latex', 'FontWeight', 'normal', 'fontsize', 18)
  ylabel('relative residual', 'fontsize', 18)

  % Save data.
  filename = sprintf('formulae_comparison_condVU.dat');
  fid = fopen(filename, 'w');
  for i = 1:nconds
    fprintf(fid, '%.5e\t%.5e\t%.5e\t%.5e\t\\\\\n',...
            condVU(i),...
            residual_direct_naive(i), residual_direct_new(i),...
            condVU(i)*eps());
  end

  %% LOCAL FUNCTIONS

  function [X ,condVU] = sqrtm_direct_naive(alpha, U, V)
  % Compute the square root of EPSILON*I + U*V' using the direct formula (1.4).
    [n,k] = size(U);
    M = V'*U;
    condVU = cond(M, 2);
    X = M;
    X(1:k+1:k^2) = X(1:k+1:k^2) + alpha;
    X = sqrtm(X);
    X(1:k+1:k^2) = X(1:k+1:k^2) - sqrt(alpha);
    X = U * (M \ X) * V';
    X(1:n+1:n^2) = X(1:n+1:n^2) + sqrt(alpha);
  end

  function [X ,condVU] = sqrtm_direct_new(alpha, U, V)
  % Compute the square root of EPSILON*I + U*V' using the direct formula (1.8).
    [n,k] = size(U);
    M = V'*U;
    condVU = cond(M, 2);
    X = M;
    X(1:k+1:k^2) = X(1:k+1:k^2) + alpha;
    X = sqrtm(X);
    X(1:k+1:k^2) = X(1:k+1:k^2) + sqrt(alpha);
    X = U * (X \  V');
    X(1:n+1:n^2) = X(1:n+1:n^2) + sqrt(alpha);
  end

end