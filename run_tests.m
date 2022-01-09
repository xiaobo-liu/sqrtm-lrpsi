%% This script runs all the tests.
warning off
addpath('data', 'include', 'methods', 'tests');
rng(0)

%% SECTION 1 - Example showing that (1.4) is unstable
n = 100;
k = 10;
alpha = 1.0;
ratios = [0.01:0.01:1];
conds = 10.^[0:16];

fprintf('Confirming instability of (1.4)...');
global_init = tic();
test_direct_formulae(n, k, alpha, ratios, conds)
fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);

%% UNUSED EXPERIMENT
% n = 1000;
% tol = 1e-15;       % Stopping criterion.
% scale = true;      % Use scaled iterations.
% rank_deficient = false;   % Use full-rank updates.
% for nonsym = [false]
%     nonneg_mat = true;
%     for alpha = [1e-6]
%         fprintf('Running full-rank accuracy test... [alpha = %.3e]', alpha);
%         global_init = tic();
%         ratios = [0.01, 0.05:0.1:1];
%         test_accuracy(n, alpha, ratios, scale, tol, rank_deficient, nonsym, nonneg_mat);
%         fprintf(' done\t[%5.2f min]\n', toc(global_init) / 60);
%     end
% end

%% SECTION 5.1 - Quality
n = 100;
tol = 1e-15;       % Stopping criterion.
scale = true;      % Use scaled iterations.
rank_deficient = false;   % Use full-rank updates.
for nonsym = [false true]
  for nonneg_mat = [false true] % nonnegative entries.
    alpha = 1.0;
    fprintf('Running full-rank accuracy test... [alpha = %.1e]', alpha);
    global_init = tic();
    ratios = [0.01, 0.05:0.01:1];
    test_accuracy(n, alpha, ratios, scale, tol, rank_deficient, nonsym, nonneg_mat);
    fprintf(' done\t[%5.2f min]\n', toc(global_init) / 60);
  end
  nonneg_mat = true;
  for alpha = [0.1 0.001]
    fprintf('Running full-rank accuracy test... [alpha = %.1e]', alpha);
    global_init = tic();
    ratios = [0.01, 0.05:0.01:1];
    test_accuracy(n, alpha, ratios, scale, tol, rank_deficient, nonsym, nonneg_mat);
    fprintf(' done\t[%5.2f min]\n', toc(global_init) / 60);
  end
end

%% SECTION 5.1 - Results are only mentioned, plot is not shown.
% n = 100;
% tol = 1e-15;       % Stopping criterion.
% scale = true;      % Use scaled iterations. on this simpler problem
% rank_deficient = false;   % Use full-rank updates.
% for nonsym = [false]
%   nonneg_mat = true;
%   for alpha = [1e-6]
%     fprintf('Running full-rank accuracy test... [alpha = %.1e]', alpha);
%     global_init = tic();
%     ratios = [0.01, 0.05:0.1:1];
%     test_accuracy(n, alpha, ratios, scale, tol, rank_deficient, nonsym, nonneg_mat);
%     fprintf(' done\t[%5.2f min]\n', toc(global_init) / 60);
%   end
% end

%% SECTION 5.2 - Timings
tol = 1e-15;          % Stopping criterion.
scale = true;         % Use scaled iterations.
rank_deficient = false; % U is full-rank.
nonneg_mat = false;
ratios = [0.03, 0.05:0.05:1];
for nonsym = [false true]
  for alpha = [0.1]
    for n = [1000, 4000, 7000, 10000]
      fprintf('Running full-rank timing test for alpha = %.1e, n = %5d... ', alpha, n);
      global_init = tic();
      test_time(n, alpha, ratios, scale, tol, rank_deficient, nonsym, nonneg_mat);
      fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);
    end
  end
end

%% SECTION 5.3 - Positive definite matrices from applications
fprintf('Running test with Lingvo matrices... ');
global_init = tic();
print_lingvo_table
test_lingvo
fprintf('done\t[%5.2f min]', toc(global_init) / 60);

exit
