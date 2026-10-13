function t01_three_way_crosscheck()
%T01_THREE_WAY_CROSSCHECK  A0/A1/A2：三条独立路径交叉验证
%
%   ── 对应验收标准 ────────────────────────────────────────────
%   A1  三层非对称平板：TE0/TE1 的 n_eff 与解析解吻合到相对误差 1e-6
%   A2  N 层 TMM 退化到三层时，与 A1 结果一致
%   A0  （本项目新增，见交接报告讨论）三条独立路径互验
%
%   ── 三条路径的独立性从哪来 ──────────────────────────────────
%   路径1 analytic3：解超越方程 κd = arctan(γ_s/κ)+arctan(γ_c/κ)+mπ，二分求根。
%   路径2 tmm      ：构造 2×2 传输矩阵连乘，求色散函数零点。数学路径完全不同，
%                    但对多层结构同样成立（路径1 只对三层有效）。
%   路径3 fdm      ：把 Helmholtz 方程直接离散成稀疏矩阵本征值问题。
%                    连"色散方程"这个中间产物都不要。
%
%   路径1 与路径2 共享同一个物理前提（Maxwell + 界面连续），但公式形式已由
%   t00 用符号运算独立验证过；路径3 与它们不共享任何一行推导。
%   三方一致的结论，可信度与两方一致不是一个量级。
%
%   ── 判据 ────────────────────────────────────────────────────
%   |analytic − tmm| / n  < 1e-9    （两者都是超越方程求根，应到机器精度量级）
%   |analytic − fdm| / n  < 1e-6    （A1 的要求；FDM 受离散化与外推限制）

fprintf('\n========== t01 三路交叉验证 ==========\n');

cases = {
    '对称    3.40 / 3.50 / 3.40,  d = 2.0', 3.40, 3.50, 3.40, 2.0
    '弱不对称 3.40 / 3.50 / 3.35,  d = 2.0', 3.40, 3.50, 3.35, 2.0
    '强不对称 3.40 / 3.50 / 3.00,  d = 2.0', 3.40, 3.50, 3.00, 2.0
    '反向不对称 3.30 / 3.50 / 3.40, d = 1.5', 3.30, 3.50, 3.40, 1.5
};

lam      = 1.0;
n_modes  = 4;
worst_12 = 0;
worst_13 = 0;
all_pass = true;

for c = 1:size(cases, 1)
    name  = cases{c,1};
    n_sub = cases{c,2};  n_core = cases{c,3};
    n_cov = cases{c,4};  d_core = cases{c,5};

    L = build_layers(layers_three(n_sub, n_core, n_cov, d_core));

    m1 = solve_slab(L, lam, struct('solver', 'analytic3', 'n_modes', n_modes));
    m2 = solve_slab(L, lam, struct('solver', 'tmm',       'n_modes', n_modes));
    m3 = solve_slab(L, lam, struct('solver', 'fdm',       'n_modes', n_modes));

    fprintf('\n── %s ──\n', name);
    fprintf('   模数：analytic=%d  tmm=%d  fdm=%d\n', numel(m1), numel(m2), numel(m3));

    if numel(m1) ~= numel(m2)
        fprintf('   [FAIL] analytic 与 tmm 的模数不一致\n');
        all_pass = false;
    end

    ncmp = min([numel(m1), numel(m2), numel(m3)]);
    for k = 1:ncmp
        n1 = real(m1(k).n_eff);
        n2 = real(m2(k).n_eff);
        n3 = real(m3(k).n_eff);

        e12 = abs(n1 - n2)/n1;
        e13 = abs(n1 - n3)/n1;
        worst_12 = max(worst_12, e12);
        worst_13 = max(worst_13, e13);

        fprintf(['   m=%d  n_eff = %.12f | tmm Δ=%.2e | fdm Δ=%.2e' ...
                 '  (fdm: %.12f)\n'], k-1, n1, e12, e13, n3);
    end

    % 模数一致性：fdm 的截断可能丢掉接近截止的高阶模，此处只告警不判失败
    if numel(m3) < numel(m1)
        fprintf('   [WARN] fdm 少了 %d 个模（通常是最接近截止的高阶模，见 A4 讨论）\n', ...
            numel(m1) - numel(m3));
    end
end

fprintf('\n── 汇总 ──\n');
ok12 = worst_12 < 1e-9;
ok13 = worst_13 < 1e-6;
report(ok12, sprintf('A2  max|analytic−tmm|/n = %.2e  （判据 < 1e-9）', worst_12));
report(ok13, sprintf('A1  max|analytic−fdm|/n = %.2e  （判据 < 1e-6）', worst_13));

if ok12 && ok13 && all_pass
    fprintf('t01 通过。\n');
else
    fprintf('t01 未通过。\n');
end
end

% ══════════════════════════════════════════════════════════════════════
function report(ok, msg)
if ok
    fprintf('  [PASS] %s\n', msg);
else
    fprintf('  [FAIL] %s\n', msg);
end
end
