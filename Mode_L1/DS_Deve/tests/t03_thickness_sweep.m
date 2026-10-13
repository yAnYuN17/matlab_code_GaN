function t03_thickness_sweep()
%T03_THICKNESS_SWEEP  A4：扫描芯层厚度，模式数按预期增加
%
%   用【真实 GaN 量级参数】（450 nm，GaN 芯 / AlGaN 包层，Δn ≈ 0.026），
%   扫描芯层厚度，验证：
%     (1) 模式数随厚度单调不减；
%     (2) 新模出现的厚度阈值与解析预测一致
%             d_c,m = [arctan(√a) + mπ] / (k0·√(n_core² − n_high²))
%
%   ── 为什么"厚度增大 ⇒ 模式数增加"值得单独验 ─────────────────
%   物理图像：固定厚度时能容纳的横向半波长数是有限的；厚度增大 ⇒ 允许更多
%   横向半波长 ⇒ 更多模式。V ∝ d 把这句话定量化：V 每跨过 arctan(√a)+mπ，
%   就多出一个模。
%
%   代码里若把 V 的定义写错（比如漏了 √ 里的平方、或把 k0 写成 k0/2），
%   曲线形状会整体变形，但【A1/A2 依然能过】——因为解析路径与 TMM 路径
%   用的是同一个 d，两边的错误会互相抵消。所以 A4 是对"厚度依赖关系"
%   的独立检查，补的是 A1/A2 覆盖不到的那一块。
%
%   ── 顺带说明：GaN LD 垂直方向为什么是少模的 ─────────────────
%   本脚本用的 Δn 只有 0.026（GaN vs Al0.06GaN 的文献实测值）。
%   在 450 nm 下这意味着要把芯层做到 ~1.4 μm 才出现 TE2，
%   而实际器件的垂直结构总厚度通常远小于此 —— 所以 GaN 基 LD 的
%   垂直方向天然是基模，多模问题全部发生在【横向】（阶段3 的内容）。
%   这个结论不是靠"感觉"，而是由本脚本的截止厚度算出来的。

fprintf('\n========== t03 厚度扫描（A4）==========\n');

% 真实 GaN 量级参数 @450 nm
lam    = 0.45;
n_core = 2.430;      % GaN 波导
n_sub  = 2.404;      % Al0.06GaN 包层（文献：与 GaN 差约 0.026）
n_cov  = 2.374;      % Al0.10GaN 包层

n_high = max(n_sub, n_cov);
n_low  = min(n_sub, n_cov);
a      = (n_high^2 - n_low^2)/(n_core^2 - n_high^2);
k0     = 2*pi/lam;
VperD  = k0*sqrt(n_core^2 - n_high^2);      % V = VperD * d

fprintf('参数：n_core=%.3f  n_sub=%.3f  n_cov=%.3f  lambda=%.3f um\n', ...
        n_core, n_sub, n_cov, lam);
fprintf('不对称度 a = %.4f   V/d = %.5f um^-1\n\n', a, VperD);

d_scan = 0.10:0.01:2.50;

D = [];  Nmodes = [];
n_eff_by_order = nan(numel(d_scan), 6);

for i = 1:numel(d_scan)
    d  = d_scan(i);
    L  = build_layers(layers_three(n_sub, n_core, n_cov, d));
    ne = solve_analytic3(L, lam, struct('pol', 'TE', 'n_modes', 6));

    D(end+1)      = d;                                   %#ok<AGROW>
    Nmodes(end+1) = numel(ne);                           %#ok<AGROW>
    for k = 1:numel(ne)
        n_eff_by_order(i, k) = real(ne(k));
    end
end

% 模式数单调性
mono = all(diff(Nmodes) >= 0);
report(mono, '模式数随厚度单调不减');

% ── 截止厚度核对 ─────────────────────────────────────────────
fprintf('\n── 各阶模出现的厚度阈值 ──\n');
all_pass = true;
d_c_theory = zeros(1, 6);

for m = 0:3
    d_th = (atan(sqrt(a)) + m*pi)/VperD;
    d_c_theory(m+1) = d_th;

    d_num = find_cutoff(n_sub, n_core, n_cov, lam, m, 0.01, 3.5);
    if isnan(d_num)
        fprintf('   m=%d  理论 d_c = %.5f um   数值：扫描范围内未出现\n', m, d_th);
        continue;
    end
    err = abs(d_num - d_th)/d_th;
    if err > 1e-4, all_pass = false; end
    fprintf('   m=%d  理论 d_c = %.5f um   数值 d_c = %.5f um   相对偏差 %.2e\n', ...
            m, d_th, d_num, err);
end

% ── 绘图 ─────────────────────────────────────────────────────
fig = figure('Visible', 'off', 'Position', [80 80 1100 480]);

subplot(1,2,1); hold on; grid on; box on;
for k = 1:size(n_eff_by_order,2)
    plot(D, n_eff_by_order(:,k), 'LineWidth', 1.4, ...
         'DisplayName', sprintf('TE%d', k-1));
end
yline(n_high, 'k--', 'DisplayName', 'n_{high} (cladding)');
xlabel('core thickness d (\mum)');
ylabel('n_{eff}');
title('t03  n_{eff} vs core thickness');
legend('Location', 'southeast');

subplot(1,2,2); hold on; grid on; box on;
stairs(D, Nmodes, 'LineWidth', 1.6, 'Color', [0.85 0.33 0.10]);
for m = 0:3
    if d_c_theory(m+1) <= max(d_scan)
        xline(d_c_theory(m+1), 'k--', 'HandleVisibility', 'off');
    end
end
ylim([0 max(Nmodes)+0.5]);
xlabel('core thickness d (\mum)');
ylabel('# of guided TE modes');
title('t03  mode count vs thickness (dashed = analytic cutoff)');

exportgraphics(fig, fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'out', 't03_thickness.png'), 'Resolution', 150);
close(fig);

fprintf('\n── 结论 ──\n');
report(all_pass, '各阶模的截止厚度与解析预测一致（相对偏差 < 1e-4）');
if mono && all_pass
    fprintf('t03 通过。\n');
else
    fprintf('t03 未通过。\n');
end
end

% ══════════════════════════════════════════════════════════════════════
function nm = count_modes(n_sub, n_core, n_cov, lam, d)
L  = build_layers(layers_three(n_sub, n_core, n_cov, d));
ne = solve_analytic3(L, lam, struct('pol', 'TE', 'n_modes', 8));
nm = numel(ne);
end

% ══════════════════════════════════════════════════════════════════════
function d_c = find_cutoff(n_sub, n_core, n_cov, lam, m, d_lo, d_hi)
need = m + 1;
if count_modes(n_sub, n_core, n_cov, lam, d_lo) >= need
    d_c = d_lo;  return;
end
if count_modes(n_sub, n_core, n_cov, lam, d_hi) < need
    d_c = NaN;   return;
end
for it = 1:45
    dm = 0.5*(d_lo + d_hi);
    if count_modes(n_sub, n_core, n_cov, lam, dm) >= need
        d_hi = dm;
    else
        d_lo = dm;
    end
end
d_c = 0.5*(d_lo + d_hi);
end

% ══════════════════════════════════════════════════════════════════════
function report(ok, msg)
if ok
    fprintf('  [PASS] %s\n', msg);
else
    fprintf('  [FAIL] %s\n', msg);
end
end
