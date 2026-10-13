function run_stage1(mode)
%RUN_STAGE1  阶段1 总驱动：一次跑完全部验证并生成交付图
%
%   用法（在 Mode_L1 目录下）：
%       startup
%       run_stage1            % 全部：验证 + 出图
%       run_stage1('figs')    % 只重新生成交付图（改了绘图参数时用）
%
%   ── 产出（对应交接报告 §7 的交付清单）───────────────────────
%     1. 可运行的求解器代码          src/
%     2. 三条路径对比图              out/fig1_three_way.png    （A1/A2）
%     3. V 数色散曲线                out/t02_vnumber.png       （A3）
%     4. 厚度扫描                    out/t03_thickness.png     （A4）
%     5. 模场剖面 + Γ                out/fig3_profile.png      （A6/A7）
%     6. 本阶段说明                  README.md
%
%   ── 运行时间 ────────────────────────────────────────────────
%   约 3–5 分钟。绝大部分时间花在 FDM 路径（t01 要构建并求解两个
%   12000 / 23999 阶的稀疏本征值问题，每个用例约 20–40 s）。
%   若只想快速冒烟测试，单独跑 t00 / t02 / t03 / t05 / t06。

if nargin < 1 || isempty(mode), mode = 'all'; end

fprintf('\n');
fprintf('==============================================================\n');
fprintf('  GaN 基 LD 光学结构设计 —— 阶段1\n');
fprintf('  非对称多层平板波导本征方程\n');
fprintf('==============================================================\n');

root   = fileparts(mfilename('fullpath'));
outdir = fullfile(root, 'out');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% ═══ 验证 ═════════════════════════════════════════════════════
if ~strcmpi(mode, 'figs')
    t00_symbolic_check();
    t01_three_way_crosscheck();
    t02_vnumber_dispersion();
    t03_thickness_sweep();
    t04_normalization();
    t05_asymmetry_physics();
    t06_conventions();
end

% ═══ 交付图 ═══════════════════════════════════════════════════
fig1_three_way(outdir);
fig3_profile(outdir);

fprintf('\n');
fprintf('==============================================================\n');
fprintf('  阶段1 完成。产出图见：%s\n', outdir);
fprintf('  说明文档见：%s\n', fullfile(root, 'README.md'));
fprintf('==============================================================\n\n');
end

% ══════════════════════════════════════════════════════════════════════
function fig1_three_way(outdir)
% 三条独立路径的 n_eff 对比（A0/A1/A2）
% analytic 与 TMM 画密集曲线（快），FDM 画散点（慢，只取几个厚度）

n_sub = 3.40;  n_core = 3.50;  n_cov = 3.30;
lam   = 1.0;

d_dense = 0.30:0.05:3.00;
d_fdm   = [0.50 1.00 1.50 2.00 2.50 3.00];

n_ord = 4;
A = nan(numel(d_dense), n_ord);
T = nan(numel(d_dense), n_ord);
F = nan(numel(d_fdm),   n_ord);

n_lo = max(n_sub, n_cov) + 1e-9;

for i = 1:numel(d_dense)
    L = build_layers(layers_three(n_sub, n_core, n_cov, d_dense(i)));

    a = solve_analytic3(L, lam, struct('pol','TE','n_modes',n_ord));
    t = find_roots_complex(@(x) dispersion_tmm(x, L, lam, 'TE'), ...
                           n_lo, n_core - 1e-9, n_ord + 2, struct());

    for k = 1:min(n_ord, numel(a)), A(i,k) = real(a(k)); end
    for k = 1:min(n_ord, numel(t)), T(i,k) = real(t(k)); end
end

for i = 1:numel(d_fdm)
    L = build_layers(layers_three(n_sub, n_core, n_cov, d_fdm(i)));
    f = solve_fdm(L, lam, n_ord, 'TE', struct());
    for k = 1:min(n_ord, numel(f)), F(i,k) = real(f(k)); end
end

fig = figure('Visible','off','Position',[80 80 980 560]);
hold on; grid on; box on;
cols = lines(n_ord);

for k = 1:n_ord
    plot(d_dense, A(:,k), '-',  'LineWidth', 2.2, 'Color', cols(k,:), ...
         'DisplayName', sprintf('TE%d  analytic', k-1));
    plot(d_dense, T(:,k), '--', 'LineWidth', 1.1, 'Color', cols(k,:), ...
         'DisplayName', sprintf('TE%d  TMM', k-1));
    plot(d_fdm,   F(:,k), 'o',  'MarkerSize', 8, 'LineWidth', 1.2, ...
         'MarkerEdgeColor', 'k', 'MarkerFaceColor', cols(k,:), ...
         'DisplayName', sprintf('TE%d  FDM', k-1));
end

yline(max(n_sub,n_cov), 'k:', 'DisplayName', 'cladding index');
xlabel('core thickness d (\mum)');
ylabel('n_{eff}');
title('fig1  Three independent paths: analytic / TMM / FDM');
legend('Location', 'eastoutside', 'FontSize', 7);
exportgraphics(fig, fullfile(outdir, 'fig1_three_way.png'), 'Resolution', 150);
close(fig);

fprintf('\n[交付图] fig1_three_way.png 已生成\n');
end

% ══════════════════════════════════════════════════════════════════════
function fig3_profile(outdir)
% GaN LD 典型层序下的模场剖面与 Γ（A6/A7）

lam = 0.45;
L   = build_layers(layers_gan_ld_typical());
m   = solve_slab(L, lam, struct('pol','TE','n_modes',3));

fig = figure('Visible','off','Position',[80 80 1050 520]);
hold on; grid on; box on;

cols = lines(numel(m));
for k = 1:numel(m)
    plot(m(k).y, abs(m(k).E).^2, '-', 'LineWidth', 1.8, 'Color', cols(k,:), ...
         'DisplayName', sprintf('TE%d  (n_{eff}=%.5f, \\Gamma=%.4f)', ...
                                k-1, real(m(k).n_eff), m(k).Gamma));
end

yb = m(1).yb;
% 窗口必须覆盖【整个外延结构】：有效芯区（InGaN 波导 + QW）在
% y ≈ 1.20–1.41 之间，峰值也在那里。窗口开小了会把峰切在边缘，
% 看起来像"模场往 p 侧发散"，实则只是没画全。
ymin = -0.15;
ymax = yb(end,1) + 0.06;         % 顶层下沿再往上一点，露出 p-clad

yl = ylim;                       % 先用曲线定下纵轴范围

for k = 1:numel(L)
    lo = max(yb(k,1), ymin);
    hi = min(yb(k,2), ymax);
    if hi <= lo, continue; end
    xline(lo, ':', 'Color', [0.65 0.65 0.65], 'HandleVisibility','off');
    if L(k).is_active
        % 有源区用底色标出。用 xregion 而不是 patch：它自动铺满当前纵轴、
        % 且不改变 ylim。用 patch 给固定高度（如 1e3）会把纵轴顶到 1000，
        % 三条曲线全被压成一条贴地直线——本项目开发中踩过这个坑。
        xregion(lo, hi, 'FaceColor', [1 0.85 0.40], 'FaceAlpha', 0.55, ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    text(0.5*(lo+hi), yl(2)*0.02, L(k).tag, 'Rotation', 90, 'FontSize', 7, ...
         'HorizontalAlignment','center', 'VerticalAlignment','bottom');
end
xline(yb(end,1), ':', 'Color', [0.65 0.65 0.65], 'HandleVisibility','off');

ylim(yl);
xlim([ymin ymax]);
xlabel('y (\mum)');
ylabel('|\phi(y)|^2   (normalized: \int w dy = 1)');
title('fig3  Mode profiles of a typical GaN LD stack  (\lambda = 450 nm, TE)');
legend('Location', 'northeast');
exportgraphics(fig, fullfile(outdir, 'fig3_profile.png'), 'Resolution', 150);
close(fig);

fprintf('\n[交付图] fig3_profile.png 已生成\n');

fprintf('\n── 典型 GaN LD 层序的模场与 Γ ──\n');
for k = 1:numel(m)
    fprintf('   TE%d  n_eff = %.8f    Γ = %.6f\n', k-1, real(m(k).n_eff), m(k).Gamma);
end
end
