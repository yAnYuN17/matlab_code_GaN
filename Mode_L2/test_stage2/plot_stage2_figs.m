% ─────────────────────────────────────────────────────
% 脚本：plot_stage2_figs（生成阶段 2 交付图到 Mode_L2\out）
% 输入：无
% 输出：out\s7_thickness.png、out\s7_substrate.png、out\s8_three_lambda.png、out\s6_tmm.png
% 调用函数：Model_Self_L2_7L_EBL('s7') / ('s8') / ('s6')
% 说明：
%   1) 不用 yyaxis（E3 坑：右轴会整体盖住左轴），改用上下两个共享 x 的子图。
%   2) 图内文字用英文，避免 -batch 下中文字体缺失时出现方框。
%   3) 衬底扫描图的纵轴必须收紧到 2.4171 附近，否则 Re(n_eff) 的 1e-6 变化看不出来。
% ─────────────────────────────────────────────────────
base = fileparts(fileparts(mfilename('fullpath')));
addpath(base);
outdir = fullfile(base, 'out');
if ~exist(outdir, 'dir'), mkdir(outdir); end

%% 1) S7：厚度扫描（α_i 与 Γ 分开画）
M  = Model_Self_L2_7L_EBL('s7');
tt = [M.thick.t];
GG = [M.thick.G];
AA = [M.thick.am];

f = figure('Visible', 'off', 'Position', [100 100 900 700]);
subplot(2,1,1);
plot(tt, AA, 'o-', 'LineWidth', 1.6); grid on;
ylabel('\alpha_i  [cm^{-1}]'); title('S7-B5: internal loss vs waveguide thickness (KIO, TE0)');
subplot(2,1,2);
plot(tt, GG, 's-', 'LineWidth', 1.6); grid on;
xlabel('t  [\mum]  (n-wave = p-wave = t)'); ylabel('\Gamma');
title('confinement factor vs thickness');
print(f, fullfile(outdir, 's7_thickness.png'), '-dpng', '-r150'); close(f);

%% 2) S7：衬底折射率扫描（纵轴收紧 + 标出无根区）
ns = [M.substrate.nsub];
ne = [M.substrate.ne];
bd = [M.substrate.bound];
f = figure('Visible', 'off', 'Position', [100 100 900 620]);
plot(ns(bd), ne(bd), 'o-', 'LineWidth', 1.6); hold on;
yl = [min(ne(bd)) - 2e-6, max(ne(bd)) + 2e-6];
patch([M.nsub_th 2.45 2.45 M.nsub_th], [yl(1) yl(1) yl(2) yl(2)], ...
      [1 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
yline(M.nsub_th, 'k--', sprintf('last bound n_{sub} = %.3f', M.nsub_th), ...
      'LabelHorizontalAlignment', 'left');
ylim(yl); xlim([1.70 2.45]); grid on;
xlabel('n_{sub}'); ylabel('Re(n_{eff})');
text(2.408, mean(yl), {'no root in window', '(bound-mode solver,', 'leaky regime -> stage 3)'}, ...
     'FontSize', 9, 'HorizontalAlignment', 'right');
title('S7-B6: substrate index scan (red band = bound-mode solver finds no root)');
print(f, fullfile(outdir, 's7_substrate.png'), '-dpng', '-r150'); close(f);

%% 3) S8：三波长（每层 P·α 堆叠）
S  = Model_Self_L2_7L_EBL('s8');
lams = [S.lambda] * 1000;
tags = {'substrate','n-clad','n-wave','QW','p-wave','EBL','p-clad'};
PA = zeros(numel(tags), numel(S));
for i = 1:numel(S)
    PA(:, i) = (S(i).Pl .* S(i).alpha(:)')';
end
f = figure('Visible', 'off', 'Position', [100 100 900 650]);
subplot(2,1,1);
bar(lams, PA', 'stacked'); grid on;
ylabel('P \cdot \alpha  [cm^{-1}]'); legend(tags, 'Location', 'eastoutside');
title('S8: layer loss contribution vs wavelength (n fixed at 450 nm values)');
subplot(2,1,2);
plot(lams, [S.am], 'o-', 'LineWidth', 1.6); grid on;
xlabel('\lambda  [nm]'); ylabel('\alpha_i  [cm^{-1}]');
title('total \alpha_i (trend only for 405/520 nm)');
print(f, fullfile(outdir, 's8_three_lambda.png'), '-dpng', '-r150'); close(f);

%% 4) S6：TMM 数值稳定性（精度 vs 可解性分开画）
R  = Model_Self_L2_7L_EBL('s6');
t1 = R.table1;
al = [t1.alpha];
ad = [t1.ad];
rel = [t1.rel];
found = [t1.found];

f = figure('Visible', 'off', 'Position', [100 100 950 620]);
okm = found & al > 0;
relp = max(rel(okm), 1e-17);          % 精确为 0 的点画在 1e-17 地板上（log 轴不能画 0）
loglog(al(okm), relp, 'o-', 'LineWidth', 1.6); hold on;
ax = gca; ax.XScale = 'log'; ax.YScale = 'log';
bad = ~found;
if any(bad)
    xb = min(al(bad));
    patch([xb 2e5 2e5 xb], [1e-18 1e-18 1e-10 1e-10], [1 0.85 0.85], ...
          'EdgeColor', 'none', 'FaceAlpha', 0.6);
    text(xb*1.1, 3e-14, {'no root found', 'Im(n_eff) ~ 0.1', '(root-finder limit, not matrix)'}, 'FontSize', 9, 'Interpreter', 'none');
end
xline(1200, 'k--', '\alpha d = 0.12');
grid on; xlabel('\alpha  [cm^{-1}]'); ylabel('|\Delta\alpha| / \alpha');
title('S6: TMM accuracy vs \alpha (d = 0.12 \mum; zeros plotted at 1e-17)');
ylim([1e-18 1e-10]);
print(f, fullfile(outdir, 's6_tmm.png'), '-dpng', '-r150'); close(f);

fprintf('FIG DONE -> %s\n', outdir);
