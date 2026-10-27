% ─────────────────────────────────────────────────────
% 脚本：plot_stage2_figs（生成阶段 2 交付图到 Mode_L2\out）
% 输入：无
% 输出：out\s7_thickness.png、out\s7_substrate.png、out\s8_three_lambda.png
% 调用函数：Model_Self_L2_7L_EBL('solve', opt)
% 说明：
%   1) 不用 yyaxis（E3 坑：右轴会整体盖住左轴），改用上下两个共享 x 的子图。
%   2) 图内文字用英文，避免 -batch 下中文字体缺失时出现方框。
%   3) 衬底扫描图的纵轴必须收紧到 2.4171 附近，否则 Re(n_eff) 的 1e-6 变化看不出来。
%   4) 2026-10-22 起图 1–3 的数据走 'solve' 接口（thick / nsub / lambda），
%      与 test_s7_channels 用同一条通路，图与判定不会分家。
% ─────────────────────────────────────────────────────
base = fileparts(fileparts(mfilename('fullpath')));
addpath(base);
outdir = fullfile(base, 'out');
if ~exist(outdir, 'dir'), mkdir(outdir); end
k0 = 2*pi/0.45;

%% 1) 厚度扫描（α_i 与 Γ 分开画）
tlist = [0.04 0.06 0.08 0.10 0.12 0.15 0.18 0.22 0.26 0.30];
GG = zeros(size(tlist)); AA = zeros(size(tlist));
for it = 1:numel(tlist)
    th = containers.Map({'n-wave','p-wave'}, {tlist(it), tlist(it)});
    S = Model_Self_L2_7L_EBL('solve', struct('scenario', 'kio', 'thick', th));
    GG(it) = S(1).Gamma;
    AA(it) = 2*k0*imag(S(1).ne)*1e4;
end

f = figure('Visible', 'off', 'Position', [100 100 900 700]);
subplot(2,1,1);
plot(tlist, AA, 'o-', 'LineWidth', 1.6); grid on;
ylabel('\alpha_i  [cm^{-1}]'); title('S7-B5: internal loss vs waveguide thickness (KIO, TE0)');
subplot(2,1,2);
plot(tlist, GG, 's-', 'LineWidth', 1.6); grid on;
xlabel('t  [\mum]  (n-wave = p-wave = t)'); ylabel('\Gamma');
title('confinement factor vs thickness');
print(f, fullfile(outdir, 's7_thickness.png'), '-dpng', '-r150'); close(f);

%% 2) 衬底折射率扫描（纵轴收紧 + 标出无根区）
nsub = [1.70:0.05:2.40, 2.400:0.002:2.440];
ne = nan(size(nsub)); bd = false(size(nsub));
for k = 1:numel(nsub)
    S = Model_Self_L2_7L_EBL('solve', struct('scenario', 'kio', 'nsub', nsub(k)));
    if S(1).ok
        ne(k) = real(S(1).ne);
        bd(k) = ne(k) > nsub(k);
    end
end
nb    = find(bd, 1, 'last');
th_v  = nsub(nb);

f = figure('Visible', 'off', 'Position', [100 100 900 620]);
plot(nsub(bd), ne(bd), 'o-', 'LineWidth', 1.6); hold on;
yl = [min(ne(bd)) - 2e-6, max(ne(bd)) + 2e-6];
patch([th_v 2.45 2.45 th_v], [yl(1) yl(1) yl(2) yl(2)], ...
      [1 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
ylim(yl); xlim([1.70 2.45]); grid on;
xlabel('n_{sub}'); ylabel('Re(n_{eff})');
% 阈值线画在 ylim 之外会看不见，改用红色带边界 + 文字说明
% 注记放左上空白区（曲线在右边陡升，左下与右侧都会被压到）
text(1.72, yl(2) - 0.02*diff(yl), {sprintf('last bound n_{sub} = %.3f', th_v), ...
     'no root in window', '(bound-mode solver,', 'leaky regime -> stage 3)'}, ...
     'FontSize', 9, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
title('S7-B6: substrate index scan (red band = bound-mode solver finds no root)');
print(f, fullfile(outdir, 's7_substrate.png'), '-dpng', '-r150'); close(f);

%% 3) 三波长（每层加权贡献 P·α·n_r/Re(n_eff) 堆叠）
%    2026-10-27 修：原来画的是未加权的 P·α，堆叠和 6.1670 与下格 α_i（6.1819）差 0.24%。
%    严格微扰式是 α_i = Σ P_j·α_j·n_rj/Re(n_eff)，所以柱段必须用 S(1).contrib ——
%    这样上图柱高之和与下图 α_i 是同一个数（模型里 contrib 就是这么定义的）。
lams = [0.405 0.450 0.520];
tags = {'substrate','n-clad','n-wave','QW','p-wave','EBL','p-clad'};
PA = zeros(numel(tags), numel(lams)); AM = zeros(1, numel(lams));
for i = 1:numel(lams)
    S = Model_Self_L2_7L_EBL('solve', struct('lambda', lams(i), 'scenario', 'kio'));
    PA(:, i) = S(1).contrib(:);
    AM(i)    = 2*(2*pi/lams(i))*imag(S(1).ne)*1e4;
end
f = figure('Visible', 'off', 'Position', [100 100 900 650]);
subplot(2,1,1);
bar(lams*1000, PA', 'stacked'); grid on;
ylabel('P \cdot \alpha \cdot n_r / Re(n_{eff})  [cm^{-1}]'); legend(tags, 'Location', 'eastoutside');
title('S8: layer contribution to \alpha_i vs wavelength (weighted; n fixed at 450 nm)');
subplot(2,1,2);
plot(lams*1000, AM, 'o-', 'LineWidth', 1.6); grid on;
xlabel('\lambda  [nm]'); ylabel('\alpha_i  [cm^{-1}]');
title('total \alpha_i (trend only for 405/520 nm)');
print(f, fullfile(outdir, 's8_three_lambda.png'), '-dpng', '-r150'); close(f);

% 4) S6：TMM 数值稳定性的图（out\s6_tmm.png）由 test_s6_tmm_stability.m 自己出，
%    S6 搬出模型后出图跟着证据走，本脚本不再重复画。

fprintf('FIG DONE -> %s\n', outdir);
