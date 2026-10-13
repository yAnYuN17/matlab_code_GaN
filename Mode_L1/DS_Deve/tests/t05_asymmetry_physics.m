function t05_asymmetry_physics()
%T05_ASYMMETRY_PHYSICS  A6：非对称结构的物理趋势
%
%   验收要求：模场峰值偏向【折射率较高的一侧】，且该侧穿透深度更大。
%
%   ⚠ 交接报告特别强调："应用到真实器件时，偏向哪一侧取决于两侧包层的
%     实际折射率与厚度，必须按真实层序判断，不要预设方向。"
%
%   所以本脚本【不】检查"偏向下方"或"偏向上方"，而是构造一对镜像结构
%   （把上下包层对调），验证：
%        高折射率包层在哪一侧，模场就偏向哪一侧，穿透深度就在哪一侧更大。
%
%   若代码里把方向写死（例如无脑假设"衬底侧折射率高"），镜像测试会立刻
%   失败，而单向测试照样通过——这就是镜像设计的价值所在。
%
%   ── 三条独立判据 ────────────────────────────────────────────
%   (1) 质心位置 y_c：偏到哪一侧；
%   (2) 两侧包层内的功率占比 P_layer：高折射率侧更大；
%   (3) 从数值模场拟合出的衰减常数 γ：高折射率侧更小（穿透 1/γ 更深），
%       并与解析值 k0·√(n_eff²−n_clad²) 对照。
%   三条判据用三种不同的量描述同一件事，互相不共享误差来源。

fprintf('\n========== t05 非对称物理（A6）==========\n');

lam  = 1.00;
n_core = 3.50;
d_core = 1.50;

% 镜像对
LA = build_layers(layers_three(3.40, n_core, 3.20, d_core));  % 高包层在下
LB = build_layers(layers_three(3.20, n_core, 3.40, d_core));  % 高包层在上

cases = {'A  高包层在下 (n_sub=3.40 > n_cov=3.20)', LA, 'down'; ...
         'B  高包层在上 (n_cov=3.40 > n_sub=3.20)', LB, 'up'};

all_pass = true;

for c = 1:size(cases,1)
    name    = cases{c,1};
    L       = cases{c,2};
    expect  = cases{c,3};
    n_sub   = real(L(1).n);
    n_cov   = real(L(end).n);

    m = solve_slab(L, lam, struct('pol','TE','n_modes',1));
    assert(~isempty(m), 't05: %s 未找到束缚模', name);

    Y = m(1).y;
    E = m(1).E;
    w = abs(E).^2;
    P = m(1).P_layer(1);      % 下包层功率占比（L(1) 是半无限底层）
    Pc = m(1).P_layer(3);     % 上包层功率占比

    yc    = trapz(Y, Y.*w) / trapz(Y, w);
    y_mid = d_core/2;

    fprintf('\n── %s ──\n', name);
    fprintf('   n_eff = %.10f\n', real(m(1).n_eff));
    fprintf('   芯层中心 y = %.4f   模场质心 y_c = %.4f\n', y_mid, yc);
    fprintf('   下包层功率占比 = %.6f   上包层功率占比 = %.6f\n', P, Pc);

    % 判据 (1) 质心偏向
    if strcmp(expect,'down')
        ok1 = yc < y_mid;
    else
        ok1 = yc > y_mid;
    end

    % 判据 (2) 功率占比
    if strcmp(expect,'down')
        ok2 = P > Pc;
    else
        ok2 = Pc > P;
    end

    % 判据 (3) 衰减常数：从数值模场拟合，与解析值对照
    g_sub_ana = 2*pi/lam*sqrt(real(m(1).n_eff)^2 - n_sub^2);
    g_cov_ana = 2*pi/lam*sqrt(real(m(1).n_eff)^2 - n_cov^2);
    g_sub_fit = fit_gamma(Y, w, -0.70, -0.10);
    g_cov_fit = fit_gamma(Y, w, d_core+0.10, d_core+0.70);

    fprintf('   下侧 γ: 解析 %.6f  拟合 %.6f （穿透深度 %.4f um）\n', ...
            g_sub_ana, g_sub_fit, 1/g_sub_ana);
    fprintf('   上侧 γ: 解析 %.6f  拟合 %.6f （穿透深度 %.4f um）\n', ...
            g_cov_ana, g_cov_fit, 1/g_cov_ana);

    e_fit = max(abs(g_sub_fit-g_sub_ana)/g_sub_ana, abs(g_cov_fit-g_cov_ana)/g_cov_ana);
    ok3_fit = e_fit < 1e-6;

    if strcmp(expect,'down')
        ok3 = g_sub_ana < g_cov_ana;      % 高折射率侧 γ 小 ⇒ 穿透深
    else
        ok3 = g_cov_ana < g_sub_ana;
    end

    report(ok1, '质心偏向高折射率包层一侧');
    report(ok2, '高折射率包层的功率占比更大');
    report(ok3, '高折射率侧的衰减常数更小（穿透深度更大）');
    report(ok3_fit, sprintf('数值拟合 γ 与解析值一致（最大相对偏差 %.2e）', e_fit));

    if ~(ok1 && ok2 && ok3 && ok3_fit), all_pass = false; end
end

% ═══ 真实 GaN LD 层序 ═════════════════════════════════════════
fprintf('\n── 真实 GaN LD 典型层序（按实际层序判断，不预设方向）──\n');
Lgd = build_layers(layers_gan_ld_typical());
mg  = solve_slab(Lgd, 0.45, struct('pol','TE','n_modes',2));

fprintf('   找到 %d 个束缚模\n', numel(mg));
for k = 1:numel(mg)
    fprintf('   TE%d  n_eff = %.8f   Γ = %.6f\n', ...
            k-1, real(mg(k).n_eff), mg(k).Gamma);
end

if ~isempty(mg)
    P = mg(1).P_layer;
    fprintf('   各层功率占比：\n');
    for k = 1:numel(Lgd)
        fprintf('     %-10s %8.5f%%   n=%.3f  d=%.4f\n', ...
                Lgd(k).tag, P(k)*100, real(Lgd(k).n), Lgd(k).d);
    end
    tags    = {Lgd.tag};
    P_nside = sum(P(strcmp(tags, 'n-clad') | strcmp(tags, 'n-wave')));
    P_pside = sum(P(strcmp(tags, 'p-wave') | strcmp(tags, 'p-clad') | strcmp(tags, 'cover')));
    fprintf('   n 侧（n-clad + n-wave）功率占比    = %.4f%%\n', P_nside*100);
    fprintf('   p 侧（p-wave + p-clad + cover）占比 = %.4f%%\n', P_pside*100);
    if P_nside > P_pside
        fprintf('   → 光偏向 n 侧（本结构折射率较高的一侧），与 A6 判据一致\n');
    else
        fprintf('   → 光偏向 p 侧\n');
    end
end

fprintf('\n── 结论 ──\n');
if all_pass
    fprintf('  [PASS] 镜像测试：模场偏向跟随折射率高的包层，不跟随"上/下"符号\n');
    fprintf('t05 通过。\n');
else
    fprintf('  [FAIL] 镜像测试未通过\n');
    fprintf('t05 未通过。\n');
end
end

% ══════════════════════════════════════════════════════════════════════
function g = fit_gamma(y, w, y_lo, y_hi)
% 在 y ∈ [y_lo, y_hi] 内对 ln(w) 做线性拟合，提取衰减常数 γ。
%   下包层（y<0）: E ∝ exp(γy)   ⇒ ln w ∝ 2γy   ⇒ γ = slope/2
%   上包层（y>d）: E ∝ exp(−γy)  ⇒ ln w ∝ −2γy  ⇒ γ = −slope/2
mask = (y >= y_lo) & (y <= y_hi);
assert(sum(mask) > 10, 'fit_gamma: 拟合区间采样点太少');

yy = y(mask);
ww = w(mask);
p  = polyfit(yy, log(ww), 1);

if y_hi <= 0
    g = p(1)/2;
else
    g = -p(1)/2;
end
end

% ══════════════════════════════════════════════════════════════════════
function report(ok, msg)
if ok
    fprintf('   [PASS] %s\n', msg);
else
    fprintf('   [FAIL] %s\n', msg);
end
end
