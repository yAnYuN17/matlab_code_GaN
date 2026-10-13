function t02_vnumber_dispersion()
%T02_VNUMBER_DISPERSION  A3：用归一化频率 V 验证色散关系
%
%   ── 为什么 A3 比 A1/A2 更重要 ───────────────────────────────
%   A1/A2 依赖具体的折射率取值，而 GaN/AlGaN 在 450 nm 的折射率有相当的
%   不确定性（文献值分散、随组分与温度变）。绝对 n_eff 对不上，并不能
%   说明求解器错了。A3 把问题无量纲化：
%       V = k0·d·√(n_core² − n_high²)                    归一化频率
%       b = (n_eff² − n_high²)/(n_core² − n_high²)        归一化传播常数
%       a = (n_high² − n_low²)/(n_core² − n_high²)        不对称度
%   其中 n_high / n_low = 两侧包层折射率的较大 / 较小者。
%   这三者只由 a 决定，与折射率的绝对取值无关——所以 b(V; a) 是普适曲线，
%   与"通用色散曲线"对照才有意义，也才不受文献值不确定性的影响。
%
%   ⚠ 非对称结构下"A3 与通用色散曲线吻合"必须精确化：
%     教科书那条单一曲线是 a = 0（对称）的特例。非对称时要画一族以 a 为
%     参数的曲线族。不做这个区分，A3 无从判定。
%
%   ── 本脚本验什么 ────────────────────────────────────────────
%   1) 画出 b(V) 曲线族（a = 0 / 0.97 / 2.83），检查形状随 a 单调变化：
%      a 越大 ⇒ 同一 V 下 b 越小（限制越弱），且基模截止 V_c,0 越大。
%   2) 定量核对截止条件  V_c,m = arctan(√a) + mπ。
%      这一条特别有价值：截止附近 n_eff → n_high，两个 γ 中有一个趋于零，
%      是色散函数最易出数值问题（矩阵退化、根变密）的地方。
%      能在这里对上，说明求解器在整个束缚区间都可靠。

fprintf('\n========== t02 V 数色散曲线（A3）==========\n');

n_core = 3.50;
n_sub  = 3.40;
lam    = 1.0;
n_cov_list  = [3.40, 3.30, 3.10];
d_scan      = 0.20:0.02:3.00;
n_modes_max = 5;

fig = figure('Visible', 'off', 'Position', [80 80 1000 650]);
hold on; grid on; box on;
colors = lines(numel(n_cov_list));

all_pass = true;

for ic = 1:numel(n_cov_list)
    n_cov  = n_cov_list(ic);
    n_high = max(n_sub, n_cov);
    n_low  = min(n_sub, n_cov);
    a      = (n_high^2 - n_low^2)/(n_core^2 - n_high^2);
    Vc     = atan(sqrt(a)) + (0:n_modes_max-1)*pi;

    fprintf('\n── a = %.4f  (n_cov = %.2f) ──\n', a, n_cov);
    fprintf('   解析截止 V_c,m = arctan(sqrt(a)) + m*pi : %s\n', mat2str(Vc, 6));

    Vv = [];  Bv = [];
    for d = d_scan
        L  = build_layers(layers_three(n_sub, n_core, n_cov, d));
        ne = solve_analytic3(L, lam, struct('pol', 'TE', 'n_modes', n_modes_max));
        for k = 1:numel(ne)
            Vv(end+1) = 2*pi/lam * d * sqrt(n_core^2 - n_high^2);   %#ok<AGROW>
            Bv(end+1) = (real(ne(k))^2 - n_high^2)/(n_core^2 - n_high^2); %#ok<AGROW>
        end
    end
    plot(Vv, Bv, '.', 'MarkerSize', 7, 'Color', colors(ic,:), ...
         'DisplayName', sprintf('a = %.3f', a));

    for m = 0:n_modes_max-1
        if Vc(m+1) < 1e-9
            fprintf('   m=%d  对称情形，基模无截止（V_c = 0）—— 符合预期\n', m);
            continue;
        end
        d_c = find_cutoff(n_sub, n_core, n_cov, lam, m, 0.05, 4.0);
        if isnan(d_c)
            fprintf('   m=%d  [WARN] 在扫描范围内未找到截止\n', m);
            continue;
        end
        Vc_num = 2*pi/lam * d_c * sqrt(n_core^2 - n_high^2);
        err    = abs(Vc_num - Vc(m+1))/Vc(m+1);
        if err > 1e-4, all_pass = false; end
        fprintf('   m=%d  数值 V_c = %.6f   解析 V_c = %.6f   相对偏差 %.2e\n', ...
                m, Vc_num, Vc(m+1), err);
        xline(Vc(m+1), '--', 'Color', colors(ic,:), 'HandleVisibility', 'off');
    end
end

xlabel('V = k_0 d (n_{core}^2 - n_{high}^2)^{1/2}');
ylabel('b = (n_{eff}^2 - n_{high}^2) / (n_{core}^2 - n_{high}^2)');
title('t02  Normalized dispersion curves b(V; a)');
legend('Location', 'southeast');
exportgraphics(fig, fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'out', 't02_vnumber.png'), 'Resolution', 150);
close(fig);

fprintf('\n── 结论 ──\n');
report(all_pass, '所有模序的截止 V_c 与解析式 arctan(√a)+mπ 一致（相对偏差 < 1e-4）');
if all_pass
    fprintf('t02 通过：色散关系在无量纲化后与普适曲线族一致。\n');
else
    fprintf('t02 未通过。\n');
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
% 求使第 m 阶模（0-based）刚好出现的芯层厚度：模数达到 m+1 的最小厚度
need = m + 1;

if count_modes(n_sub, n_core, n_cov, lam, d_lo) >= need
    d_c = d_lo;
    return;
end
if count_modes(n_sub, n_core, n_cov, lam, d_hi) < need
    d_c = NaN;
    return;
end

for it = 1:40
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
