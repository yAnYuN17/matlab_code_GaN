function t04_normalization()
%T04_NORMALIZATION  A5：模场归一化 ∫w(y)dy = 1
%
%   ── 归一化的定义（写在 normalize_mode.m 首行，此处复述）────────
%       TE : w(y) = |E_y(y)|²
%       TM : w(y) = |H_y(y)|² / n(y)²
%       目标：∫_{-∞}^{+∞} w(y) dy = 1
%   权重正比于坡印廷矢量的 z 分量（功率密度），所以 TE 与 TM 的权重不同。
%
%   ── 本脚本验三件事 ──────────────────────────────────────────
%   (1) 归一化后积分确实为 1 —— 多种结构、TE/TM 各来一遍；
%   (2) Γ 与归一化【无关】：用未归一化模场与已归一化模场算出的 Γ 必须一致。
%       这条是 Γ 定义自洽性的直接证据（分子分母同阶，归一化因子相消）。
%       若哪天有人把 Γ 的定义改成"除以有源区厚度"之类的非齐次形式，
%       这条会立刻失败。
%   (3) 半无限层的采样延伸足够长：把延伸长度加倍后 Γ 不应有可察觉变化。
%       若这条不成立，说明场尾被截断，Γ 会系统性偏小——
%       而 Γ 偏小会直接传导到阶段4 的模式增益，是那种"看起来正常但数值
%       长期偏了"的错误。

fprintf('\n========== t04 归一化（A5）==========\n');

L3a = build_layers(layers_three(3.40, 3.50, 3.40, 2.00));
L3b = build_layers(layers_three(3.40, 3.50, 3.20, 1.50));
Lgd = build_layers(layers_gan_ld_typical());

names = {'三层对称     (TE)', '三层非对称   (TE)', '三层非对称   (TM)', ...
         'GaN LD 典型   (TE)', 'GaN LD 典型   (TM)'};
Ls    = {L3a,   L3b,   L3b,   Lgd,   Lgd};
lams  = [1.00,  1.00,  1.00,  0.45,  0.45];
pols  = {'TE',  'TE',  'TM',  'TE',  'TM'};

all_pass = true;
worst_norm = 0;
worst_gam  = 0;

for c = 1:numel(names)
    L   = Ls{c};
    lam = lams(c);
    pol = pols{c};

    modes = solve_slab(L, lam, struct('pol', pol, 'n_modes', 3));
    fprintf('\n── %s   (%d 个模) ──\n', names{c}, numel(modes));

    if isempty(modes)
        fprintf('   [WARN] 未找到束缚模，跳过\n');
        continue;
    end

    for k = 1:numel(modes)
        nv = modes(k).norm_val;
        worst_norm = max(worst_norm, abs(nv - 1));
        ok_n = abs(nv - 1) < 1e-10;

        % Γ 与归一化无关
        [y, Eraw, n_y] = mode_profile(real(modes(k).n_eff), L, lam, pol, struct());
        G_raw = confinement_factor(y, Eraw, n_y, modes(k).yb, L, pol);
        G_nrm = modes(k).Gamma;

        if isnan(G_nrm)
            gam_txt = '(结构无有源区，Γ = NaN)';
            ok_g = true;
        else
            dG = abs(G_raw - G_nrm);
            worst_gam = max(worst_gam, dG);
            ok_g = dG < 1e-12;
            gam_txt = sprintf('Γ=%.10f  |Γ_raw−Γ_norm|=%.2e', G_nrm, dG);
        end

        if ~ok_n || ~ok_g, all_pass = false; end
        fprintf('   m=%d  n_eff=%.10f  ∫w dy = %.15f  %s  [%s]\n', ...
                k-1, real(modes(k).n_eff), nv, gam_txt, ternary(ok_n && ok_g));
    end
end

% ── (3) 采样延伸长度：截断充分性 + Γ 稳定性 ───────────────────
fprintf('\n── 半无限层延伸长度（GaN LD 典型层序，TE0）──\n');
modes1 = solve_slab(Lgd, 0.45, struct('pol','TE','n_modes',1,'tail_dec',20));
modes2 = solve_slab(Lgd, 0.45, struct('pol','TE','n_modes',1,'tail_dec',40));

% 判据 A：截断是否真的充分 —— 直接看采样区间端点处的功率密度相对峰值。
% 这是"场尾有没有被切掉"的直接证据，与积分精度无关。
w1   = abs(modes1(1).E).^2;
edge = max([w1(1), w1(end)]) / max(w1);
ok_edge = edge < 1e-15;

% 判据 B：Γ 在延伸长度加倍时的稳定性。
% ⚠ 注意残差的来源：tail_dec 加倍而每层采样点数不变 ⇒ 半无限层的网格步长
%   加倍 ⇒ 复合梯形公式的误差（−(h²/12)∫w''dy）约放大 4 倍。所以这个残差
%   是【积分采样密度】造成的，不是物理截断。观测值 ~3e-10，
%   对应 Γ 的相对变化 ~1e-8 量级。判据因此设在 1e-8 而不是机器精度——
%   把判据设在机器精度只会逼人去调采样点数，而那个数本身没有物理含义。
G1 = modes1(1).Gamma;
G2 = modes2(1).Gamma;
dG = abs(G2 - G1);
ok_ext = dG < 1e-8;

fprintf('   tail_dec=20 → Γ = %.12f\n', G1);
fprintf('   tail_dec=40 → Γ = %.12f\n', G2);
fprintf('   Γ 差值 = %.3e\n', dG);
fprintf('   端点功率密度 / 峰值 = %.3e（判据 < 1e-15，即场尾确已衰减到零）\n', edge);

fprintf('\n── 结论 ──\n');
report(worst_norm < 1e-10, sprintf('A5  归一化后积分与 1 的最大偏差 = %.2e（判据 < 1e-10）', worst_norm));
report(worst_gam  < 1e-12, sprintf('Γ 与归一化无关的最大偏差 = %.2e（判据 < 1e-12）', worst_gam));
report(ok_edge,            sprintf('场尾截断充分：端点功率密度/峰值 = %.2e（判据 < 1e-15）', edge));
report(ok_ext,             sprintf('延伸加倍后 Γ 变化 = %.2e（判据 < 1e-8，残差来自梯形积分采样密度）', dG));

if all_pass && ok_ext
    fprintf('t04 通过。\n');
else
    fprintf('t04 未通过。\n');
end
end

% ══════════════════════════════════════════════════════════════════════
function s = ternary(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

% ══════════════════════════════════════════════════════════════════════
function report(ok, msg)
if ok
    fprintf('  [PASS] %s\n', msg);
else
    fprintf('  [FAIL] %s\n', msg);
end
end
