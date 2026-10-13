function crosscheck_v1_layers()
%CROSSCHECK_V1_LAYERS  对 Mode_L1_V1.0/ 的两份层结构做独立核查
%
%   2026-09-26 新增。用途：补交接文件 §4.1.1 第 3 项的缺口 ——
%   Mode_L1_V1.0/ 的两份层结构（layers_gan_ld_ebl / layers_gan_ld_mqw9）
%   从未跑过 DS_Deve 的验证，而阶段 2 要在它们上面做。
%
%   本脚本【只读】现有文件，不修改任何交付物、不写任何图。
%
%   ── 核查项 ──────────────────────────────────────────────────
%     P1  7 层结构：tmm 与 fdm 两条独立路径互验 + 自洽量
%     P2  9 层结构：同上
%     P3  副本同步：Model_S1_7L_EBL.m 与 layers_gan_ld_ebl.m + solve_slab
%         ⚠ 二者是【同一份代码的两个副本】（V1.0/README §5.6 自述），
%           故本项属"同步性检查"，【不是】独立验证。
%     P4  三层退化：analytic3 / tmm / fdm 三路径互验
%         analytic3 只对严格三层有效，故 P1/P2 拿不到第三路 —— 这是
%         代码的结构性限制，不是本脚本的取舍。
%     P5  层分割不变性：把 n-wave 一分为二（同折射率），物理结构不变，
%         tmm 结果必须逐位不变。这是对 TMM 连乘在真实多层结构上的检查。
%     P6  fdm 的 Richardson 收敛比（三档网格），应 ≈ 4
%
%   用法（在 Mode_L1 目录下）：
%       addpath('Mode_L1_V1.0'); crosscheck_v1_layers

fprintf('\n===== crosscheck_v1_layers =====\n');
fprintf('MATLAB %s\n', version);

here = fileparts(mfilename('fullpath'));
m1   = fileparts(here);
addpath(fullfile(m1, 'DS_Deve', 'src'));
addpath(fullfile(m1, 'DS_Deve', 'layers'));
addpath(fullfile(here, 'layers'));
addpath(here);

lam = 0.45;          % 450 nm
NM  = 3;
optT = struct('solver', 'tmm',       'n_modes', NM);
optF = struct('solver', 'fdm',       'n_modes', NM);
optA = struct('solver', 'analytic3', 'n_modes', NM);

np = 0; nf = 0;

% ══════════════════════════════════════════════════════════════
%  P1 / P2  两份层结构：tmm vs fdm
% ══════════════════════════════════════════════════════════════
cname = {'7-layer  layers_gan_ld_ebl', '9-layer  layers_gan_ld_mqw9'};
cget  = {@layers_gan_ld_ebl, @layers_gan_ld_mqw9};

for c = 1:2
    fprintf('\n--- P%d  %s ---\n', c, cname{c});
    L = build_layers(cget{c}());
    n_clad = max(real(L(1).n), real(L(end).n));
    n_max  = max(real([L.n]));
    fprintf('    layers = %d   n_eff window = (%.6f, %.6f)\n', ...
        numel(L), n_clad, n_max);

    A = solve_slab(L, lam, optT);
    B = solve_slab(L, lam, optF);

    fprintf('    %-6s %-16s %-16s %-10s %-12s\n', ...
        'mode', 'tmm', 'fdm', 'd_n/n', 'Gamma(tmm)');
    worst = 0;
    kmax  = min(numel(A), numel(B));
    for k = 1:kmax
        d = abs(real(A(k).n_eff) - real(B(k).n_eff)) / real(A(k).n_eff);
        worst = max(worst, d);
        fprintf('    TE%-4d %-16.10f %-16.10f %-10.2e %-12.6f\n', ...
            k-1, real(A(k).n_eff), real(B(k).n_eff), d, A(k).Gamma);
    end
    [np, nf] = chk(np, nf, worst < 1e-6, ...
        sprintf('P%d tmm-vs-fdm  max|dn|/n = %.2e  (< 1e-6)', c, worst));
    [np, nf] = chk(np, nf, numel(A) == numel(B), ...
        sprintf('P%d mode count: tmm=%d  fdm=%d', c, numel(A), numel(B)));

    for k = 1:numel(A)
        tag = sprintf('P%d TE%d', c, k-1);
        [np, nf] = chk(np, nf, abs(A(k).norm_val - 1) < 1e-12, ...
            sprintf('%s  norm_val = 1  (got %.6e)', tag, A(k).norm_val));
        s = sum(A(k).P_layer);
        [np, nf] = chk(np, nf, abs(s - 1) < 1e-12, ...
            sprintf('%s  sum(P_layer) = 1  (got %.16f)', tag, s));
        ga = sum(A(k).P_layer([L.is_active]));
        [np, nf] = chk(np, nf, abs(ga - A(k).Gamma) < 1e-12, ...
            sprintf('%s  Gamma == sum(P_active)  (diff %.2e)', tag, abs(ga - A(k).Gamma)));
        [np, nf] = chk(np, nf, ...
            real(A(k).n_eff) > n_clad && real(A(k).n_eff) < n_max, ...
            sprintf('%s  n_eff inside window', tag));
    end

    if c == 1, A7 = A; L7 = L; end
end

% ══════════════════════════════════════════════════════════════
%  P3  副本同步：Model_S1_7L_EBL (单文件) vs solve_slab (多文件)
% ══════════════════════════════════════════════════════════════
fprintf('\n--- P3  copy-sync: Model_S1_7L_EBL vs solve_slab(tmm) ---\n');
fprintf('    NOTE: same code, two copies -> sync check, NOT independent\n');

oldvis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');    % 不弹窗、不写盘
C = Model_S1_7L_EBL();
set(0, 'DefaultFigureVisible', oldvis);
close all;

dmax = 0; gmax = 0;
for k = 1:min(numel(A7), numel(C))
    d = abs(real(A7(k).n_eff) - real(C(k).n_eff));
    g = abs(A7(k).Gamma     - C(k).Gamma);
    dmax = max(dmax, d); gmax = max(gmax, g);
    fprintf('    TE%d  multi=%.10f  single=%.10f  dn=%.2e  dGamma=%.2e\n', ...
        k-1, real(A7(k).n_eff), real(C(k).n_eff), d, g);
end
[np, nf] = chk(np, nf, dmax == 0, ...
    sprintf('P3 n_eff bit-identical (max dn = %.2e)', dmax));
[np, nf] = chk(np, nf, gmax == 0, ...
    sprintf('P3 Gamma bit-identical (max dGamma = %.2e)', gmax));

% ══════════════════════════════════════════════════════════════
%  P4  三层退化：analytic3 / tmm / fdm 三路径
% ══════════════════════════════════════════════════════════════
fprintf('\n--- P4  3-layer degenerate: analytic3 / tmm / fdm ---\n');
% 一个合法的三层结构：衬底 / 1.429 um 芯 / p-clad 半无限
% （1.429 = 7 层结构全部有限层厚度之和；芯折射率取波导值 2.455）
L3 = build_layers(struct( ...
    'd',         { Inf,   1.429,  Inf   }, ...
    'n',         { 1.770, 2.455,  2.374 }, ...
    'tag',       {'sub', 'core', 'cov'  }, ...
    'is_active', { false, true,  false  }));

P = solve_slab(L3, lam, optA);
Q = solve_slab(L3, lam, optT);
R = solve_slab(L3, lam, optF);

w12 = 0; w13 = 0;
kmax = min([numel(P), numel(Q), numel(R)]);
for k = 1:kmax
    e12 = abs(real(P(k).n_eff) - real(Q(k).n_eff)) / real(P(k).n_eff);
    e13 = abs(real(P(k).n_eff) - real(R(k).n_eff)) / real(P(k).n_eff);
    w12 = max(w12, e12); w13 = max(w13, e13);
    fprintf('    TE%d  analytic=%.10f  tmm=%.10f  fdm=%.10f  (d12=%.1e d13=%.1e)\n', ...
        k-1, real(P(k).n_eff), real(Q(k).n_eff), real(R(k).n_eff), e12, e13);
end
[np, nf] = chk(np, nf, w12 < 1e-9, ...
    sprintf('P4 max|analytic-tmm|/n = %.2e  (< 1e-9)', w12));
[np, nf] = chk(np, nf, w13 < 1e-6, ...
    sprintf('P4 max|analytic-fdm|/n = %.2e  (< 1e-6)', w13));

% ══════════════════════════════════════════════════════════════
%  P5  层分割不变性：n-wave 0.120 -> 0.060 + 0.060
% ══════════════════════════════════════════════════════════════
fprintf('\n--- P5  layer-split invariance (n-wave 0.120 -> 0.060+0.060) ---\n');
Lsp = build_layers(struct( ...
    'd',         { Inf, 1.200, 0.060,      0.060,      0.009, 0.080, 0.020,  Inf   }, ...
    'n',         { 1.770,2.404, 2.455,      2.455,      2.550, 2.455, 2.343,  2.374 }, ...
    'tag',       {'substrate','n-clad','n-wave1','n-wave2','QW','p-wave','EBL','p-clad' }, ...
    'is_active', { false,false, false,      false,      true,  false, false,  false }));
S = solve_slab(Lsp, lam, optT);

dmax = 0;
for k = 1:min(numel(A7), numel(S))
    d = abs(real(A7(k).n_eff) - real(S(k).n_eff));
    dmax = max(dmax, d);
    fprintf('    TE%d  %d-layer=%.12f  %d-layer=%.12f  dn=%.2e\n', ...
        k-1, numel(L7), real(A7(k).n_eff), numel(Lsp), real(S(k).n_eff), d);
end
[np, nf] = chk(np, nf, dmax < 1e-12, ...
    sprintf('P5 split-invariant n_eff (max dn = %.2e, < 1e-12)', dmax));

% ══════════════════════════════════════════════════════════════
%  P6  fdm Richardson 收敛比（⚠ 必须用严格嵌套的网格序列）
% ══════════════════════════════════════════════════════════════
fprintf('\n--- P6  fdm grid convergence ratio, NESTED grids (expect ~4) ---\n');
% ⚠ 本项第一版写成"比较 fdm(N).coarse 与 fdm(2N).coarse"，得到 ratio = 2.72
%   并误判为 FAIL。根因是【测试设计错，不是求解器错】：
%       solve_fdm(N)   内部网格：N 点 -> h = R/(N-1)
%                      2N-1 点 -> h = R/(2N-2) = 上一个的一半  ← 严格嵌套
%   但 fdm(3000).coarse 是 R/2999，fdm(6000).coarse 是 R/5999，
%   2999 与 5999 之间没有 2 倍关系 —— 不嵌套，O(h^2) 的干净收敛比不成立。
%   这正是《开发坑记录》B3「Richardson 外推要求网格严格嵌套」的同一条。
%   正确序列（来自 solve_fdm 返回的 info，全部严格嵌套）：
%       fdm(3000).n_coarse  ->  R/2999
%       fdm(3000).n_fine    ->  R/5998
%       fdm(6000).n_fine    ->  R/11998
[~,~,~,ia] = solve_fdm(L7, lam, NM, 'TE', struct('N', 3000));
[~,~,~,ib] = solve_fdm(L7, lam, NM, 'TE', struct('N', 6000));
n1 = real(ia.n_coarse); n2 = real(ia.n_fine); n3 = real(ib.n_fine);
mm = min([numel(n1), numel(n2), numel(n3)]);
for k = 1:mm
    den = n2(k) - n3(k);
    if abs(den) > 0
        fprintf('    TE%d  R/2999: %.10f  R/5998: %.10f  R/11998: %.10f  ratio=%.3f\n', ...
            k-1, n1(k), n2(k), n3(k), (n1(k)-n2(k))/den);
    end
end

% ⚠ 观察项，不设 PASS/FAIL。两次运行的结论：
%   1) 【成立】fdm 相对 tmm 的偏差随网格加密而下降（TE0：9.0e-7 -> 7.6e-8 -> 3.0e-8），
%      收敛趋势明确；P1/P2 的 tmm-vs-fdm 一致性（<1e-8）已足够支撑阶段 1 的结论。
%   2) 【不成立】相邻网格的比值非单调（TE0 为 -9.2 / +7.2 / +2.7），
%      干净的 O(h^2) 收敛比【无法】从本数据判出。故"ratio ≈ 4"不能当验收判据。
%      〔推断，未验证〕可能与界面相对网格位置变化导致的误差系数振荡有关。
%      原 A1/A2 用"绝对偏差 < 1e-6"而非收敛比，是正确的判据设计。
fprintf('    [NOTE] convergence-ratio check is OBSERVATIONAL, not a pass/fail criterion\n');
fprintf('           (absolute deviation shrinks 9.0e-7 -> 7.6e-8 -> 3.0e-8, but the\n');
fprintf('            ratio is non-monotonic; do NOT use ratio~4 as acceptance)\n');

% ══════════════════════════════════════════════════════════════
%  P7  各层功率占比全表（供与 V1.0/README §2.3 的表格逐项对照）
% ══════════════════════════════════════════════════════════════
fprintf('\n--- P7  P_layer of 7-layer stack (compare with V1.0/README 2.3) ---\n');
fprintf('    %-12s %-12s %-12s %-12s\n', 'layer', 'TE0[%]', 'TE1[%]', 'TE2[%]');
for k = 1:numel(L7)
    fprintf('    %-12s %-12.4f %-12.4f %-12.4f\n', L7(k).tag, ...
        A7(1).P_layer(k)*100, A7(2).P_layer(k)*100, A7(3).P_layer(k)*100);
end
ss = sum([A7(1).P_layer; A7(2).P_layer; A7(3).P_layer], 2);
fprintf('    %-12s %-12.4f %-12.4f %-12.4f\n', 'SUM', ...
    ss(1)*100, ss(2)*100, ss(3)*100);

% ══════════════════════════════════════════════════════════════
%  P8  TM 分支   ★★ 本项的 FAIL 是真 bug，不是测试缺陷 ★★
% ══════════════════════════════════════════════════════════════
% 结论（已用两条独立证据确认）：
%   src/solve_fdm.m 的签名接收 pol，但【函数体内从未使用】它 ——
%   `grep -n "pol" src/solve_fdm.m` 只命中第 1、4 行（函数签名与其文档行）。
%   它的离散矩阵用的是 TE 的 Helmholtz 形式
%         d²E/dy² + k0²[n²(y) − n_eff²]·E = 0
%   而 TM 应当是
%         d/dy[ (1/n²)·dH_y/dy ] + k0²[1 − n_eff²/n²]·H_y = 0
%   两者的离散算子不同（TM 需 (1/n²) 加权与 n 的导数项）。
%   ⇒ solver='fdm' + pol='TM' 会【静默返回 TE 的结果】。不报错、不警告。
%
% 直接证据（最小对比：layers_three(3.40,3.50,3.20,2.0), lam = 1.0）：
%     fdm TE = 3.493280527698 3.473254411654 3.440568315963
%     fdm TM = 3.493280527698 3.473254411654 3.440568315963   <- 逐位相同
%     tmm TE = 3.493280530768 ...
%     tmm TM = 3.493110306428 ...                             <- TM 与 TE 差 ~1.7e-4
%     max|fdm"TM" − tmm"TE"| = 2.79e-08   ⇒ fdm 的 "TM" 实为 TE
%
% 影响范围：
%   · A1–A7 全部只测 TE；唯一涉及 TM 的 t04_normalization 走的是
%     solve_slab 的默认路径 tmm ⇒ 本 bug 从未被任何测试覆盖，是静默错误。
%   · 当前阶段（TE 主线 + tmm 生产路径）实际影响 = 0。
%   · 危险场景：后续阶段若用 fdm 交叉验证 TM 的 tmm，会看到 ~1e-4 的偏差。
%     那不是 tmm 错，而是 fdm 错；若据此去"修" tmm，会把对的改成错的
%     （与《开发坑记录》A4「两条'独立'路径可能是伪独立」同类）。
%   · 附带缺口：t00_symbolic_check 只推导了 TE 的色散方程，TM 的解析形式
%     从未经符号独立核对。修好 fdm 的 TM 正好补上这条第三路径。
%
% ⚠ 本项不修改任何文件，仅报告。是否修 solve_fdm.m 由用户决定。
fprintf('\n--- P8  TM branch  [FAIL here = REAL BUG in solve_fdm.m] ---\n');
oTm = struct('pol', 'TM', 'solver', 'analytic3', 'n_modes', NM);
uTm = struct('pol', 'TM', 'solver', 'tmm',       'n_modes', NM);
vTm = struct('pol', 'TM', 'solver', 'fdm',       'n_modes', NM);

Pm = solve_slab(L3, lam, oTm);
Qm = solve_slab(L3, lam, uTm);
Rm = solve_slab(L3, lam, vTm);
w12 = 0; w13 = 0;
for k = 1:min([numel(Pm), numel(Qm), numel(Rm)])
    e12 = abs(real(Pm(k).n_eff) - real(Qm(k).n_eff)) / real(Pm(k).n_eff);
    e13 = abs(real(Pm(k).n_eff) - real(Rm(k).n_eff)) / real(Pm(k).n_eff);
    w12 = max(w12, e12); w13 = max(w13, e13);
    fprintf('    3L TM%d  analytic=%.10f  tmm=%.10f  fdm=%.10f  (d12=%.1e d13=%.1e)\n', ...
        k-1, real(Pm(k).n_eff), real(Qm(k).n_eff), real(Rm(k).n_eff), e12, e13);
end
[np, nf] = chk(np, nf, w12 < 1e-9, ...
    sprintf('P8 TM 3-layer max|analytic-tmm|/n = %.2e  (< 1e-9)', w12));
[np, nf] = chk(np, nf, w13 < 1e-6, ...
    sprintf('P8 TM 3-layer max|analytic-fdm|/n = %.2e  (< 1e-6)', w13));

T7 = solve_slab(L7, lam, uTm);
F7 = solve_slab(L7, lam, vTm);
w = 0;
for k = 1:min(numel(T7), numel(F7))
    e = abs(real(T7(k).n_eff) - real(F7(k).n_eff)) / real(T7(k).n_eff);
    w = max(w, e);
    fprintf('    7L TM%d  tmm=%.10f  fdm=%.10f  Gamma=%.6f  (d=%.1e)\n', ...
        k-1, real(T7(k).n_eff), real(F7(k).n_eff), T7(k).Gamma, e);
end
[np, nf] = chk(np, nf, w < 1e-6, ...
    sprintf('P8 TM 7-layer max|tmm-fdm|/n = %.2e  (< 1e-6)', w));
% TE 与 TM 必须不同（若相同，说明 pol 参数被忽略 —— 静默失效）
[np, nf] = chk(np, nf, abs(real(Qm(1).n_eff) - real(Q(1).n_eff)) > 1e-6, ...
    sprintf('P8 TM != TE  (TM0=%.8f  TE0=%.8f)', real(Qm(1).n_eff), real(Q(1).n_eff)));

% ══════════════════════════════════════════════════════════════
%  P9  P_layer 高精度打印（README §2.3 的 EBL/TE0 写 2.5647）
% ══════════════════════════════════════════════════════════════
fprintf('\n--- P9  high-precision P_layer, 7-layer TE0 [%%] ---\n');
for k = 1:numel(L7)
    fprintf('    %-12s %.8f\n', L7(k).tag, A7(1).P_layer(k)*100);
end

fprintf('\n===== summary: %d PASS, %d FAIL =====\n', np, nf);
if nf == 0
    fprintf('crosscheck_v1_layers: OK\n');
else
    fprintf('crosscheck_v1_layers: NOT OK\n');
end
end

% ══════════════════════════════════════════════════════════════
function [np, nf] = chk(np, nf, ok, msg)
if ok
    np = np + 1;
    fprintf('    [PASS] %s\n', msg);
else
    nf = nf + 1;
    fprintf('    [FAIL] %s\n', msg);
end
end
