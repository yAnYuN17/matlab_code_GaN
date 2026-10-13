function t00_symbolic_check()
%T00_SYMBOLIC_CHECK  用符号运算独立推导三层色散方程（防"自证"）
%
%   ── 这个脚本存在的原因 ──────────────────────────────────────
%   验收标准 A1 要求"三层解析解的 n_eff 与解析解吻合 1e-6"。但注意：
%   【解析解本身也是数值求根得来的】。如果手写的超越方程推导有错
%   （比如把 TM 的界面条件写成 TE 的、或把衰减解写成了增长解），
%   路径1 与路径2 会一起错，A1/A2 双双通过而结论是错的。
%   这就是"自证"漏洞——两条路径再"独立"，只要共享同一条错误推导就没用。
%
%   本脚本从 Maxwell 方程 + 界面连续条件出发，用符号运算【重新推导】
%   一遍色散方程，再与 dispersion_analytic3.m 里手写的那条逐项比对。
%   它验证的是"公式的形式"，与被验证的数值求解器完全独立。
%
%   ── 推导（TE 模，中心坐标）──────────────────────────────────
%   core (−d/2 < y < d/2)： E = A·cos(κy) + B·sin(κy),  κ = k0√(n_core²−n_eff²)
%   sub  (y < −d/2)      ： E = C·exp(γ_s·(y + d/2)),   γ_s = k0√(n_eff²−n_sub²)
%   cov  (y > +d/2)      ： E = D·exp(−γ_c·(y − d/2)),  γ_c = k0√(n_eff²−n_cov²)
%   两侧包层都取【衰减解】（Re γ > 0），这是束缚模的前提。
%
%   在 y = −d/2 与 y = +d/2 分别令 E 与 dE/dy 连续，消去 C、D，
%   得到关于 (A, B) 的两个齐次方程；有非零解要求系数行列式为零。
%
%   记号：S = sin(κd/2), C = cos(κd/2)
%       y = −d/2 消去 C 后：  A·(κS − γ_s·C) + B·(κC + γ_s·S) = 0
%       y = +d/2 消去 D 后：  A·(−κS + γ_c·C) + B·(κC + γ_c·S) = 0

fprintf('\n========== t00 符号推导核对 ==========\n');

syms kap gs gc dthk real

S = sin(kap*dthk/2);
C = cos(kap*dthk/2);

% 2×2 系数行列式
DET = (kap*S - gs*C)*(kap*C + gc*S) ...
    - (kap*C + gs*S)*(-kap*S + gc*C);
DET = simplify(expand(DET));

expected = sin(kap*dthk)*(kap^2 - gs*gc) - kap*(gs + gc)*cos(kap*dthk);
chk = simplify(DET - expected);

fprintf('符号推导行列式      : %s\n', char(DET));
fprintf('预期等价形式        : %s\n', char(expected));
fprintf('两者之差（应为 0）  : %s\n', char(chk));

pass_sym = isequal(chk, sym(0));

% 由 DET = 0 解出 tan(κd) 的形式
%   sin(κd)(κ² − γ_s γ_c) = κ(γ_s + γ_c) cos(κd)
%   ⇒ tan(κd) = κ(γ_s+γ_c)/(κ² − γ_s γ_c)
% 而 arctan(γ_s/κ) + arctan(γ_c/κ) 的正切恰好是同一个式子：
%   (γ_s/κ + γ_c/κ)/(1 − γ_s γ_c/κ²) = κ(γ_s+γ_c)/(κ² − γ_s γ_c)
fprintf('\n由此得 tan(κd) = κ(γ_s+γ_c)/(κ² − γ_s·γ_c)，\n');
fprintf('而 arctan(γ_s/κ)+arctan(γ_c/κ) 的正切是同一个式子，\n');
fprintf('故 κd = arctan(γ_s/κ) + arctan(γ_c/κ) + mπ 成立。\n');

% ── 数值核对：把手写解析解算出的根代回符号行列式，应得 0 ──────
fprintf('\n── 数值核对：解析根代回符号行列式 ──\n');

n_sub = 3.40; n_core = 3.50; n_cov = 3.30; d_core = 2.0; lam = 1.0;
L = build_layers(layers_three(n_sub, n_core, n_cov, d_core));
modes = solve_slab(L, lam, struct('solver', 'analytic3', 'n_modes', 4));
k0 = 2*pi/lam;

max_det = 0;
for k = 1:numel(modes)
    ne    = real(modes(k).n_eff);
    kap_v = k0*sqrt(n_core^2 - ne^2);
    gs_v  = k0*sqrt(ne^2 - n_sub^2);
    gc_v  = k0*sqrt(ne^2 - n_cov^2);

    dv = double(subs(DET, {kap, gs, gc, dthk}, {kap_v, gs_v, gc_v, d_core}));
    max_det = max(max_det, abs(dv));

    fprintf('  m=%d  n_eff=%.12f   行列式=%.3e\n', k-1, ne, dv);
end

% 行列式值的量纲 ~ κ³，用 κ³ 归一化后再判定
kap_ref = k0*sqrt(n_core^2 - n_sub^2);
scale   = kap_ref^3;
pass_num = (max_det/scale) < 1e-10;

fprintf('\n  最大 |行列式| / κ³ = %.3e  （判据 < 1e-10）\n', max_det/scale);

fprintf('\n── 结论 ──\n');
report(pass_sym, '符号推导与预期形式一致');
report(pass_num, '解析根使符号行列式为零');
if pass_sym && pass_num
    fprintf('t00 通过：手写色散方程的形式已被独立验证，A1 不构成自证。\n');
else
    fprintf('t00 失败。\n');
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
