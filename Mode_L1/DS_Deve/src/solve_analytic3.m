function neff_list = solve_analytic3(L, lambda, opts)
%SOLVE_ANALYTIC3  路径1：逐个模序解三层超越方程
%
%   neff_list = solve_analytic3(L, lambda, opts)
%
%   opts : .pol      'TE' | 'TM'   默认 'TE'
%          .n_modes  最多求多少阶     默认 5
%
%   只对【严格三层】结构有效（L(1) 下半无限包层 / L(2) 芯层 / L(3) 上半无限包层）。
%
%   ── 为什么不需要扫描 ────────────────────────────────────────
%   色散函数关于 n_eff 严格单调递减（证明见 dispersion_analytic3.m 头注释），
%   所以给定模序 m，束缚区间内至多一个根。只需检查区间两端符号再二分：
%   不会漏根、不会跳根、不需要扫描。这是三条路径里最快的一条，
%   因此适合被厚度扫描之类的"外层循环"反复调用。
%
%   ── 模序与 n_eff 的对应 ─────────────────────────────────────
%   m = 0 给出最大的 n_eff（基模），m 越大 n_eff 越小。返回顺序即为降序，
%   与 solve_slab 的输出顺序一致。
%
%   ⚠ 与 TMM 路径的差异：本函数【不做】模场重建，只返回 n_eff。
%     需要模场与 Γ 时请用 solve_slab。

if nargin < 3, opts = struct(); end
if ~isfield(opts, 'pol'),     opts.pol     = 'TE'; end
if ~isfield(opts, 'n_modes'), opts.n_modes = 5;    end

assert(numel(L) == 3, ...
    ['solve_analytic3: 只支持严格三层结构（当前 %d 层）。\n' ...
     '  多层结构请用 solver = ''tmm''。'], numel(L));

k0     = 2*pi/lambda;
n_core = real(L(2).n);
n1     = real(L(1).n);
n3     = real(L(3).n);

if n1 >= n3
    n_high = n1;  n_low = n3;
else
    n_high = n3;  n_low = n1;
end

% ── 精确截止判据（不要换成"在 n_high+δ 处取函数值"）──────────
% 模序 m 存在的条件是色散函数在 n_eff → n_high⁺ 处的极限值 g_m > 0。
% 该极限有闭式解：
%     TE : g_m = κ_max·d − arctan(γ_low/κ_max) − mπ
%     TM : g_m = κ_max·d − arctan((n_core/n_low)²·γ_low/κ_max) − mπ
% 其中 κ_max = k0·√(n_core²−n_high²)，γ_low = k0·√(n_high²−n_low²)。
% TE 情形等价于 V > arctan(√a) + mπ，即经典的截止条件。
%
% ⚠ 不要图省事改成"看 f(n_high + δ) 的符号"。那样得到的截止点会系统性
%   偏高 arctan(γ(δ)/κ) ≈ k0·√(2·n_high·δ)/κ：δ = 1e-9 时偏差约 2e-4。
%   这个量级足以让 t02/t03 的截止核对卡在 1e-4 附近过不去，而它完全是
%   【判据本身引入的假偏差】，与求解器精度无关。本项目开发中确实撞过一次，
%   排查花的时间远超写下这段闭式的时间。
kappa_max = k0*sqrt(n_core^2 - n_high^2);
gam_low   = k0*sqrt(complex(n_high^2 - n_low^2));

if strcmpi(opts.pol, 'TE')
    cut_term = atan(gam_low/kappa_max);
else
    cut_term = atan((n_core/n_low)^2 * gam_low/kappa_max);
end

delta = 1e-13 * max(1, n_high);     % 靠近 n_high 时的最小可分辨间隙
n_hi  = n_core - 1e-12*max(1, n_core);

neff_list = [];
m_max = opts.n_modes + 2;           % 多试几阶，够用即止

for m = 0:m_max
    g = kappa_max*L(2).d - cut_term - m*pi;
    if real(g) <= 0
        break;                       % 精确截止判据：该阶及以上都不存在
    end

    fh = @(x) dispersion_analytic3(x, L, lambda, opts.pol, m);

    % 二分的左端点取 n_high 本身，并用【极限值 g】作为该点的函数值。
    % 不能改成"在 n_high+δ 处取 f 的值"：当 g 本身只有 1e-5 量级时
    % （模刚刚超过截止），f(n_high+δ) 与 g 的差会盖过 g，导致误判。
    a  = n_high;
    fa = g;
    b  = n_hi;

    for it = 1:200
        c = 0.5*(a + b);
        if (c - n_high) <= delta || (b - a) <= 1e-15*max(1, abs(c))
            break;
        end
        fc = fh(c);
        if real(fc) == 0
            break;
        end
        if sign(real(fc)) == sign(real(fa))
            a = c; fa = fc;
        else
            b = c;
        end
    end

    % 根若贴到 n_high 上（模式刚过截止），绝对误差 < delta，相对误差 < 1e-13，
    % 对任何下游用途都足够；此处不再额外报错。
    neff_list(end+1) = 0.5*(a + b);   %#ok<AGROW>
    if numel(neff_list) >= opts.n_modes
        break;
    end
end
end
