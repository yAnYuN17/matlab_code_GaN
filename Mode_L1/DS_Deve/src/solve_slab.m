function modes = solve_slab(layers, lambda, opts)
%SOLVE_SLAB  非对称多层平板波导本征模求解器（阶段1 主入口）
%
%   modes = solve_slab(layers, lambda, opts)
%
%   ── 输入 ────────────────────────────────────────────────────
%   layers : struct 数组，字段
%       .d          厚度 [μm]，【首尾两层必须为 Inf】表示半无限
%       .n          复折射率 n' + 1i·n''（阶段1 取 n'' = 0）
%       .tag        层名（可省略）
%       .is_active  是否属于有源区，决定 Γ 分子的积分范围（可省略）
%                   层序：L(1) 在底（衬底侧），L(end) 在顶（覆盖层侧）
%
%   lambda : 真空波长 [μm]
%
%   opts   : （全部可省略，有默认值）
%       .pol      'TE' | 'TM'                     默认 'TE'
%       .n_modes  最多返回多少阶模                  默认 5
%       .solver   'analytic3' | 'tmm' | 'fdm'     默认 'tmm'
%       .n_scan   根搜索扫描点数（tmm）             默认 20000
%       .npts     模场重建每层采样点数               默认 800
%
%   ── 输出 ────────────────────────────────────────────────────
%   modes : struct 数组，按 n_eff 降序（modes(1) 是基模），字段
%       .n_eff     复有效折射率            ← 从第一天就是复数
%       .beta      复传播常数 [μm^-1] = k0·n_eff
%       .order     模序 m（0 = 基模）
%       .y, .E     模场剖面（E 已归一化，∫w dy = 1）
%       .n_y       与 y 对应的折射率
%       .Gamma     限制因子，定义见 confinement_factor.m 首行
%       .norm_val  归一化积分值（应恒为 1，用于自检）
%       .pol, .lambda, .solver, .yb
%
%   ── 为后续阶段预留的接口（现在不实现，但结构已支持）─────────
%   (1) 复折射率：layers.n 与 modes.n_eff 全程是 complex。阶段2 只需给
%       Mg 掺杂层与衬底赋虚部，本文件的接口一行都不用改。届时
%       find_roots_complex 会自动切换到复根分支。
%   (2) 开放边界 / PML：首尾层已允许"半无限 + 复折射率"。阶段2 要算衬底
%       泄漏损耗时，在 layers 两端各插一层 PML，并在 dispersion_tmm 里
%       换成散射矩阵（S 矩阵）法即可——本文件通过 opts.solver 分派，
%       加一个 'tmm_smatrix' 分支就能接上。
%   (3) 有源区标记：.is_active 已经把 Γ 的分子范围数据化了。阶段4 引入
%       增益后，Γ·g 直接对 is_active 层求和即可。
%
%   ── 三条求解路径（A0 三路交叉验证）─────────────────────────
%   'analytic3'  三层解析超越方程 + 二分（路径1）。只对严格三层有效。
%   'tmm'        N 层传输矩阵 det 式零点 + 扫描（路径2）。任意层数。
%   'fdm'        1D 有限差分 + eigs（路径3）。不经过色散方程。
%
%   生产计算用 'tmm'；'analytic3' 与 'fdm' 是它的两把标尺。

if nargin < 3, opts = struct(); end
if ~isfield(opts, 'pol'),     opts.pol     = 'TE';  end
if ~isfield(opts, 'n_modes'), opts.n_modes = 5;     end
if ~isfield(opts, 'solver'),  opts.solver  = 'tmm'; end

L   = build_layers(layers);
k0  = 2*pi/lambda;
N   = numel(L);
yb  = layer_bounds(L);

n_all      = real([L.n]);
n_clad_max = max(real(L(1).n), real(L(N).n));
n_max      = max(n_all);

% 束缚模的搜索区间：两侧包层中折射率较高者 < n_eff < 结构内最大折射率
tiny = 1e-9 * max(1, n_max);
n_lo = n_clad_max + tiny;
n_hi = n_max      - tiny;

assert(n_hi > n_lo, ...
    ['solve_slab: 搜索区间为空（n_lo = %g ≥ n_hi = %g）。\n' ...
     '  这通常意味着结构里没有任何一层能当波导芯：\n' ...
     '  束缚模要求存在折射率高于【两侧半无限包层】的有限层。'], n_lo, n_hi);

switch lower(opts.solver)
    case 'analytic3'
        neff_list = solve_analytic3(L, lambda, opts);   % 见 src/solve_analytic3.m

    case 'tmm'
        fh = @(x) dispersion_tmm(x, L, lambda, opts.pol);
        % 多找 2 个根再截断：多层结构的色散函数不保证单调，偶尔会在
        % 真实模式附近出现数值伪根。多留余量可以避免伪根把真模挤掉。
        neff_list = find_roots_complex(fh, n_lo, n_hi, opts.n_modes + 2, opts);
        if numel(neff_list) > opts.n_modes
            neff_list = neff_list(1:opts.n_modes);
        end

    case 'fdm'
        neff_list = solve_fdm(L, lambda, opts.n_modes, opts.pol, opts);

    otherwise
        error('solve_slab: 未知的 solver "%s"（可选 analytic3 / tmm / fdm）', opts.solver);
end

neff_list = neff_list(:).';

% 按 n_eff 降序重排（TE0 有最大的 n_eff）
[~, ord]  = sort(real(neff_list), 'descend');
neff_list = neff_list(ord);

if isempty(neff_list)
    warning('solve_slab: 未找到任何束缚模（搜索区间 [%g, %g]）', n_lo, n_hi);
end

% ═══ 逐模重建模场、归一化、算 Γ ═══════════════════════════════
modes = struct([]);
for j = 1:numel(neff_list)
    ne = neff_list(j);

    [y, E, n_y] = mode_profile(ne, L, lambda, opts.pol, opts);
    [E, norm_val] = normalize_mode(y, E, n_y, opts.pol);
    [G, det] = confinement_factor(y, E, n_y, yb, L, opts.pol);

    modes(j).n_eff    = ne;
    modes(j).beta     = k0 * ne;
    modes(j).order    = j - 1;
    modes(j).y        = y;
    modes(j).E        = E;
    modes(j).n_y      = n_y;
    modes(j).Gamma    = G;
    modes(j).norm_val = norm_val;
    modes(j).pol      = opts.pol;
    modes(j).lambda   = lambda;
    modes(j).solver   = opts.solver;
    modes(j).yb       = yb;
    modes(j).P_layer  = det.P_layer;
end
end
