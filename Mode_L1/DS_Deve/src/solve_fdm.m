function [neff_list, ys, Es, info] = solve_fdm(L, lambda, n_modes, pol, opts)
%SOLVE_FDM  1D 有限差分 + 稀疏本征值求解（路径 3，独立于色散方程）
%
%   [neff_list, ys, Es, info] = solve_fdm(L, lambda, n_modes, pol, opts)
%
%   ── 这条路为什么必须存在 ────────────────────────────────────
%   路径1（解析超越方程）与路径2（TMM）虽然数学路径不同，但都由同一套
%   "分层 + 界面连续 → 色散方程"的推导而来。若推导本身写错（比如 TM 的
%   界面条件写成 TE 的、或 γ 取了增长解分支），两条路径会【一起错】，
%   A1/A2 双双通过而结果是错的。
%
%   路径3 连"色散方程"这个中间产物都不要：直接把 Helmholtz 方程
%        d²E/dy² + k0²[n²(y) − n_eff²]·E = 0
%   在 y 上离散成矩阵本征值问题。它唯一继承的物理前提是 Maxwell 方程本身，
%   与路径1/2 不共享任何一行推导。三方一致的结论，可信度与两方一致
%   不是一个量级。
%
%   ── 离散格式与"亚格子平均"──────────────────────────────────
%   TE 模的弱形式为
%        ∫ (dE/dy)(dv/dy) dy − k0²∫ n²(y)·E·v dy + k0²n_eff²∫ E·v dy = 0
%   在均匀网格上得到三对角矩阵
%        A(j,j) = −2/h² + k0²·⟨n²⟩_j ,  A(j,j±1) = 1/h²
%   其中 ⟨n²⟩_j 是第 j 个控制体 [y_j−h/2, y_j+h/2] 内 n² 的【精确平均值】。
%
%   ⚠ 这个平均不是可选的优化。若图省事直接取网格点所在层的 n²，折射率
%     阶跃界面在网格移动时会引入 O(h) 量级的噪声，收敛阶从 O(h²) 掉下来，
%     Richardson 外推失效，路径3 就永远对不到 1e-6。控制体平均把界面
%     信息解析地积进矩阵，恢复了 O(h²) 与干净的外推。
%
%   ── Richardson 外推 ────────────────────────────────────────
%   二阶中心差分的本征值误差 ~ C·h²。用 N 与 2N−1 两个网格（后者严格包含
%   前者的节点）算两次，外推
%        n_eff ≈ (4·n_eff(h/2) − n_eff(h)) / 3
%   把误差降到 O(h⁴)。info 里保留了外推前后的值，可以现场检查收敛比
%   (n1−n2)/(n2−n3) 是否接近 4——接近就说明 O(h²) 成立、外推可信；
%   不接近说明网格还不够细或边界截断不够远，此时不要相信外推值。
%
%   ── 为什么必须用 shift-invert 而不能用 'lr' / 'la' ──────────
%   矩阵的尺度是 1/h²。取 h ≈ 1.2 nm 时 1/h² ≈ 6.9e5，而待求本征值
%   λ = k0²n_eff² ≈ 484。也就是说：我们想找的是谱的【最高端】，但那一端
%   只占整个谱宽的 ~1e-3，且相邻模式的本征值间隔远小于谱宽。
%   在这种"目标在谱顶、且极度密集"的情形下，ARPACK 在 'lr' / 'la' 模式下
%   的 Krylov 空间分辨不出相邻模式，会直接报 NotAllEigsConverged。
%
%   改用 shift-invert：给一个略高于谱顶的位移 σ，eigs 转而求解
%   (A − σI)⁻¹ 的最大本征值。这一步把谱顶的几个模式映射成"最占优"的几个，
%   同时把它们的间隔大幅拉开，收敛又快又稳。位移只需一个，对三对角矩阵
%   的 LU 分解是 O(N)，代价可忽略。
%
%   ── 边界处理与其代价 ────────────────────────────────────────
%   两端取 Dirichlet (E = 0)。真实结构是半无限的，截断会引入两个误差：
%     (1) 边界反射扰动本征值，量级 ~ exp(−2γ·L_ext)。取 L_ext 使 γL_ext ≳ 30，
%         该项 ~ 1e-26，可忽略。
%     (2) 把连续谱的辐射模离散成一批"假离散模"。它们混在返回结果里，
%         靠"n_eff 必须大于两侧包层折射率"这一物理判据剔除。
%
%   输入 opts
%       .N          粗网格点数，默认 12000
%       .L_ext      两侧延伸长度 [μm]，默认自动（30 个衰减长度，下限 3 μm）
%       .richardson 是否做外推，默认 true
%       .tol        eigs 收敛容差，默认 1e-12

if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'N'),          opts.N = 12000;  end
if ~isfield(opts, 'richardson'), opts.richardson = true; end
if ~isfield(opts, 'tol'),        opts.tol = 1e-12; end

k0 = 2*pi/lambda;
yb = layer_bounds(L);
d_total = yb(end, 1);

n_all     = real([L.n]);
n_clad_lo = real(L(1).n);
n_clad_hi = real(L(end).n);
n_max     = max(n_all);

% 自动定延伸长度：用"包层与结构最大折射率的中点"估一个 γ
if isfield(opts, 'L_ext') && ~isempty(opts.L_ext)
    L_ext_sub = opts.L_ext(1);
    L_ext_cov = opts.L_ext(min(2, numel(opts.L_ext)));
else
    n_est = 0.5*(max(n_clad_lo, n_clad_hi) + n_max);
    L_ext_sub = est_extent(k0, n_est, n_clad_lo);
    L_ext_cov = est_extent(k0, n_est, n_clad_hi);
end

y_min = -L_ext_sub;
y_max =  d_total + L_ext_cov;

% ═══ 两个网格各解一次 ════════════════════════════════════════
[lam1, ~,  ~  ] = fdm_core(opts.N);
[lam2, V2, ys]  = fdm_core(2*opts.N - 1);

n1 = sqrt(lam1)/k0;
n2 = sqrt(lam2)/k0;

if opts.richardson
    m     = min(numel(n1), numel(n2));
    n_use = (4*n2(1:m) - n1(1:m))/3;
    V_use = V2(:, 1:m);
else
    m     = numel(n2);
    n_use = n2;
    V_use = V2;
end

% ═══ 剔除假模（辐射模被 Dirichlet 边界离散化出来的伪解）═══════
n_clad_max = max(n_clad_lo, n_clad_hi);
keep = real(n_use) > n_clad_max + 1e-6;

n_keep    = n_use(keep);
V_keep    = V_use(:, keep);

m_out = min(numel(n_keep), n_modes);
neff_list = n_keep(1:m_out).';
V_out     = V_keep(:, 1:m_out);

% ═══ 组装模场（补回两端 Dirichlet 零点）══════════════════════
Npts   = numel(ys);
Es     = zeros(Npts, m_out);
for j = 1:m_out
    v = zeros(Npts, 1);
    v(2:Npts-1) = V_out(:, j);
    Es(:, j) = real(v) / max(abs(real(v)));
end

info = struct();
info.n_coarse = n1;
info.n_fine   = n2;
info.n_extrap = neff_list;
info.N_coarse = opts.N;
info.N_fine   = 2*opts.N - 1;
info.L_ext    = [L_ext_sub, L_ext_cov];
info.y_range  = [y_min, y_max];
info.richardson = opts.richardson;

% ══════════════════════════════════════════════════════════════
    function [lam, V, ys_] = fdm_core(Npts_req)
        Npts = Npts_req;
        ys_  = linspace(y_min, y_max, Npts);
        h    = ys_(2) - ys_(1);

        % 控制体 [y−h/2, y+h/2] 内 n² 的精确平均
        n2 = zeros(1, Npts);
        ylo = ys_ - h/2;
        yhi = ys_ + h/2;
        for kk = 1:numel(L)
            lo = max(yb(kk,1), y_min);
            hi = min(yb(kk,2), y_max);
            if hi <= lo, continue; end
            ov = max(0, min(yhi, hi) - max(ylo, lo));
            n2 = n2 + (L(kk).n^2) * (ov / h);
        end

        e = ones(Npts, 1);
        A = spdiags([e, -2*e, e], -1:1, Npts, Npts) / h^2 ...
            + spdiags((k0^2 * n2(:)), 0, Npts, Npts);

        A_int = A(2:Npts-1, 2:Npts-1);

        n_req = min(max(n_modes + 4, 8), size(A_int,1) - 1);

        % shift-invert：位移取在谱顶之上，见函数头注释
        sigma = k0^2 * (n_max + 0.1)^2;

        % 注意：不要设 eo.issym。使用数值 sigma（shift-invert）时 eigs 会
        % 忽略该字段并发出 MATLAB:eigs:IgnoredOptionIssym 警告——
        % 警告本身无害，但会掩盖真正需要关注的收敛警告。这里显式不设。
        eo = struct('tol', opts.tol, 'maxit', 1000);

        [V, D] = eigs(A_int, n_req, sigma, eo);
        lam = diag(D);

        [~, o] = sort(real(lam), 'descend');
        lam = lam(o);
        V   = V(:, o);
    end
end

% ══════════════════════════════════════════════════════════════════════
function Lext = est_extent(k0, n_est, n_clad)
% 按 30 个衰减长度估计需要的延伸量，下限 3 μm
g = real(k0 * sqrt(complex(n_est^2 - n_clad^2)));
if g > 1e-9
    Lext = max(3, 30/g);
else
    Lext = 3;
end
end
