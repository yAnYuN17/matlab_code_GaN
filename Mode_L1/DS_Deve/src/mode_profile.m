function [y, E, n_y, yb] = mode_profile(n_eff, L, lambda, pol, opts)
%MODE_PROFILE  由 n_eff 重建模场剖面 φ(y)
%
%   [y, E, n_y, yb] = mode_profile(n_eff, L, lambda, pol, opts)
%
%   输入
%       n_eff   : 有效折射率（由任一路径求出，标量）
%       L       : build_layers 规范化后的层结构
%       lambda  : 真空波长 [μm]
%       pol     : 'TE' | 'TM'
%       opts    : .npts       每层采样点数，默认 800
%                 .tail_dec   半无限层向外延伸多少个衰减长度，默认 20
%                 .tail_max   延伸长度上限 [μm]，默认 50
%
%   输出
%       y    : 1×M 横坐标 [μm]，单调递增，L(2) 下沿为 y = 0
%       E    : 1×M 模场，TE 为 E_y，TM 为 H_y（未归一化）
%       n_y  : 1×M 对应位置的折射率（用于 Γ 与 TM 归一化的权重）
%       yb   : N×2 各层的 [y_lo, y_hi]
%
%   ── 重建方法 ────────────────────────────────────────────────
%   从底部半无限层的边界状态出发，逐层用与 dispersion_tmm 【完全相同】的
%   层矩阵向上传输，在每个采样点求值。所以重建出的 φ(y) 一定满足
%   界面连续条件——这是把它当作"真解"而不是"拟合曲线"的前提。
%
%   ── 半无限层要延伸多深 ──────────────────────────────────────
%   默认向外延伸 20 个衰减长度：场幅衰减到 e^-20 ≈ 2e-9，对归一化积分
%   （∝ 场幅平方）的贡献约 4e-18，完全可忽略。Γ 的分子分母都用到这个积分，
%   截断不足会让 Γ 系统性偏小。若某侧 γ 很小（模式接近截止），延伸长度会
%   自动放大——此时务必检查 tail_max 有没有被撞到，撞到就说明该模式
%   已接近"几乎不束缚"，Γ 的数值不可信。

if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'npts'),     opts.npts     = 800; end
if ~isfield(opts, 'tail_dec'), opts.tail_dec = 20;  end
if ~isfield(opts, 'tail_max'), opts.tail_max = 50;  end

k0 = 2*pi/lambda;
N  = numel(L);
yb = layer_bounds(L);
d_total = yb(N, 1);

gam_sub = k0 * sqrt(complex(n_eff^2 - L(1).n^2));
gam_cov = k0 * sqrt(complex(n_eff^2 - L(N).n^2));

depth_sub = tail_depth(gam_sub, opts);
depth_cov = tail_depth(gam_cov, opts);

is_TE = strcmpi(pol, 'TE');

y_cell = cell(1, N);
E_cell = cell(1, N);
n_cell = cell(1, N);

% ═══ 第 1 层：底部半无限层 ═══════════════════════════════════
yy = linspace(-depth_sub, 0, opts.npts);
y_cell{1} = yy(1:end-1);                       % 去掉 y=0，避免与下一层重复
E_cell{1} = exp(gam_sub * y_cell{1});
n_cell{1} = L(1).n * ones(size(y_cell{1}));

% 在 y = 0 处的状态向量（与 dispersion_tmm 的底部边界条件一致）
if is_TE
    v = [1; gam_sub];
else
    v = [1; gam_sub / L(1).n^2];
end
cur_y = 0;

% ═══ 第 2 … N-1 层：有限层，逐层传输 ════════════════════════
for k = 2:N-1
    dk = L(k).d;
    nk = L(k).n;
    kz = k0 * sqrt(complex(nk^2 - n_eff^2));

    yy = linspace(cur_y, cur_y + dk, opts.npts);
    dy = yy - cur_y;

    if abs(kz) < 1e-14
        S = dy;
    else
        S = sin(kz*dy) / kz;
    end
    C  = cos(kz*dy);
    SN = sin(kz*dy);

    if is_TE
        a12 = S;
        a21 = -kz * SN;
    else
        a12 = (nk^2) * S;
        a21 = -(kz / (nk^2)) * SN;
    end
    a11 = C;
    a22 = C;

    vv = a11 * v(1) + a12 * v(2);

    % 该层上沿的新状态（供下一层使用）
    v = [a11(end)*v(1) + a12(end)*v(2); ...
         a21(end)*v(1) + a22(end)*v(2)];

    if k < N-1
        y_cell{k} = yy(1:end-1);
        E_cell{k} = vv(1:end-1);
    else
        y_cell{k} = yy;
        E_cell{k} = vv;
    end
    n_cell{k} = nk * ones(size(y_cell{k}));
    cur_y = cur_y + dk;
end

% ═══ 第 N 层：顶部半无限层 ═══════════════════════════════════
yy = linspace(d_total, d_total + depth_cov, opts.npts);
yy = yy(2:end);                                % 去掉 y=d_total，避免重复
y_cell{N} = yy;
E_cell{N} = exp(-gam_cov * (yy - d_total)) * v(1);
n_cell{N} = L(N).n * ones(size(yy));

% ═══ 拼接 ════════════════════════════════════════════════════
y   = [y_cell{:}];
E   = [E_cell{:}];
n_y = [n_cell{:}];

% 当折射率全为实数且 n_eff 为实数时，清掉运算过程中产生的 0 虚部
if imag(n_eff) == 0 && all(imag([L.n]) == 0)
    E   = real(E);
    n_y = real(n_y);
end

assert(numel(y) == numel(E) && numel(y) == numel(n_y), ...
    'mode_profile: 拼接后长度不一致');
end

% ══════════════════════════════════════════════════════════════════════
function dep = tail_depth(gam, opts)
% 半无限层向外延伸的深度：20 个衰减长度，但上下都有夹逼
g = real(gam);
if g > 1e-9
    dep = min(opts.tail_max, opts.tail_dec / g);
    dep = max(dep, 0.1);
else
    % γ 近零 ⇒ 该侧几乎不衰减（模式接近截止），用上限并靠上层断言兜底
    dep = opts.tail_max;
end
end
