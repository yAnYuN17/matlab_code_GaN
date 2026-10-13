function [f, aux] = dispersion_tmm(n_eff, L, lambda, pol)
%DISPERSION_TMM  N 层平板波导的传输矩阵色散函数（路径 2）
%
%   f = dispersion_tmm(n_eff, L, lambda, 'TE')
%
%   对任意层数成立（首尾为半无限）。f(n_eff) = 0 的零点即本征模。
%   n_eff 可为标量或数组，逐元素计算（内部完全向量化，扫描上万点无压力）。
%
%   ── 状态向量与层矩阵 ────────────────────────────────────────
%   TE 模取 v = [E_y ; dE_y/dy]，界面连续条件为 E_y 与 dE_y/dy 连续。
%   TM 模取 v = [H_y ; (1/n²)·dH_y/dy]，界面连续条件为 H_y 与 (1/n²)dH_y/dy 连续。
%   每层的传输矩阵（把状态向量从层下沿传到上沿，层厚 d、横向波数 kz）：
%       TE:  M = [ cos(kz·d)        sin(kz·d)/kz     ]
%                [ -kz·sin(kz·d)    cos(kz·d)        ]
%       TM:  M = [ cos(kz·d)        (n²/kz)·sin(kz·d)]
%                [ -(kz/n²)·sin(kz·d)  cos(kz·d)     ]
%
%   ── 色散方程怎么来的 ──────────────────────────────────────
%   把中间所有有限层乘成一个总矩阵，作用在底部半无限层的边界状态上：
%       底部（离开波导向 -y 衰减）: v_bot ∝ [1 ; s_sub]
%       顶部（离开波导向 +y 衰减）: v_top(2)/v_top(1) = -s_cov
%   其中对 TE 取 s = γ，对 TM 取 s = γ/n²，γ 为该半无限层的衰减常数。
%   两式联立即得
%       f(n_eff) = (M·[1;s_sub])(2) + s_cov·(M·[1;s_sub])(1) = 0
%
%   ⚠ 这个 f 的构造【没有经过任何近似或离散化】。它与解析色散方程共享的
%     只有"Maxwell 方程 + 界面连续"这一条物理前提，但数学路径完全不同：
%     路径1 解超越方程，路径2 求 2×2 矩阵连乘的行列式式零点。
%     这正是 A0 三路交叉验证里"两条路径相互独立"的含义。
%
%   ── 数值陷阱（阶段1 已处理，阶段2 需重新审视）───────────────
%   (1) kz → 0 时 sin(kz·d)/kz 是 0/0。本函数在 |kz| 极小时改用极限值 d
%       （即 sin(kz·d)/kz → d），避免产生 NaN 假根。
%       触发场景：n_eff 逼近某一层的折射率。层数越多、折射率台阶越密，
%       这种情况越频繁——交接报告里"根会变密"讲的就是它。
%   (2) 半无限层的 γ 必须是衰减解（Re γ > 0）。MATLAB 的 sqrt 对复数返回
%       实部非负的分支，恰好满足；但阶段2 一旦引入复折射率与 PML，
%       必须显式检查 Re γ 的符号，否则会解出物理上不存在的"增长模"。
%   (3) 本实现是直接矩阵连乘（TMM）。当层内有吸收、且层厚 × 虚部很大时，
%       矩阵元会指数增长，导致精度丢失。阶段2 算衬底泄漏损耗时应切换到
%       散射矩阵（S 矩阵）法。切换点已经隔离好：本函数是唯一需要改的地方，
%       它对外的签名和零点语义保持不变。
%
%   输入
%       n_eff   : 试探有效折射率（标量或数组）
%       L       : build_layers 规范化后的层结构
%       lambda  : 真空波长 [μm]
%       pol     : 'TE' | 'TM'
%
%   输出
%       f       : 色散函数值，形状与 n_eff 相同
%       aux     : struct，含 gamma_sub / gamma_cov / M（标量输入时）

if nargin < 4 || isempty(pol), pol = 'TE'; end

is_TE = strcmpi(pol, 'TE');
N = numel(L);

in_size = size(n_eff);
x = n_eff(:).';                       % 1×M 行向量

k0 = 2*pi/lambda;

% ── 两侧半无限层的衰减常数 ──────────────────────────────────
gam_sub = k0 * sqrt(x.^2 - L(1).n^2);
gam_cov = k0 * sqrt(x.^2 - L(N).n^2);

if is_TE
    s_sub = gam_sub;
    s_cov = gam_cov;
else
    s_sub = gam_sub / L(1).n^2;
    s_cov = gam_cov / L(N).n^2;
end

% ── 中间层传输矩阵连乘（从下往上，逐次左乘）─────────────────
M11 = ones(size(x));
M12 = zeros(size(x));
M21 = zeros(size(x));
M22 = ones(size(x));

for k = 2:N-1
    nk = L(k).n;
    dk = L(k).d;

    kz = k0 * sqrt(nk^2 - x.^2);
    [a11, a12, a21, a22] = layer_matrix(kz, dk, nk, pol);

    % M_acc ← M_layer · M_acc
    nM11 = a11.*M11 + a12.*M21;
    nM12 = a11.*M12 + a12.*M22;
    nM21 = a21.*M11 + a22.*M21;
    nM22 = a21.*M12 + a22.*M22;
    M11 = nM11; M12 = nM12; M21 = nM21; M22 = nM22;
end

% ── 色散函数 ────────────────────────────────────────────────
% v_bot = [1; s_sub]；要求 v_top(2) + s_cov·v_top(1) = 0
v1 = M11 + M12.*s_sub;
v2 = M21 + M22.*s_sub;
f  = v2 + s_cov.*v1;

f = reshape(f, in_size);

if nargout > 1
    aux = struct('gam_sub', gam_sub, 'gam_cov', gam_cov, ...
                 'M', [M11; M12; M21; M22]);
end
end
