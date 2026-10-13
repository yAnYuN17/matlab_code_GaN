function [f, aux] = dispersion_analytic3(n_eff, L, lambda, pol, m)
%DISPERSION_ANALYTIC3  三层非对称平板波导的解析色散函数（路径 1）
%
%   f = dispersion_analytic3(n_eff, L, lambda, 'TE', m)
%
%   仅对【严格三层】结构有效，即 size(L,2) == 3：
%       L(1) 下包层（半无限）
%       L(2) 芯层
%       L(3) 上包层（半无限）
%
%   f(n_eff) = 0 即本征方程的解。n_eff 可以是标量或任意形状数组，逐元素计算。
%
%   ── TE 模色散方程 ──────────────────────────────────────────
%       κ·d = arctan(γ_sub/κ) + arctan(γ_cov/κ) + mπ
%   其中
%       κ     = k0·√(n_core² − n_eff²)      芯层横向波数
%       γ_sub = k0·√(n_eff² − n_sub²)      下包层衰减常数
%       γ_cov = k0·√(n_eff² − n_cov²)      上包层衰减常数
%
%   ── TM 模色散方程 ──────────────────────────────────────────
%       κ·d = arctan((n_core/n_sub)²·γ_sub/κ)
%           + arctan((n_core/n_cov)²·γ_cov/κ) + mπ
%
%   ── 为什么这个函数能当"基准" ────────────────────────────────
%   它不经过任何矩阵构造或离散化，方程形式可用 Symbolic Math Toolbox 从
%   Maxwell 方程 + 界面连续条件重新推导核对（见 tests/t00_symbolic_check.m）。
%   这是本项目唯一一处"公式本身可被独立验证"的地方。
%
%   ── 单调性（决定了求根策略）────────────────────────────────
%   在束缚模区间 n_eff ∈ (n_high, n_core) 内，κ 随 n_eff 增大而减小、
%   γ 随 n_eff 增大而增大，故 arctan(γ/κ) 单调增、κ·d 单调减，
%   合起来 f 关于 n_eff 【严格单调递减】。
%   ⇒ 给定模序 m，区间内至多一个根，可直接用二分法，不会漏根、不会跳根。
%   （这条性质只对三层成立；多层结构的 f 不保证单调，故 TMM 路径必须扫描。）
%
%   ── 截止条件 ──────────────────────────────────────────────
%   TE_m 模存在所需的最小归一化频率：
%       V_c,m = arctan(√a) + mπ
%   其中  V = k0·d·√(n_core² − n_high²)， n_high = max(n_sub, n_cov)，
%         a = (n_high² − n_low²)/(n_core² − n_high²) 为不对称度，n_low = min(n_sub, n_cov)。
%   注意 a = 0（对称）时 V_c,0 = 0，即【对称波导基模永不截止】；
%   而 a > 0 时 V_c,0 > 0，【非对称波导连基模都有截止】——这是 A4 要验的点。
%
%   输入
%       n_eff   : 试探有效折射率（标量或数组）
%       L       : build_layers 规范化后的三层结构
%       lambda  : 真空波长 [μm]
%       pol     : 'TE' | 'TM'
%       m       : 模序（非负整数）
%
%   输出
%       f       : 色散函数值，形状与 n_eff 相同
%       aux     : struct，含 κ / γ_sub / γ_cov / V / a，供诊断用

if nargin < 4 || isempty(pol), pol = 'TE'; end
if nargin < 5 || isempty(m),   m   = 0;    end

assert(numel(L) == 3, ...
    'dispersion_analytic3: 只支持严格三层结构，当前 %d 层。多层请用 dispersion_tmm。', ...
    numel(L));

in_size = size(n_eff);
x = n_eff(:).';                       % 统一成 1×M 行向量，便于向量化

k0      = 2*pi/lambda;
n_sub   = L(1).n;
n_core  = L(2).n;
n_cov   = L(3).n;
d       = L(2).d;

kappa   = k0 * sqrt(n_core^2 - x.^2);     % 芯层横向波数
gam_sub = k0 * sqrt(x.^2 - n_sub^2);      % 下包层衰减常数
gam_cov = k0 * sqrt(x.^2 - n_cov^2);      % 上包层衰减常数

if strcmpi(pol, 'TE')
    f = kappa*d - atan(gam_sub./kappa) - atan(gam_cov./kappa) - m*pi;
else
    f = kappa*d ...
        - atan((n_core/n_sub)^2 * gam_sub./kappa) ...
        - atan((n_core/n_cov)^2 * gam_cov./kappa) - m*pi;
end

f = reshape(f, in_size);

if nargout > 1
    n_high = max(real(n_sub), real(n_cov));
    n_low  = min(real(n_sub), real(n_cov));
    aux = struct( ...
        'kappa',   reshape(kappa,   in_size), ...
        'gam_sub', reshape(gam_sub, in_size), ...
        'gam_cov', reshape(gam_cov, in_size), ...
        'V',       k0*d*sqrt(complex(n_core^2 - n_high^2)), ...
        'a',       (n_high^2 - n_low^2)/real(n_core^2 - n_high^2), ...
        'n_high',  n_high, ...
        'n_low',   n_low);
end
end
