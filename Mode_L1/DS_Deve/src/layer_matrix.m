function [a11, a12, a21, a22] = layer_matrix(kz, d, n, pol)
%LAYER_MATRIX  单层的 2×2 传输矩阵（TE / TM 通用）
%
%   [a11,a12,a21,a22] = layer_matrix(kz, d, n, 'TE')
%
%   把状态向量从层下沿传到上沿：  v_out = [a11 a12; a21 a22] · v_in
%       TE : v = [E_y ; dE_y/dy]           界面连续：E_y 与 dE_y/dy
%       TM : v = [H_y ; (1/n²)·dH_y/dy]    界面连续：H_y 与 (1/n²)dH_y/dy
%
%   其中层厚 d（标量，单位 μm）、折射率 n（可为复数）、
%   横向波数 kz = k0·√(n² − n_eff²)（可为复数）。
%
%   kz 支持数组（逐元素），d 必须是标量——扫描时对同一层批量求值用。
%
%   ── kz → 0 的处理（必须做，否则会静默出错）──────────────────
%   矩阵元素含 sin(kz·d)/kz，当 kz = 0 时是 0/0。此处取极限值 d。
%   触发场景：试探的 n_eff 恰好落入某一层的折射率。层数越多、折射率台阶
%   越密，撞上的概率越大——这正是交接报告所记"根会变密"的数值根源。
%   若不处理会产生 NaN，而 NaN 在符号比较中既不大于零也不小于零，
%   会让根搜索静默漏掉整个模式的根。TE 与 TM 两个分支都要替换。

assert(isscalar(d) && d > 0, 'layer_matrix: 层厚 d 必须为正标量');

is_TE = strcmpi(pol, 'TE');

kd = kz .* d;
S  = sin(kd) ./ kz;
C  = cos(kd);

bad = (abs(kz) < 1e-14);
if any(bad(:))
    S(bad) = d;                       % sin(kz·d)/kz → d
end

if is_TE
    a11 = C;
    a12 = S;
    a21 = -kz .* sin(kd);
    a22 = C;
else
    a11 = C;
    a12 = (n^2) .* S;
    a21 = -(kz / (n^2)) .* sin(kd);
    a22 = C;
end
end
