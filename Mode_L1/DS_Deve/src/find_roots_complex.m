function [roots_n, info] = find_roots_complex(fun, n_lo, n_hi, n_roots, opts)
%FIND_ROOTS_COMPLEX  在一维区间内搜索色散函数的零点（实根与复根统一接口）
%
%   [roots_n, info] = find_roots_complex(fun, n_lo, n_hi, n_roots, opts)
%
%   输入
%       fun     : 函数句柄，fun(x) 接受标量或数组，返回同尺寸的色散函数值
%       n_lo    : 搜索下界（通常是较高的包层折射率，束缚模要求 n_eff > n_lo）
%       n_hi    : 搜索上界（通常是结构内最大折射率，n_eff < n_hi）
%       n_roots : 最多返回多少个根（默认 5）
%       opts    : .n_scan   扫描点数，默认 20000
%                 .tol      收敛容差（相对），默认 1e-13
%                 .maxit    最大迭代次数，默认 200
%                 .merge    根的合并距离（相对），默认 1e-6
%
%   输出
%       roots_n : 【按 n_eff 降序】排列的前 n_roots 个根
%                 （阶段1 为实数，阶段2 起为复数）
%       info    : struct，含 .mode ('real'|'complex') / .n_found / .scanned
%
%   ── 两条分支，自动选择 ──────────────────────────────────────
%   阶段1（全部折射率为实）：色散函数在实轴上取实值。此时用"扫描定位
%   符号变化 + 二分"，收敛有保证、不会漏根。二分的可靠性来自闭区间
%   符号相反的充分条件，不依赖函数的单调性或光滑性。
%
%   阶段2 起（引入复折射率）：色散函数取复值，"符号变化"失去意义。
%   此时自动切换到"|f| 局部极小定位 + 复数牛顿迭代"。
%
%   ⚠ 复数牛顿是有可能收敛到伪根的（|f| 的局部极小但不是零点）。阶段2
%     启用后必须加一道后验检查：把根代回 |f| 与相邻点比较，确认它是真零点
%     而不是鞍点或极小值平台。阶段1 用不到，但接口先留在这里。
%
%   ── 为什么扫描要 20000 点 ───────────────────────────────────
%   n_eff 逼近某一层折射率时，色散函数的零点会变密（多层结构尤其明显）。
%   按典型区间宽度 0.2 估算，20000 点给出 1e-5 的采样间隔，足以分辨
%   绝大多数密集根。若某次扫描发现可疑的窄区间，就把 opts.n_scan 加大
%   重扫，而不是改用更松的容差。

if nargin < 4 || isempty(n_roots), n_roots = 5;    end
if nargin < 5 || isempty(opts),    opts    = struct(); end
if ~isfield(opts, 'n_scan'), opts.n_scan = 20000; end
if ~isfield(opts, 'tol'),    opts.tol    = 1e-13; end
if ~isfield(opts, 'maxit'),  opts.maxit  = 200;   end
if ~isfield(opts, 'merge'),  opts.merge  = 1e-6;  end

assert(n_hi > n_lo, 'find_roots_complex: 要求 n_hi > n_lo');

xs = linspace(n_lo, n_hi, opts.n_scan);
fs = fun(xs);
fs = fs(:).';

fscale = max(abs(fs));
if fscale == 0, fscale = 1; end

% 判定函数是否在实轴上取实值
is_real = all(abs(imag(fs)) <= 1e-10 * fscale);

if is_real
    roots_n = scan_bisect(fun, xs, real(fs), opts);
    info.mode = 'real';
else
    roots_n = scan_newton(fun, xs, fs, opts);
    info.mode = 'complex';
end

% 去重 + 降序 + 截断
%
% ⚠ 必须保留【最大的】n_roots 个根，而不是最小的。
%   波导问题里 n_eff 越大 ⇒ 模序越低 ⇒ 物理上越重要（基模的 n_eff 最大）。
%   若按升序截断，得到的会是"最靠近截止的一批高阶模"，基模被直接挤掉。
%   这个错误不会报错、不会警告，只会安静地返回错的模式集——
%   本项目开发中已经踩过一次，务必保留这条注释。
if ~isempty(roots_n)
    roots_n = sort(roots_n(:).', 'descend');
    keep = [true, abs(diff(roots_n)) > opts.merge * max(1, abs(roots_n(1:end-1)))];
    roots_n = roots_n(keep);
    if numel(roots_n) > n_roots
        roots_n = roots_n(1:n_roots);
    end
end

info.n_found = numel(roots_n);
info.scanned = opts.n_scan;
end

% ══════════════════════════════════════════════════════════════════════
function roots_n = scan_bisect(fun, xs, fr, opts)
% 实值分支：扫描符号变化区间，再二分细化
roots_n = [];
s = sign(fr);
s(s == 0) = 1;                       % 恰好为零的采样点，并入正侧
idx = find(s(1:end-1) .* s(2:end) < 0);

for j = 1:numel(idx)
    a  = xs(idx(j));
    b  = xs(idx(j)+1);
    fa = fun(a);

    for it = 1:opts.maxit
        c = 0.5*(a + b);
        if (b - a) <= opts.tol * max(1, abs(c))
            break;
        end
        fc = fun(c);
        if real(fc) == 0
            a = c; b = c;
            break;
        end
        if sign(real(fc)) == sign(real(fa))
            a = c; fa = fc;
        else
            b = c;
        end
    end

    roots_n(end+1) = 0.5*(a + b);   %#ok<AGROW>
end
end

% ══════════════════════════════════════════════════════════════════════
function roots_n = scan_newton(fun, xs, fs, opts)
% 复值分支：定位 |f| 的局部极小，再用复数牛顿细化
mag = abs(fs);
cand = [];
for j = 2:numel(xs)-1
    if mag(j) < mag(j-1) && mag(j) <= mag(j+1)
        cand(end+1) = xs(j);        %#ok<AGROW>
    end
end

roots_n = [];
for x0 = cand
    r = newton_complex(fun, x0, opts);
    if ~isempty(r)
        roots_n(end+1) = r;         %#ok<AGROW>
    end
end
end

% ══════════════════════════════════════════════════════════════════════
function r = newton_complex(fun, x0, opts)
% 复数牛顿迭代，带步长阻尼。无法收敛时返回空。
h      = 1e-7;                       % 数值导数步长（对 n_eff ~ O(1) 合适）
r      = x0;
n_iter = opts.maxit;

for it = 1:n_iter
    fx = fun(r);
    f1 = fun(r + h);
    f2 = fun(r - h);
    df = (f1 - f2) / (2*h);

    if abs(df) < 1e-300
        r = [];
        return;
    end

    step = fx / df;
    if abs(step) > 0.5                 % 阻尼：单步不超过 0.5
        step = 0.5 * step / abs(step);
    end
    r = r - step;

    if abs(step) < 1e-14 * max(1, abs(r))
        break;
    end
    if it == n_iter
        r = [];                        % 未收敛，丢弃
        return;
    end
end

% 后验检查：|f| 必须显著小于邻域
fr  = abs(fun(r));
eps_ = 1e-3 * max(1, abs(r));
nb  = max(abs(fun(r + eps_)), abs(fun(r - eps_)));
if ~(fr < 0.5 * nb)
    r = [];
end
end
