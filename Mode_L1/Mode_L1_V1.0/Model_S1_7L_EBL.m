function modes = Model_S1_7L_EBL(action)
%Model_S1_7L_EBL  7 层 GaN 基蓝光 LD 平板波导 —— 单文件求解与绘图
%
%   用法（把本文件放到任意目录，cd 过去即可，【不需要任何其他文件】）：
%
%       Model_S1_7L_EBL              % 求解 + 画图（图留在屏幕上，可缩放、可取点）
%       Model_S1_7L_EBL('save')      % 求解 + 导出 PNG 到 out/ebl_modes.png
%       Model_S1_7L_EBL('9')         % 换用 9 层 MQW 结构
%       Model_S1_7L_EBL('9save')     % 9 层 + 导出
%       m = Model_S1_7L_EBL(...)     % 同时把求解结果返回，接着在命令行算别的
%
%   ── 为什么是"单文件自包含" ──────────────────────────────────
%   全部依赖都写在本文件【下方】的局部函数里，不引用 src/ 下的任何东西。
%   把一个 .m 拷到任何地方就能跑，从参数定义一直到出图。
%
%   ── 命名规则 ────────────────────────────────────────────────
%       Model_S<阶段号>_<后缀>     后缀 = 层数 + 特征
%   本文件 Model_S1_7L_EBL = 阶段 1 · 7 层 · 含 EBL。
%   换结构时请【另存为新文件并改后缀】（例如 Model_S1_9L_MQW），
%   不要在本文件上改数字——保留它作为对照组。
%
%   ⚠ 文件名 / 函数名必须是 ASCII 标识符。中文名（如"一条龙"）虽然
%     exist('一条龙.m') 会返回 2（文件确实找得到），但调用时会报
%     "文本字符无效...请检查是否为 ASCII 字符"。文件【内容】里的中文
%     注释与字符串不受影响，可以正常显示。
%
%   ── 文件结构（从上往下读）──────────────────────────────────
%       §1  参数定义    ← 要改结构、波长、偏振，只改这一段
%       §2  求解        色散函数求根 → n_eff
%       §3  后处理      由 n_eff 重建模场 → 归一化 → 算 Γ
%       §4  结果输出    把关键数字打到命令行
%       §5  绘图        折射率 + 模场，两张子图
%       §6  局部函数    上面各步用到的全部函数（原来散在 src/ 里）
%
%   ── 物理问题 ────────────────────────────────────────────────
%   光沿 z 传播，沿外延方向 y 被折射率台阶限制。要解的是
%       d²E/dy² + k₀²·[n²(y) − n_eff²]·E = 0
%   的束缚模。每层内部 E 是正弦/指数函数，界面处 E 与 dE/dy 连续
%   （TE；TM 换成 H_y 与 (1/n²)dH_y/dy）。把这些条件凑成一个只含 n_eff
%   的方程（色散方程），它的零点就是本征值。
%
%   ── 单位与符号约定（全项目统一，不得另立门户）────────────────
%       长度 μm ｜ k₀ = 2π/λ [μm⁻¹] ｜ β = k₀·n_eff [μm⁻¹]
%       时间因子 exp(−iωt)，沿 +z 取 exp(+iβz)
%       ⇒ α = 2·k₀·Im(n_eff)  [cm⁻¹] 乘 1e4
%       ⇒ Im(n_eff) > 0 表示【损耗】，< 0 表示【净增益】
%       Γ ≡ ∫_有源区 w dy / ∫_{-∞}^{+∞} w dy
%           TE: w = |E_y|²      TM: w = |H_y|²/n²
%       分子取【全部 is_active 层的并集，不含势垒】；分母取无穷积分。

% ══════════════════════════════════════════════════════════════════
%  §1  参数定义
% ══════════════════════════════════════════════════════════════════

if nargin < 1 || isempty(action), action = 'show'; end
do_save  = contains(lower(string(action)), 'save');
use9     = contains(lower(string(action)), '9');

lambda  = 0.45;        % 真空波长 [μm] = 450 nm
pol     = 'TE';        % 'TE' | 'TM'
n_modes = 3;           % 最多求几阶模

% ── 层结构 ────────────────────────────────────────────────────
%    默认 7 层：衬底 / n包层 / n波导 / QW / p波导 / EBL / p包层
%    换成 9 层（阱-垒-阱 MQW）：带 '9' 参数运行，或把下面两行换成
%    layers = stack_9layer_mqw();        两个函数的定义都在 §6 末尾。
if use9
    layers = stack_9layer_mqw();
else
    layers = stack_7layer();
end

% 求解选项
opts = struct( ...
    'n_scan',   20000, ...   % 色散函数扫描点数（根搜索用）
    'tol',      1e-13, ...   % 求根收敛容差（相对）
    'maxit',    200,   ...   % 最大迭代次数
    'merge',    1e-6,  ...   % 去重时两个根的最小相对间距
    'npts',     800,   ...   % 模场重建时每层采样点数
    'tail_dec', 20,    ...   % 半无限层向外延伸多少个衰减长度
    'tail_max', 50);         % 延伸长度上限 [μm]

% ══════════════════════════════════════════════════════════════════
%  §2  求解：色散函数求根 → n_eff
% ══════════════════════════════════════════════════════════════════

L  = build_layers(layers);      % 校验层结构，折射率强制转 complex
k0 = 2*pi/lambda;
N  = numel(L);
yb = layer_bounds(L);           % 各层的 [y_lo, y_hi]，L(2) 下沿定为 y = 0

% 束缚模的搜索区间：n_eff 必须【高于两侧包层中折射率较高的那个】，
% 且【低于结构内最大的折射率】。低于前者就不是束缚模（会泄漏出去）。
n_clad_max = max(real(L(1).n), real(L(N).n));
n_max      = max(real([L.n]));
tiny       = 1e-9 * max(1, n_max);
n_lo       = n_clad_max + tiny;
n_hi       = n_max      - tiny;

assert(n_hi > n_lo, ...
    '搜索区间为空（n_lo = %g ≥ n_hi = %g）：结构里必须存在折射率高于两侧包层的有限层', ...
    n_lo, n_hi);

% 多找 2 个根再截断：多层结构的色散函数不保证单调，偶尔会在真模附近
% 出现数值伪根，留余量避免伪根把真模挤掉。
fh = @(x) dispersion_tmm(x, L, lambda, pol);
neff_list = find_roots_complex(fh, n_lo, n_hi, n_modes + 2, opts);
if numel(neff_list) > n_modes
    neff_list = neff_list(1:n_modes);
end
neff_list = neff_list(:).';

% 按 n_eff 降序（n_eff 越大 ⇒ 模序越低 ⇒ modes(1) 是基模）
[~, ord]  = sort(real(neff_list), 'descend');
neff_list = neff_list(ord);

% ══════════════════════════════════════════════════════════════════
%  §3  后处理：逐模重建模场 → 归一化 → 算 Γ
% ══════════════════════════════════════════════════════════════════
%  ⚠ 顺序不能换：先重建场，再归一化，最后算 Γ。
%    且 normalize_mode 与 confinement_factor 必须传同一个 pol——
%    TE 与 TM 的功率权重不同（TE 用 |E_y|²，TM 用 |H_y|²/n²）。

modes = struct([]);
for j = 1:numel(neff_list)
    ne = neff_list(j);

    [y, E, n_y]   = mode_profile(ne, L, lambda, pol, opts);
    [E, norm_val] = normalize_mode(y, E, n_y, pol);
    [G, det]      = confinement_factor(y, E, n_y, yb, L, pol);

    modes(j).n_eff    = ne;              % 有效折射率（本征值）
    modes(j).beta     = k0 * ne;         % 传播常数 [μm^-1]
    modes(j).order    = j - 1;           % 模序，0 = 基模
    modes(j).y        = y;               % 横坐标 [μm]
    modes(j).E        = E;               % 归一化后的模场 φ(y)
    modes(j).n_y      = n_y;             % 对应位置的折射率
    modes(j).Gamma    = G;               % 限制因子
    modes(j).norm_val = norm_val;        % 归一化积分值，应恒为 1
    modes(j).P_layer  = det.P_layer;     % 各层功率占比，总和 = 1
    modes(j).yb       = yb;
end

% ══════════════════════════════════════════════════════════════════
%  §4  结果输出
% ══════════════════════════════════════════════════════════════════

fprintf('\n');
fprintf('  ── 层结构（%d 层，λ = %g nm，%s 模）──\n', N, lambda*1000, pol);
fprintf('   #  tag          d[um]     n        P_layer(TE0)   P_layer(TE1)\n');
fprintf('  --------------------------------------------------------------------\n');
for k = 1:N
    if isinf(L(k).d), dtxt = 'Inf'; else, dtxt = sprintf('%.3f', L(k).d); end
    p0 = modes(1).P_layer(k)*100;
    if numel(modes) >= 2, p1 = modes(2).P_layer(k)*100; else, p1 = NaN; end
    fprintf('   %d  %-11s  %-8s %.3f    %9.4f%%    %9.4f%%\n', ...
            k, L(k).tag, dtxt, real(L(k).n), p0, p1);
end

fprintf('\n  ── 各阶模 ──\n');
for k = 1:numel(modes)
    fprintf('   TE%d   n_eff = %.8f   Γ = %.6f   (∫w dy = %.15f)\n', ...
            k-1, real(modes(k).n_eff), modes(k).Gamma, modes(k).norm_val);
end

% 自检：Γ 必须等于各 is_active 层 P_layer 之和。
% ⚠ 这条检查不是形式。Γ 曾经写成 trapz(y(mask), w(mask))，其中 mask 是
%   各 is_active 层的并集；当这些层【不连续】时（MQW 的阱之间隔着垒），
%   y(mask) 把几段拼成一个连续数组，trapz 于是在上一段末点与下一段首点
%   之间拉了一条直线，把中间的势垒整段积了进去——不报错、不警告，
%   实测让双层阱结构（3/6/3 nm）的 Γ 从 0.01704 被抬到 0.03409（约 2 倍）。
%   现在改为逐层单独积分再求和，两者严格相等。
act = find([L.is_active]);
if ~isempty(act)
    dG = max(abs(arrayfun(@(k) sum(k.P_layer(act)) - k.Gamma, modes)));
    fprintf('\n  ── 自检 ──\n');
    fprintf('   Γ 与各 is_active 层 P_layer 之和的最大偏差 = %.3e  (应 = 0)\n', dG);
    fprintf('   各层 P_layer 之和的最大偏差 = %.3e  (应 = 0)\n', ...
            max(abs(arrayfun(@(k) sum(k.P_layer) - 1, modes))));
end
fprintf('\n');

% ══════════════════════════════════════════════════════════════════
%  §5  绘图
% ══════════════════════════════════════════════════════════════════
%  两张子图，各看各的东西：
%    上图——整个外延结构：看模式被限制在哪个区域、两侧各穿透多深。
%    下图——芯区放大：看 QW 的折射率尖峰与 EBL 的凹陷。
%  两张图横轴差一个量级，这正是处理薄层时必须单独放大看的原因——
%  20 nm 的 EBL 在上图里不到一个像素宽。

if do_save, vis = 'off'; else, vis = 'on'; end
fig = figure('Visible', vis, 'Position', [60 40 1040 760], ...
             'Name', 'GaN LD slab waveguide', 'Color', 'w');

ax1 = subplot(2,1,1);
draw_stack(ax1, modes(1).y, real(modes(1).n_y), modes, yb, L, ...
           [-0.6 2.6], 'northwest');
title(ax1, sprintf('Whole epitaxial stack   (\\lambda = %g nm, %s modes)', ...
                   lambda*1000, pol));

ax2 = subplot(2,1,2);
draw_stack(ax2, modes(1).y, real(modes(1).n_y), modes, yb, L, ...
           [1.10 1.62], 'northeast');
title(ax2, 'Zoom on the core:  n-wave / QW / p-wave / EBL');

if do_save
    outdir = fullfile(fileparts(mfilename('fullpath')), 'out');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    exportgraphics(fig, fullfile(outdir, 'ebl_modes.png'), 'Resolution', 160);
    close(fig);
    fprintf('[saved] out/ebl_modes.png\n\n');
end
end

% ══════════════════════════════════════════════════════════════════
%  §6  局部函数
% ══════════════════════════════════════════════════════════════════
%  以下就是原来散在 src/ 里的全部代码。写成同一文件的局部函数后，
%  本文件不再依赖任何外部 .m。

% ──────────────────────────────────────────────────────────────────
function layers = stack_7layer()
%STACK_7LAYER  7 层结构：衬底 / n包层 / n波导 / QW / p波导 / EBL / p包层
%
%   ⚠ 四个字段数组 d / n / tag / is_active 【必须等长】，第 i 个元素
%     描述同一层。加层就是在四个数组的同一位置各插一个元素。
%     插错位置不会报错，只会安静地解出一个物理上不存在的结构。
%
%   ⚠ 首尾两层厚度必须是 Inf（半无限）。截断最外层等于凭空插进两个界面，
%     会解出物理上不存在的假离散模。build_layers 对这条直接报错。
%
%   折射率取值锚在 GaN ≈ 2.43 @450 nm，以及"GaN 与 Al0.06GaN 折射率差
%   ≈ 0.026"这条文献对比上，按 Al 组分线性外推（斜率 ≈ 0.433/单位组分）：
%       Al0.06GaN → 2.404 ｜ Al0.10GaN → 2.374 ｜ Al0.20GaN → 2.343
%   属于【量级可靠、末位不可靠】的典型值，不适用于任何定量对标。
%
%   ⚠ EBL 的 2.343 比上下两层的 2.455 / 2.374 都低，形成折射率凹陷——
%     这是整个结构里唯一的"反向"台阶。它把光往 n 侧推（光学作用），
%     同时阻挡电子（电学作用，本阶段不涉及）。

layers = struct( ...
    'd', ...
        { Inf,      1.200,      0.120,      0.009,      0.080,      0.020,      Inf      }, ...
    'n', ...
        { 1.770,    2.404,      2.455,      2.550,      2.455,      2.343,      2.374    }, ...
    'tag', ...
        {'substrate','n-clad',  'n-wave',   'QW',       'p-wave',   'EBL',      'p-clad' }, ...
    'is_active', ...
        { false,    false,      false,      true,       false,      false,      false    });
end

% ──────────────────────────────────────────────────────────────────
function layers = stack_9layer_mqw()
%STACK_9LAYER_MQW  9 层结构：把单层等效 QW 展开成 阱-垒-阱
%
%   其余 6 层与 stack_7layer 完全相同。改动只有：原来第 4 层那个 9 nm 的
%   QW，换成了 3 nm 阱 / 6 nm 垒 / 3 nm 阱 三层，层数 7 → 9。
%
%   ⚠ is_active 只标在两个阱上，【不标垒】。这正是本项目 Γ 的口径：
%        Γ = (阱1 + 阱2 的功率) / 总功率
%     若把垒也标上，Γ 会变大但含义变了（变成"MQW 区平均"，两者不可混用）。

layers = struct( ...
    'd', ...
        { Inf,   1.200,   0.120,   0.003,   0.006,   0.003,   0.080,   0.020,   Inf    }, ...
    'n', ...
        { 1.770, 2.404,   2.455,   2.550,   2.430,   2.550,   2.455,   2.343,   2.374  }, ...
    'tag', ...
        {'substrate','n-clad','n-wave','QW1','barrier','QW2','p-wave','EBL','p-clad'}, ...
    'is_active', ...
        { false, false,   false,   true,    false,   true,    false,   false,   false  });
end

% ──────────────────────────────────────────────────────────────────
function L = build_layers(layers)
%BUILD_LAYERS  校验并规范化层结构
%
%   输出 L 与输入同构，但 n 已强制转为 complex，缺省字段已补齐。
%   这一步是硬约束的集中处：任何不合法输入在这里报错，而不是留到求解器
%   里变成 NaN 或假模。

assert(isstruct(layers) && ~isempty(layers), 'build_layers: layers 必须是 struct 数组');
N = numel(layers);
assert(N >= 3, ...
    'build_layers: 至少需要 3 层（两侧半无限包层 + 至少一个有限层），当前 %d 层', N);

L = layers;
for k = 1:N
    assert(isfield(L(k), 'd') && ~isempty(L(k).d), 'build_layers: 第 %d 层缺少字段 d', k);
    assert(isfield(L(k), 'n') && ~isempty(L(k).n), 'build_layers: 第 %d 层缺少字段 n', k);

    if ~isfield(L(k), 'tag') || isempty(L(k).tag)
        L(k).tag = sprintf('layer%d', k);
    end
    if ~isfield(L(k), 'is_active') || isempty(L(k).is_active)
        L(k).is_active = false;
    end

    L(k).n = complex(L(k).n);        % 复折射率接口：从第一天起就是复数

    assert(real(L(k).n) > 0, ...
        'build_layers: 第 %d 层 (%s) 折射率实部必须为正，当前 %g', ...
        k, L(k).tag, real(L(k).n));
    assert(imag(L(k).n) >= 0, ...
        ['build_layers: 第 %d 层 (%s) 折射率虚部为负 (%g)。\n' ...
         '  按本项目约定 Im(n) > 0 表示【损耗】；若该层是增益介质，\n' ...
         '  请在阶段4 的增益模块中处理，不要在此处给负虚部。'], ...
        k, L(k).tag, imag(L(k).n));
end

assert(isinf(L(1).d) && L(1).d > 0, 'build_layers: 第 1 层（衬底侧）厚度必须为 +Inf');
assert(isinf(L(N).d) && L(N).d > 0, 'build_layers: 第 %d 层（覆盖层侧）厚度必须为 +Inf', N);

if N > 2
    dmid = [L(2:N-1).d];
    assert(all(dmid > 0) && all(isfinite(dmid)), 'build_layers: 中间层厚度必须为有限正值');
end

act = [L.is_active];
assert(~act(1) && ~act(N), 'build_layers: 半无限层不能标记为有源区');
end

% ──────────────────────────────────────────────────────────────────
function yb = layer_bounds(L)
%LAYER_BOUNDS  由层结构算出每层的 y 区间（N×2，第 k 行是 [y_lo, y_hi]）
%
%   约定：把 L(2) 的下沿放在 y = 0，y 随层号递增。首尾半无限层分别记为
%   [-Inf, 0] 与 [d_total, +Inf]。这个原点只影响绘图横轴，不影响物理量。

N = numel(L);
yb = zeros(N, 2);
yb(1, :) = [-Inf, 0];

y = 0;
for k = 2:N-1
    yb(k, :) = [y, y + L(k).d];
    y = y + L(k).d;
end

yb(N, :) = [y, Inf];
end

% ──────────────────────────────────────────────────────────────────
function [a11, a12, a21, a22] = layer_matrix(kz, d, n, pol)
%LAYER_MATRIX  单层的 2×2 传输矩阵（TE / TM 通用）
%
%   把状态向量从层下沿传到上沿： v_out = [a11 a12; a21 a22] · v_in
%       TE : v = [E_y ; dE_y/dy]           界面连续：E_y 与 dE_y/dy
%       TM : v = [H_y ; (1/n²)·dH_y/dy]    界面连续：H_y 与 (1/n²)dH_y/dy
%
%   ⚠ kz → 0 时矩阵元含 sin(kz·d)/kz 是 0/0，此处取极限值 d。
%     触发场景：试探的 n_eff 恰好落入某一层的折射率。层数越多、折射率
%     台阶越密，撞上的概率越大。不处理会产生 NaN，而 NaN 在符号比较中
%     既不大于零也不小于零，会让根搜索【静默漏掉整个模式】。

assert(isscalar(d) && d > 0, 'layer_matrix: 层厚 d 必须为正标量');
is_TE = strcmpi(pol, 'TE');

kd = kz .* d;
S  = sin(kd) ./ kz;
C  = cos(kd);

bad = (abs(kz) < 1e-14);
if any(bad(:))
    S(bad) = d;                      % sin(kz·d)/kz → d
end

if is_TE
    a11 = C;      a12 = S;
    a21 = -kz .* sin(kd);           a22 = C;
else
    a11 = C;      a12 = (n^2) .* S;
    a21 = -(kz / (n^2)) .* sin(kd); a22 = C;
end
end

% ──────────────────────────────────────────────────────────────────
function [f, aux] = dispersion_tmm(n_eff, L, lambda, pol)
%DISPERSION_TMM  N 层平板波导的色散函数 f(n_eff)，零点即本征模
%
%   任意层数成立（首尾半无限）。n_eff 可为标量或数组，逐元素计算。
%
%   ── 色散方程怎么来的 ────────────────────────────────────────
%   把中间所有有限层乘成一个总矩阵 M，作用在底部半无限层的边界状态上：
%       底部（离开波导向 -y 衰减）: v_bot ∝ [1 ; s_sub]
%       顶部（离开波导向 +y 衰减）: v_top(2)/v_top(1) = -s_cov
%   其中对 TE 取 s = γ，对 TM 取 s = γ/n²。两式联立即得
%       f(n_eff) = (M·[1;s_sub])(2) + s_cov·(M·[1;s_sub])(1) = 0
%
%   这个构造【没有经过任何近似或离散化】，只用到"Maxwell 方程 + 界面连续"。

if nargin < 4 || isempty(pol), pol = 'TE'; end
is_TE = strcmpi(pol, 'TE');
N = numel(L);

in_size = size(n_eff);
x = n_eff(:).';                       % 1×M 行向量
k0 = 2*pi/lambda;

% 两侧半无限层的衰减常数（必须取衰减解，Re γ > 0）
gam_sub = k0 * sqrt(x.^2 - L(1).n^2);
gam_cov = k0 * sqrt(x.^2 - L(N).n^2);

if is_TE
    s_sub = gam_sub;              s_cov = gam_cov;
else
    s_sub = gam_sub / L(1).n^2;   s_cov = gam_cov / L(N).n^2;
end

% 中间层传输矩阵连乘（从下往上，逐次左乘）
M11 = ones(size(x));  M12 = zeros(size(x));
M21 = zeros(size(x)); M22 = ones(size(x));

for k = 2:N-1
    nk = L(k).n;
    dk = L(k).d;
    kz = k0 * sqrt(nk^2 - x.^2);
    [a11, a12, a21, a22] = layer_matrix(kz, dk, nk, pol);

    nM11 = a11.*M11 + a12.*M21;
    nM12 = a11.*M12 + a12.*M22;
    nM21 = a21.*M11 + a22.*M21;
    nM22 = a21.*M12 + a22.*M22;
    M11 = nM11; M12 = nM12; M21 = nM21; M22 = nM22;
end

v1 = M11 + M12.*s_sub;
v2 = M21 + M22.*s_sub;
f  = v2 + s_cov.*v1;
f  = reshape(f, in_size);

if nargout > 1
    aux = struct('gam_sub', gam_sub, 'gam_cov', gam_cov, 'M', [M11; M12; M21; M22]);
end
end

% ──────────────────────────────────────────────────────────────────
function [roots_n, info] = find_roots_complex(fun, n_lo, n_hi, n_roots, opts)
%FIND_ROOTS_COMPLEX  在一维区间内搜索零点（实根与复根统一接口）
%
%   阶段1（折射率全实）色散函数在实轴上取实值，用"扫描符号变化 + 二分"，
%   收敛有保证、不会漏根。阶段2 引入复折射率后自动切到"|f| 局部极小 +
%   复数牛顿"分支。返回的 roots_n 按 n_eff 【降序】排列。

if nargin < 4 || isempty(n_roots), n_roots = 5; end
if nargin < 5 || isempty(opts), opts = struct(); end
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

is_real = all(abs(imag(fs)) <= 1e-10 * fscale);

if is_real
    roots_n = scan_bisect(fun, xs, real(fs), opts);
    info.mode = 'real';
else
    roots_n = scan_newton(fun, xs, fs, opts);
    info.mode = 'complex';
end

% 去重 + 降序 + 截断
% ⚠ 必须保留【最大的】n_roots 个根。波导问题里 n_eff 越大 ⇒ 模序越低
%   ⇒ 物理上越重要（基模的 n_eff 最大）。按升序截断留下的是"最靠近截止的
%   一批高阶模"，基模被直接挤掉——不报错、不警告，只给出错的模式集。
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

% ──────────────────────────────────────────────────────────────────
function roots_n = scan_bisect(fun, xs, fr, opts)
%SCAN_BISECT  实值分支：扫描符号变化区间，再二分细化
%   二分只依赖"闭区间两端符号相反"这个充分条件，不要求函数单调或光滑。
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

% ──────────────────────────────────────────────────────────────────
function roots_n = scan_newton(fun, xs, fs, opts)
%SCAN_NEWTON  复值分支：定位 |f| 的局部极小，再用复数牛顿细化
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

% ──────────────────────────────────────────────────────────────────
function r = newton_complex(fun, x0, opts)
%NEWTON_COMPLEX  复数牛顿迭代，带步长阻尼。无法收敛时返回空。
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

% 后验检查：|f| 必须显著小于邻域，剔除鞍点/极小值平台造成的伪根
fr  = abs(fun(r));
eps_ = 1e-3 * max(1, abs(r));
nb  = max(abs(fun(r + eps_)), abs(fun(r - eps_)));
if ~(fr < 0.5 * nb)
    r = [];
end
end

% ──────────────────────────────────────────────────────────────────
function [y, E, n_y, yb] = mode_profile(n_eff, L, lambda, pol, opts)
%MODE_PROFILE  由 n_eff 重建模场剖面 φ(y)
%
%   从底部半无限层的边界状态出发，逐层用与 dispersion_tmm 【完全相同】的
%   层矩阵向上传输，在每个采样点求值。所以重建出的 φ(y) 一定满足界面连续
%   条件——这是把它当作"真解"而不是"拟合曲线"的前提。
%
%   输出 E 未归一化；TE 为 E_y，TM 为 H_y。

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

% ── 第 1 层：底部半无限层 ─────────────────────────────────────
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

% ── 第 2 … N-1 层：有限层，逐层传输 ──────────────────────────
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
        a12 = S;                 a21 = -kz * SN;
    else
        a12 = (nk^2) * S;        a21 = -(kz / (nk^2)) * SN;
    end
    a11 = C;  a22 = C;

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

% ── 第 N 层：顶部半无限层 ────────────────────────────────────
yy = linspace(d_total, d_total + depth_cov, opts.npts);
yy = yy(2:end);                                % 去掉 y=d_total，避免重复
y_cell{N} = yy;
E_cell{N} = exp(-gam_cov * (yy - d_total)) * v(1);
n_cell{N} = L(N).n * ones(size(yy));

% ── 拼接 ────────────────────────────────────────────────────
y   = [y_cell{:}];
E   = [E_cell{:}];
n_y = [n_cell{:}];

% 折射率全实且 n_eff 为实时，清掉运算中产生的 0 虚部
if imag(n_eff) == 0 && all(imag([L.n]) == 0)
    E   = real(E);
    n_y = real(n_y);
end

assert(numel(y) == numel(E) && numel(y) == numel(n_y), ...
    'mode_profile: 拼接后长度不一致');
end

% ──────────────────────────────────────────────────────────────────
function dep = tail_depth(gam, opts)
%TAIL_DEPTH  半无限层向外延伸的深度：默认 20 个衰减长度（场幅 e^-20 ≈ 2e-9）
%   延伸不足会让 Γ 系统性偏小（分子分母都用这个积分）。
%   若某侧 γ 很小（模式接近截止），延伸长度会自动放大到上限 tail_max。
g = real(gam);
if g > 1e-9
    dep = min(opts.tail_max, opts.tail_dec / g);
    dep = max(dep, 0.1);
else
    dep = opts.tail_max;
end
end

% ──────────────────────────────────────────────────────────────────
function [E, norm_val] = normalize_mode(y, E, n_y, pol)
%NORMALIZE_MODE  归一化模场，使 ∫ w(y) dy = 1
%
%   w(y) 正比于时间平均坡印廷矢量的 z 分量（功率密度）。按本项目约定
%   E(y,z,t) = φ(y)·exp(i(βz − ωt)) 可得
%       TE:  S_z ∝ β·|E_y|²            （与 n 无关）
%       TM:  S_z ∝ (β/n²)·|H_y|²
%   所以 TM 的权重【必须】带 1/n²。两者混用会让 Γ 与后续的功率、增益
%   计算系统性偏离。
%
%   返回的 norm_val 是归一化【后】复算的积分，应恒为 1。

if strcmpi(pol, 'TE')
    w = abs(E).^2;
else
    w = abs(E).^2 ./ max(real(n_y).^2, eps);
end

norm_before = trapz(y, w);
assert(isfinite(norm_before) && norm_before > 0, ...
    'normalize_mode: 归一化积分非正或非有限 (%g)，模场重建有问题', norm_before);

E = E / sqrt(norm_before);

if strcmpi(pol, 'TE')
    w2 = abs(E).^2;
else
    w2 = abs(E).^2 ./ max(real(n_y).^2, eps);
end
norm_val = trapz(y, w2);
end

% ──────────────────────────────────────────────────────────────────
function [G, detail] = confinement_factor(y, E, n_y, yb, L, pol)
%CONFINEMENT_FACTOR  限制因子 Γ（本项目统一口径，全阶段复用）
%
%        Γ ≡  ∫_{有源区}  w(y) dy  /  ∫_{-∞}^{+∞} w(y) dy
%
%        TE 模 :  w(y) = |E_y(y)|²
%        TM 模 :  w(y) = |H_y(y)|² / n(y)²
%
%   物理含义：【有源区内的净功率占总功率的比例】。无量纲，0 ≤ Γ ≤ 1。
%
%   ── 为什么文献里的 Γ 彼此不可比（两处独立分歧，光钉死一处没用）──
%     (1) 分子范围：单个阱 / 全部阱层的并集 / MQW 区（含势垒）。
%         阱垒厚度相当时，这三者能差 2–3 倍。
%     (2) 分母范围：无穷积分，还是截断到外延结构。
%         对漏模只有无穷积分是自洽的。
%   本项目取：(1) = 全部 is_active 层的并集；(2) = 无穷积分。
%
%   配套口径：Γ(全部阱) × g(单阱材料增益)。若改用 MQW 区平均材料增益，
%   必须另算一个 Γ_MQW，两者【不可混用】。
%
%   ── detail.P_layer ──────────────────────────────────────────
%   各层功率占比（长度 = 层数，总和 = 1）。它不依赖 is_active，任何结构
%   都能算，用于非对称性判据、损耗归属分析与衬底泄漏量化。

if strcmpi(pol, 'TE')
    w = abs(E).^2;
else
    w = abs(E).^2 ./ max(real(n_y).^2, eps);
end

P_tot = trapz(y, w);
assert(P_tot > 0, 'confinement_factor: 总功率积分非正');

% ── 各层功率占比（与 is_active 无关）─────────────────────────
N = numel(L);
P_layer = zeros(1, N);
for k = 1:N
    if isinf(yb(k,1))
        mk = (y <= yb(k,2));
    elseif isinf(yb(k,2))
        mk = (y >= yb(k,1));
    else
        mk = (y >= yb(k,1) & y <= yb(k,2));
    end
    if any(mk)
        P_layer(k) = trapz(y(mk), w(mk)) / P_tot;
    end
end

act = [L.is_active];
if ~any(act)
    % 没有有源区 ⇒ Γ 在数学上无定义。返回 NaN 而不是 0：
    % 0 会被下游误读成"模式与有源区完全不重叠"这个【物理结论】，
    % 而 NaN 会显式传播，逼使用者面对"结构没标有源区"这个【输入问题】。
    G = NaN;
    detail = struct('P_act', NaN, 'P_tot', P_tot, 'P_layer', P_layer);
    return;
end

% ── 有源区功率：必须【逐层单独积分再求和】─────────────────────
% ⚠ 不能写成 P_act = trapz(y(mask), w(mask))，其中 mask 是各 is_active
%   层的并集。当这些层【不连续】时（多层 MQW，阱与阱之间隔着势垒），
%   y(mask) 会把几段拼接成一个连续数组，trapz 于是在上一段的末点与下一段
%   的首点之间拉一条直线，把中间的势垒整段积了进去。
%   这个错误不报错、不警告，只让 Γ 偏大：实测双层阱结构（3/6/3 nm）
%   Γ 由应有的 0.01704 被抬到 0.03409，约 2 倍。
%   逐层积分天然按层分段，不会跨越任何间隙，且与上面 P_layer 的口径一致。
P_act = 0;
for k = 1:N
    if ~L(k).is_active, continue; end
    if isinf(yb(k,1))
        mk = (y <= yb(k,2));
    elseif isinf(yb(k,2))
        mk = (y >= yb(k,1));
    else
        mk = (y >= yb(k,1) & y <= yb(k,2));
    end
    if any(mk)
        P_act = P_act + trapz(y(mk), w(mk));
    end
end
G = P_act / P_tot;

detail = struct('P_act', P_act, 'P_tot', P_tot, 'P_layer', P_layer);
end

% ──────────────────────────────────────────────────────────────────
function draw_stack(ax, y, n_y, m, yb, L, xr, legloc)
%DRAW_STACK  在一个坐标轴上同时画折射率分布与各阶模场（双纵轴）
%
%   ⚠ 纵轴范围一律【由当前窗口内的数据算出】。任何写死的纵轴上限都会把
%     峰值切掉，而切掉峰值在视觉上极容易被误读成"模场形状不对"：
%     |φ|² 归一化后峰值 = 1/(有效宽度)，对本结构约 3.4，不是 1。
%
%   ⚠ 层底色必须画在【左轴】，与折射率线同轴。yyaxis 的两根轴是分层
%     渲染的（右轴整体盖在左轴之上，与绘制先后无关）；若把底色画在右轴，
%     它会把折射率线整条盖成灰色。

axes(ax); %#ok<LAXES>
cla(ax);
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

% ── 只统计绘图窗口内的数据，范围才不会被窗口外的尾场拉偏 ──────
sel  = (y >= xr(1)) & (y <= xr(2));
nsel = n_y(sel);

nmin = min(nsel);
nmax = max(nsel);
pad_n = 0.08*(nmax - nmin) + 0.01;

pk = 0;
for k = 1:numel(m)
    w = abs(m(k).E).^2;
    pk = max(pk, max(w(sel)));
end
ytop = pk * 1.10;                    % 右轴上限：留 10% 余量给图例和峰

% ── 左轴：层底色 + 折射率阶梯 ────────────────────────────────
yyaxis(ax, 'left');
ylim(ax, [nmin - pad_n, nmax + pad_n]);
yl = ylim(ax);

for k = 1:numel(L)
    lo = max(yb(k,1), xr(1));
    hi = min(yb(k,2), xr(2));
    if hi <= lo, continue; end
    if L(k).is_active
        patch(ax, [lo hi hi lo], [yl(1) yl(1) yl(2) yl(2)], [1 0.85 0.35], ...
              'FaceAlpha', 0.50, 'EdgeColor','none', 'HandleVisibility','off');
    else
        patch(ax, [lo hi hi lo], [yl(1) yl(1) yl(2) yl(2)], [0.92 0.92 0.92], ...
              'FaceAlpha', 0.50, 'EdgeColor','none', 'HandleVisibility','off');
    end
end
plot(ax, y, n_y, 'k-', 'LineWidth', 1.6, 'HandleVisibility','off');
ylabel(ax, 'refractive index n');

% ── 右轴：各阶模 |φ|² + 层边界标注 ───────────────────────────
yyaxis(ax, 'right');
ylim(ax, [0 ytop]);

cols = lines(numel(m));
for k = 1:numel(m)
    plot(ax, y, abs(m(k).E).^2, '-', 'LineWidth', 1.7, 'Color', cols(k,:), ...
         'DisplayName', sprintf('TE%d  (n_{eff}=%.5f, \\Gamma=%.4f)', ...
                                k-1, real(m(k).n_eff), m(k).Gamma));
end
ylabel(ax, '|\phi(y)|^2    (normalized: \int w dy = 1)');

for k = 2:numel(L)
    xb = yb(k,1);
    if xb > xr(1) && xb < xr(2)
        plot(ax, [xb xb], [0 ytop], ':', 'Color',[0.55 0.55 0.55], ...
             'HandleVisibility','off');
    end
end
for k = 1:numel(L)
    lo = max(yb(k,1), xr(1));
    hi = min(yb(k,2), xr(2));
    if hi <= lo, continue; end
    text(ax, 0.5*(lo+hi), 0.015*ytop, L(k).tag, 'Rotation', 90, 'FontSize', 7.5, ...
         'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
         'Color', [0.25 0.25 0.25]);
end

% ── 收尾：xlim 必须在所有 yyaxis 切换之后设，否则会被重置 ────
xlim(ax, xr);
legend(ax, 'Location', legloc, 'FontSize', 8);
hold(ax, 'off');
end
