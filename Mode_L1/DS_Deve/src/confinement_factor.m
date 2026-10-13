function [G, detail] = confinement_factor(y, E, n_y, yb, L, pol)
%CONFINEMENT_FACTOR  限制因子 Γ（本项目统一口径，阶段1–6 全程复用）
%
%   [G, detail] = confinement_factor(y, E, n_y, yb, L, pol)
%
%   ── 定义（第一行写死，不得在别处另立门户）──────────────────
%
%         Γ ≡  ∫_{有源区}  w(y) dy  /  ∫_{-∞}^{+∞} w(y) dy
%
%   其中功率密度权重
%         TE 模 :  w(y) = |E_y(y)|²
%         TM 模 :  w(y) = |H_y(y)|² / n(y)²
%
%   物理含义：【有源区内的净功率占总功率的比例】。无量纲，0 ≤ Γ ≤ 1。
%
%   ── 这个口径配套什么，不配套什么 ────────────────────────────
%   配套的是"用阱的材料增益 g 写模式增益"：  g_mode = Γ · g。
%   若某处改用"MQW 区平均材料增益"，必须另算一个 Γ_MQW（分子取整个 MQW 区
%   含势垒），两者【不可混用】——混用会让模式增益偏离 2–3 倍。
%
%   ── 为什么文献里的 Γ 数值彼此不可比 ─────────────────────────
%   有两处独立的分歧，光钉死一处没用：
%     (1) 分子范围：单个阱 / 全部阱层的并集 / MQW 区（含势垒）。
%         阱垒厚度相当时，这三者能差 2–3 倍。
%     (2) 分母范围：无穷积分，还是截断到外延结构。
%         对漏模，只有无穷积分是自洽的；截断会给出与截断位置相关的伪 Γ。
%   本项目取：(1) = 全部 is_active 层的并集；(2) = 无穷积分（数值上由
%   mode_profile 的 tail_dec 保证截断误差 < 1e-17）。
%
%   ── 输出的 detail.P_layer ───────────────────────────────────
%   各层的功率占比（长度 = 层数，顺序与 L 一致，总和为 1）。
%   它不依赖 is_active，任何结构都能算，用于：
%     · A6 的非对称性判据（光偏向哪一侧）
%     · 阶段2 的损耗归属分析（哪一层贡献了多少吸收）
%     · 阶段3 的衬底泄漏量化
%   ⇒ 这是本阶段为后续阶段预留的最有用的一个输出，不要删。

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
        m = (y <= yb(k,2));
    elseif isinf(yb(k,2))
        m = (y >= yb(k,1));
    else
        m = (y >= yb(k,1) & y <= yb(k,2));
    end
    if any(m)
        P_layer(k) = trapz(y(m), w(m)) / P_tot;
    end
end

act = [L.is_active];
if ~any(act)
    % 没有有源区 ⇒ Γ 在数学上无定义。此处返回 NaN 而不是 0：
    % 0 会被下游误读成"模式与有源区完全不重叠"这个【物理结论】，
    % 而 NaN 在任何后续运算中都会显式传播，逼使用者面对"结构没标有源区"
    % 这个【输入问题】。两者性质完全不同，不能用一个数糊过去。
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
