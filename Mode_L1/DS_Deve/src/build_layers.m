function L = build_layers(layers)
%BUILD_LAYERS  校验并规范化层结构，供 solve_slab 使用
%
%   L = build_layers(layers)
%
%   输入 layers : struct 数组（1×N），字段
%       .d          厚度 [μm]，【首尾两层必须为 Inf】表示半无限
%       .n          复折射率 n' + 1i·n''（阶段1 取 n'' = 0）
%       .tag        层名字符串，供 Γ 与后续阶段引用（可省略）
%       .is_active  逻辑值，标记该层是否属于有源区，即 Γ 分子的积分范围
%                   （可省略，默认 false）
%
%   输出 L : 规范化后的 struct 数组，n 已强制为 complex 类型
%
%   ── 层序方向约定 ───────────────────────────────────────────
%       L(1) 是最底层（衬底侧），L(end) 是最顶层（覆盖层侧）。
%       y 坐标随层号递增：L(1) → L(N) 对应 y 从小到大。
%       这一点会贯穿模场剖面、Γ 积分与所有绘图，不得中途反转。
%
%   ── 为什么首尾必须半无限 ───────────────────────────────────
%       把最外层截断成有限厚度，等于人为插入一个额外的介质界面，
%       从而解出物理上不存在的"假离散模"。真实器件的包层外侧在光学
%       意义上确实是半无限介质（衬底、空气或金属），必须如实建模。
%       这条是阶段1 的硬约束：本函数会直接报错而不是静默容忍。
%
%   ── 为后续阶段预留的接口 ───────────────────────────────────
%       (1) 复折射率：.n 从第一天起就是 complex。阶段2 只需给 n'' 赋值，
%           本函数与求解器一行都不用改。
%       (2) 开放边界：首尾层允许"半无限 + 复折射率"。阶段2 要算衬底泄漏
%           损耗时，在 L(1) 或 L(end) 外侧再插一层 PML / 复坐标变换层即可，
%           层结构的表示能力已经支持。

assert(isstruct(layers) && ~isempty(layers), ...
    'build_layers: layers 必须是 struct 数组');
N = numel(layers);
assert(N >= 3, ...
    'build_layers: 至少需要 3 层（两侧半无限包层 + 至少一个有限层），当前 %d 层', N);

L = layers;

for k = 1:N
    assert(isfield(L(k), 'd') && ~isempty(L(k).d), ...
        'build_layers: 第 %d 层缺少字段 d', k);
    assert(isfield(L(k), 'n') && ~isempty(L(k).n), ...
        'build_layers: 第 %d 层缺少字段 n', k);

    if ~isfield(L(k), 'tag') || isempty(L(k).tag)
        L(k).tag = sprintf('layer%d', k);
    end
    if ~isfield(L(k), 'is_active') || isempty(L(k).is_active)
        L(k).is_active = false;
    end

    % 复折射率接口：从第一天起就是复数
    L(k).n = complex(L(k).n);

    assert(real(L(k).n) > 0, ...
        'build_layers: 第 %d 层 (%s) 折射率实部必须为正，当前 %g', ...
        k, L(k).tag, real(L(k).n));
    assert(imag(L(k).n) >= 0, ...
        ['build_layers: 第 %d 层 (%s) 折射率虚部为负 (%g)。\n' ...
         '  按本项目约定 Im(n) > 0 表示【损耗】；若该层是增益介质，\n' ...
         '  请在阶段4 的增益模块中处理，不要在此处给负虚部。'], ...
        k, L(k).tag, imag(L(k).n));
end

% 首尾半无限
assert(isinf(L(1).d) && L(1).d > 0, ...
    'build_layers: 第 1 层（最底层衬底侧）厚度必须为 +Inf');
assert(isinf(L(N).d) && L(N).d > 0, ...
    'build_layers: 第 %d 层（最顶层覆盖层侧）厚度必须为 +Inf', N);

% 中间层厚度
if N > 2
    dmid = [L(2:N-1).d];
    assert(all(dmid > 0) && all(isfinite(dmid)), ...
        'build_layers: 中间层厚度必须为有限正值');
end

% 有源区必须存在且位于内部（不能是半无限层）
act = [L.is_active];
assert(~act(1) && ~act(N), 'build_layers: 半无限层不能标记为有源区');
end
