function yb = layer_bounds(L)
%LAYER_BOUNDS  由层结构算出每层的 y 区间
%
%   yb = layer_bounds(L)   返回 N×2 数组，第 k 行为 [y_lo, y_hi]
%
%   约定：把 L(2) 的下沿放在 y = 0，y 随层号递增（与 build_layers 的
%   "L(1) 在底、L(end) 在顶"一致）。首尾半无限层分别记为 [-Inf, 0] 与 [d, +Inf]。
%
%   这个 y 原点是人为选定的（定在第一个有限层的下沿），只影响绘图横轴，
%   不影响任何物理量。

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
