function [E, norm_val] = normalize_mode(y, E, n_y, pol)
%NORMALIZE_MODE  归一化模场，使 ∫ w(y) dy = 1
%
%   [E, norm_val] = normalize_mode(y, E, n_y, pol)
%
%   norm_val 是归一化【后】的积分值，应恒为 1，用于自检。
%
%   ── 归一化的定义（A5 的判据，必须在第一行写死）───────────────
%       TE 模 : w(y) = |E_y(y)|²
%       TM 模 : w(y) = |H_y(y)|² / n(y)²
%       归一化目标：  ∫_{-∞}^{+∞} w(y) dy = 1
%
%   ── 为什么 TE 与 TM 的权重不同 ──────────────────────────────
%   w(y) 正比于时间平均坡印廷矢量的 z 分量 S_z(y)，即【功率密度】。
%   按本项目约定 E(y,z,t) = φ(y)·exp(i(βz - ωt))，可以算出
%       TE:  S_z ∝ β·|E_y|²            （与 n 无关）
%       TM:  S_z ∝ (β/n²)·|H_y|²
%   所以 TM 的权重里必须带 1/n²。把两者混用会让 Γ 与后续的功率、增益
%   计算系统性偏离——这是"Γ 文献值不可比"之外，另一个容易踩的坑。
%
%   ⚠ 归一化必须在【模场重建之后、Γ 计算之前】做，且 Γ 与归一化必须用
%     同一个 pol。solve_slab 内部已按这个顺序调用。

if strcmpi(pol, 'TE')
    w = abs(E).^2;
else
    w = abs(E).^2 ./ max(real(n_y).^2, eps);
end

norm_before = trapz(y, w);

assert(isfinite(norm_before) && norm_before > 0, ...
    'normalize_mode: 归一化积分非正或非有限 (%g)，模场重建有问题', norm_before);

E = E / sqrt(norm_before);

% 复算归一化【后】的积分作为返回的自检值——它应恒为 1。
% 返回这个值而不是归一化前的值：调用方拿到的数字要么是 1（正常），
% 要么暴露归一化没做对，不存在第三种解释，比"某个不知量纲的大数"有用。
if strcmpi(pol, 'TE')
    w2 = abs(E).^2;
else
    w2 = abs(E).^2 ./ max(real(n_y).^2, eps);
end
norm_val = trapz(y, w2);
end
