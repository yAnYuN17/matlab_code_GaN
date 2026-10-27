% ─────────────────────────────────────────────────────
% 脚本：test_s7_selftest（用原文 Table II 反查 Table I 系数抄录是否正确）
% 输入：无
% 输出：命令行逐行对照与末行 "S7SELF TEST pass=1/0"
% 调用函数：Model_Self_L2_7L_EBL('sigma', opt)
% 说明：
%   1) Table II 给的是 Γ_i、N_i、Γ_iα_i（520 nm 绿光 LD），由 Γ_iα_i/Γ_i 反推 α_i，
%      再与 N_i·σ(λ,pol,x) 比。判据直接比 Γα（原文只印到 0.1 cm^-1，半个末位即 0.05），
%      容差取 0.06；不用 α 的相对偏差——小 Γα 行会被舍入放大成假偏差。
%   2) QW 行用**双载流子**：电子与空穴各 2e19 cm^-3，σ(α) = σ_e + σ_h。
%      掺杂层只有一种载流子，QW 里两种都在，这是本表唯一一行要相加的。
%      2026-10-22 之前该行只取空穴，预测 49.14 vs 原文反推 90.0（差 1.83 倍）；
%      加上电子通道后 90.38，Γα 预测 2.711 vs 原文 2.7（|Δ|=0.011），故改为参与判据。
%   3) 行里的 x 是原文那一层的组分（p-guiding 0.10、QW 0.30），不是本结构的 0.04/0.15：
%      这一步验的是"系数抄对了没有"，不是"本结构算得对不对"。
% ─────────────────────────────────────────────────────
base = fileparts(fileparts(mfilename('fullpath')));
addpath(base);

% 列：层名 / Γ_i / N[cm^-3] / 载流子（'e' | 'h' | 'he' 表示两种各 N）/ x / Γα_原文
rows = { ...
    'p-cladding', 0.23, 3e19, 'h',  0.00, 9.5; ...
    'p-guiding',  0.16, 1e19, 'h',  0.10, 2.9; ...
    'EBL',        0.03, 2e19, 'h',  0.00, 0.8; ...
    'QW',         0.03, 2e19, 'he', 0.30, 2.7; ...
    'n-guiding',  0.19, 2e18, 'e',  0.10, 0.4; ...
    'n-cladding', 0.27, 2e18, 'e',  0.00, 0.2 };

fprintf('\n===== Table I 参数化 vs 原文 Table II（520 nm）=====\n');
fprintf('%-12s %6s %10s %9s %9s %9s %9s %9s\n', ...
    '层', 'Γ', 'N[cm^-3]', 'α_原文', 'α_预测', 'Γα_原文', 'Γα_预测', '|ΔΓα|');
dev = zeros(1, size(rows, 1));
for k = 1:size(rows, 1)
    tag = rows{k,1}; G = rows{k,2}; Nc = rows{k,3}; cr = rows{k,4};
    xl  = rows{k,5}; GA = rows{k,6};
    switch cr
        case 'he'
            ce = Model_Self_L2_7L_EBL('sigma', struct('lambda', 0.52, 'pol', 'TE', ...
                                                      'carrier', 'e', 'x', xl));
            ch = Model_Self_L2_7L_EBL('sigma', struct('lambda', 0.52, 'pol', 'TE', ...
                                                      'carrier', 'h', 'x', xl));
            sig = ce.total + ch.total;
        otherwise
            c   = Model_Self_L2_7L_EBL('sigma', struct('lambda', 0.52, 'pol', 'TE', ...
                                                      'carrier', cr, 'x', xl));
            sig = c.total;
    end
    a_pred = Nc * sig;
    dGA    = abs(G*a_pred - GA);
    dev(k) = dGA;
    fprintf('%-12s %6.2f %10.2e %9.2f %9.2f %9.2f %9.2f %9.3f\n', ...
        tag, G, Nc, GA/G, a_pred, GA, G*a_pred, dGA);
end
mx = max(dev);
fprintf('全部六行 max|ΔΓα| = %.3f cm^-1（判据 0.06 = 原文一个末位的舍入）\n', mx);
fprintf('\nS7SELF TEST pass=%d\n', mx <= 0.06);
