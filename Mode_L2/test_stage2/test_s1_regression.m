% ─────────────────────────────────────────────────────
% 脚本：test_s1_regression（S1 回归：阶段2 副本与阶段1 原件逐位一致）
% 输入：无
% 输出：命令行两行 DIFF + 末行 "S1 TEST pass=1/0"；四个差值全为 0 才通过
% 调用函数：Model_Self_L1_7L_EBL('show') / ('9')，Model_Self_L2_7L_EBL('show') / ('9')
% 说明：
%   1) 比对 7 层与 9 层两组，逐项看 n_eff、Γ、P_layer、远场 FWHM 与 th_pk。
%   2) 判据是全 0（逐位相同），不是"差很小"；α=0 时虚部精确为 0 才是这条判据成立的前提。
%   3) 每次改动求解器、层构建、根搜索之后都要重跑；阶段3 动 find_roots_complex 的
%      merge 阈值（遗留缺陷 L3）时，这是第一条要过的线。
% ─────────────────────────────────────────────────────
here = fileparts(mfilename('fullpath'));
base = fileparts(here);
addpath(fullfile(fileparts(base), 'Mode_L1', 'Mode_L1_V1.0'));
addpath(base);

a7 = Model_Self_L1_7L_EBL('show'); b7 = Model_Self_L2_7L_EBL('show');
a9 = Model_Self_L1_7L_EBL('9');    b9 = Model_Self_L2_7L_EBL('9');
d7 = cmp(a7, b7, '7L');
d9 = cmp(a9, b9, '9L');

worst = max([d7, d9]);
fprintf('S1 TEST pass=%d  (max|diff|=%.3g)\n', worst == 0, worst);

% ── 小工具 ──
function d = cmp(x, y, tag)
    dn = max(abs([x.n_eff] - [y.n_eff]));
    dG = max(abs([x.Gamma] - [y.Gamma]));
    dP = 0; dF = 0;
    for k = 1:numel(x)
        dP = max(dP, max(abs(x(k).P_layer - y(k).P_layer)));
        dF = max(dF, abs(x(k).ff.fwhm  - y(k).ff.fwhm));
        dF = max(dF, abs(x(k).ff.th_pk - y(k).ff.th_pk));
    end
    fprintf('%s DIFF dn=%.18g dG=%.18g dP=%.18g dF=%.18g\n', tag, dn, dG, dP, dF);
    d = max([dn, dG, dP, dF]);
end
