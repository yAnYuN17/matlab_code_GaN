% ─────────────────────────────────────────────────────
% 测试：test_s6_tmm_stability（阶段2 S6 验收入口：TMM 失真起点）
% 输入：无
% 输出：α 扫描表 + 失真起点，末行 S6 TEST 打印 onset
% 调用函数：Model_Self_L2_7L_EBL('s6')
% 说明：
%   1.三层对称结构、三层同赋吸收，α 从 0 扫到 1e5 cm⁻¹；
%   2.参考解为三层解析色散（复数牛顿独立求根），误差 = TMM − 解析；
%   3.判据：目标量级（α ≤ 100 cm⁻¹）内 dα/α < 1e-6 则可跳过 S 矩阵。
% ─────────────────────────────────────────────────────
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
out = Model_Self_L2_7L_EBL('s6');
fprintf('S6 TEST onset=%g  max_rel(alpha<=100)=%.3e  max_rel_all=%.3e\n', ...
    out.onset, out.max_rel_upto_100, out.max_rel_all);
