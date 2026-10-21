% ─────────────────────────────────────────────────────
% 测试：test_s5_perturbation（阶段2 S5 验收入口：B1/B2/B3 微扰标尺）
% 输入：无
% 输出：命令行打印，末行 S5 TEST pass=1 为通过
% 调用函数：Model_Self_L2_7L_EBL('s5')
% 说明：
%   1.B1 弱吸收极限：逐层单独赋 α ∈ {1,3,10} cm⁻¹，与 P_layer 比较；
%   2.B2 根后验：每个根的 |f|/邻域最小值 < 0.5；
%   3.B3 多通道叠加：6 层同时赋 α（各 ≤10），与 Σ P_j α_j 及严格式比较。
% ─────────────────────────────────────────────────────
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
out = Model_Self_L2_7L_EBL('s5');
fprintf('S5 TEST pass=%d  b1_rel_naive_max=%.3e b1_rel_strict_max=%.3e b3_rel=%.3e\n', ...
    out.pass, out.b1_rel_naive_max, out.b1_rel_strict_max, out.b3_rel_sum);
