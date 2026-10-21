% ─────────────────────────────────────────────────────
% 测试：test_s4_rootcheck（阶段2 S4 验收入口：复根后验检查）
% 输入：无
% 输出：命令行打印，末行 S4 TEST pass=1 为通过
% 调用函数：Model_Self_L2_7L_EBL('rootcheck')
% 说明：
%   1.模型内动作 'rootcheck' 跑三组用例：真实 7 层结构、构造真根、构造伪根；
%   2.判据：真实结构每根 ratio < 0.5，构造真根留 2 个，构造伪根留 0 个。
% ─────────────────────────────────────────────────────
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
out = Model_Self_L2_7L_EBL('rootcheck');
fprintf('S4 TEST pass=%d n_real=%d n_true=%d n_pseudo=%d ratioMax=%.3e\n', ...
    out.pass, out.n_real, out.n_true, out.n_pseudo, max(out.ratio_real));
