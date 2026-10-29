% ─────────────────────────────────────────────────────
% 启动器：run_all_stage3（阶段 3 的唯一入口，只做转发与退出码）
% 输入：无
% 输出：Mode_L3/out/ 下的日志/七图/七csv；退出码非零 = 有 FAIL
% 调用函数：Model_Self_L3_7L_EBL（无参数 → 阶段 3 全流程）
% 说明：
%   1.按施工单 §4，入口名与 matlab -batch 的命令行不变；
%   2.全部物理、判据、出图都在 Mode_L3/Model_Self_L3_7L_EBL.m 的 e3_* 段里，
%     这里只负责把它调起来、按 FAIL 数决定退出码；
%   3.日志由主文件用 UTF-8 的 fopen 写（Mode_L3/out/stage3_tests.log），不用 diary。
% ─────────────────────────────────────────────────────

addpath(fileparts(fileparts(mfilename('fullpath'))));
S = Model_Self_L3_7L_EBL();
if isstruct(S) && isfield(S, 'fail') && S.fail > 0
    error('阶段 3 有 %d 项未通过，见 Mode_L3/out/stage3_tests.log', S.fail);
end
