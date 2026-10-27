% ───────────────────────────────────────────────────────────────────────────────────
% 脚本：test_s5_perturbation（B1 弱吸收极限 + B3 多通道叠加自洽）
% 输入：无
% 输出：命令行两张表 + B1/B3 判定（末行 S5 TEST pass=…）
% 调用函数：Model_Self_L2_7L_EBL('solve', opt)
% 说明：
%   1.B 步后本脚本自己算判定，不再调模型里的 s5 入口；全部物理经 'solve' 取数，
%     每次调用内部重建层结构，逐层赋 α 不会互相污染（坑 H5）；
%   2.判据口径 A（2026-10-20 用户定）：与 P_j·n_rj/Re(n_eff) 比，相对偏差 <1e-2。
%     'solve' 返回的 contrib(j) 就是 P_j·α_j·n_rj/Re(n_eff)，除以 α_j 即该参考值；
%   3.字面 ΣP_jα_j 的偏差另列：它等于 |n_rj/n_eff − 1|，是物理不是数值误差；
%   4.线性度：α_j 从 1 到 10 时 α_mode/α_j 的变化 <1e-2（弱吸收极限的适用边界）。
% ───────────────────────────────────────────────────────────────────────────────────
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
addpath('C:/Users/L/Desktop/UG/GaN_LD/waveguide/Mode_L1/Mode_L1_V1.0');
addpath('C:/Users/L/Desktop/UG/GaN_LD/waveguide/Mode_L1/Mode_L1_V1.0/layers');

k0  = 2*pi/0.45;
S0  = Model_Self_L2_7L_EBL('solve', struct());       % 参考态：全层无吸收
if numel(S0) < 1 || ~S0(1).ok
    fprintf('S5 TEST pass=0  (参考态没解出来)\n');
    return
end
tags = S0(1).tag;
P    = S0(1).P_layer;
ne0  = real(S0(1).ne);
fprintf('参考态 Re(n_eff)=%.12f  Γ=%.9f  ΣP−1=%.2e\n', ne0, S0(1).Gamma, sum(P)-1);

fprintf('\n===== B1 逐层弱吸收（TE0，450 nm，α ∈ {1,3,10} cm^-1）=====\n');
fprintf('%-10s %9s %14s %13s %11s %11s\n', 'layer', 'P(%)', 'alpha_i/alpha', 'P*nr/neff', 'naive_rel', 'strict_rel');
alist   = [1 3 10];
r1      = zeros(1, numel(tags));
rs      = zeros(1, numel(tags));
lin     = zeros(1, numel(tags));
ratio1  = zeros(1, numel(tags));
ratio10 = zeros(1, numel(tags));
for j = 1:numel(tags)
    rel_naive = zeros(1, numel(alist));
    rel_stric = zeros(1, numel(alist));
    ratio     = zeros(1, numel(alist));
    for ia = 1:numel(alist)
        % 层名含 '-'（n-clad），不能做 struct 字段名 ⇒ 用 containers.Map 传 α 表
        Mj = containers.Map(); Mj(tags{j}) = alist(ia);
        S  = Model_Self_L2_7L_EBL('solve', struct('alpha', Mj));
        ai = 2*k0*imag(S(1).ne)*1e4;
        ref   = S(1).contrib(j) / alist(ia);          % = P_j·n_rj/Re(ne)
        ratio(ia)     = ai / alist(ia);
        rel_naive(ia) = ai/alist(ia)/P(j) - 1;
        rel_stric(ia) = ai/alist(ia)/ref - 1;
    end
    r1(j)  = ratio(1);  ratio1(j)  = ratio(1);
    ratio10(j) = ratio(3);
    rs(j)  = max(abs(rel_stric));
    lin(j) = abs(ratio(3)/ratio(1) - 1);
    fprintf('%-10s %9.4f %14.6f %13.6f %11.2e %11.2e\n', tags{j}, 100*P(j), ...
        ratio(1), ratio(1)/P(j), max(abs(rel_naive)), rs(j));
end

fprintf('\n===== B3 六层同时赋 α（各 ≤10 cm^-1）=====\n');
comb = containers.Map();
comb('n-clad') = 3; comb('n-wave') = 5; comb('QW')  = 10;
comb('p-wave') = 4; comb('EBL')    = 2; comb('p-clad') = 6;
S  = Model_Self_L2_7L_EBL('solve', struct('alpha', comb));
ai = 2*k0*imag(S(1).ne)*1e4;
naive = sum(P .* S(1).alpha);
stric = sum(S(1).contrib);
fprintf('alpha_mode = %.6f   ΣP·α = %.6f（偏 %.2e）   ΣP·α·n_r/n = %.6f（偏 %.2e）\n', ...
    ai, naive, abs(ai/naive-1), stric, abs(ai/stric-1));

B1 = max(rs) < 1e-2 && max(lin) < 1e-2;
B3 = abs(ai/stric - 1) < 2e-2;

fprintf('\n── 判定 ──\n');
fprintf('B1 弱吸收极限（严格式 <1e-2，线性度 <1e-2）: %d\n', B1);
fprintf('    严格偏差 max = %.2e    线性度 max = %.2e\n', max(rs), max(lin));
fprintf('B3 多通道叠加（严格式 <2e-2）              : %d   （%.2e）\n', B3, abs(ai/stric-1));
fprintf('S5 TEST pass=%d  (pass=%d fail=%d)\n', B1 && B3, B1+B3, 2-(B1+B3));
