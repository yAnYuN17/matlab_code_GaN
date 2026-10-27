% ─────────────────────────────────────────────────────
% 脚本：test_s7_channels（S7 验收入口：四通道分解 + 厚度扫描 + 衬底扫描）
% 输入：无
% 输出：命令行 PASS/FAIL 与末行 "S7 TEST pass=1/0"
% 调用函数：Model_Self_L2_7L_EBL('solve', opt)
% 说明：
%   1) 判据全部写成绝对偏差或与原文末位对齐的量，不用"看起来差不多"。
%   2) 2026-10-22 起全部数据走 'solve' 接口（thick / nsub 改写口），
%      不再依赖模型内的 's7' 动作——模型只负责解，判定留在这里。
%   3) Table I 系数抄录的反查单独放在 test_s7_selftest.m。
% ─────────────────────────────────────────────────────
base = fileparts(fileparts(mfilename('fullpath')));
addpath(base);

npass = 0; nfail = 0;
k0 = 2*pi/0.45;

% ── T1 基准：零吸收 ──
S0 = Model_Self_L2_7L_EBL('solve', struct('scenario', 'zero'));
d1 = abs(real(S0(1).ne) - 2.417128456985);
d2 = abs(S0(1).Gamma - 0.026775013);
d3 = abs(sum(S0(1).P_layer) - 1);
c1 = d1 < 1e-9 && d2 < 1e-9 && d3 < 1e-12;
fprintf('%-52s %s   dn=%.1e dG=%.1e dSumP=%.1e\n', ...
    'T1  基准 Re(n_eff)/Γ/ΣP（阶段 1 归档值）', pf(c1), d1, d2, d3);
npass = npass + c1; nfail = nfail + ~c1;

% ── T2 两场景 α_i（KIO 全层 / MIN 只给 EBL+p-clad）──
Sk = Model_Self_L2_7L_EBL('solve', struct('scenario', 'kio'));
Sm = Model_Self_L2_7L_EBL('solve', struct('scenario', 'min'));
ak = 2*k0*imag(Sk(1).ne)*1e4;
am = 2*k0*imag(Sm(1).ne)*1e4;
d4 = abs(ak - 6.1816); d5 = abs(am - 2.4578);
c2 = d4 < 5e-3 && d5 < 5e-3;
fprintf('%-52s %s   KIO=%.4f MIN=%.4f\n', ...
    'T2  450 nm TE0 α_i（KIO 6.1816 / MIN 2.4578）', pf(c2), ak, am);
npass = npass + c2; nfail = nfail + ~c2;

% ── T3 严格式 α_i = Σ P·α·n_r/n_eff（一阶微扰的残差 <1e-3）──
r1 = abs(ak/sum(Sk(1).contrib) - 1);
r2 = abs(am/sum(Sm(1).contrib) - 1);
c3 = r1 < 1e-3 && r2 < 1e-3;
fprintf('%-52s %s   min=%.2e max=%.2e\n', ...
    'T3  严格式 ΣP·α·n_r/n_eff（<1e-3）', pf(c3), min(r1, r2), max(r1, r2));
npass = npass + c3; nfail = nfail + ~c3;

% ── T4 通道求和 = 逐层贡献之和 ──
d6 = 0;
for k = [1 2]
    ch = Sk(k).ch;
    d6 = max(d6, abs(ch.sub + ch.n + ch.act + ch.p + ch.contact - sum(Sk(k).contrib)));
end
c4 = d6 < 1e-12;
fprintf('%-52s %s   max|Sum_ch - Sum_contrib|=%.1e\n', 'T4  通道求和自洽（<1e-12）', pf(c4), d6);
npass = npass + c4; nfail = nfail + ~c4;

% ── T5 p 侧通道主导（TE0，KIO）──
c5 = Sk(1).ch.p > Sk(1).ch.n && Sk(1).ch.p > 4*Sk(1).ch.act;
fprintf('%-52s %s   p=%.4f n=%.4f act=%.4f\n', 'T5  p 侧通道最大（TE0，KIO）', ...
    pf(c5), Sk(1).ch.p, Sk(1).ch.n, Sk(1).ch.act);
npass = npass + c5; nfail = nfail + ~c5;

% ── T6 厚度扫描（B5）：n-wave = p-wave = t ──
tlist = [0.04 0.06 0.08 0.10 0.12 0.15 0.18 0.22 0.26 0.30];
tt = []; GG = []; AA = []; G5 = [];
for it = 1:numel(tlist)
    th = containers.Map({'n-wave','p-wave'}, {tlist(it), tlist(it)});
    S = Model_Self_L2_7L_EBL('solve', struct('scenario', 'kio', 'thick', th));
    tt(end+1) = tlist(it);                       %#ok<SAGROW>
    GG(end+1) = S(1).Gamma;                      %#ok<SAGROW>
    AA(end+1) = 2*k0*imag(S(1).ne)*1e4;          %#ok<SAGROW>
    G5(end+1) = S(1).Gamma*500 - AA(end);        %#ok<SAGROW>
end
[~, ip] = max(GG);
[ap, ia] = max(AA);
c6 = abs(tt(ip) - 0.12) < 1e-12 && abs(ap - 6.6607) < 5e-3 && ...
     AA(end) < ap && all(G5 > 0);
fprintf('%-52s %s   Γ峰 t=%.2f；α峰 t=%.2f %.4f；α(0.30)=%.4f\n', ...
    'T6  厚度扫描（Γ峰 0.12、α峰 0.10、α 回落）', pf(c6), tt(ip), tt(ia), ap, AA(end));
npass = npass + c6; nfail = nfail + ~c6;

% ── T7 衬底扫描（B6）：束缚截止 ──
nsub = [1.70:0.05:2.40, 2.400:0.002:2.440];
ne = nan(1, numel(nsub)); bound = false(1, numel(nsub));
for k = 1:numel(nsub)
    S = Model_Self_L2_7L_EBL('solve', struct('scenario', 'kio', 'nsub', nsub(k)));
    if S(1).ok
        ne(k) = real(S(1).ne);
        bound(k) = ne(k) > nsub(k);
    end
end
nb = find(bound, 1, 'last');
th_v = nsub(nb);
c7 = th_v >= 2.410 && th_v < 2.420 && nb < numel(nsub) && ~S(1).ok;
fprintf('%-52s %s   n_sub_th=%.3f（下一个点无根）；Re(n_eff)=%.9f\n', ...
    'T7  束缚模截止阈值（2.410~2.420，且下一点无根）', pf(c7), th_v, ne(nb));
npass = npass + c7; nfail = nfail + ~c7;

fprintf('\nS7 TEST pass=%d  (pass=%d fail=%d)\n', nfail == 0, npass, nfail);

% ── 小工具 ──
function s = pf(c), if c, s = 'PASS'; else, s = 'FAIL'; end, end
