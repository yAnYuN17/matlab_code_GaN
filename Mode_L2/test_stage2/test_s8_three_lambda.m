% ───────────────────────────────────────────────────────────────────────────────────
% 脚本：test_s8_three_lambda（B8：405/450/520 nm 同结构对比）
% 输入：无
% 输出：命令行表格 + B8 判定 T1–T4（末行 S8 TEST pass=…）
% 调用函数：Model_Self_L2_7L_EBL('solve', opt)
% 说明：
%   1.三档波长共用 450 nm 的折射率（没有可引用的 n(λ) 来源），随 λ 变的只有
%     k0 与 Kioupakis 的 α(λ) ⇒ 405/520 只能报趋势，不能报绝对 α_i 结论；
%   2.模型文件已不带 s8 入口（B 步），本脚本全部经 'solve' 取数，物理仍在模型里；
%   3.判据：T1 450 nm 锚定 A 步的 α_i（相对 <1e-3）；T2 α_i 随 λ 单调增；
%     T3 每层贡献之和 = α_i（相对 <1e-3）；T4 Γ 随 λ 单调减。
% ───────────────────────────────────────────────────────────────────────────────────
addpath('C:/Users/L/Desktop/UG/GaN_LD/waveguide/Mode_L1/Mode_L1_V1.0');
addpath('C:/Users/L/Desktop/UG/GaN_LD/waveguide/Mode_L1/Mode_L1_V1.0/layers');
addpath('C:/Users/L/Desktop/UG/GaN_LD/waveguide/Mode_L2');

lam = [0.405, 0.450, 0.520];
k0s = 2*pi ./ lam;

fprintf('声明：n(λ) 固定为 450 nm 的常数，405/520 nm 的 α 只作趋势，不作定量结论。\n');

S  = cell(1, numel(lam));
for i = 1:numel(lam)
    S{i} = Model_Self_L2_7L_EBL('solve', struct('lambda', lam(i), 'scenario', 'kio'));
    s    = S{i}(1);                                  % 只取 TE0
    ai   = 2*k0s(i)*imag(s.ne)*1e4;
    sc   = sum(s.contrib);
    fprintf('\nλ=%3.0f nm  Re(n_eff)=%.9f  Γ=%.6f  α_i=%.4f cm^-1（Σ贡献 %.4f，偏 %.1e）\n', ...
        lam(i)*1000, real(s.ne), s.Gamma, ai, sc, abs(sc/ai - 1));
    fprintf('  通道：n 侧 %.4f / 有源区 %.4f / p 侧 %.4f  [cm^-1]\n', s.ch.n, s.ch.act, s.ch.p);
    tg = s.tag;
    for j = 1:numel(tg)
        fprintf('    %-10s P=%8.4f%%   α=%9.3f   P·α·n_r/n=%8.4f\n', ...
            tg{j}, 100*s.P_layer(j), s.alpha(j), s.contrib(j));
    end
end

% ── B8 判定 ───────────────────────────────────────────────────────────────────────
ai_all = zeros(1, numel(lam));
for i = 1:numel(lam)
    ai_all(i) = 2*k0s(i)*imag(S{i}(1).ne)*1e4;
end

T1 = abs(ai_all(2)/6.181587 - 1) < 1e-3;             % 450 nm 锚 A 步（S7 首跑 6.181587）
T2 = all(diff(ai_all) > 0);
T3 = true;
for i = 1:numel(lam)
    T3 = T3 && abs(sum(S{i}(1).contrib)/ai_all(i) - 1) < 1e-3;
end
T4 = all(diff([S{1}(1).Gamma, S{2}(1).Gamma, S{3}(1).Gamma]) < 0);

fprintf('\n── B8 ──\n');
fprintf('T1 450 nm 锚定 A 步        : %d  (α_i=%.4f vs 6.181587)\n', T1, ai_all(2));
fprintf('T2 α_i 随 λ 单调增          : %d  (%.4f → %.4f → %.4f)\n', T2, ai_all);
fprintf('T3 Σ贡献 = α_i（相对<1e-3） : %d\n', T3);
fprintf('T4 Γ 随 λ 单调减            : %d  (%.6f → %.6f → %.6f)\n', T4, ...
    S{1}(1).Gamma, S{2}(1).Gamma, S{3}(1).Gamma);
fprintf('S8 TEST pass=%d  (pass=%d fail=%d)\n', all([T1 T2 T3 T4]), ...
    sum([T1 T2 T3 T4]), sum(~[T1 T2 T3 T4]));
