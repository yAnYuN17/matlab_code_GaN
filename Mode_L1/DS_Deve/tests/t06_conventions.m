function t06_conventions()
%T06_CONVENTIONS  A7：归一化定义与单位约定的一致性自检
%
%   验收要求：代码里有明确的归一化定义注释 + 单位约定注释。
%
%   这一条通常靠人眼检查，但"人眼检查过"和"每次改动后仍然成立"是两回事。
%   本脚本把它变成可重复执行的检查，并在最后加了一个真正有技术含量的部分：
%   【复折射率接口冒烟测试】—— 阶段1 全程 n'' = 0，接口是否真的为阶段2 留好了，
%   只有拿一个带虚部的折射率真跑一遍才知道。
%
%   ⚠ 复折射率测试【不判失败】：阶段1 的设计目标不含复根求解，
%     这里只报告"已经能跑"还是"阶段2 需要补什么"。若已能跑出复数 n_eff，
%     正确性判据是 Im(n_eff) > 0（吸收介质 ⇒ 损耗 ⇒ 按本项目约定为正）。

fprintf('\n========== t06 约定自检（A7）==========\n');

root = fileparts(fileparts(mfilename('fullpath')));
ok_all = true;

% ── (1) 关键注释是否在位 ─────────────────────────────────────
fprintf('\n── 注释与约定 ──\n');

req = {
    fullfile(root,'startup.m'),                  '单位约定',      '单位约定（全项目唯一定义处）'
    fullfile(root,'src','normalize_mode.m'),     '∫',             '归一化积分的定义'
    fullfile(root,'src','normalize_mode.m'),     'TE 模',         'TE/TM 权重区分'
    fullfile(root,'src','confinement_factor.m'), 'Γ ≡',           'Γ 的定义式'
    fullfile(root,'src','confinement_factor.m'), '有源区',         'Γ 分子的范围约定'
    fullfile(root,'src','build_layers.m'),       '半无限',         '最外层半无限的要求'
    fullfile(root,'src','solve_slab.m'),         '复折射率',       '复折射率接口说明'
    fullfile(root,'src','dispersion_analytic3.m'),'截止',          '截止条件'
};

for i = 1:size(req,1)
    f = req{i,1};
    if ~isfile(f)
        fprintf('   [FAIL] 文件不存在：%s\n', f);
        ok_all = false;
        continue;
    end
    txt = fileread(f);
    has = contains(txt, req{i,2});
    if ~has, ok_all = false; end
    fprintf('   [%s] %-42s ← %s\n', ternary(has), ...
            sprintf('%s : "%s"', relpath(f, root), req{i,2}), req{i,3});
end

% ── (2) 数值约定自洽性 ───────────────────────────────────────
fprintf('\n── 数值约定 ──\n');

lam = 0.45;
L   = build_layers(layers_gan_ld_typical());
m   = solve_slab(L, lam, struct('pol','TE','n_modes',2));
k0  = 2*pi/lam;

maxb = 0;
for k = 1:numel(m)
    maxb = max(maxb, abs(m(k).beta - k0*m(k).n_eff));
end
ok_beta = maxb < 1e-12;
report(ok_beta, sprintf('β = k0·n_eff 一致（最大偏差 %.2e）', maxb));

% 单位换算：α [cm^-1] = 2*k0*Im(n_eff)*1e4
% 阶段1 全实数，Im(n_eff) = 0 ⇒ α = 0
maxIm = max(abs(imag([m.n_eff])));
ok_real = maxIm == 0;
report(ok_real, sprintf('阶段1 的 n_eff 为实数（max|Im| = %.1e）', maxIm));

% 归一化自检值
maxnv = max(abs([m.norm_val] - 1));
report(maxnv < 1e-10, sprintf('归一化自检值恒为 1（最大偏差 %.2e）', maxnv));

% ── (3) 复折射率接口冒烟测试（阶段2 的预备检查）─────────────
fprintf('\n── 复折射率接口冒烟测试（阶段2 预备）──\n');

loss = 1e-4;                       % 给两侧包层加一点吸收
Lc = build_layers(layers_three(3.40 + 1i*loss, 3.50, 3.20 + 1i*loss, 1.50));

fprintf('   输入：n_sub = 3.4000 + %.1ei,  n_cov = 3.2000 + %.1ei\n', loss, loss);

try
    mc = solve_slab(Lc, 1.0, struct('solver','tmm','n_modes',1));
    if isempty(mc)
        fprintf('   [INFO] 复数分支未找到根 —— 阶段2 需要实现复数根搜索\n');
    else
        ne = mc(1).n_eff;
        fprintf('   输出：n_eff = %.10f %+.10fi\n', real(ne), imag(ne));
        if imag(ne) > 0
            fprintf('   [PASS] Im(n_eff) > 0，符号与"吸收介质 → 损耗"的约定一致\n');
        else
            fprintf('   [WARN] Im(n_eff) ≤ 0，符号约定或根搜索分支需在阶段2 复核\n');
        end
    end
catch ME
    fprintf('   [INFO] 复数分支尚未接通（%s）—— 属预期，阶段2 实现\n', ME.identifier);
end

% ── 结论 ─────────────────────────────────────────────────────
fprintf('\n── 结论 ──\n');
if ok_all && ok_beta && ok_real && maxnv < 1e-10
    fprintf('t06 通过。\n');
else
    fprintf('t06 未通过。\n');
end
end

% ══════════════════════════════════════════════════════════════════════
function s = ternary(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

% ══════════════════════════════════════════════════════════════════════
function r = relpath(f, root)
r = strrep(f, [root filesep], '');
end

% ══════════════════════════════════════════════════════════════════════
function report(ok, msg)
if ok
    fprintf('   [PASS] %s\n', msg);
else
    fprintf('   [FAIL] %s\n', msg);
end
end
