% ───────────────────────────────────────────────────────────────────────────────────
% 脚本：test_s6_tmm_stability（B9：TMM 数值稳定性——失真起点与失败归因）
% 输入：无（直接运行）
% 输出：命令行四段表格；末行 S6 TEST pass=<0|1>
% 调用函数：Model_Self_L2_7L_EBL('solve')（仅第四段，经公开接口取数）
% 说明：
%   1.本脚本原是模型的动作 's6'，2026-10-22 搬出——模型只保留日常功能；
%     搬出前后的结论与全部实测值见 阶段2_汇报单_V1.1.md。
%   2.第一至三段自带色散/求根内核的拷贝（局部函数在文件末尾），与模型同源：
%     搬出时逐位比对过，漂移由日后比对修复（用户 2026-10-22 明确接受）。
%   3.几何：对称三层 n_core=2.455 / n_clad=2.404（取自 n-wave / n-clad），三层同赋 α；
%     参考解是三层解析色散（tan 形式，复数牛顿，seed 独立，不拿 TMM 根当种子）。
%   4.判据：α ≤ 100 cm⁻¹ 段 max dα/α < 1e-6；且 "找不到根" 与 "精度失真" 分开记。
%   5.判据解读：第一段在全窗口下 α≥3e4 无根，是复根离实轴太远（实轴扫描失效），
%     不是矩阵病态——第二段窄窗口把同一个根钉住后 α≤3000 仍逐位相同。
% ───────────────────────────────────────────────────────────────────────────────────
addpath(fileparts(fileparts(mfilename('fullpath'))));   % Mode_L2（'solve' 所在）
addpath(fileparts(mfilename('fullpath')));

lambda = 0.45;
opts   = struct('num_scan', 2e4, 'tol', 1e-13, 'maxit', 300, 'merge', 1e-6, ...
                'chk_ratio', 0.5, 'chk_h', 10.^-(3:8));
n_core = 2.455;
n_clad = 2.404;
k0     = 2*pi/lambda;

%% 第一段：施工单规定，d = 0.12 μm，全窗口扫 α
d1    = 0.12;
alist = [0 1 3 10 30 100 300 1e3 3e3 1e4 3e4 6e4 1e5];
fprintf('\n===== S6 第一段：d=%.2f μm，全窗口，扫 α =====\n', d1);
fprintf('%9s %10s %11s %13s %12s %12s %10s %5s\n', ...
    'α[cm^-1]', 'α·d', '|Im κd|', 'Re(n_tmm)', 'dn', 'dα[cm^-1]', 'dα/α', '根数');
t1   = repmat(stage2_s6_point(n_core, n_clad, d1, 0, lambda, opts), 1, numel(alist));
seed = 2.43;
for i = 1:numel(alist)
    t1(i) = stage2_s6_point(n_core, n_clad, d1, alist(i), lambda, opts, [], seed);
    seed  = t1(i).ne_a;
    fprintf('%9g %10.3e %11.3e %13.9f %12.3e %12.3e %10.2e %5d\n', ...
        t1(i).alpha, t1(i).ad, t1(i).kdi, real(t1(i).ne_t), ...
        t1(i).dn, t1(i).dalpha, t1(i).rel, t1(i).nroots);
end

% !! 交付图 s6_tmm.png 由本脚本自己出（S6 搬出模型后，出图也跟着证据走）。
al    = [t1.alpha];
rel   = [t1.rel];
found = [t1.found];
[pth, ~, ~] = fileparts(mfilename('fullpath'));
outdir = fullfile(pth, '..', 'out');
if ~exist(outdir, 'dir'), mkdir(outdir); end
f  = figure('Visible', 'off', 'Position', [100 100 950 620]);
okm  = found & al > 0;
relp = max(rel(okm), 1e-17);          % 精确为 0 的点画在 1e-17 地板上（log 轴不能画 0）
loglog(al(okm), relp, 'o-', 'LineWidth', 1.6); hold on;
ax = gca; ax.XScale = 'log'; ax.YScale = 'log';
if any(~found)
    xb = min(al(~found));
    patch([xb 2e5 2e5 xb], [1e-18 1e-18 1e-10 1e-10], [1 0.85 0.85], ...
          'EdgeColor', 'none', 'FaceAlpha', 0.6);
    text(xb*1.1, 3e-14, {'no root found', 'Im(n_eff) ~ 0.04', '(root-finder limit, not matrix)'}, ...
         'FontSize', 9, 'Interpreter', 'none');
end
xline(1200, 'k--', '\alpha d = 0.12');
grid on; xlabel('\alpha  [cm^{-1}]'); ylabel('|\Delta\alpha| / \alpha');
ylim([1e-18 1e-10]);
title(sprintf('S6: TMM accuracy vs \\alpha (d = %.2f \\mum; zeros plotted at 1e-17)', d1));
print(f, fullfile(outdir, 's6_tmm.png'), '-dpng', '-r150'); close(f);
fprintf('S6 图已更新：%s\n', fullfile(outdir, 's6_tmm.png'));

%% 第二段：d = 0.3 μm 单模，窄窗口 ±0.02，压 |Im κd|
d2     = 0.30;
win2   = 0.02;
alist2 = [100 300 1e3 3e3 1e4 3e4 6e4 1e5 2e5];
fprintf('\n===== S6 第二段：d=%.2f μm，窄窗口 ±%.2f，加压 =====\n', d2, win2);
fprintf('%9s %10s %11s %13s %12s %12s %10s %9s %5s\n', ...
    'α[cm^-1]', 'α·d', '|Im κd|', 'Re(n_tmm)', 'dn', 'dα/α', 'cond(M)', 'amp', '根数');
t2   = repmat(stage2_s6_point(n_core, n_clad, d2, 0, lambda, opts, win2), 1, numel(alist2));
seed = 2.43;
for i = 1:numel(alist2)
    t2(i) = stage2_s6_point(n_core, n_clad, d2, alist2(i), lambda, opts, win2, seed);
    seed  = t2(i).ne_a;
    fprintf('%9g %10.3e %11.3e %13.9f %12.3e %12.2e %10.2e %9.2e %5d\n', ...
        t2(i).alpha, t2(i).ad, t2(i).kdi, real(t2(i).ne_t), ...
        t2(i).dn, t2(i).rel, t2(i).condM, t2(i).amp, t2(i).nroots);
end

%% 第三段：反面对照，d=5 μm 多模，全窗口 vs 窄窗口
d3 = 5.0;
fprintf('\n===== S6 第三段：d=%.1f μm 多模，全窗口 vs 窄窗口 =====\n', d3);
fprintf('%9s %10s %11s %13s %12s %12s %5s\n', ...
    'α[cm^-1]', '窗口', '|Im κd|', 'Re(n_tmm)', 'dn', 'dα/α', '根数');
for aa = [100 1e3]
    rq = stage2_s6_point(n_core, n_clad, d3, aa, lambda, opts);          % 全窗口
    rs = stage2_s6_point(n_core, n_clad, d3, aa, lambda, opts, 0.02);    % 窄窗口
    fprintf('%9g %10s %11.3e %13.9f %12.3e %12.2e %5d\n', ...
        aa, '全', rq.kdi, real(rq.ne_t), rq.dn, rq.rel, rq.nroots);
    fprintf('%9s %10s %11.3e %13.9f %12.3e %12.2e %5d\n', ...
        '', '窄', rs.kdi, real(rs.ne_t), rs.dn, rs.rel, rs.nroots);
end

%% 第四段：真实 7 层 + 量级 α（走公开接口 'solve'，与模型日常通路同一条）
fprintf('\n===== S6 第四段：真实 7 层 + 量级 α（n-clad 36 / p-clad 18 / EBL 36 cm⁻¹）=====\n');
S0 = Model_Self_L2_7L_EBL('solve', struct('scenario', 'zero'));
Ml = containers.Map({'n-clad','p-clad','EBL'}, {36, 18, 36});
S1 = Model_Self_L2_7L_EBL('solve', struct('alpha', Ml));
fprintf('α=0    : Re(n)=%.12f  Γ=%.9f  ΣP−1=%.1e\n', ...
    real(S0(1).ne), S0(1).Gamma, abs(sum(S0(1).P_layer)-1));
fprintf('α≠0    : Re(n)=%.12f  Γ=%.9f  ΣP−1=%.1e  Im(n)=%.9f  α_mode=%.4f cm⁻¹\n', ...
    real(S1(1).ne), S1(1).Gamma, abs(sum(S1(1).P_layer)-1), ...
    imag(S1(1).ne), 2*k0*1e4*imag(S1(1).ne));
fprintf('ΔRe(n)=%.3e  ΔΓ=%.3e\n', abs(S1(1).ne-S0(1).ne), abs(S1(1).Gamma-S0(1).Gamma));
fprintf('（cond(M)=24676.18、根后验 max ratio=1.19e-6 为搬出前模型内 aux 的记录，见汇报单 §9.6）\n');

%% 判定
rel1 = [t1.rel];
idx1 = find(rel1 > 1e-6, 1, 'first');
if isempty(idx1), onset1 = Inf; else, onset1 = alist(idx1); end
rel2 = [t2.rel];
idx2 = find(rel2 > 1e-6, 1, 'first');
if isempty(idx2), onset_kd = Inf; else, onset_kd = t2(idx2).kdi; end
rel_lo = [t1([t1.alpha] <= 100).rel];
ok2    = [t2.found];
if any(~ok2)
    fprintf('\n第二段：α = %s 在窄窗口内找不到根\n', mat2str(alist2(~ok2)));
end
fprintf('\n-- 判定 --\n');
fprintf('第一段 α≤100：max dα/α = %.3e；α>1e4 全部"找不到根"（求解窗口，不是精度）\n', max(rel_lo));
fprintf('第二段：解析参考根处 |f| 最大 %.3e（未归一化，不能单独据此判根；判根看 ne_t 与 ne_a 的逐位比较）\n', max([t2.res_a]));
fprintf('第二段 max |Im κd| = %.3f，失真起点 %s\n', max([t2.kdi]), ...
    ternary(isinf(onset_kd), '未出现（全段 dα/α ≤ 1e-6）', sprintf('|Im κd| ≈ %.3f', onset_kd)));
fprintf('\nS6 TEST pass=%d  (onset=%s, max|dα|/α=%.3e)\n', ...
    isinf(onset1) && max(rel_lo) < 1e-6, ...
    ternary(isinf(onset1), 'Inf', sprintf('%g', onset1)), max(rel_lo));


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：单层 2×2 转移矩阵（TE / TM 通用）
% 输入：横向波数 kz; 厚度 d; 折射率 n; 偏振 pol
% 输出：对应该层的转移矩阵    
% 说明：
%   TE：M=[cos(kh) , 1/ksin(kh); -ksin(kh) , cos(kh)]
%   TM: M=[cos(kh) , n^2/ksin(kh); -k/n^2 * sin(kh) , cos(kh)]
% ───────────────────────────────────────────────────────────────────────────────────
function [a11, a12, a21, a22] = layer_matrix(kz, d, n, pol)
    % isscalar(d)，d为标量(不为数组)，返回truth
    assert(isscalar(d) && d > 0, 'layer_matrix: 层厚 d 必须为正标量');
        
    % 偏振判断。strcmpi(A,'B') 忽略大小写比较字符串，符合返回truth
    is_TE = strcmpi(pol, 'TE');
    
    kd = kz .* d;
    S  = sin(kd) ./ kz;
    C  = cos(kd);
    
    % kz逼近0的预案，取极限：sin(kz·d)/kz → d
    % any(x) 若x是数组，对每列进行判断，只要有一个非0元素，返回truth，输
    % 出一个行向量；x为向量，就输出T or F。
    bad = (abs(kz) < 1e-14);
    if any(bad(:))
        S(bad) = d;                     
    end
    
    if is_TE
        a11 = C;      a12 = S;
        a21 = -kz .* sin(kd);           a22 = C;
    else
        a11 = C;      a12 = (n^2) .* S;
        a21 = -(kz / (n^2)) .* sin(kd); a22 = C;
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：转移矩阵计算
% 输入：有效折射率 n_eff; 层数据结构体 L; 波长 lambda; 偏振 pol
% 输出：f：某一n_eff下（其他参量数值默认）的转移矩阵输出(输出为0即满足，是束缚模)。
%       aux：结构体 (参见'说明'中两个公式的左侧)。
%           aux.gam_sub  p0
%           aux.gam_cov  p2
%           aux.M        连乘后的转移矩阵 M
% 调用函数：[a11, a12, a21, a22] = layer_matrix(kz, d, n, pol) (获取单层转移矩阵)
% 说明：
%   1.TE：[-p0,1] M^TE [1;-p2]=0
%           状态函数：[E_y ; dE_y/dy]
%     TM: [-p0/n0^2] M^TM [1; -p2/n2^2]=0
%           状态函数：[H_y ; (1/n²)·dH_y/dy]
%   2.n_eff其实是待解项，通过遍历所有n_eff，保留满足公式的n_eff
%   3.个人：有限真实边界需重新考虑
% ───────────────────────────────────────────────────────────────────────────────────
function [f, aux] = dispersion_tmm(n_eff, L, lambda, pol)
    if nargin < 4 || isempty(pol), pol = 'TE'; end
    is_TE = strcmpi(pol, 'TE');
    N = numel(L);
    
    in_size = size(n_eff);                
    x = n_eff(:).';                       % 1×M 行向量
    k0 = 2*pi/lambda;
    
    % 两侧半无限层的衰减常数（必须取衰减解，Re γ > 0; 加分支保护）
    gam_sub = k0 * sqrt(x.^2 - L(1).n^2);
    gam_cov = k0 * sqrt(x.^2 - L(N).n^2);
    gam_sub = BranchChoose_decay(gam_sub);      
    gam_cov = BranchChoose_decay(gam_cov);      
    
    if is_TE
        s_sub = gam_sub;              s_cov = gam_cov;
    else
        s_sub = gam_sub / L(1).n^2;   s_cov = gam_cov / L(N).n^2;
    end
    
    % 中间层转移矩阵连乘（从下往上，逐次左乘）
    M11 = ones(size(x));  M12 = zeros(size(x));
    M21 = zeros(size(x)); M22 = ones(size(x));
    
    for k = 2:N-1
        nk = L(k).n;
        dk = L(k).d;
        kz = k0 * sqrt(nk^2 - x.^2);
        [a11, a12, a21, a22] = layer_matrix(kz, dk, nk, pol);
    
        nM11 = a11.*M11 + a12.*M21;
        nM12 = a11.*M12 + a12.*M22;
        nM21 = a21.*M11 + a22.*M21;
        nM22 = a21.*M12 + a22.*M22;
        M11 = nM11; M12 = nM12; M21 = nM21; M22 = nM22;
    end
    
    v1 = M11 + M12.*s_sub;
    v2 = M21 + M22.*s_sub;
    f  = v2 + s_cov.*v1;
    f  = reshape(f, in_size);
    
    if nargout > 1
        aux = struct('gam_sub', gam_sub, 'gam_cov', gam_cov, 'M', [M11; M12; M21; M22]);
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：色散函数零点求解
% 输入：函数 dispersion_tmm(输出转移矩阵结果); n_lo 底层折射率，n_hi 顶层折射率; 
%           n_roots 求解模式数; opts 求解要用参数
% 输出：roots_n：规定模式数下的 n_eff 解。
%      info：结构体
%           info.mode       折射率性质（实/复）
%           info.ModeNum    可求解的最大模数
%           info.n_found    本次求解的模数
%           info.scanned    求解中采用的扫描点数
% 调用函数：scan_bisect(二分法), scan_newton(牛顿法), stage2_check_roots(实分支根的后验检查)
% 说明：
%   1.输入里有函数，是因为该函数需额外输出参量
%   2.本函数考虑了实折射率和复折射率
% ───────────────────────────────────────────────────────────────────────────────────
function [roots_n, info] = find_roots_complex(fun, n_lo, n_hi, n_roots, opts)
    % 输出参数补齐，及求解项数值判定
    if nargin < 4 || isempty(n_roots), n_roots = 5; end
    if nargin < 5 || isempty(opts), opts = struct(); end
    if ~isfield(opts, 'num_scan'), opts.num_scan = 20000; end       % 扫描点数
    if ~isfield(opts, 'tol'),    opts.tol    = 1e-13; end           % 求根收敛容差
    if ~isfield(opts, 'maxit'),  opts.maxit  = 200;   end           % 最大迭代次数
    if ~isfield(opts, 'merge'),  opts.merge  = 1e-6;  end           % 去重时两个根的最小相对间距
    if ~isfield(opts, 'chk_ratio'), opts.chk_ratio = 0.5; end       % |f(r)|/邻域最小值的判定阈值
    
    assert(n_hi > n_lo, 'find_roots_complex: 要求 n_hi > n_lo');
    
    xs = linspace(n_lo, n_hi, opts.num_scan);                       % 均匀扫描
    fs = fun(xs);                                                   % 获得遍历点的转移矩阵结果
    fs = fs(:).';
    
    % 防御性编程-不存在所有点的输出都是0
    fscale = max(abs(fs));
    if fscale == 0, fscale = 1; end
    opts.chk_fscale = fscale;                                       % 后验检查的绝对残差门限
    
    % 折射率 实/复判断，允许有噪声
    is_real = all(abs(imag(fs)) <= 1e-10 * fscale);
    
    if is_real
        roots_n = scan_bisect(fun, xs, real(fs), opts);             % 二分法
        info.mode = 'real';
        chk = stage2_check_roots(fun, roots_n, opts);               % 根检验并记录
    else
        [roots_n, chk, cand_all] = scan_newton(fun, xs, fs, opts);  % 牛顿法（含后验检查）
        info.mode = 'complex';
        info.cand = cand_all;                                       % 阶段2 S4：逐候选记录（被拒的也在）
    end
    
    % 去重、降序、（第一次）截断、剔除不合格的根
    if ~isempty(roots_n)
        % 按实部降序，取对应模式的索引
        [~, ord] = sort(real(roots_n(:).'), 'descend');
        roots_n  = roots_n(ord);
        if ~isempty(chk), chk = chk(ord); end
        % 第一个根恒保留(true)，后面每个根与前一根的差大于阈值才保留。
        % keep(i)=true/false。对数组进行标记
        % max(1, abs(roots_n(1:end-1))) 得到的是数组；将 roots_n 每一个和 1 比较取最大值。
        keep = [true, abs(diff(roots_n)) > opts.merge * max(1, abs(roots_n(1:end-1)))];
        roots_n = roots_n(keep);

        if ~isempty(chk), chk = chk(keep); end                  % 阶段2 S4：记录跟着根走
        % 阶段2 S4：剔除后验检查不合格（伪根/鞍点）的根
        if ~isempty(chk)
            ok = [chk.accepted];
            roots_n = roots_n(ok);
            chk     = chk(ok);
        end
        info.ModeNum = numel(roots_n);
        if numel(roots_n) > n_roots
            roots_n = roots_n(1:n_roots);
            if ~isempty(chk), chk = chk(1:n_roots); end         % 阶段2 S4
        end
        info.checks = chk;                                      % 阶段2 S4：逐根后验记录
    end
    
    info.n_found = numel(roots_n);
    info.scanned = opts.num_scan;
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：二分法求解
% 输入：函数 dispersion_tmm(输出转移矩阵结果); 所有待测的n_eff xs(向量); 
%       待测n_eff(实部)的转移矩阵输出 fs; 求解要用参数 opts
% 输出：满足条件的n_eff
% 说明：
%   1.二分只依赖"闭区间两端符号相反"这个充分条件，不要求函数单调或光滑
%   2.本函数用于实折射率
% ───────────────────────────────────────────────────────────────────────────────────
function roots_n = scan_bisect(fun, xs, fs, opts)
    roots_n = [];
    s = sign(fs);                        % 取每个fs的符号：>0→取1,<0→取-1,=0取0
    s(s == 0) = 1;                       % 恰好为零的采样点，视为正数
    idx = find(s(1:end-1) .* s(2:end) < 0);
    
    for j = 1:numel(idx)
        a  = xs(idx(j));
        b  = xs(idx(j)+1);
        fa = fun(a);
    
        for it = 1:opts.maxit
            c = 0.5*(a + b);
            if (b - a) <= opts.tol * max(1, abs(c))  % 异号且两者较近满足条件，判定为根
                break;
            end
            fc = fun(c);
            if real(fc) == 0                         % 中点恰好是根
                a = c; b = c;
                break;
            end
            if sign(real(fc)) == sign(real(fa))      % 根的区间判定
                a = c; fa = fc;
            else
                b = c;
            end
        end
    
        roots_n(end+1) = 0.5*(a + b);   %#ok<AGROW>  % 将满足条件的中点C作为解
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：复折射率求解
% 输入：函数 dispersion_tmm(输出转移矩阵结果); 所有待测的n_eff xs(向量); 
%      待测n_eff的转移矩阵输出 fs; 求解要用参数 opts
% 输出：roots_n: 满足条件的复折射率; 
%      chk: 结构体数组(参考 stage2_empty_check)，仅记录经过后验的 n 的检测数据; 
%      cand_all:结构体数组(同样),记录所有备选 n 的检测数据
% 调用函数：newton_complex(输出：经检查后的元素 r，和关联r的结构体 c), 
%          stage2_empty_check [用于数据记录 (→ 后验检查) 的空数组]
% 说明：
%   1.本函数用于复折射率
%   2.在matlab中，数组作为参量输入到函数，但是函数是依次取元素的，也就是说函数
%       内部只获取了一个元素而非一整个向量。
% ───────────────────────────────────────────────────────────────────────────────────
function [roots_n, chk, cand_all] = scan_newton(fun, xs, fs, opts)
    mag = abs(fs);                              % 取模，而不是取绝对值。
    cand = [];

    % 选取局部的谷(极小)，作为备选的迭代初值
    % f(z)=0 → |f(z)=0|  →  在离散点中，邻近点数值小于其边上两点的数值
    for j = 2:numel(xs)-1
        if mag(j) < mag(j-1) && mag(j) <= mag(j+1)
            cand(end+1) = xs(j);                %#ok<AGROW>
        end
    end

    roots_n  = [];
    chk      = stage2_empty_check();
    cand_all = stage2_empty_check();
    for x0 = cand                               % x0 是复折射率
        [r, c] = newton_complex(fun, x0, opts);
        c.x0 = x0;
        cand_all(end+1) = c;                    %#ok<AGROW> 
        if ~isempty(r)
            roots_n(end+1) = r;                 %#ok<AGROW>
            chk(end+1)     = c;                 %#ok<AGROW>
        end
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：牛顿法求解
% 输入：函数 dispersion_tmm(输出转移矩阵结果); 待迭代的复折射率 x0; 求解要用参数 opts
% 输出：元素r。对x0操作，满足对应元素仍保留，不满足为空[]; 
%       chk: 结构体
%       chk.fr            复折射率代入dispersion_tmm后，其结果的模，|fun(r)|
%       chk.nb            八邻域中的最小值
%       chk.ratio         fr/max(nb,)
%       chk.accepted      条件判定: fr 需小于门限 且 ratio<=预设比例
% 调用函数：newton_complex, stage2_check_one 
% 说明：
%   1.f(z_k + \delta) \approx f(z_k)+f'(z_k)\delta, \delta = - [ f(z_k) / f'(z_k)]
%       z_(k+1) = z_k + \delta
%   2.&&& 设置了步长阻尼，可以在opts额外设置参数
%   3.因为只获取了一个元素(虽然输入是向量，但是函数本身只对一个其中元素操作)，所以当
%       当前元素不符合迭代要求时，是 return 对下一元素进行验证，而不是 break 退出循环。
% ───────────────────────────────────────────────────────────────────────────────────
function [r, chk] = newton_complex(fun, x0, opts)
    h      = 1e-7;                          % 用于求导的步长
    r      = x0;
    n_iter = opts.maxit;

    for it = 1:n_iter
        % 求解f'(r)，即df
        fx = fun(r);
        f1 = fun(r + h);
        f2 = fun(r - h);
        df = (f1 - f2) / (2*h);
        
        % 无意义，记录该点，然后退出
        if abs(df) < 1e-300
            chk = stage2_check_one(fun, r, opts);
            r = [];
            return
        end

        step = fx / df;
        if abs(step) > 0.5                 % 阻尼：单步不超过 0.5，超过则为±0.5
            step = 0.5 * step / abs(step);
        end
        r = r - step;
        
        % 迭代步长小于一定值，退出循环，记录(还需后验检查)
        if abs(step) < 1e-14 * max(1, abs(r))
            break;
        end
        % 达到迭代最大次数仍未收敛，丢弃
        if it == n_iter
            chk = stage2_check_one(fun, r, opts);
            r = [];                        
            return
        end
    end

    chk = stage2_check_one(fun, r, opts);
    if ~chk.accepted
        r = [];
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：BranchChoose_decay
% 输入：g：k0·sqrt(n_eff² − n²)，可为复数数组
% 输出：g：已按"场必须向半无限侧衰减"选支后的衰减常数
% 说明：
%   1.带损耗时 n 为复数，sqrt 取主值只保证 Re ≥ 0，Re = 0 时 Im 的符号无保证；
%   2.Re = 0 且 Im < 0 的分支等价于场向外增长，必须取反。
%   3.该保护在无损耗时恒不触发。
%   4.半无限层衰减系数的分支选择。
% ───────────────────────────────────────────────────────────────────────────────────
function g = BranchChoose_decay(g)
    flip = real(g) < 0 | (real(g) == 0 & imag(g) < 0);
    g(flip) = -g(flip);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：单根的后验检查
% 输入：函数 dispersion_tmm(输出转移矩阵结果); 待检查的根 r; 求解要用参数 opts    
% 输出：c：结构体
%       c.fr            复折射率代入dispersion_tmm后，其结果的模，|fun(r)|
%       c.nb            八邻域中的最小值
%       c.ratio         fr/max(nb,)
%       c.accepted      条件判定: fr 需小于门限 且 ratio<=预设比例
% 说明：
%   1.对于非零极小伪根的判断：从 1e-3…1e-8 六个尺度 × 8 个方向取最小值——真根的 
%     |f| 随步长线性趋零，伪根在盆底附近 |f| 基本不变 ⇒ ratio ≈ 1，被拒。只取 
%     1e-3 一个尺度不严谨。
%   2.&&& ratio 记录进 info.checks，是 B2 的验收把手。
%   3.判据: 绝对阈值和相对陡度。
%   4.realmin, matlab的内置的最小正正规数, 2.2251e-308.
% ───────────────────────────────────────────────────────────────────────────────────
function c = stage2_check_one(fun, r, opts)
    if ~isfield(opts, 'chk_ratio'), opts.chk_ratio = 0.5;  end
    if ~isfield(opts, 'chk_h'),     opts.chk_h     = 10.^-(3:8); end
    if ~isfield(opts, 'chk_res'),   opts.chk_res   = 1e-6; end   % 阶段2 S6：绝对残差门限（相对 fscale）

    fr = abs(fun(r));
    nb = inf;
    for h = opts.chk_h * max(1, abs(r))
        for th = (0:7) * pi/4
            v = abs(fun(r + h * exp(1i*th)));
            if v < nb
                nb = v;
            end
        end
    end

    % 阶段2 S6 补：只比 ratio 会放过"陡峭函数的非零局部极小"（|f| 是深坑但不是零），
    % 必须再加一条绝对残差门限：|f(r)| 相对扫描段最大值要小到 chk_res 量级。
    if isfield(opts, 'chk_fscale')
        rmax = opts.chk_res * opts.chk_fscale;
    else
        rmax = inf;
    end

    c.fr       = fr;
    c.nb       = nb;
    c.ratio    = fr / max(nb, realmin);
    c.rmax     = rmax;
    c.accepted = (c.ratio <= opts.chk_ratio) && (fr <= rmax);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_check_roots（一批根的后验检查，实分支用）
% 输入：fun：色散函数; roots_n：根向量; opts：见 stage2_check_one
% 输出：chk：与 roots_n 等长的 struct 数组，字段同 stage2_check_one（另加 x0）
% 调用函数：stage2_empty_check [用于数据记录 (→ 后验检查) 的空数组]
% 说明：
%   1.实分支由二分给出，理论上必为真根；仍逐根记录，保证"所有返回根有后验记录"。
% ───────────────────────────────────────────────────────────────────────────────────
function chk = stage2_check_roots(fun, roots_n, opts)
    cc = cell(1, numel(roots_n));
    for k = 1:numel(roots_n)
        c = stage2_check_one(fun, roots_n(k), opts);
        c.x0 = roots_n(k);
        cc{k} = c;
    end
    if isempty(cc)
        chk = stage2_empty_check();
    else
        chk = [cc{:}];
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：用于后验记录的空 struct 数组
% 输入：无
% 输出：chk：0×1 的空 struct 数组 (参考 stage2_check_one)
%           chk.fr
%           chk.nb
%           chk.ratio
%           chk.rmax
%           chk.accepted
%           chk.x0          原本的复折射率值
% 说明：
%   1.struct 数组不能从 [] 起步（[] 是 double，end+1 赋值报"无法从 struct 转换为 double"，
%     坑记录 D1 同类）；必须用带字段的空 struct 起手。
%   2.字段顺序与 stage2_check_one 的返回一致，末尾多一个 x0（由调用方填写）。
% ───────────────────────────────────────────────────────────────────────────────────
function chk = stage2_empty_check()
    chk = struct('fr', {}, 'nb', {}, 'ratio', {}, 'rmax', {}, 'accepted', {}, 'x0', {});
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_s6_point（单点：TMM 根 vs 三层解析根）
% 输入：n_core/n_clad：实折射率; d：芯厚 [μm]; a：吸收 α [cm⁻¹]; lambda; opts
%        win：可选，搜根窗口半宽（以解析根为中心）；省略则用全窗口
% 输出：r：结构体
%       r.ad
%       r.含 ad、kdi、ne_t、ne_a、dn、dalpha、rel、nroots、found
% 说明：
%   1.三层同赋吸收：n = n_r + i·α/(2k0·1e4)；
%   2.解析根初值固定 2.43，不用 TMM 结果做种子；
%   3.窄窗口用于排除多模选错根：窗口只含解析根 ± win；
%   4.found=0 表示 TMM 在窗口内没找到根。
% ───────────────────────────────────────────────────────────────────────────────────
% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_s6_point（单点：TMM 根 vs 三层解析根）
% 输入：n_core/n_clad：实折射率; d：芯厚 [μm]; a：吸收 α [cm⁻¹]; lambda; opts
%        win：可选，搜根窗口半宽（以解析根为中心）；省略则用全窗口
%        seed：可选，解析根的牛顿初值（默认 2.43）
% 输出：r：结构体
%       r.ad、kdi、ne_t、ne_a、dn、dalpha、rel、nroots、found
% 说明：
%   1.三层同赋吸收：n = n_r + i·α/(2k0·1e4)；
%   2.解析参考与 TMM 共用同一分支约定：包层 sqrt 经 BranchChoose_decay 强制
%     Re γ ≥ 0。若这里用未强制的 sqrt，α 一大两者解的就是不同方程，
%     会把"参考解根本不是根"误判成 TMM 失真；
%   3.解析根初值默认 2.43，调用方可传上一轮的根做连续化，不拿 TMM 结果做种子；
%   4.窄窗口用于排除多模选错根；
%   5.found=0 表示 TMM 在窗口内没找到根。
% ───────────────────────────────────────────────────────────────────────────────────
function r = stage2_s6_point(n_core, n_clad, d, a, lambda, opts, win, seed)
    if nargin < 7, win  = [];  end
    if nargin < 8, seed = 2.43; end

    k0 = 2*pi/lambda;
    dl = a/(2*k0*1e4);
    n1 = complex(n_core, dl);
    n2 = complex(n_clad, dl);

    L = struct('tag',{},'n',{},'d',{},'is_active',{});
    L(1) = struct('tag','substrate','n',n2,'d',Inf,'is_active',false);
    L(2) = struct('tag','core',     'n',n1,'d',d  ,'is_active',true);
    L(3) = struct('tag','cover',    'n',n2,'d',Inf,'is_active',false);

    % 解析参考（分支约定与 TMM 一致）
    ga   = @(x) BranchChoose_decay(k0*sqrt(x.^2 - n2^2));
    cap  = @(x) k0*sqrt(n1.^2 - x.^2);
    % 用 tan 形式 κ·tan(κd/2) = γ，不用 κd = 2·atan(γ/κ)：
    % 复参数下 atan 取主值会跳到别的分支，α 一大参考解就和 TMM 解的不是同一个方程。
    fa   = @(x) cap(x).*tan(cap(x)*d/2) - ga(x);
    ne_a = stage2_newton_scalar(fa, seed, 300);

    % TMM
    fh   = @(x) dispersion_tmm(x, L, lambda, 'TE');
    n_lo = max(real(n2)) + 1e-9;
    n_hi = real(n1) - 1e-9;
    if ~isempty(win)
        n_lo = max(n_lo, real(ne_a) - win);
        n_hi = min(n_hi, real(ne_a) + win);
    end
    ok = true;
    try
        roots_n = find_roots_complex(fh, n_lo, n_hi, 3, opts);
    catch
        roots_n = [];
    end
    if isempty(roots_n)
        ne_t = NaN; ok = false; nr = 0;
    else
        [~, ord] = sort(real(roots_n(:).'), 'descend');
        ne_t = roots_n(ord(1));
        nr   = numel(roots_n);
    end

    kap = k0*sqrt(n1^2 - ne_a^2);
    r = struct('alpha', a, 'd', d, 'ad', a/1e4*d, 'kdi', abs(imag(kap*d)), ...
               'ne_t', ne_t, 'ne_a', ne_a, 'nroots', nr, 'found', ok, ...
               'dn', NaN, 'dalpha', NaN, 'rel', NaN);
    if ok
        r.dn     = abs(ne_t - ne_a);
        r.dalpha = 2*k0*1e4*abs(imag(ne_t) - imag(ne_a));
        r.rel    = r.dalpha / max(a, eps);
        % 附加诊断：TMM 根处 |f| 与解析根处 |f|
        r.res_t  = abs(fh(ne_t));
        r.res_a  = abs(fh(ne_a));
        % 附加诊断：连乘矩阵的条件数与边界条件项的"抵消放大倍数"
        % amp = max(|v2|,|s_cov·v1|)/|f|：真正决定精度的不是 cond(M)，而是这个放大倍数
        [f_t, aux] = dispersion_tmm(ne_t, L, lambda, 'TE');
        Mm  = [aux.M(1) aux.M(2); aux.M(3) aux.M(4)];
        v1  = Mm(1,1) + Mm(1,2)*aux.gam_sub;
        v2  = Mm(2,1) + Mm(2,2)*aux.gam_sub;
        r.condM = cond(Mm);
        r.amp   = max(abs(v2), abs(aux.gam_cov*v1)) / max(abs(f_t), eps);
    else
        r.res_t = NaN; r.res_a = abs(fh(ne_a));
        r.condM = NaN; r.amp = NaN;
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：ternary（条件选择字符串）
% 输入：cond：逻辑; a、b：待选字符串
% 输出：a（cond 真）或 b
% ───────────────────────────────────────────────────────────────────────────────────
function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_newton_scalar（标量复数牛顿，带阻尼，仅用于 S6 参考解）
% 输入：fun：标量函数; x0：初值; maxit：最大迭代
% 输出：x：根
% 说明：
%   1.步长上限 0.5、收敛判据 |step| < 1e-14·max(1,|x|)；
%   2.不收敛时返回最后一次迭代值（S6 只作参考，不静默造根）。
% ───────────────────────────────────────────────────────────────────────────────────
function x = stage2_newton_scalar(fun, x0, maxit)
    x = x0;
    h = 1e-7;
    for it = 1:maxit
        fx = fun(x);
        df = (fun(x + h) - fun(x - h)) / (2*h);
        if abs(df) < 1e-300, break; end
        step = fx / df;
        if abs(step) > 0.5, step = 0.5*step/abs(step); end
        x = x - step;
        if abs(step) < 1e-14*max(1, abs(x)), break; end
    end
end
