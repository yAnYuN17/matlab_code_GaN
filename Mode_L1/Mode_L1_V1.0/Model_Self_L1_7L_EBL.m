function modes = Model_Self_L1_7L_EBL(action)
%   ── 命名规则 ────────────────────────────────────────────────
%      函数文件，函数名必须与文件名<.m>一致
%      直接"运行"，以0输入调用该函数
%      可以在命令行调用，如"modes = Model_Self_L1_7L_EBL('test'); "实现有输入调用'save'
%      也可以在另一个文件中调用该函数
%
%   ── 单位与符号约定（全项目统一，不得另立门户）────────────────
%       长度单位 μm 
%       时间因子 exp(−iωt)，沿 +z 取 exp(+iβz)
%       ⇒ α = 2·k₀·Im(n_eff)  [cm⁻¹] 乘 1e4
%       ⇒ Im(n_eff)  >0 [损耗]，<0 [净增益]
%       Γ ≡ ∫_有源区 w dy / ∫_{-∞}^{+∞} w dy
%           TE: w = |E_y|²      TM: w = |H_y|²/n²
%
%   ── 说明 ──────────────────────────────────────────────────
%       远场发散角，未实现
    %% 初始化
    clc;                                                        %清空命令行窗口
    close all;                                                  %关闭所有已打开的图窗
    % clear;                                                      %清除变量
    
    
    % 模型选择（已有，或者自己再建；下述代码可自行添加或修改）
    % 没有输入或者输入为空，默认action为"show"。nargin为函数实际收到的输入参数数量。
    % 检查字符串中是否含有"save"和"9"，有则为真。转字符串转小写。
    if nargin < 1 || isempty(action)
        action = 'show'; 
    end
    do_save  = contains(lower(string(action)), 'save');
    use9     = contains(lower(string(action)), '9');
    usetest  = contains(lower(string(action)), 'test');
    
    if use9
        layers = stack_9layer_mqw();
    elseif usetest
        layers = test_layer();
    else
        layers = stack_7layer();
    end

    % 赋值
    lambda  = 0.45;                 % 真空波长 [μm] = 450 nm
    k0 = 2*pi/lambda;
    pol     = 'TE';                 % 'TE' | 'TM'
    n_modes = 3;                    % 最多求几阶模  
    L  = build_layers(layers);      % 校验层结构，折射率强制转 complex
    N  = numel(L);
    yb = layer_bounds(L);           % 定义各层的y区间
    savedir = 'out';                % 图片保存地址-后缀
    
    
    % 求解项赋值
    opts = struct( ...
        'n_scan',   20000, ...   % 色散函数扫描点数（根搜索用）
        'tol',      1e-13, ...   % 求根收敛容差（相对）
        'maxit',    200,   ...   % 最大迭代次数
        'merge',    1e-6,  ...   % 去重时两个根的最小相对间距
        'npts',     800,   ...   % 模场重建时每层采样点数
        'tail_dec', 20,    ...   % 半无限层向外延伸多少个衰减长度
        'tail_max', 50);         % 延伸长度上限 [μm]
    
    %% 色散函数求根 → n_eff
    % 束缚模存在条件判定
    n_clad_max = max(real(L(1).n), real(L(N).n));
    n_max      = max(real([L.n]));
    tiny       = 1e-9 * max(1, n_max);
    n_lo       = n_clad_max + tiny;
    n_hi       = n_max      - tiny;
    
    assert(n_hi > n_lo, ...
        '搜索区间为空（n_lo = %g ≥ n_hi = %g）：结构里必须存在折射率高于两侧包层的有限层', ...
        n_lo, n_hi);
    
    % 转移矩阵求解
    fh = @(x) dispersion_tmm(x, L, lambda, pol);
    [neff_list, find_date] = find_roots_complex(fh, n_lo, n_hi, n_modes + 2, opts);
    % 做了两次截断，因为求根时多求了几个解
    if numel(neff_list) > n_modes
        neff_list = neff_list(1:n_modes); 
    end
    neff_list = neff_list(:).';
    
    % 模式排序。n_eff 越大 ⇒ 模序越低。
    [~, ord]  = sort(real(neff_list), 'descend');
    neff_list = neff_list(ord);

    %% 重建模场、归一化、Γ求解   
    modes = struct([]);                      % 结构体数组
    for j = 1:numel(neff_list)
        ne = neff_list(j);
    
        [y, E, n_y]   = mode_profile(ne, L, lambda, pol, opts);
        [E, norm_val] = normalize_mode(y, E, n_y, pol);
        [G, det]      = confinement_factor(y, E, n_y, yb, L, pol);
    
        modes(j).n_eff    = ne;              % 有效折射率（本征值）
        modes(j).beta     = k0 * ne;         % 传播常数 [μm^-1]
        modes(j).order    = j - 1;           % 模序，0 = 基模
        modes(j).y        = y;               % 快轴坐标 [μm]
        modes(j).E        = E;               % 归一化后的模场 φ(y)
        modes(j).n_y      = n_y;             % 对应位置的折射率
        modes(j).Gamma    = G;               % 限制因子
        modes(j).norm_val = norm_val;        % 归一化积分值，应恒为 1
        modes(j).P_layer  = det.P_layer;     % 各层功率占比，总和 = 1
        modes(j).yb       = yb;
    end
    
    %% 输出 删除原命令行输出，保留自检。后续可调整为 UI 或者图片        

    fprintf('── 自检 ──\n');
    % 模式模场归一化自检
    fprintf('──模式归一化自检(应为1)──\n');
    for j = 1:numel(neff_list)
        fprintf('模式 %s%d: 自检数 = %.3e \n',pol, modes(j).order ,modes(j).norm_val);
    end

    % 模式限制因子自检
    % 限制因子 Gamma = 有源区功率占比之和 sum(P_layer(act))
    % arrayfun(@(k) 表达式, modes) 把 modes 的每个元素依次取出来，赋给变量 k，
    % 计算"表达式”，把每个结果按顺序拼成一个数组返回
    act = find([L.is_active]);    
    if ~isempty(act)
        dG = max(abs(arrayfun(@(k) sum(k.P_layer(act)) - k.Gamma, modes)));        
        fprintf('所有模式中，Γ 与有源层功率占比之和的最大偏差 = %.3e  (应 = 0)\n', dG);
        fprintf('各层功率占比之和的最大偏差 = %.3e  (应 = 0)\n', ...
                max(abs(arrayfun(@(k) sum(k.P_layer) - 1, modes))));
    end

    fprintf('\n── 说明 ──\n');
    fprintf('可求解的模式个数: %d \n', find_date.ModeNum);

    %% 绘图（制表）
    % 表一： 纵轴：层； 横轴：厚度，折射率，各个模式的功率占比，损耗贡献
    % 图一： 纵轴：层厚度：横轴：折射率变化，各模式的模场归一化分布，附限制因子数值
    % 绘图函数并没对TE、TM进行区分，默认为TE，有需要可以后期修改。
    % 横坐标区间的标定？。？

    % 根据 do_save 决定图窗是否显示
    if do_save, vis = 'off'; else, vis = 'on'; end
    fig = figure('Visible', vis, 'Position', [60 40 1040 760], ...
                 'Name', 'GaN LD slab waveguide', 'Color', 'w');
    
    ax1 = subplot(2,1,1);
    draw_stack(ax1, modes(1).y, real(modes(1).n_y), modes, yb, L, ...
               [-0.6 2.6], 'northwest');
    title(ax1, sprintf('Whole epitaxial stack   (\\lambda = %g nm, %s modes)', ...
                       lambda*1000, pol));
    
    ax2 = subplot(2,1,2);
    draw_stack(ax2, modes(1).y, real(modes(1).n_y), modes, yb, L, ...
               [1.10 1.62], 'northeast');
    title(ax2, 'Zoom on the core:  n-wave / QW / p-wave / EBL');

    %% 保存 
    if do_save
        % 在运行函数的文件目录里生成一个out文件夹，用于保存数据
        [pathstr, name, ~]=fileparts(mfilename('fullpath'));
        outdir = fullfile(pathstr, [name, '_', savedir]);
        if ~exist(outdir, 'dir')
            mkdir(outdir); 
        end
        % 图形导出函数(图形窗口句柄，输出路径，分辨率)
        exportgraphics(fig, fullfile(outdir, 'ebl_modes.png'), 'Resolution', 160);
        close(fig);
        fprintf('\n图片已保存，地址：%s \n\n',outdir);
    end
end

%% 函数

% ───────────────────────────────────────────────────────────────────────────────────
% 测试结构：衬底 / …… / 包层
% 首尾两层厚度必须是 Inf（半无限）
% ───────────────────────────────────────────────────────────────────────────────────
function layers = test_layer()
    layers = struct( ...
        'd', ...
            { Inf,     1.5,     Inf      }, ...
        'n', ...
            { 1.770,       2.550,      2.374    }, ...
        'tag', ...
            {'substrate',   'QW',      'clad' }, ...
        'is_active', ...
            { false,     true,      false    });
end

% ───────────────────────────────────────────────────────────────────────────────────
% 7 层结构：衬底 / n包层 / n波导 / QW / p波导 / EBL / p包层
% 首尾两层厚度必须是 Inf（半无限）
% ───────────────────────────────────────────────────────────────────────────────────
function layers = stack_7layer()
    layers = struct( ...
        'd', ...
            { Inf,      1.200,      0.120,      0.009,      0.080,      0.020,      Inf      }, ...
        'n', ...
            { 1.770,    2.404,      2.455,      2.550,      2.455,      2.343,      2.374    }, ...
        'tag', ...
            {'substrate','n-clad',  'n-wave',   'QW',       'p-wave',   'EBL',      'p-clad' }, ...
        'is_active', ...
            { false,    false,      false,      true,       false,      false,      false    });
end


% ───────────────────────────────────────────────────────────────────────────────────
% 9 层结构：把单层等效 QW 展开成 阱-垒-阱
% Γ = (阱1 + 阱2 的功率) / 总功率
% ───────────────────────────────────────────────────────────────────────────────────
function layers = stack_9layer_mqw()
    layers = struct( ...
        'd', ...
            { Inf,   1.200,   0.120,   0.003,   0.006,   0.003,   0.080,   0.020,   Inf    }, ...
        'n', ...
            { 1.770, 2.404,   2.455,   2.550,   2.430,   2.550,   2.455,   2.343,   2.374  }, ...
        'tag', ...
            {'substrate','n-clad','n-wave','QW1','barrier','QW2','p-wave','EBL','p-clad'}, ...
        'is_active', ...
            { false, false,   false,   true,    false,   true,    false,   false,   false  });
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：层结构校验、规范化
% 输入：层结构体
% 输出：L 与输入同构，但 n 强制转为 complex，缺省字段已补齐；任何不合法输入在这里报错
% 说明：field-字段；任一字段与数值都要有
% ───────────────────────────────────────────────────────────────────────────────────
function L = build_layers(layers)
    assert(isstruct(layers) && ~isempty(layers), 'build_layers: layers 必须是 struct 数组');
    N = numel(layers);
    assert(N >= 3, ...
        'build_layers: 至少需要 3 层（两侧半无限包层 + 至少一个有限层），当前 %d 层', N);
    %layers为至少3层的结构体数组
    %assert(条件A，输出'B'),如果A不满足，输出B
    
    L = layers;
    for k = 1:N
        assert(isfield(L(k), 'd') && ~isempty(L(k).d), 'build_layers: 第 %d 层缺少字段 d', k);
        assert(isfield(L(k), 'n') && ~isempty(L(k).n), 'build_layers: 第 %d 层缺少字段 n', k);
    
        if ~isfield(L(k), 'tag') || isempty(L(k).tag)
            L(k).tag = sprintf('layer%d', k);
        end
        if ~isfield(L(k), 'is_active') || isempty(L(k).is_active)
            L(k).is_active = false;
        end
    
        % 复折射率定义+复折射率数值判断
        L(k).n = complex(L(k).n);            
        assert(real(L(k).n) > 0, ...
            'build_layers: 第 %d 层 (%s) 折射率实部必须为正，当前 %g', ...
            k, L(k).tag, real(L(k).n));
        assert(imag(L(k).n) >= 0, ...
            ['build_layers: 第 %d 层 (%s) 折射率虚部为负 (%g)。\n' ...
             '  按本项目约定 Im(n) > 0 表示【损耗】；若该层是增益介质，\n' ...
             '  请在阶段4 的增益模块中处理，不要在此处给负虚部。'], ...
            k, L(k).tag, imag(L(k).n));
    end
    
    % 厚度数值判读
    assert(isinf(L(1).d) && L(1).d > 0, 'build_layers: 第 1 层（衬底侧）厚度必须为 +Inf');
    assert(isinf(L(N).d) && L(N).d > 0, 'build_layers: 第 %d 层（覆盖层侧）厚度必须为 +Inf', N);
    if N > 2
        dmid = [L(2:N-1).d];
        assert(all(dmid > 0) && all(isfinite(dmid)), 'build_layers: 中间层厚度必须为有限正值');
    end
    
    act = [L.is_active];
    assert(~act(1) && ~act(N), 'build_layers: 半无限层不能标记为有源区');
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：y轴标定
% 输入：规范化后的层结构体
% 输出：二维数组。每层对应y值区间的上、下限。
% 说明：
%   依照各层厚度，明确各层的y区间。
%   L(2) 的下沿位于 y = 0，递增
% ───────────────────────────────────────────────────────────────────────────────────
function yb = layer_bounds(L)
    N = numel(L);
    yb = zeros(N, 2);
    yb(1, :) = [-Inf, 0];
    
    y = 0;
    for k = 2:N-1
        yb(k, :) = [y, y + L(k).d];
        y = y + L(k).d;
    end
    
    yb(N, :) = [y, Inf];
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：单层 2×2 传输矩阵（TE / TM 通用）
% 输入：横向波数 kz; 厚度 d; 折射率 n; 偏振 pol
% 输出：对应该层的传输矩阵    
% 说明：
%   TE：M=[cos(kh) , 1/ksin(kh); -ksin(kh) , cos(kh)]
%   TM: M=[cos(kh) , n^2/ksin(kh); -k/n^2 * sin(kh) , cos(kh)]
% ───────────────────────────────────────────────────────────────────────────────────
function [a11, a12, a21, a22] = layer_matrix(kz, d, n, pol)
    % 能否转移到校验函数？
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
% 函数：传输矩阵计算
% 输入：有效折射率 n_eff; 层数据结构体 L; 波长 lambda; 偏振 pol
% 输出：f：某一n_eff下（其他参量数值默认）的传输矩阵输出(输出为0即满足，是束缚模)。
%       aux：结构体，包含连乘后的传输矩阵、p0、p2(参见'说明'中两个公式的左侧)。
% 调用函数：[a11, a12, a21, a22] = layer_matrix(kz, d, n, pol)
% 说明：
%   TE：[-p0,1] M^TE [1;-p2]=0
%   TM: [-p0/n0^2] M^TM [1; -p2/n2^2]=0
%   n_eff其实是待解项，通过遍历所有n_eff，保留满足公式的n_eff
%   有限真实边界需重新考虑（个人觉得）
%   复折射率也是这个条件？？
% ───────────────────────────────────────────────────────────────────────────────────
%       TE : v = [E_y ; dE_y/dy]           界面连续：E_y 与 dE_y/dy
%       TM : v = [H_y ; (1/n²)·dH_y/dy]    界面连续：H_y 与 (1/n²)dH_y/dy
% ──────────────────────────────────────────────────────────────────
function [f, aux] = dispersion_tmm(n_eff, L, lambda, pol)
    if nargin < 4 || isempty(pol), pol = 'TE'; end
    is_TE = strcmpi(pol, 'TE');
    N = numel(L);
    
    in_size = size(n_eff);                
    x = n_eff(:).';                       % 1×M 行向量
    k0 = 2*pi/lambda;
    
    % 两侧半无限层的衰减常数（必须取衰减解，Re γ > 0）
    gam_sub = k0 * sqrt(x.^2 - L(1).n^2);
    gam_cov = k0 * sqrt(x.^2 - L(N).n^2);
    
    if is_TE
        s_sub = gam_sub;              s_cov = gam_cov;
    else
        s_sub = gam_sub / L(1).n^2;   s_cov = gam_cov / L(N).n^2;
    end
    
    % 中间层传输矩阵连乘（从下往上，逐次左乘）
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
% 输入：函数 dispersion_tmm; n_lo 底层折射率，n_hi 顶层折射率; 
%           n_roots 求解模式数; opts 求解要用参数
% 输出：roots_n：规定模式数下的 n_eff 解。
%      info：结构体，记录本次求解的模数, 对应扫描点数, 最大求解模式数
% 调用函数：scan_bisect，scan_newton
% 说明：
%   输入里有函数，是因为该函数需额外输出参量
%   本函数考虑了实折射率和复折射率
% ───────────────────────────────────────────────────────────────────────────────────
function [roots_n, info] = find_roots_complex(fun, n_lo, n_hi, n_roots, opts)
    % 输出参数补齐，及求解项数值判定
    if nargin < 4 || isempty(n_roots), n_roots = 5; end
    if nargin < 5 || isempty(opts), opts = struct(); end
    if ~isfield(opts, 'n_scan'), opts.n_scan = 20000; end %扫描点数
    if ~isfield(opts, 'tol'),    opts.tol    = 1e-13; end %求根收敛容差
    if ~isfield(opts, 'maxit'),  opts.maxit  = 200;   end %最大迭代次数
    if ~isfield(opts, 'merge'),  opts.merge  = 1e-6;  end %去重时两个根的最小相对间距
    
    assert(n_hi > n_lo, 'find_roots_complex: 要求 n_hi > n_lo');
    
    xs = linspace(n_lo, n_hi, opts.n_scan);               %均匀扫描
    fs = fun(xs);
    fs = fs(:).';
    
    % 防御性编程-不存在所有点的输出都是0
    fscale = max(abs(fs));
    if fscale == 0, fscale = 1; end
    
    % 折射率 实/复判断，允许有噪声
    is_real = all(abs(imag(fs)) <= 1e-10 * fscale);
    
    if is_real
        roots_n = scan_bisect(fun, xs, real(fs), opts);   %二分法
        info.mode = 'real';
    else
        roots_n = scan_newton(fun, xs, fs, opts);         %牛顿法
        info.mode = 'complex';
    end
    
    % 去重、降序、截断
    if ~isempty(roots_n)
        roots_n = sort(roots_n(:).', 'descend');
        % 第一个根恒保留(true)，后面每个根与前一根的差大于阈值才保留。
        % keep(i)=true/false。对数组进行标记
        % max(1, abs(roots_n(1:end-1))) 得到的是数组；将roots_n每一个和1比较取最大值。
        keep = [true, abs(diff(roots_n)) > opts.merge * max(1, abs(roots_n(1:end-1)))];
        roots_n = roots_n(keep);
        info.ModeNum = numel(roots_n);
        if numel(roots_n) > n_roots
            roots_n = roots_n(1:n_roots);
        end
    end
    
    info.n_found = numel(roots_n);
    info.scanned = opts.n_scan;
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：二分法求解
% 输入：函数 dispersion_tmm; 所有待测的n_eff xs(向量); 
%       待测n_eff的转移矩阵输出 fr(向量); 求解要用参数 opts
% 输出：满足条件的n_eff
% 说明：
%   二分只依赖"闭区间两端符号相反"这个充分条件，不要求函数单调或光滑
%   用于实折射率
% ───────────────────────────────────────────────────────────────────────────────────
function roots_n = scan_bisect(fun, xs, fr, opts)
    roots_n = [];
    s = sign(fr);                        % 取每个fr的符号：>0→取1,<0→取-1,=0取0
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
% 函数：牛顿法求解
% 输入：函数 dispersion_tmm; 所有待测的n_eff xs(向量); 
%       待测n_eff的转移矩阵输出 fr(向量); 求解要用参数 opts
% 输出：满足条件的n_eff
% 调用函数：newton_complex
% 说明：
%   用于复折射率
% ───────────────────────────────────────────────────────────────────────────────────
function roots_n = scan_newton(fun, xs, fs, opts)
    %SCAN_NEWTON  复值分支：定位 |f| 的局部极小，再用复数牛顿细化
    mag = abs(fs);
    cand = [];
    for j = 2:numel(xs)-1
        if mag(j) < mag(j-1) && mag(j) <= mag(j+1)
            cand(end+1) = xs(j);        %#ok<AGROW>
        end
    end
    
    roots_n = [];
    for x0 = cand
        r = newton_complex(fun, x0, opts);
        if ~isempty(r)
            roots_n(end+1) = r;         %#ok<AGROW>
        end
    end
end


% ──────────────────────────────────────────────────────────────────
function r = newton_complex(fun, x0, opts)
    %NEWTON_COMPLEX  复数牛顿迭代，带步长阻尼。无法收敛时返回空。
    h      = 1e-7;                       % 数值导数步长（对 n_eff ~ O(1) 合适）
    r      = x0;
    n_iter = opts.maxit;
    
    for it = 1:n_iter
        fx = fun(r);
        f1 = fun(r + h);
        f2 = fun(r - h);
        df = (f1 - f2) / (2*h);
    
        if abs(df) < 1e-300
            r = [];
            return;
        end
    
        step = fx / df;
        if abs(step) > 0.5                 % 阻尼：单步不超过 0.5
            step = 0.5 * step / abs(step);
        end
        r = r - step;
    
        if abs(step) < 1e-14 * max(1, abs(r))
            break;
        end
        if it == n_iter
            r = [];                        % 未收敛，丢弃
            return;
        end
    end
    
    % 后验检查：|f| 必须显著小于邻域，剔除鞍点/极小值平台造成的伪根
    fr  = abs(fun(r));
    eps_ = 1e-3 * max(1, abs(r));
    nb  = max(abs(fun(r + eps_)), abs(fun(r - eps_)));
    if ~(fr < 0.5 * nb)
        r = [];
    end
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：模场剖面 φ(y) 重建
% 输入：有效折射率 n_eff; 层数据结构体 L; 波长 lambda; 
%       偏振 pol; 求解要用参数 opts
% 输出：对应n_eff的模场分布[横轴y，纵轴E(模场)、n(折射率)，每层的y轴上下界]
% 调用函数：tail_depth
% 说明：
%   从底部半无限层的边界状态出发，用 TMM 逐层向上传输
%   输出 E 未归一化；TE 为 E_y，TM 为 H_y。
%   层-厚度y-模场E-折射率n 一一对应。一个元胞对应一层，一层对应含有npts个结点的厚度，每
%                                 个厚度的值(E or n)不一样的。
% ───────────────────────────────────────────────────────────────────────────────────
function [y, E, n_y, yb] = mode_profile(n_eff, L, lambda, pol, opts)
    if nargin < 5 || isempty(opts), opts = struct(); end
    if ~isfield(opts, 'npts'),     opts.npts     = 800; end  % 模场重建时每层采样点数
    if ~isfield(opts, 'tail_dec'), opts.tail_dec = 20;  end  % 半无限层向外延伸的衰减长度
    if ~isfield(opts, 'tail_max'), opts.tail_max = 50;  end  % 延伸长度上限 [μm]
    
    k0 = 2*pi/lambda;
    N  = numel(L);
    yb = layer_bounds(L);       % yb 还要输出吗？
    d_total = yb(N, 1);         % 总厚度（不考虑第1、N层）
    
    gam_sub = k0 * sqrt(complex(n_eff^2 - L(1).n^2));
    gam_cov = k0 * sqrt(complex(n_eff^2 - L(N).n^2));
    
    depth_sub = tail_depth(gam_sub, opts);
    depth_cov = tail_depth(gam_cov, opts);
    
    is_TE = strcmpi(pol, 'TE');
    
    y_cell = cell(1, N);
    E_cell = cell(1, N);
    n_cell = cell(1, N);
    
    % ── 第 1 层：底部半无限层 ─────────────────────────────────────
    yy = linspace(-depth_sub, 0, opts.npts);
    y_cell{1} = yy(1:end-1);                       % 去掉 y=0，避免与下一层重复
    E_cell{1} = exp(gam_sub * y_cell{1});
    n_cell{1} = L(1).n * ones(size(y_cell{1}));
    
    % 在 y = 0 处的状态向量（与 dispersion_tmm 的底部边界条件一致）
    if is_TE
        v = [1; gam_sub];
    else
        v = [1; gam_sub / L(1).n^2];
    end
    cur_y = 0;
    
    % ── 第 2 … N-1 层：和传输矩阵一致 ──────────────────────────
    for k = 2:N-1
        dk = L(k).d;
        nk = L(k).n;
        kz = k0 * sqrt(complex(nk^2 - n_eff^2));
    
        yy = linspace(cur_y, cur_y + dk, opts.npts);
        dy = yy - cur_y;        % 有什么意义？
        
        %取极限
        if abs(kz) < 1e-14
            S = dy;
        else
            S = sin(kz*dy) / kz;
        end
        C  = cos(kz*dy);
        SN = sin(kz*dy);
    
        if is_TE
            a12 = S;                 a21 = -kz * SN;
        else
            a12 = (nk^2) * S;        a21 = -(kz / (nk^2)) * SN;
        end
        a11 = C;  a22 = C;
    
        vv = a11 * v(1) + a12 * v(2);
    
        % 该层上沿的新状态（供下一层使用）
        v = [a11(end)*v(1) + a12(end)*v(2); ...
             a21(end)*v(1) + a22(end)*v(2)];
    
        if k < N-1
            y_cell{k} = yy(1:end-1);
            E_cell{k} = vv(1:end-1);
        else
            y_cell{k} = yy;
            E_cell{k} = vv;
        end
        n_cell{k} = nk * ones(size(y_cell{k}));
        cur_y = cur_y + dk;
    end
    
    % ── 第 N 层：顶部半无限层 ────────────────────────────────────
    yy = linspace(d_total, d_total + depth_cov, opts.npts);
    yy = yy(2:end);                                % 去掉 y=d_total，避免重复
    y_cell{N} = yy;
    E_cell{N} = exp(-gam_cov * (yy - d_total)) * v(1);
    n_cell{N} = L(N).n * ones(size(yy));
    
    % ── 拼接 ────────────────────────────────────────────────────
    y   = [y_cell{:}];
    E   = [E_cell{:}];
    n_y = [n_cell{:}];
    
    % 折射率全实且 n_eff 为实时，清掉运算中产生的 0 虚部    如果折射率全实，哪里来的虚部？
    if imag(n_eff) == 0 && all(imag([L.n]) == 0)
        E   = real(E);
        n_y = real(n_y);
    end
    
    assert(numel(y) == numel(E) && numel(y) == numel(n_y), ...
        'mode_profile: 拼接后长度不一致');
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：模场衰减厚度
% 输入：横向衰减常数 gam; 求解要用参数 opts
% 输出：衰减厚度
% 说明：
%   默认 20 个衰减长度
%   延伸不足会让 Γ 系统性偏小
% ───────────────────────────────────────────────────────────────────────────────────
function dep = tail_depth(gam, opts)
    g = real(gam);
    if g > 1e-9
        dep = min(opts.tail_max, opts.tail_dec / g);
        dep = max(dep, 0.1);
    else
        dep = opts.tail_max;
    end
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：模场归一化
% 输入：坐标系 E-y,n-y; 偏振 pol; 模式序数 ord
% 输出：E：归一化后的模场分布; norm_val：自检数
% 说明：
%   E(y,z,t) = φ(y)·exp(i(βz − ωt))
%   TE:  S_z ∝ β·|E_y|²
%   TM:  S_z ∝ (β/n²)·|H_y|²
%   norm_val 用于归一化后复算，应恒为 1
% ───────────────────────────────────────────────────────────────────────────────────
function [E, norm_val] = normalize_mode(y, E, n_y, pol)
    if strcmpi(pol, 'TE')
        w = abs(E).^2;
    else
        w = abs(E).^2 ./ max(real(n_y).^2, eps);  % eps 两个相邻浮点数的差 2.22e-16
    end
    
    norm_before = trapz(y, w);                    % 梯形法数值积分
    assert(isfinite(norm_before) && norm_before > 0, ...
        'normalize_mode: 归一化积分非正或非有限 (%g)，模场重建有问题', norm_before);
    
    E = E / sqrt(norm_before);
    
    if strcmpi(pol, 'TE')
        w2 = abs(E).^2;
    else
        w2 = abs(E).^2 ./ max(real(n_y).^2, eps);
    end
    norm_val = trapz(y, w2);                       % 用于自检
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：限制因子
% 输入：坐标系 E-y,E-n; 厚度 yb; 层数据结构体 L; 偏振 pol
% 输出：G：限制因子; 
%       detail：结构体，含有源区总功率 P_act，Ld总功率 P_tot，各层功率(向量) P_layer
% 说明：
%   口径需要统一。
%   本函数：全部 is_active 层的并集; 无穷积分
%   MQWs 是全部阱层的并集。考虑势垒。
%   分母范围需要明确：无穷积分，还是截断到外延结构。对漏模只有无穷积分是自洽的。
%   此函数可以考虑和模场归一化合并。
% ───────────────────────────────────────────────────────────────────────────────────
function [G, detail] = confinement_factor(y, E, n_y, yb, L, pol)
    if strcmpi(pol, 'TE')
        w = abs(E).^2;
    else
        w = abs(E).^2 ./ max(real(n_y).^2, eps);
    end
    
    P_tot = trapz(y, w);
    assert(P_tot > 0, 'confinement_factor: 总功率积分非正');
    
    % ── 各层功率占比（与 is_active 无关）─────────────────────────
    N = numel(L);
    P_layer = zeros(1, N);
    % 将顶层、底层区分出来；获取各层的区间
    for k = 1:N
        if isinf(yb(k,1))           % 底层
            mk = (y <= yb(k,2));
        elseif isinf(yb(k,2))       % 顶层
            mk = (y >= yb(k,1));
        else
            mk = (y >= yb(k,1) & y <= yb(k,2));
        end
        if any(mk)
            P_layer(k) = trapz(y(mk), w(mk)) / P_tot;
        end
    end
    
    % 判断是否有有源区，无返回NaN
    % 不返回 0，避免误认为"模式与有源区完全不重叠"
    act = [L.is_active];
    if ~any(act)
        G = NaN;
        detail = struct('P_act', NaN, 'P_tot', P_tot, 'P_layer', P_layer);
        return;
    end
    
    % 逐层单独积分再求和     % 这个和上面是不是部分重复了，可以合并的；或者获取索引再求和？
    P_act = 0;
    for k = 1:N
        if ~L(k).is_active, continue; end
        if isinf(yb(k,1))
            mk = (y <= yb(k,2));
        elseif isinf(yb(k,2))
            mk = (y >= yb(k,1));
        else
            mk = (y >= yb(k,1) & y <= yb(k,2));
        end
        if any(mk)
            P_act = P_act + trapz(y(mk), w(mk));
        end
    end
    G = P_act / P_tot;
    
    detail = struct('P_act', P_act, 'P_tot', P_tot, 'P_layer', P_layer);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：绘图、制表
% 输入：坐标系 ax; 快轴坐标 y; 折射率(对应快轴坐标) n_y; 模式数据结构体 m;
%       各层上下限 yb; 层数据结构体 L; 横轴显示范围 xr; 图例位置(字符串) legloc
% 输出：可视化图片，数据表格
% 说明：
%   在一个坐标轴上同时画折射率分布(左轴)与各阶模场(右轴)
%   纵轴范围由当前窗口内的数据算出。
%   层底色必须画在左轴，与折射率线同轴。yyaxis 的两根轴是分层渲染（右轴整体盖在左轴之上）;
%   若把底色画在右轴，它会把折射率线整条盖成灰色。
% ───────────────────────────────────────────────────────────────────────────────────
function draw_stack(ax, y, n_y, m, yb, L, xr, legloc)
    % 绘图
    %     |φ|² 归一化后峰值 = 1/(有效宽度)，对本结构约 3.4，不是 1。
    %
    %   ⚠ 层底色必须画在【左轴】，与折射率线同轴。yyaxis 的两根轴是分层
    %     渲染的（右轴整体盖在左轴之上，与绘制先后无关）；若把底色画在右轴，
    %     它会把折射率线整条盖成灰色。
    
    axes(ax); 
    cla(ax);                                % 清空该轴内容
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    
    % 选择对应横坐标范围 xr 的数据区间
    sel  = (y >= xr(1)) & (y <= xr(2));
    nsel = n_y(sel);
    
    % 绘图-折射率的步进值
    n_min = min(nsel);
    n_max = max(nsel);
    pad_n = 0.08*(n_max - n_min) + 0.01;    
    
    % 模场-纵轴的区间范围
    pk = 0;
    for k = 1:numel(m)                      % numel(m) 总模式数
        w = abs(m(k).E).^2;                 
        pk = max(pk, max(w(sel)));
    end
    ytop = pk * 1.10;                       % 右轴上限：留 10% 余量给图例和峰
    

    % 左轴：层底色 + 折射率阶梯 
    yyaxis(ax, 'left');
    ax.YColor = 'k';
    ylim(ax, [n_min - pad_n, n_max + pad_n]);
    yl = ylim(ax);
    
    % 高亮 QW，画折射率阶梯
    for k = 1:numel(L)
        lo = max(yb(k,1), xr(1));
        hi = min(yb(k,2), xr(2));
        if hi <= lo, continue; end
        if L(k).is_active
            patch(ax, [lo hi hi lo], [yl(1) yl(1) yl(2) yl(2)], [1 0.85 0.35], ...
                  'FaceAlpha', 0.50, 'EdgeColor','none', 'HandleVisibility','off');
        else
            patch(ax, [lo hi hi lo], [yl(1) yl(1) yl(2) yl(2)], [0.92 0.92 0.92], ...
                  'FaceAlpha', 0.50, 'EdgeColor','none', 'HandleVisibility','off');
        end
    end
    plot(ax, y, n_y, 'k-', 'LineWidth', 1.6, 'HandleVisibility','off');
    ylabel(ax, 'refractive index n');
    

    % 右轴：各阶模 |φ|² + 层边界 + 层标注
    yyaxis(ax, 'right');
    ylim(ax, [0 ytop]);
    
    % 模场
    cols = lines(numel(m));                 % 生成一组区分度好的颜色
    % is_TE = strcmpi(pol, 'TE');
    % if is_TE
        for k = 1:numel(m)
            plot(ax, y, abs(m(k).E).^2, '-', 'LineWidth', 1.7,  ...
                 'Color', cols(k,:),'DisplayName',  ...
                 sprintf('TE%d  (n_{eff}=%.5f, \\Gamma=%.4f)',k-1, real(m(k).n_eff), m(k).Gamma));
        end
        ylabel(ax, '|\phi(y)|^2');
    % else
    % end
    
    % 层边界
    for k = 2:numel(L)
        xb = yb(k,1);
        if xb > xr(1) && xb < xr(2)
            plot(ax, [xb xb], [0 ytop], ':', 'Color',[0.55 0.55 0.55], ...
                 'HandleVisibility','off');
        end
    end

    % 层标注
    for k = 1:numel(L)
        lo = max(yb(k,1), xr(1));
        hi = min(yb(k,2), xr(2));
        if hi <= lo, continue; end
        text(ax, 0.5*(lo+hi), 0.015*ytop, L(k).tag, 'Rotation', 90, 'FontSize', 7.5, ...
             'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
             'Color', [0.25 0.25 0.25]);
    end
    
    % xlim 必须在所有 yyaxis 切换之后设，否则会被重置
    xlim(ax, xr);
    legend(ax, 'Location', legloc, 'FontSize', 8);  % 生成图例
    hold(ax, 'off');
end
