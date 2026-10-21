function modes = Model_Self_L1_7L_EBL_test(action)
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
    % !! 验证分支：action 含 'verify' 时只跑验证电池，不画图、不求解主结构
    % !! 判据用与 TMM 无关的参照物（有限差分 / 解析行列式 / FFT），见文件末尾 run_verify
    if contains(lower(string(action)), 'verify')
        run_verify();
        return;
    end
    
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
        'num_scan',   20000, ...   % 色散函数扫描点数（根搜索用）
        'tol',      1e-13, ...   % 求根收敛容差（相对）
        'maxit',    200,   ...   % 最大迭代次数
        'merge',    1e-6,  ...   % 去重时两个根的最小相对间距
        'npts',     800,   ...   % 模场重建时每层采样点数
        'tail_dec', 20,    ...   % 半无限层向外延伸多少个衰减长度
        'tail_max', 50,    ...   % 延伸长度上限 [μm]
        'num_ff',     2001);       % 远场角分布采样点数
    
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
        ff            = far_field(y, E, k0, opts);
    
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
        modes(j).ff       = ff;              % 快轴远场（结构体—角轴/分布/FWHM/峰值角）
    end
    
    %% 输出 删除原命令行输出，调整为表格输出；保留自检。后续可添加 UI。        

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

    %% 绘图
    % 画一个三维的脊形波导示意图 （→ UI ?）
    % 画一个面向光束的x-y的模场分布

    % 图一： 纵轴：层厚度：横轴：折射率变化，各模式的模场归一化分布，附限制因子数值
    % 横坐标区间的标定
    pad = 0.6;                                    % 全图两侧留白 [μm]
    xr_full = [yb(2,1) - pad, yb(N,1) + pad];

    iact = find([L.is_active]);                   % 有源区层号
    pad_core = 0.3;                               % 芯区两侧留白 [μm]
    if isempty(iact)
        xr_zoom = xr_full;                        % 无有源区则退回全图
    else
        xr_zoom = [min(yb(iact,1)) - pad_core, ...
                   max(yb(iact,2)) + pad_core];
    end
    xr_full = [-0.6 2.6];
    xr_zoom = [1.10 1.62];

    % 根据 do_save 决定图窗是否显示
    if do_save, vis = 'off'; else, vis = 'on'; end
    fig1 = figure('Visible', vis, 'Position', [60 40 1040 760], ...
                 'Name', 'GaN LD slab waveguide', 'Color', 'w');
    
    ax1 = subplot(2,1,1);
    draw_stack(ax1, modes(1).y, real(modes(1).n_y), modes, yb, L, ...
               xr_full, 'northwest', pol);
    title(ax1, sprintf('Whole epitaxial stack   (\\lambda = %g nm, %s modes)', ...
                       lambda*1000, pol));
    
    % 有源区部分局部放大图
    ax2 = subplot(2,1,2);
    draw_stack(ax2, modes(1).y, real(modes(1).n_y), modes, yb, L, ...
               xr_zoom, 'northeast', pol);
    title(ax2, 'Zoom on the core:  n-wave / QW / p-wave / EBL');



    % 图二：快轴远场角分布
    fig2 = figure('Visible', vis, 'Position', [120 80 900 420], ...
                  'Name', 'GaN LD far field (fast axis)', 'Color', 'w');
    axf = axes(fig2);
    hold(axf, 'on'); grid(axf, 'on'); box(axf, 'on');

    thMax = 0;
    cols  = lines(numel(modes));
    for k = 1:numel(modes)
        fk = modes(k).ff;
        plot(axf, fk.theta*180/pi, fk.I, '-', 'LineWidth', 1.6, 'Color', cols(k,:), ...
             'DisplayName', sprintf('%s%d  (FWHM = %.2f°, peak = %.2f°)', ...
                                    pol, k-1, fk.fwhm*180/pi, fk.th_pk*180/pi));
        thMax = max(thMax, abs(fk.th_pk*180/pi) + fk.fwhm*180/pi);
    end
    thMax = min(90, max(15, thMax*1.2));

    % 半高线：FWHM 的判据线。画出来一眼能看出那个数是"量到哪一段" ——
    % 峰值角偏离 0 的模（多瓣），FWHM 只量到峰值那一瓣。
    yline(axf, 0.5, ':', 'half', 'Color', [0.4 0.4 0.4], ...
          'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');

    xlim(axf, [-thMax, thMax]);
    xlabel(axf, '\theta  [deg]   (fast axis, relative to z)');
    ylabel(axf, '|\Phi(\theta)|^2   (normalized to peak)');
    title(axf, sprintf('Fast-axis far field   (\\lambda = %g nm, %s)', lambda*1000, pol));
    legend(axf, 'Location', 'northeast', 'FontSize', 8);
    hold(axf, 'off');

    %% 制表
    % uitable 是 UI 控件，exportgraphics 导不出它 —— 'save' 模式只导图。
    % 两个关键点：
    %   1.ColumnWidth 默认 'auto' 只按内容算宽，控件拉多大都不会铺满 —— 必须显式设列宽；
    %   2.传统 figure 的 uitable 没有 RowHeight 属性、行高不可控，
    %      所以高度按行数贴合内容（宁贴紧、不留白），窗口高度也随内容。

   
    scr  = get(groot, 'ScreenSize');                   % [1 1 屏宽 屏高]，像素
    w3   = min(1500, scr(3) - 60);                     % 窗口宽：上限 1500 px，且不出屏
    rowH = 32;                                         % 单行高 [px]（FontSize 16 约 28，留余量）
    gap  = 26;                                         % 两块表的间隙 [px]
    mar  = 14;                                         % 窗口内边距 [px]
    hT1  = (numel(modes) + 1) * rowH + 10;             % 表一高：数据行 + 表头
    hT2  = (N + 1)           * rowH + 10;              % 表二高：数据行 + 表头

    fig3 = figure('Visible', vis, 'Name', 'GaN LD - data tables', 'Color', 'w');
    fig3.InnerPosition = [30 60 w3, hT1 + hT2 + gap + 2*mar];    % 直接定客户区，免估标题栏

    % 表一：纵轴：模式;  横轴：n_eff, beta, Γ, theta_FWHM, 峰值角, 归一化自检数
    mCols = {'mode', 'n_eff', 'beta [1/um]', 'Gamma', ...
             'FWHM [deg]', 'peak [deg]', 'norm check'};
    mData = cell(numel(modes), numel(mCols));
    for k = 1:numel(modes)
        mData(k, :) = { sprintf('%s%d', pol, k-1), ...
                        sprintf('%.8f', real(modes(k).n_eff)), ...
                        sprintf('%.6f', real(modes(k).beta)), ...
                        sprintf('%.6f', modes(k).Gamma), ...
                        sprintf('%.2f',   modes(k).ff.fwhm*180/pi), ...
                        sprintf('%.2f',   modes(k).ff.th_pk*180/pi), ...
                        sprintf('%.3e',   modes(k).norm_val) };
    end
    t1 = uitable(fig3, 'Data', mData, 'ColumnName', mCols, 'FontSize', 16, ...
                 'Units', 'pixels');


    % 表二： 纵轴：层;  横轴：厚度，折射率，各个模式的功率占比，损耗贡献
    lCols = [{'#', 'tag', 'd [um]', 'n'}, ...
             arrayfun(@(j) sprintf('P(%s%d) [%%]', pol, j-1), ...
                      1:numel(modes), 'UniformOutput', false)];
    lData = cell(N, numel(lCols));
    for k = 1:N
        dstr = 'Inf';
        if isfinite(L(k).d), dstr = sprintf('%.4f', L(k).d); end
        lData(k, 1:4) = { sprintf('%d', k), L(k).tag, dstr, sprintf('%.4f', real(L(k).n)) };
        for j = 1:numel(modes)
            lData{k, 4+j} = sprintf('%.4f', 100*modes(j).P_layer(k));
        end
    end
    t2 = uitable(fig3, 'Data', lData, 'ColumnName', lCols, 'FontSize', 16, ...
                 'Units', 'pixels');

    % 首次布局：列宽按窗口宽等分铺满、两块表上下排好；拉伸窗口时自动重排
    relayout_tables(fig3, t1, t2, mCols, lCols);
    fig3.SizeChangedFcn = @(s, ~) relayout_tables(s, t1, t2, mCols, lCols);


    %% 保存 
    if do_save
        % 在运行函数的文件目录里生成一个out文件夹，用于保存数据
        [pathstr, name, ~]=fileparts(mfilename('fullpath'));
        outdir = fullfile(pathstr, [name, '_', savedir]);
        if ~exist(outdir, 'dir')
            mkdir(outdir); 
        end

        % 图形导出函数(图形窗口句柄，输出路径，分辨率)
        % 表格 figure 是 uitable，导不出也不关
        exportgraphics(fig1,  fullfile(outdir, 'ebl_modes.png'),  'Resolution', 160);
        exportgraphics(fig2, fullfile(outdir, 'fast_field.png'), 'Resolution', 160);
        close([fig1, fig2]);        
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
    %layers为至少3层的结构体数组
    %assert(条件A，输出'B'),如果A不满足，输出B
    assert(isstruct(layers) && ~isempty(layers), 'build_layers: layers 必须是 struct 数组');
    N = numel(layers);
    assert(N >= 3, ...
        'build_layers: 至少需要 3 层（两侧半无限包层 + 至少一个有限层），当前 %d 层', N);
   
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
%   1.依照各层厚度，明确各层的y区间。
%   2.L(2) 的下沿位于 y = 0，递增
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
%       aux：结构体 (参见'说明'中两个公式的左侧)。
%           aux.gam_sub  p0
%           aux.gam_cov  p2
%           aux.M        连乘后的传输矩阵 M
% 调用函数：[a11, a12, a21, a22] = layer_matrix(kz, d, n, pol)
% 说明：
%   1.TE：[-p0,1] M^TE [1;-p2]=0
%     TM: [-p0/n0^2] M^TM [1; -p2/n2^2]=0
%   2.n_eff其实是待解项，通过遍历所有n_eff，保留满足公式的n_eff
%   3.个人：有限真实边界需重新考虑；复折射率也是这个条件？
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
%      info：结构体
%           info.mode       折射率性质（实/复）
%           info.ModeNum    可求解的最大模数
%           info.n_found    本次求解的模数
%           info.scanned    求解中采用的扫描点数
% 调用函数：scan_bisect，scan_newton
% 说明：
%   1.输入里有函数，是因为该函数需额外输出参量
%   2.本函数考虑了实折射率和复折射率
% ───────────────────────────────────────────────────────────────────────────────────
function [roots_n, info] = find_roots_complex(fun, n_lo, n_hi, n_roots, opts)
    % 输出参数补齐，及求解项数值判定
    if nargin < 4 || isempty(n_roots), n_roots = 5; end
    if nargin < 5 || isempty(opts), opts = struct(); end
    if ~isfield(opts, 'num_scan'), opts.num_scan = 20000; end   %扫描点数
    if ~isfield(opts, 'tol'),    opts.tol    = 1e-13; end       %求根收敛容差
    if ~isfield(opts, 'maxit'),  opts.maxit  = 200;   end       %最大迭代次数
    if ~isfield(opts, 'merge'),  opts.merge  = 1e-6;  end       %去重时两个根的最小相对间距
    
    assert(n_hi > n_lo, 'find_roots_complex: 要求 n_hi > n_lo');
    
    xs = linspace(n_lo, n_hi, opts.num_scan);                   %均匀扫描
    fs = fun(xs);
    fs = fs(:).';
    
    % 防御性编程-不存在所有点的输出都是0
    fscale = max(abs(fs));
    if fscale == 0, fscale = 1; end
    
    % 折射率 实/复判断，允许有噪声
    is_real = all(abs(imag(fs)) <= 1e-10 * fscale);
    
    if is_real
        roots_n = scan_bisect(fun, xs, real(fs), opts);         %二分法
        info.mode = 'real';
    else
        roots_n = scan_newton(fun, xs, fs, opts);               %牛顿法
        info.mode = 'complex';
    end
    
    % 去重、降序、（第一次）截断
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
    info.scanned = opts.num_scan;
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：二分法求解
% 输入：函数 dispersion_tmm; 所有待测的n_eff xs(向量); 
%       待测n_eff的转移矩阵输出 fr(向量); 求解要用参数 opts
% 输出：满足条件的n_eff
% 说明：
%   1.二分只依赖"闭区间两端符号相反"这个充分条件，不要求函数单调或光滑
%   2.本函数用于实折射率
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
%   本函数用于复折射率
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
%   1.从底部半无限层的边界状态出发，用 TMM 逐层向上传输
%   2.输出 E 未归一化；TE 为 E_y，TM 为 H_y。
%   3.层-厚度y-模场E-折射率n 一一对应。一个元胞对应一层，一层对应含有npts个结点的厚度，每
%                                 个厚度的值(E or n)不一样的。
%   4.个人：yb还要输出吗
% ───────────────────────────────────────────────────────────────────────────────────
function [y, E, n_y, yb] = mode_profile(n_eff, L, lambda, pol, opts)
    if nargin < 5 || isempty(opts), opts = struct(); end
    if ~isfield(opts, 'npts'),     opts.npts     = 800; end  % 模场重建时每层采样点数
    if ~isfield(opts, 'tail_dec'), opts.tail_dec = 20;  end  % 半无限层向外延伸的衰减长度
    if ~isfield(opts, 'tail_max'), opts.tail_max = 50;  end  % 延伸长度上限 [μm]
    
    k0 = 2*pi/lambda;
    N  = numel(L);
    yb = layer_bounds(L);       
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
        dy = yy - cur_y;        
        
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
    
    % 折射率全实且 n_eff 为实时，清掉运算中产生的 0 虚部; 
    % 虚部来自 complex() 的类型强制，需处理避免绘图警告
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
%   1.默认 20 个衰减长度
%   2.延伸不足会让 Γ 系统性偏小
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
%   1.E(y,z,t) = φ(y)·exp(i(βz − ωt))
%   2.TE:  S_z ∝ β·|E_y|²
%     TM:  S_z ∝ (β/n²)·|H_y|²
%   3.norm_val 用于归一化后复算自检，应恒为 1
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
%       detail：结构体
%           detail.P_act    有源区总功率
%           detail.P_tot    激光器总功率
%           detail.P_layer  各层功率(向量)
% 说明：
%   1.口径需要统一。
%   2.本函数：全部 is_active 层的并集; 无穷积分
%   3.MQWs 是全部阱层的并集。需考虑把势垒排除。
%   4.分母范围需要明确：无穷积分，还是截断到外延结构。对漏模只有无穷积分是自洽的。
%   5.此函数可以考虑和模场归一化合并。
% ───────────────────────────────────────────────────────────────────────────────────
function [G, detail] = confinement_factor(y, E, n_y, yb, L, pol)
    if strcmpi(pol, 'TE')
        w = abs(E).^2;
    else
        w = abs(E).^2 ./ max(real(n_y).^2, eps);
    end
    
    P_tot = trapz(y, w);
    assert(P_tot > 0, 'confinement_factor: 总功率积分非正');
    
    % ── 各层功率（与 is_active 无关）与有源区功率 ────────────────
    N = numel(L);
    P_layer = zeros(1, N);
    P_act   = 0;
    % 将顶层、底层区分出来；获取各层的区间
    for k = 1:N
        if isinf(yb(k,1))           % 底层
            mk = (y <= yb(k,2));
        elseif isinf(yb(k,2))       % 顶层
            mk = (y >= yb(k,1));
        else
            mk = (y >= yb(k,1) & y <= yb(k,2));
        end
        if ~any(mk), continue; end

        Pk = trapz(y(mk), w(mk));           % 本层积分只算一次
        P_layer(k) = Pk / P_tot;
        if L(k).is_active
            P_act = P_act + Pk;
        end
    end
    
    % 判断是否有有源区，无则返回NaN
    % 不返回 0，避免误认为"模式与有源区完全不重叠"
    act = [L.is_active];
    if ~any(act)
        G = NaN;
        detail = struct('P_act', NaN, 'P_tot', P_tot, 'P_layer', P_layer);
        return;
    end
    
    G = P_act / P_tot;
    
    detail = struct('P_act', P_act, 'P_tot', P_tot, 'P_layer', P_layer);
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：快轴远场角分布与 FWHM 发散角
% 输入：快轴坐标 y; 归一化模场 E; 真空波数 k0; 求解要用参数 opts
% 输出：ff：结构体
%         ff.theta   远场角轴 [rad]，范围 ±π/2（半空间）
%         ff.I       角分布 |Φ(θ)|²，已归一化到峰值 1
%         ff.fwhm    半高全宽（FWHM）[rad]
%         ff.th_pk   峰值所在的角度 [rad]
% 调用函数：fwhm_angle
% 说明：
%   1.口径取 FWHM（半高全宽）
%   2.远场角分布 = 近场的空间傅里叶变换（Fraunhofer 近似）：
%       Φ(θ) ∝ ∫ φ(y)·exp(−i·k₀·sinθ·y) dy
%   3.θ 是相对 z 轴的夹角。
%   4.近场需要截断，截断是因为：模场理论上延伸到无穷，只不过在顶层、底层衰减
%     至趋于 0，但是仍是有数值的。因此需要截断，相当于 *一个方窗，即远场多一
%     次 sinc 卷积，影响到 FWHM。存疑时的查法：把 opts.tail_dec由 20 改成 
%     40 重跑，FWHM 不变即证明截断无影响。
%   5.多瓣模（th_pk ≠ 0）的 fwhm 必须连 th_pk 一起读。同一个数字在不同模上
%     含义不同，报告时必须带上 th_pk。逻辑：FWHM 从峰值角出发，寻找两端 =half 
%     的位置（角轴上为对应角度）。存在两种情况：1.该峰与另一个峰都 >half，FWHM
%     囊括多个峰；2.该峰的两端凹陷 <half,FWHM 仅包含这一个峰。  
%   6. θ = 0 不对应任何一个 y，远场 |Φ|² 与 y 的原点选取无关。
% ───────────────────────────────────────────────────────────────────────────────────
function ff = far_field(y, E, k0, opts)
    if nargin < 4 || isempty(opts), opts = struct(); end
    if ~isfield(opts, 'num_ff'), opts.num_ff = 2001; end

    y = reshape(y, 1, []);                             % 统一成行向量
    E = reshape(E, 1, []);

    theta = linspace(-pi/2, pi/2, opts.num_ff);        % 角轴 [rad]
    ky    = k0 * sin(theta);                           % 横向波数 [μm^-1]

    % 傅里叶积分。被积函数 φ(y)·exp(−i·ky·y) 在层界面上连续，且 y 单调递增，
    % 所以可以直接 trapz —— 与 Γ 那处的"分段积分"不是一回事：那里跳的是 mask，
    % 这里积的是连续函数。
    % 2 指沿矩阵的第 2 维积分
    Phi = trapz(y, E .* exp(-1i * ky(:) * y), 2);      % num_ff×1

    I = abs(Phi).^2;
    I = I / max(I);                                    % 峰值归一

    [dw, th_pk] = fwhm_angle(theta, I);

    ff = struct('theta', theta, 'I', I, 'fwhm', dw, 'th_pk', th_pk);
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：角分布的半高全宽
% 输入：角轴 theta [rad]（单调递增）; 角分布 I_theta（已按峰值归一）
% 输出：dw：半高全宽 [rad]; th_pk：峰值所在的角度 [rad]
% 说明：
%   1.自峰值向两侧找第一个低于半高的采样点，在跨半高的相邻两点间线性插值。
%   2.任一侧在半空间内都没降到半高时返回 NaN —— NaN 表示"无定义"，
%   3.不要用一个大数代替：那会被下游读成"发散角很大"这个物理结论。
% ───────────────────────────────────────────────────────────────────────────────────
function [dw, th_pk] = fwhm_angle(theta, I_theta)
    half = 0.5;
    n = numel(theta);
    [~, i0] = max(I_theta);
    th_pk = theta(i0);

    % 峰值落在角轴端点上时，无法向该侧找半高
    if i0 <= 1 || i0 >= n
        dw = NaN;
        return;
    end

    iR = find(I_theta(i0:end)  < half, 1, 'first');
    iL = find(I_theta(i0:-1:1) < half, 1, 'first');
    if isempty(iR) || isempty(iL)
        dw = NaN;                                      % 半空间内未降到半高
        return;
    end
    iR = i0 + iR - 1;
    iL = i0 - iL + 1;

    % 跨半高的相邻两点间线性插值（interp1 要求 X 单调递增）
    thR = interp1([I_theta(iR) I_theta(iR-1)], [theta(iR) theta(iR-1)], half);
    thL = interp1([I_theta(iL) I_theta(iL+1)], [theta(iL) theta(iL+1)], half);

    dw = abs(thR - thL);
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：模场分布绘图
% 输入：坐标系 ax; 快轴坐标 y; 折射率(对应快轴坐标) n_y; 模式数据结构体 m;
%       各层上下限 yb; 层数据结构体 L; 横轴显示范围 xr; 图例位置(字符串) legloc;
%       偏振 pol
% 输出：模场分布图
% 说明：
%   1.在一个坐标轴上同时画折射率分布(左轴)与各阶模场(右轴)
%   2.纵轴范围由当前窗口内的数据算出。
%   3.|φ|² 归一化后峰值 = 1/(有效宽度)
%   4.TE 画功率密度 |E_y|²，TM 画 |H_y|²/n²，两者量级不同，纵轴上限跟着切换。
%   5.层底色必须画在左轴，与折射率线同轴。yyaxis 的两根轴是分层渲染（右轴整体盖在
%       左轴之上，与绘制先后无关）；若把底色画在右轴，它会把折射率线整条盖成灰色。
% ───────────────────────────────────────────────────────────────────────────────────
function draw_stack(ax, y, n_y, m, yb, L, xr, legloc, pol)
    if nargin < 9 || isempty(pol), pol = 'TE'; end
    is_TE = strcmpi(pol, 'TE');

    % 功率权重与纵轴名
    if is_TE
        wfun = @(E) abs(E).^2;
        ylab = '|\phi(y)|^2';
    else
        wfun = @(E) abs(E).^2 ./ max(real(n_y).^2, eps);
        ylab = '|\phi(y)|^2 / n(y)^2';
    end
    
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
        w = wfun(m(k).E);                 
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
    for k = 1:numel(m)
        plot(ax, y, wfun(m(k).E), '-', 'LineWidth', 1.7,  ...
             'Color', cols(k,:),'DisplayName',  ...
             sprintf('%s%d  (n_{eff}=%.5f, \\Gamma=%.4f)', pol, k-1, ...
                     real(m(k).n_eff), m(k).Gamma));
    end
    ylabel(ax, ylab);
    
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


% ───────────────────────────────────────────────────────────────────────────────────
%% ═══════════════════════════════════════════════════════════════════════════════════
% 验证电池（阶段 1）
% 目的：用与 TMM 求解器无关的参照物判定"模型算得对不对"，而不是让模型自证。
%   V0 回归   —— 与阶段2 已公布的无吸收基准对齐
%   V1 有限差 —— 把 Helmholtz 方程直接离散求本征值（方法完全不同的独立解）
%   V2 解析   —— 三层平板界面匹配行列式（闭式，可到 1e-12）
%   V3 不变量 —— 归一化 / ΣP=1 / Γ=ΣP_active / 模序与窗口
%   V4 远场   —— 用 FFT 独立算一遍角分布宽度
%   V5 收敛   —— num_scan / npts / tail_dec / num_ff
%   V6 物理   —— Γ 随波导厚度、n_eff 随波长的单调性
%   V7 量级   —— 405/450/520 的远场与 Nichia 2021 的 40~50°(1/e²) 对照
% 注意：V7 用的各层折射率是"典型量级"，不是 Nichia 原文数据（原文未公布层结构）。
% ═══════════════════════════════════════════════════════════════════════════════════
function run_verify()
    fprintf('\n========== 阶段1 模型验证（Model_Self_L1_7L_EBL_test）==========\n');
    res = [];
    lam    = 0.45;
    k0     = 2*pi/lam;
    opts_v = struct('num_scan',20000,'tol',1e-13,'maxit',200,'merge',1e-6, ...
                    'npts',800,'tail_dec',20,'tail_max',50,'num_ff',2001);

    %% V0 回归锚点
    fprintf('\n-- V0 回归：7 层堆 @450 nm --\n');
    L7 = build_layers(stack_7layer());
    [ne0, G0, fw0, e20] = solve_first(L7, lam, 'TE', opts_v);
    fprintf('   TE0  n_eff = %.12f   Gam = %.9f   FWHM = %.2f deg   1/e2 = %.2f deg\n', ...
            real(ne0), G0, fw0, e20);
    fprintf('   阶段2 公布的无吸收基准 n_eff = 2.417128457\n');
    res(end+1) = chk('V0  n_eff 与阶段2 基准一致 (|d| < 2e-9)', ...
                     abs(real(ne0)-2.417128457) < 2e-9, ...
                     sprintf('d = %.3e', real(ne0)-2.417128457));

    %% V1 有限差分独立解（三对角反迭代，不用 eigs）
    fprintf('\n-- V1 有限差分本征值（网格离散 + Thomas 反迭代，与 TMM 无共享代码）--\n');
    cs = { '对称   2.550/2.374  d=1.5',  2.374, 2.550, 2.374, 1.500 ; ...
           '对称   2.550/2.374  d=0.6',  2.374, 2.550, 2.374, 0.600 ; ...
           '非对称 1.770/2.550/2.374 d=1.5', 1.770, 2.550, 2.374, 1.500 };
    w1 = 0; w1c = 0;
    for c = 1:size(cs,1)
        n0 = cs{c,2}; n1 = cs{c,3}; n2 = cs{c,4}; d = cs{c,5};
        Lt = build_layers(layers3(n0,n1,n2,d));
        ne_t = solve_all(Lt, lam, 'TE', 2, opts_v);
        [ys, ns] = fin_layers(Lt, 3.0);
        ne_a = fd_first(ys, ns, lam, 1.0e-3);
        ne_b = fd_first(ys, ns, lam, 5.0e-4);
        da = abs(real(ne_t(1)) - ne_a);
        db = abs(real(ne_t(1)) - ne_b);
        fprintf('   %s\n      TMM %.9f | FD(h) %.9f  d=%.1e | FD(h/2) %.9f  d=%.1e | 误差比 = %.1f\n', ...
                cs{c,1}, real(ne_t(1)), ne_a, da, ne_b, db, da/max(db,eps));
        w1 = max(w1, db); w1c = max(w1c, da);
    end
    res(end+1) = chk('V1  三层：TMM 与 FD 一致（|d|<1e-4 且加密后误差降 >=2 倍）', ...
                     w1 < 1e-4 && w1c/max(w1,eps) >= 2, ...
                     sprintf('max|d| = %.2e, 误差比 = %.1f', w1, w1c/max(w1,eps)));

    % V1b：7 层堆无解析参照，用 FD 独立复核基模
    [ys7, ns7] = fin_layers(L7, 3.0);
    [neFD, itFD, rsFD] = fd_first(ys7, ns7, lam, 5.0e-4);
    fprintf('   7L TE0  TMM %.12f | FD(h/2) %.12f  d = %.1e  [反迭代 %d 步, 残差 %.1e]\n', ...
            real(ne0), neFD, abs(neFD - real(ne0)), itFD, rsFD);
    res(end+1) = chk('V1b 7 层堆：TMM 与 FD 一致（|d| < 1e-4）', ...
                     abs(neFD - real(ne0)) < 1e-4, ...
                     sprintf('d = %.1e', abs(neFD - real(ne0))));

    %% V2 三层解析
    fprintf('\n-- V2 三层平板解析色散（界面匹配行列式，独立推导）--\n');
    w2 = 0; bad2 = 0;
    for ip = 1:2
        pols = {'TE','TM'}; pol = pols{ip};
        for c = 1:size(cs,1)
            n0 = cs{c,2}; n1 = cs{c,3}; n2 = cs{c,4}; d = cs{c,5};
            Lt = build_layers(layers3(n0,n1,n2,d));
            ne_t = solve_all(Lt, lam, pol, 8, opts_v);
            ne_r = slab_asym_ref(n0,n1,n2,d,lam,pol);
            ne_t = sort(real(ne_t(:)), 'descend');
            ne_r = sort(real(ne_r(:)), 'descend');
            m = min(numel(ne_t), numel(ne_r));
            if m == 0
                dd = Inf;
            else
                dd = max(abs(ne_t(1:m) - ne_r(1:m)));
            end
            fprintf('   %s %-32s : TMM %d 根 / 解析 %d 根  max|d| = %.2e\n', ...
                    pol, cs{c,1}, numel(ne_t), numel(ne_r), dd);
            w2 = max(w2, dd);
            if numel(ne_t) ~= numel(ne_r), bad2 = bad2 + 1; end
        end
    end
    res(end+1) = chk('V2  TMM 与解析式一致（|d| < 1e-9，根数相同）', ...
                     w2 < 1e-9 && bad2 == 0, sprintf('max|d| = %.2e, 根数不符 %d 组', w2, bad2));

    %% V3 不变量
    fprintf('\n-- V3 不变量（三套结构）--\n');
    stacks = {stack_7layer(), stack_9layer_mqw(), layers3(1.770,2.550,2.374,1.5)};
    nmx    = {'7L','9L','test'};
    e3 = 0;
    for s = 1:3
        L = build_layers(stacks{s});
        yb = layer_bounds(L);
        nel = solve_all(L, lam, 'TE', 3, opts_v);
        act = find([L.is_active]);
        for j = 1:numel(nel)
            [y, E, n_y] = mode_profile(nel(j), L, lam, 'TE', opts_v);
            [E, nv]    = normalize_mode(y, E, n_y, 'TE');
            [G, det]   = confinement_factor(y, E, n_y, yb, L, 'TE');
            y = y(:); E = E(:);
            w = abs(E).^2;
            Iall = trapz(y, w);
            mask = false(size(y));
            for k = act(:).'
                mask = mask | (y >= yb(k,1) & y <= yb(k,2));
            end
            Gdir = trapz(y(mask), w(mask)) / Iall;
            e3 = max([e3, abs(nv-1), abs(sum(det.P_layer)-1), abs(sum(det.P_layer(act))-G)]);
            % !! 9 层的 QW 只有 3 nm，会落在 npts=800 的几何网格点之间，trapz 参考值不可信
            if s ~= 2, e3 = max(e3, abs(G - Gdir)); end
            fprintf('   %-5s TE%d  n_eff=%.9f  Gamma=%.6f  |Gamma-trapz|=%.1e  norm-1=%.1e  sumP-1=%.1e\n', ...
                    nmx{s}, j-1, real(nel(j)), G, abs(G-Gdir), abs(nv-1), abs(sum(det.P_layer)-1));
        end
    end
    res(end+1) = chk('V3  不变量（norm=1, sumP=1, Gamma=sumP_act；7L/test 另与 trapz 对拍）', ...
                     e3 < 1e-4, sprintf('max 偏差 = %.2e', e3));

    %% V4 远场 FFT 复核
    fprintf('\n-- V4 远场：far_field(trapz) vs FFT（独立实现）--\n');
    w4 = 0;
    for s = [1 3]
        L = build_layers(stacks{s});
        [ne, G, fw, e2, y, E] = solve_first(L, lam, 'TE', opts_v);
        [fw_ff, e2_ff] = ff_fft_check(y, E, lam);
        d1 = abs(fw - fw_ff)/max(fw,eps);
        fprintf('   %-5s FWHM trapz %.2f deg / FFT %.2f deg (相对差 %.2f%%)   1/e2 %.2f deg\n', ...
                nmx{s}, fw, fw_ff, 100*d1, e2);
        w4 = max(w4, d1);
    end
    res(end+1) = chk('V4  远场宽度与 FFT 一致（相对差 < 1%）', w4 < 0.01, ...
                     sprintf('max 相对差 = %.2f%%', 100*w4));

    %% V5 收敛性
    fprintf('\n-- V5 收敛性（以 7 层堆 @450 nm 的 TE0 为基准）--\n');
    L = build_layers(stack_7layer());
    [ne_r, G_r, fw_r] = solve_first(L, lam, 'TE', opts_v);
    tests = { 'num_scan=2e3',  struct('num_scan',2e3); ...
              'num_scan=1e5',  struct('num_scan',1e5); ...
              'npts=200',      struct('npts',200);   ...
              'npts=3200',     struct('npts',3200);  ...
              'tail_dec=8',    struct('tail_dec',8); ...
              'tail_dec=40',   struct('tail_dec',40);...
              'num_ff=801',    struct('num_ff',801) };
    w5 = 0;
    for t = 1:size(tests,1)
        o = opts_v; f = fieldnames(tests{t,2});
        o.(f{1}) = tests{t,2}.(f{1});
        [ne_t, G_t, fw_t] = solve_first(L, lam, 'TE', o);
        fprintf('   %-12s d n_eff = %+.2e   d Gamma = %+.2e   d FWHM = %+.3f deg\n', ...
                tests{t,1}, real(ne_t-ne_r), G_t-G_r, fw_t-fw_r);
        w5 = max([w5, abs(real(ne_t-ne_r))/max(abs(real(ne_r)),eps), abs(G_t-G_r)/max(G_r,eps)]);
    end
    res(end+1) = chk('V5  参数扰动后 n_eff 与 Gamma 相对变化 < 1e-3（npts=200 属粗网格）', ...
                     w5 < 1e-3, sprintf('max 相对变化 = %.2e', w5));

    %% V6 物理趋势（非单调：应有内部极值）
    fprintf('\n-- V6 物理趋势：芯层加厚 -> Gamma 先升后降（峰在模式被 QW 拉住的区间）--\n');
    dw = [0.04 0.06 0.09 0.12 0.18 0.24 0.32];
    nG = zeros(size(dw)); nF = zeros(size(dw)); nN = zeros(size(dw));
    for i = 1:numel(dw)
        L = build_layers(stack_7layer_var(dw(i), dw(i)*2/3));
        [ne, G, fw] = solve_first(L, lam, 'TE', opts_v);
        nG(i) = G; nF(i) = fw; nN(i) = real(ne);
        fprintf('   dw = %.3f um : n_eff = %.6f  Gamma = %.6f  FWHM = %.2f deg\n', dw(i), nN(i), G, fw);
    end
    [~, iG] = max(nG);  [~, iF] = max(nF);
    okG = iG > 1 && iG < numel(dw);              % Gamma 的极大值在内部，不在端点
    okF = iF > 1 && nF(end) < nF(iF);            % FWHM 先升后转
    res(end+1) = chk('V6  Gamma 在 dw 内部出现极大值，FWHM 先升后转', okG && okF, ...
                     sprintf('Gamma 峰位 %.2f um, FWHM 峰位 %.2f um', dw(iG), dw(iF)));

    %% V7 三个波长的量级对照
    fprintf('\n-- V7 405/450/520 nm 量级对照（折射率为典型值，非 Nichia 数据）--\n');
    fprintf('   Nichia 2021（p.458）：芯片远场椭圆，宽轴 1/e2 强度处 40~50 deg\n');
    fprintf('   注：原文只给角度与脊宽，未给层结构；下表 n(λ) 为典型量级，只能看量级不能判对错。\n');
    for nm = [405 450 520]
        L = build_layers(stack_7layer_at(nm));
        [ne, G, fw, e2] = solve_first(L, nm/1000, 'TE', opts_v);
        fprintf('   %d nm : n_eff = %.6f  Gamma = %.6f  FWHM = %.2f deg  1/e2 全角 = %.2f deg\n', ...
                nm, real(ne), G, fw, e2);
    end
    res(end+1) = chk('V7  仅记录，不判对错（见报告）', true, '量级对照');

    %% V8 用 Nichia 的远场角反推芯包折射率差
    fprintf('\n-- V8 用 Nichia 的 1/e2 远场角反推芯包折射率差 dn --\n');
    fprintf('   (其余层固定，只改 n-clad；原文只给角度，这里反过来问"什么 dn 能复现")\n');
    fprintf('   n-clad    dn     n_eff      Gamma    FWHM    1/e2全角\n');
    best = [];
    for nc = 2.400:0.004:2.430
        L = build_layers(struct( ...
            'd',   {Inf,    1.200, 0.120,  0.009, 0.080, 0.020, Inf}, ...
            'n',   {1.770,  nc,    2.455,  2.550, 2.455, 2.343, 2.374}, ...
            'tag', {'substrate','n-clad','n-wave','QW','p-wave','EBL','p-clad'}, ...
            'is_active', {false, false, false, true, false, false, false}));
        [ne, G, fw, e2] = solve_first(L, lam, 'TE', opts_v);
        fprintf('   %.3f   %.3f   %.6f   %.4f   %5.2f   %5.2f\n', ...
                nc, 2.455-nc, real(ne), G, fw, e2);
        if e2 >= 40 && e2 <= 50, best(end+1) = 2.455 - nc; end
    end
    if isempty(best)
        fprintf('   -> 该区间内无解\n');
        res(end+1) = chk('V8  存在 dn 使 1/e2 全角落在 Nichia 的 40~50 deg 内', false, '无解');
    else
        fprintf('   -> dn 在 [%.3f, %.3f] 时远场角落在 Nichia 的 40~50 deg 内\n', min(best), max(best));
        res(end+1) = chk('V8  存在 dn 使 1/e2 全角落在 Nichia 的 40~50 deg 内', true, ...
                         sprintf('dn 区间 [%.3f, %.3f]', min(best), max(best)));
    end

    fprintf('\n========== 汇总：%d / %d 项通过 ==========\n', sum(res), numel(res));
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：chk（打印一条检查结果）
% ───────────────────────────────────────────────────────────────────────────────────
function ok = chk(name, cond, detail)
    ok = logical(cond);
    if ok
        fprintf('   [PASS] %s   (%s)\n', name, detail);
    else
        fprintf('   [FAIL] %s   (%s)\n', name, detail);
    end
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：layers3（构造三层平板结构）
% ───────────────────────────────────────────────────────────────────────────────────
function layers = layers3(n0, n1, n2, d)
    layers = struct('d', {Inf, d, Inf}, 'n', {n0, n1, n2}, ...
                    'tag', {'bot','core','top'}, ...
                    'is_active', {false, true, false});
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：solve_all（调本文件的求解器，返回降序排列的全部根）
% ───────────────────────────────────────────────────────────────────────────────────
function nel = solve_all(L, lam, pol, nmodes, opts)
    N    = numel(L);
    nmax = max(real([L.n]));
    tiny = 1e-9 * max(1, nmax);
    n_lo = max(real(L(1).n), real(L(N).n)) + tiny;
    n_hi = nmax - tiny;
    if n_hi <= n_lo, nel = []; return; end
    [nel, ~] = find_roots_complex(@(x) dispersion_tmm(x, L, lam, pol), ...
                                  n_lo, n_hi, nmodes + 2, opts);
    nel = sort(nel(:).', 'descend');
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：solve_first（基模的 n_eff / Gamma / 远场 FWHM / 1/e2 全角）
% ───────────────────────────────────────────────────────────────────────────────────
function [ne, Gam, fw, e2, y, E] = solve_first(L, lam, pol, opts)
    k0   = 2*pi/lam;
    nel  = solve_all(L, lam, pol, 3, opts);
    assert(~isempty(nel), 'solve_first: 窗口内无根');
    ne   = nel(1);
    yb   = layer_bounds(L);
    [y, E, n_y] = mode_profile(ne, L, lam, pol, opts);
    [E, ~]      = normalize_mode(y, E, n_y, pol);
    [Gam, ~]    = confinement_factor(y, E, n_y, yb, L, pol);
    ff          = far_field(y, E, k0, opts);
    [fw, e2]    = width_deg(ff.theta, ff.I);
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：width_deg（由角分布直接读 FWHM 与 1/e2 全角，单位 deg）
% 说明：不假设高斯形状，直接在 I(theta) 上插值半高 / e^-2 位置
% ───────────────────────────────────────────────────────────────────────────────────
function [fw, e2] = width_deg(th, I)
    I  = I(:).'; th = th(:).';
    I  = I / max(I);
    [~, i0] = max(I);
    fw = cross_width(th, I, i0, 0.5);
    e2 = cross_width(th, I, i0, exp(-2));
end

function dw = cross_width(th, I, i0, tgt)
    dw = NaN;
    iR = find(I(i0:end)   < tgt, 1, 'first');
    iL = find(I(i0:-1:1)  < tgt, 1, 'first');
    if isempty(iR) || isempty(iL), return; end
    iR = i0 + iR - 1;  iL = i0 - iL + 1;
    a  = interp1([I(iR) I(iR-1)], [th(iR) th(iR-1)], tgt);
    b  = interp1([I(iL) I(iL+1)], [th(iL) th(iL+1)], tgt);
    dw = abs(a - b) * 180/pi;
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：ff_fft_check（用 FFT 独立算远场宽度，核验 far_field 的 trapz 实现）
% 说明：模场先重采样到均匀网格，补零后 FFT，角轴 theta = asin(ky/k0)
% ───────────────────────────────────────────────────────────────────────────────────
function [fw, e2] = ff_fft_check(y, E, lam)
    k0 = 2*pi/lam;
    y  = double(y(:).');  E = real(E(:).');
    N1 = 4096;
    yy = linspace(y(1), y(end), N1);
    EE = interp1(y, E, yy, 'linear', 0);
    h  = yy(2) - yy(1);
    N2 = 2^18;
    F  = fft(EE, N2) * h;
    ky = (0:N2-1)';
    ky(ky > N2/2) = ky(ky > N2/2) - N2;
    ky = 2*pi*ky / (N2*h);
    I  = abs(F).^2;
    sel = abs(ky) < k0;
    th  = asin(ky(sel)/k0);
    Is  = I(sel);
    [th, o] = sort(th); Is = Is(o);
    [fw, e2] = width_deg(th, Is);
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：fd_slab（有限差分直接离散 Helmholtz 方程求本征值）
% 输入：n0/n1/n2 三层折射率; d 芯厚 [um]; lam 波长 [um]; h 网格步长 [um]; k 取前 k 个模
% 输出：ne 降序的本征有效折射率
% 说明：
%   1.方程 E'' + k0^2 (n^2 - ne^2) E = 0 离散为 A E = k0^2 ne^2 E，
%     A = D2 + diag(k0^2 n^2)，两端 Dirichlet。
%   2.这是与转移矩阵完全不同的解法，用来给 TMM 提供独立参照。
%   3.误差 O(h^2)，加密一倍误差应降约 4 倍。
% ───────────────────────────────────────────────────────────────────────────────────
function [ne, itc, res] = fd_first(ys, ns, lam, h)
% 说明：把 (d2/dy2 + k0^2 n^2)E = k0^2 ne^2 E 在非均匀网格上直接离散，三对角 Thomas
%       反迭代 + Rayleigh 商收敛到基模。全程不调用 eigs（本机 R2023b 的 eigs 直接返回 NaN），
%       且与 TMM 不共享任何代码，是真正独立的解法。
    k0   = 2*pi/lam;
    ymax = sum(ys);
    N    = round(ymax/h) + 1;
    y    = linspace(0, ymax, N).';
    e    = [0 cumsum(ys(:).')];
    nn   = ns(1)*ones(N,1);
    for s = 1:numel(ys)
        nn(y > e(s) & y <= e(s+1)) = ns(s);
    end
    ni    = N - 2;
    % !! 按 h^2 无量纲化：对角 -2+k0^2 n^2 h^2、次对角 1；否则 Thomas 主元在 1e6 量级上抵消
    a     = h^2*k0^2*nn(2:N-1).^2 - 2;      % 对角元（已乘 h^2）
    b     = 1;                              % 次对角元（已乘 h^2）
    % !! Sturm 计数 + 二分求最大本征值：Gershgorin 给区间，比反迭代稳，
    %    不会因为初值不当而收敛到某个无关的盒模式
    lo = min(a) - 2*abs(b);
    hi = max(a) + 2*abs(b);
    res = hi - lo;
    for itc = 1:80
        mid = (lo + hi)/2;
        if count_neg(a - mid, b) >= ni
            hi = mid;      % mid 已超过最大本征值
        else
            lo = mid;
        end
        res = hi - lo;
    end
    mu = (lo + hi)/2;
    % 本征值 mu = h^2*k0^2*ne^2
    ne = sqrt(max(mu,0))/(h*k0);
end


function c = count_neg(d, b)
% 说明：对称三对角阵 LDL^T 不选主元分解，负主元个数 = 小于 mu 的本征值个数（Sturm 序列）。
    n  = numel(d);
    c  = 0;
    dp = d(1);
    if dp == 0, dp = 1e-300; end
    if dp < 0, c = c + 1; end
    for i = 2:n
        dp = d(i) - b^2/dp;
        if abs(dp) < 1e-300, dp = sign(dp + double(dp == 0))*1e-300; end
        if dp < 0, c = c + 1; end
    end
end



function [ys, ns] = fin_layers(L, pad)
% 说明：把 Inf 厚度的衬底/包层截成有限厚度 pad，供有限差分离散使用。
    ys = [L.d];
    ns = [L.n];
    ys(~isfinite(ys)) = pad;
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：slab_asym_ref（三层平板解析色散，界面匹配直推，TE/TM 通用）
% 输入：n0 下包层 / n1 芯 / n2 上包层; d 芯厚 [um]; lam [um]; pol
% 输出：ne 降序的全部根
% 说明：
%   1.芯内 E = cos(uy) + (G0/U) sin(uy)，两侧按 exp(-gamma y) 衰减；
%      y = d 处匹配得 (G0 G2 - U^2) sin(u d) + U (G0 + G2) cos(u d) = 0。
%   2.TE: U = u, Gi = gammai；TM: U = u/n1^2, Gi = gammai/ni^2（标准替换）。
%   3.该式对奇偶模统一，无需分奇偶。
% ───────────────────────────────────────────────────────────────────────────────────
function ne = slab_asym_ref(n0, n1, n2, d, lam, pol)
    k0 = 2*pi/lam;
    lo = max(n0, n2) + 1e-10;
    hi = n1 - 1e-10;
    g  = @(x) slab_det(x, n0, n1, n2, d, k0, pol);
    xs = linspace(lo, hi, 30001);
    fs = arrayfun(g, xs);
    ne = [];
    for i = 1:numel(xs)-1
        if fs(i)*fs(i+1) < 0
            try
                ne(end+1) = fzero(g, [xs(i) xs(i+1)]);
            catch
            end
        end
    end
    ne = sort(unique(round(ne*1e10)/1e10), 'descend');
end

function f = slab_det(x, n0, n1, n2, d, k0, pol)
    u  = k0*sqrt(max(n1^2 - x^2, 0));
    g0 = k0*sqrt(max(x^2 - n0^2, 0));
    g2 = k0*sqrt(max(x^2 - n2^2, 0));
    if strcmpi(pol, 'TM')
        U  = u/n1^2;  G0 = g0/n0^2;  G2 = g2/n2^2;
    else
        U  = u;       G0 = g0;       G2 = g2;
    end
    f = (G0*G2 - U^2)*sin(u*d) + U*(G0 + G2)*cos(u*d);
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stack_7layer_var（7 层堆，波导厚度可变）
% ───────────────────────────────────────────────────────────────────────────────────
function layers = stack_7layer_var(dw, dp)
    layers = struct( ...
        'd',   {Inf,    1.200,  dw,     0.009,  dp,     0.020,  Inf}, ...
        'n',   {1.770,  2.404,  2.455,  2.550,  2.455,  2.343,  2.374}, ...
        'tag', {'substrate','n-clad','n-wave','QW','p-wave','EBL','p-clad'}, ...
        'is_active', {false, false, false, true, false, false, false});
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stack_7layer_at（7 层堆，按波长给典型折射率）
% 说明：数值是"典型量级"（GaN 色散 + Al0.06GaN 包层约 -0.026），不是 Nichia 原文数据
% ───────────────────────────────────────────────────────────────────────────────────
function layers = stack_7layer_at(nm)
    switch nm
        case 405
            nn = [1.790 2.474 2.500 2.560 2.500 2.430 2.474];
        case 450
            nn = [1.770 2.414 2.440 2.520 2.440 2.370 2.414];
        case 520
            nn = [1.760 2.354 2.380 2.520 2.380 2.310 2.354];
        otherwise
            error('stack_7layer_at: 只支持 405/450/520');
    end
    layers = struct( ...
        'd',   {Inf, 1.200, 0.120, 0.009, 0.080, 0.020, Inf}, ...
        'n',   num2cell(nn), ...
        'tag', {'substrate','n-clad','n-wave','QW','p-wave','EBL','p-clad'}, ...
        'is_active', {false, false, false, true, false, false, false});
end


% 函数：relayout_tables（把两块表格铺满窗口宽度并上下排好）
% 输入：fig（figure 句柄）; t1（模式表句柄）; t2（层表句柄）; mCols; lCols（两表的列名）
% 输出：无（就地更新 t1、t2 的 ColumnWidth 与 Position）
% 说明：
%   1.传统 figure 的 uitable 没有 RowHeight 属性、行高不可控，故高度按行数固定；
%     而 ColumnWidth 默认 'auto' 不会自己铺满控件，必须按窗口宽度显式等分。
%   2.本函数既用于首次布局，也挂在 fig 的 SizeChangedFcn 上做窗口拉伸重排。
%   3.常量 rowH / gap / mar 须与制表段保持一致。
% ───────────────────────────────────────────────────────────────────────────────────
function relayout_tables(fig, t1, t2, mCols, lCols)
    rowH = 32;                                         % 单行高 [px]，与制表段同值
    gap  = 26;                                         % 两表间隙 [px]
    mar  = 14;                                         % 窗口内边距 [px]

    wp  = get(fig, 'InnerPosition');                   % 客户区当前几何
    w   = wp(3);
    h1  = (size(get(t1, 'Data'), 1) + 1) * rowH + 10;  % 表一高：数据行 + 表头
    h2  = (size(get(t2, 'Data'), 1) + 1) * rowH + 10;
    cw1 = floor((w - 2*mar) / numel(mCols));           % 每列宽 = 可用宽 / 列数
    cw2 = floor((w - 2*mar) / numel(lCols));

    set(t1, 'ColumnWidth', num2cell(repmat(cw1, 1, numel(mCols))), ...
            'Position',    [mar, 2*mar + h2 + gap, w - 2*mar, h1]);
    set(t2, 'ColumnWidth', num2cell(repmat(cw2, 1, numel(lCols))), ...
            'Position',    [mar, mar,               w - 2*mar, h2]);
end
