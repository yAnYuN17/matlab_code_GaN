function modes = Model_Self_L2_7L_EBL(action)
%   ── 命名规则 ──────────────────────────────────────────────
%      函数文件，函数名必须与文件名<.m>一致
%      直接"运行"，以0输入调用该函数
%      可以在命令行调用，如"modes = Model_Self_L2_7L_EBL('test'); "实现有输入调用'save'
%      也可以在另一个文件中调用该函数
%
%   ── 单位与符号约定（全项目统一，不得另立门户）────────────────
%       长度单位 μm 
%       时间因子 exp(−iωt)，沿 +z 取 exp(+iβz)
%       ⇒ α = 2·k₀·Im(n_eff)  [cm⁻¹] 乘 1e4
%       ⇒ Im(n_eff)  >0 [损耗]，<0 [净增益]
%       Γ ≡ ∫_有源区 w dy / ∫_{-∞}^{+∞} w dy
%           TE: w = |E_y|²      TM: w = |H_y|²/n²
%  ── 符号说明────────────────────────────────────────────────
%       !! 仍有疑惑的地方，或者需要后续注意的内容（严重）
%       && 需要再想想，当前无法解决或者难以弄清。
%       % ## 上一版代码内容，在更新后被注释
%       % !! 代码更新后，追加的注释。
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

    % ## 复根的后验检查（阶段2） S4：新增动作 'rootcheck'——根后验检查自测（构造函数 + 真实 7 层结构）
    if contains(lower(string(action)), 'rootcheck')
        modes = stage2_rootcheck_selftest();
        return
    end

    % ## 微扰标尺 B1/B2/B3（不画图，只算）
    if contains(lower(string(action)), 's5')
        modes = stage2_s5_perturbation();
        return
    end

    % ## TMM 数值稳定性（三层均匀吸收，α 1→1e5 cm⁻¹）
    if contains(lower(string(action)), 's6')
        modes = stage2_s6_tmm_stability();
        return
    end

    % ## 阶段2 S7：新增动作 's7self'——用原文 Table II 反查 Table I 系数抄录是否正确
    if contains(lower(string(action)), 's7self')
        modes = stage2_s7_selftest();
        return
    end

    % ## 四通道分解 + 厚度扫描 + 衬底折射率扫描
    if contains(lower(string(action)), 's7')
        modes = stage2_s7_channels();
        return
    end

    % ## 同结构三波长（405/450/520 nm）对比
    if contains(lower(string(action)), 's8')
        modes = stage2_s8_three_lambda();
        return
    end
    
    % λ 提前到层结构之前；三个 stack 调用都传入 λ
    lambda = stage2_lambda();                % [μm]
    if use9
        layers = stack_9layer_mqw(lambda);   
    elseif usetest
        layers = test_layer(lambda);         
    else
        layers = stack_7layer(lambda);       
    end

    % 赋值
    k0 = 2*pi/lambda;
    pol     = 'TE';                 % 'TE' | 'TM'
    n_modes = 3;                    % 最多求几阶模  
    L  = build_layers(layers);      % 校验层结构，折射率强制转 complex
    N  = numel(L);
    yb = layer_bounds(L);           % 定义各层的y区间
    savedir = 'out';                % 图片保存地址-后缀
    
    
    % 求解项赋值
    opts = struct( ...
        'num_scan',   20000, ...    % 色散函数扫描点数（根搜索用）
        'tol',      1e-13, ...      % 求根收敛容差（相对）
        'maxit',    200,   ...      % 最大迭代次数
        'merge',    1e-6,  ...      % 去重时两个根的最小相对间距
        'npts',     800,   ...      % 模场重建时每层采样点数
        'tail_dec', 20,    ...      % 半无限层向外延伸多少个衰减长度
        'tail_max', 50,    ...      % 延伸长度上限 [μm]
        'num_ff',   2001, ...       % 远场角分布采样点数
        'chk_ratio', 0.5);          % |f(r)|/邻域最小值的判定阈值

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
% 调用函数：n_of_lambda，获取含虚部的n
% ───────────────────────────────────────────────────────────────────────────────────
function layers = test_layer(lambda)
    tags = {'substrate','QW','clad'};       
    layers = struct( ...
        'd', ...
            { Inf,     1.5,     Inf      }, ...
        'n', ...
            num2cell(n_of_lambda(lambda, tags)), ...
        'tag', ...
            tags, ...
        'is_active', ...
            { false,     true,      false    });
end

% ───────────────────────────────────────────────────────────────────────────────────
% 7 层结构：衬底 / n包层 / n波导 / QW / p波导 / EBL / p包层
% 首尾两层厚度必须是 Inf（半无限）
% 调用函数：n_of_lambda，获取含虚部的n
% ───────────────────────────────────────────────────────────────────────────────────
function layers = stack_7layer(lambda)
    tags = {'substrate','n-clad','n-wave','QW','p-wave','EBL','p-clad'};   
    layers = struct( ...
        'd', ...
            { Inf,      1.200,      0.120,      0.009,      0.080,      0.020,      Inf      }, ...
        'n', ...
            num2cell(n_of_lambda(lambda, tags)), ...
        'tag', ...
            tags, ...
        'is_active', ...
            { false,    false,      false,      true,       false,      false,      false    });
end


% ───────────────────────────────────────────────────────────────────────────────────
% 9 层结构：把单层等效 QW 展开成 阱-垒-阱
% Γ = (阱1 + 阱2 的功率) / 总功率
% 调用函数：n_of_lambda，获取含虚部的n
% ───────────────────────────────────────────────────────────────────────────────────
function layers = stack_9layer_mqw(lambda)
    tags = {'substrate','n-clad','n-wave','QW1','barrier','QW2','p-wave','EBL','p-clad'};   
    layers = struct( ...
        'd', ...
            { Inf,   1.200,   0.120,   0.003,   0.006,   0.003,   0.080,   0.020,   Inf    }, ...
        'n', ...
            num2cell(n_of_lambda(lambda, tags)), ...
        'tag', ...
            tags, ...
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
% 输出：roots_n: 满足条件的n_eff; chk: ; cand_all:
% 调用函数：newton_complex, stage2_empty_check [用于数据记录 (→ 后验检查) 的空数组]
% 说明：
%   本函数用于复折射率
% ───────────────────────────────────────────────────────────────────────────────────
function [roots_n, chk, cand_all] = scan_newton(fun, xs, fs, opts)
    %SCAN_NEWTON  复值分支：定位 |f| 的局部极小，再用复数牛顿细化
    mag = abs(fs);                      % 取模，而不是取绝对值。
    cand = [];

    % 选取局部的谷，作为备选的迭代初值
    % f(z)=0 → |f(z)=0|  →  在离散点中，邻近点数值小于其边上两点的数值
    for j = 2:numel(xs)-1
        if mag(j) < mag(j-1) && mag(j) <= mag(j+1)
            cand(end+1) = xs(j);        %#ok<AGROW>
        end
    end

    roots_n  = [];
    chk      = stage2_empty_check();
    cand_all = stage2_empty_check();
    for x0 = cand
        [r, c] = newton_complex(fun, x0, opts);
        c.x0 = x0;
        cand_all(end+1) = c;                    %#ok<AGROW>  % 全部候选都留记录
        if ~isempty(r)
            roots_n(end+1) = r;                 %#ok<AGROW>
            chk(end+1)     = c;                 %#ok<AGROW>
        end
    end
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：牛顿法求解
% 输入：函数 dispersion_tmm(输出转移矩阵结果); 待迭代的复折射率 x0(向量); 求解要用参数 opts
% 输出：roots_n: 满足条件的n_eff; chk: ; cand_all:
% 调用函数：newton_complex, stage2_check_one 
% 说明：
%   1.f(z_k + \delta) \approx f(z_k)+f'(z_k)\delta, \delta = - [ f(z_k) / f'(z_k)]
%       z_(k+1) = z_k + \delta
%   2.&& 设置了步长阻尼，可以在opts额外设置参数
% ───────────────────────────────────────────────────────────────────────────────────
function [r, chk] = newton_complex(fun, x0, opts)
    %NEWTON_COMPLEX  复数牛顿迭代，带步长阻尼。无法收敛时返回空。
    h      = 1e-7;                       % 用于求导的步长
    r      = x0;
    n_iter = opts.maxit;

    for it = 1:n_iter
        % 求解f'(r)，即df
        fx = fun(r);
        f1 = fun(r + h);
        f2 = fun(r - h);
        df = (f1 - f2) / (2*h);

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
        
        % 迭代步长小于一定值，退出记录(还需后验检查)
        if abs(step) < 1e-14 * max(1, abs(r))
            break;
        end
        if it == n_iter
            chk = stage2_check_one(fun, r, opts);
            r = [];                        % 未收敛，丢弃
            return
        end
    end

    % 阶段2 S4：后验检查（多尺度邻域，逐根记录 |f(r)|、邻域最小值、比值、判定）
    chk = stage2_check_one(fun, r, opts);
    if ~chk.accepted
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
    
    % ── 第 2 … N-1 层：和转移矩阵一致 ──────────────────────────
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

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：波长赋值
% 输入：无
% 输出：lambda：真空波长 [μm]
% 说明：
%   1.阶段2 规定 λ 只在这里改；改这一行即换波长，别处不许再出现 0.45。
%   2.!! 405/520 nm 的 n(λ) 与 α(λ) 在《阶段2 施工单》§九 里明确"无可引用来源"，
%     换波长时实部仍取 450 nm 表（n450_of_tag），结论只能写"仅限几何/趋势"。
%   3.本函数是为了减少后续波长的连续调用和赋值（到损耗那段再看看 &&）
% ───────────────────────────────────────────────────────────────────────────────────
function lambda = stage2_lambda()
    lambda = 0.45;             
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：复折射率赋值
% 输入：lambda：真空波长 [μm]; tags：层名 cell（长度 = 层数）
% 输出：nv：与 tags 等长的复折射率行向量
% 调用函数：n450_of_tag，alpha_of_layer
% 说明：
%   1.实部+虚部
%   2.Im(n) > 0 表示【损耗】，与 build_layers 的 imag(n) >= 0 断言同向。
%   3.α [cm⁻¹] 与 Im(n) 的换算：Im(n) = α / (2·k0·1e4)，k0 = 2π/λ [μm⁻¹]；
%     450 nm 时 α = 10 cm⁻¹ ↔ Im(n) = 3.581e-5（施工单 B1 的把手）。
%   4.!! 实部只有 450 nm 一组有来源
%   5.!! 无 n(λ) 色散来源。
% ───────────────────────────────────────────────────────────────────────────────────
function nv = n_of_lambda(lambda, tags)
    k0 = 2*pi/lambda;
    nv = zeros(1, numel(tags));
    for k = 1:numel(tags)
        nr = n450_of_tag(tags{k});
        a  = alpha_of_layer(tags{k}, lambda);       % [cm^-1]
        nv(k) = complex(nr, a / (2*k0*1e4));
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：450 nm 逐层折射率实部表
% 输入：tag：层名
% 输出：nr：折射率实部（无量纲）
% 说明：
%   1.数值来自阶段1, 可后续更改。
%   2.确保已有模型中的'tag'和下面是一一对应的。
%   3.test_layer 用 'clad' 这个层名，值与 p-clad 相同（2.374）。
% ───────────────────────────────────────────────────────────────────────────────────
function nr = n450_of_tag(tag)
    switch tag
        case 'substrate', nr = 1.770;
        case 'n-clad',    nr = 2.404;
        case 'n-wave',    nr = 2.455;
        case 'QW',        nr = 2.550;
        case 'QW1',       nr = 2.550;
        case 'barrier',   nr = 2.430;
        case 'QW2',       nr = 2.550;
        case 'p-wave',    nr = 2.455;
        case 'EBL',       nr = 2.343;
        case 'p-clad',    nr = 2.374;
        case 'clad',      nr = 2.374;
        otherwise
            error('n450_of_tag: 未知层名 %s', tag);
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：逐层吸收系数 α [cm⁻¹]
% 输入：tag：层名(字符串); lambda：真空波长 [μm]
% 输出：a：该层吸收系数 [cm⁻¹]
% 说明：
%   1.由全局 containers.Map 进行赋值。
%   2.使用 Map，是因为 Struct 不允许使用'-'作为字段名。
%       STAGE2_ALPHA = containers.Map({'QW'}, {10});
%   3.containers.Map(keys,value)：key=层名，value=α [cm^-1]；空表示默认 0。键区分大小写
%   4.未使用 lambda，作为预留接口。
%   5.&& 加层时，需要对 known 同步修改。
% ───────────────────────────────────────────────────────────────────────────────────
function a = alpha_of_layer(tag, lambda) %#ok<INUSD>
    global STAGE2_ALPHA
    known = {'substrate','n-clad','n-wave','QW','QW1','barrier','QW2','p-wave','EBL','p-clad','clad'};
    if ~any(strcmp(tag, known))
        error('alpha_of_layer: 未知层名 %s', tag);
    end

    % isa(A,className)，变量类型判断;  
    % isKey，只适用于 containers.Map，判断 Map 里是否有该键。
    if ~isempty(STAGE2_ALPHA) && isa(STAGE2_ALPHA, 'containers.Map') && isKey(STAGE2_ALPHA, tag)
        a = STAGE2_ALPHA(tag);          % 阶段2 S3–S7：测试/扫描注入
    else
        a = 0;                          % 默认无吸收
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
%   4.半无限层衰减常数的分支选择
% ───────────────────────────────────────────────────────────────────────────────────
function g = BranchChoose_decay(g)
    flip = real(g) < 0 | (real(g) == 0 & imag(g) < 0);
    g(flip) = -g(flip);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：单根后验检查
% 输入：fun：色散函数; r：待检查的根; opts：需 chk_ratio（判定阈值）、chk_h（尺度序列）
% 输出：c：结构体
%       c.fr=|f(r)|、
%       c.nb=邻域最小值、
%       c.ratio=fr/nb、
%       c.accepted=判定
% 说明：
%   1.原检查只取 1e-3 一个尺度，且与邻域最大值比较：盆底很平的非零极小（伪根）
%     在 1e-3 外 |f| 已经长大，"看起来像根"而被放过。
%   2.改为 1e-3…1e-8 六个尺度 × 8 个方向取最小值：真根的 |f| 随步长线性趋零，
%     伪根在盆底附近 |f| 基本不变 ⇒ ratio ≈ 1，被拒。
%   3.ratio 记录进 info.checks，是 B2 的验收把手。
%   4.阶段2 S6 补绝对残差门限 fr <= chk_res·fscale：单靠 ratio 会放过
%     "陡峭函数的非零局部极小"（S6 第三段 d=5 的反例：|f|=33 仍被判为根）。
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
            if v < nb, nb = v; end
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
% 函数：stage2_rootcheck_selftest（S4 自测：真实结构 + 构造真根 + 构造伪根）
% 输入：无
% 输出：out：结构体
%           out.n_real   
%           out.n_true   
%           out.n_pseudo 
%           out.ratio_real 
%           out.ratio_pseudo 
%           out.ncand_pseudo 
%           out.pass
% 调用函数：stage2_lambda，stack_7layer，dispersion_tmm，find_roots_complex
% 说明：
%   1.真实结构用 7 层、TE、450 nm，搜索区间与本模型主流程同口径；
%   2.构造真根 f=(x−2.40)(x−2.35)：应返回 2 个、全部通过；
%   3.构造伪根 f=|x−2.40|²+1e-9：|f| 在 2.40 取非零极小、牛顿会收敛到谷底，
%     但 ratio≈1 ⇒ 必须被后验检查拒绝（返回 0 根）；
%   4.判据：真实结构每根 ratio < chk_ratio，且 2 真根全留、伪根全拒。
% ───────────────────────────────────────────────────────────────────────────────────
function out = stage2_rootcheck_selftest()
    opts = struct('num_scan', 20000, 'tol', 1e-13, 'maxit', 200, ...
                  'merge', 1e-6, 'chk_ratio', 0.5, 'chk_h', 10.^-(3:8));

    lambda = stage2_lambda();
    Lt     = stack_7layer(lambda);
    fh     = @(x) dispersion_tmm(x, Lt, lambda, 'TE');
    n_lo   = max(real(Lt(1).n), real(Lt(end).n)) + 1e-9;
    n_hi   = max(real([Lt.n])) - 1e-9;

    [roots_real, info_real] = find_roots_complex(fh, n_lo, n_hi, 5, opts);

    % 真根：全纯二次式，两个根应全部留下
    f_true   = @(x) (x - 2.40).*(x - 2.35) + 1i*1e-11;
    [roots_true, info_true] = find_roots_complex(f_true, 2.30, 2.46, 5, opts);

    % 伪根：给全纯函数加一个复常数并不能造出伪根——它在复平面上照样有零点。
    % 真正要复现的形态是"|f| 在某条方向上降不下去"（色散函数分支切换处就是这样）。
    % 取 f=(Re x−2.40)²+1e-9：沿虚轴恒为 1e-9，任何步长都不减小 ⇒ 必须被拒。
    f_pseudo = @(x) (1 + 1i)*abs(x - 2.40).^2 + 1e-9;   % 非全纯 ⇒ 复平面上无零点
    [roots_pseudo, info_pseudo] = find_roots_complex(f_pseudo, 2.30, 2.46, 5, opts);

    fprintf('\n===== S4 根后验检查（|f(root)| / 邻域最小值）=====\n');
    stage2_print_checks('真实结构7L', roots_real,   info_real,   opts);
    stage2_print_checks('构造真根  ', roots_true,   info_true,   opts);
    stage2_print_checks('构造伪根  ', roots_pseudo, info_pseudo, opts);

    out = struct();
    out.n_real   = numel(roots_real);
    out.n_true   = numel(roots_true);
    out.n_pseudo = numel(roots_pseudo);
    out.ratio_real = [info_real.checks.ratio];
    out.ratio_pseudo = inf;
    out.ncand_pseudo = 0;
    if isfield(info_pseudo, 'cand') && ~isempty(info_pseudo.cand)
        out.ncand_pseudo = numel(info_pseudo.cand);
        out.ratio_pseudo = min([info_pseudo.cand.ratio]);
    end
    out.pass = (out.n_real == 3) && all(out.ratio_real < opts.chk_ratio) && ...
               (out.n_true == 2) && (out.n_pseudo == 0) && ...
               (out.ncand_pseudo >= 1) && (out.ratio_pseudo > opts.chk_ratio);
    fprintf('S4 SELFTEST pass=%d n_real=%d n_true=%d n_pseudo=%d\n', ...
        out.pass, out.n_real, out.n_true, out.n_pseudo);
    fprintf('S4 ratios_real = %s\n', mat2str(out.ratio_real, 4));
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_print_checks（打印一组根的后验记录）
% 输入：tag：标签; roots_n：根; info：find_roots_complex 的 info; opts：参数
% 输出：无（命令行输出）
% 说明：
%   1.每条根给 root、|f(root)|、邻域最小值、ratio、是否接受；
%   2.info.checks 为空时只打印候选数，便于判断"根被拒"还是"根本没找到候选"。
% ───────────────────────────────────────────────────────────────────────────────────
function stage2_print_checks(tag, roots_n, info, opts)
    if isfield(info, 'checks') && ~isempty(info.checks)
        for k = 1:numel(roots_n)
            c = info.checks(k);
            fprintf('%s root=%.10f%+.3e  |f|=%.3e  nb=%.3e  ratio=%.3e  ok=%d\n', ...
                tag, real(roots_n(k)), imag(roots_n(k)), c.fr, c.nb, c.ratio, c.accepted);
        end
    else
        ncand = 0;
        rmin  = NaN;
        if isfield(info, 'cand') && ~isempty(info.cand)
            ncand = numel(info.cand);
            rmin  = min([info.cand.ratio]);
        end
        fprintf('%s 无被接受的根（候选 %d 个，最小 ratio=%.3e，阈值 %.1f）\n', ...
            tag, ncand, rmin, opts.chk_ratio);
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：用于后验记录的空 struct 数组
% 输入：无
% 输出：chk：0×1 的空 struct 数组
%           chk.fr
%           chk.nb
%           chk.ratio
%           chk.rmax
%           chk.accepted
%           chk.x0
% 说明：
%   1.struct 数组不能从 [] 起步（[] 是 double，end+1 赋值报"无法从 struct 转换为 double"，
%     坑记录 D1 同类）；必须用带字段的空 struct 起手。
%   2.字段顺序与 stage2_check_one 的返回一致，末尾多一个 x0（由调用方填写）。
% ───────────────────────────────────────────────────────────────────────────────────
function chk = stage2_empty_check()
    chk = struct('fr', {}, 'nb', {}, 'ratio', {}, 'rmax', {}, 'accepted', {}, 'x0', {});
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_s5_perturbation（S5 微扰标尺：B1 弱吸收极限 / B2 根后验 / B3 多通道叠加）
% 输入：无
% 输出：out：结构体，含 B1/B2/B3 的判据值与 pass 标志
% 调用函数：stack_7layer，layer_bounds，dispersion_tmm，find_roots_complex，
%           mode_profile，normalize_mode，confinement_factor，alpha_of_layer
% 说明：
%   1.B1：逐层单独赋 α_j ∈ {1,3,10} cm⁻¹，α_mode/α_j 与 P_layer(j) 比较；
%   2.严格一阶微扰给的是 α_mode = α_j·P_j·n_rj/Re(n_eff)（S3 全层同 α 实测已证实），
%     故同时给"朴素判据"（对 P_j）与"严格判据"（对 P_j·n_rj/n_eff）；
%   3.B3：多层同时赋 α（各 ≤10），与 Σ_j P_j α_j 及严格式比较；
%   4.全程 TE、7 层、450 nm；α 通过全局表 STAGE2_ALPHA 注入，函数退出时清空。
% ───────────────────────────────────────────────────────────────────────────────────
function out = stage2_s5_perturbation()
    global STAGE2_ALPHA

    lambda = stage2_lambda();
    k0     = 2*pi/lambda;
    pol    = 'TE';
    opts   = struct('num_scan', 2e4, 'tol', 1e-13, 'maxit', 200, 'merge', 1e-6, ...
                    'npts', 800, 'tail_dec', 20, 'tail_max', 50, 'num_ff', 2001, ...
                    'chk_ratio', 0.5, 'chk_h', 10.^-(3:8));

    L  = stack_7layer(lambda);
    N  = numel(L);
    yb = layer_bounds(L);
    n_max = max(real([L.n]));
    tiny  = 1e-9 * max(1, n_max);
    n_lo  = max(real(L(1).n), real(L(N).n)) + tiny;
    n_hi  = n_max - tiny;
    fh    = @(x) dispersion_tmm(x, L, lambda, pol);
    tags  = {L.tag};

    % ── 基准：全层无吸收 ──
    STAGE2_ALPHA = [];
    [ne0, G0, P0, ck0] = stage2_s5_solve_te0(fh, L, lambda, pol, opts, yb, n_lo, n_hi);
    nr0 = real(ne0);
    a0  = 2*k0*imag(ne0)*1e4;

    fprintf('\n===== S5 基准（全层无吸收）=====\n');
    fprintf('Re(n_eff)=%.12f  Im(n_eff)=%.3e  Γ=%.9f  α_mode=%.3e cm^-1\n', ...
        nr0, imag(ne0), G0, a0);
    fprintf('ΣP_layer-1 = %.3e   根后验 max ratio = %.3e\n', sum(P0)-1, max([ck0.ratio]));

    % ── B1：逐层单独赋 α ──
    alist = [1 3 10];
    nb1   = numel(tags) * numel(alist);
    b1 = struct('tag', {}, 'nr', {}, 'P', {}, 'alpha', {}, 'a_mode', {}, ...
                'ratio', {}, 'rel_naive', {}, 'rel_strict', {}, 'ok_root', {});
    k = 0;
    for j = 1:numel(tags)
        for a = alist
            k = k + 1;
            STAGE2_ALPHA = containers.Map({tags{j}}, {a});
            % 注意：层折射率是在 stack_* 里一次性算好的，改 α 后必须重建层结构，
            % 否则拿到的还是上一轮的 n（S5 首跑 α_mode 恒 0 就是这个原因）。
            L  = stack_7layer(lambda);
            fh = @(x) dispersion_tmm(x, L, lambda, pol);
            [ne, ~, Pl, ck] = stage2_s5_solve_te0(fh, L, lambda, pol, opts, yb, n_lo, n_hi);
            am = 2*k0*imag(ne)*1e4;
            b1(k).tag        = tags{j};
            b1(k).nr         = real(n450_of_tag(tags{j}));
            b1(k).P          = P0(j);
            b1(k).alpha      = a;
            b1(k).a_mode     = am;
            b1(k).ratio      = am / a;
            b1(k).rel_naive  = abs(am/a / max(P0(j), realmin) - 1);
            b1(k).rel_strict = abs(am/a / max(P0(j)*b1(k).nr/nr0, realmin) - 1);
            b1(k).ok_root    = all([ck.ratio] < opts.chk_ratio);
        end
    end
    STAGE2_ALPHA = [];

    fprintf('\n===== B1 逐层弱吸收（TE0，450 nm）=====\n');
    fprintf('%-10s %4s %12s %12s %12s %10s %10s %10s\n', ...
        '层', 'α', 'α_mode', 'α_mode/α', 'n_r/n_eff', 'rel_朴素', 'rel_严格', '根ok');
    for k = 1:nb1
        fprintf('%-10s %4g %12.4e %12.6e %12.6f %10.2e %10.2e %10d\n', ...
            b1(k).tag, b1(k).alpha, b1(k).a_mode, b1(k).ratio, ...
            b1(k).nr/nr0, b1(k).rel_naive, b1(k).rel_strict, b1(k).ok_root);
    end

    % 单层内非线性度：α=10 与 α=1 的单位 α 衰减之比
    fprintf('\n-- B1 线性度（同层 α=10 与 α=1 的单位 α 衰减之比，理想 = 1）--\n');
    nl_max = 0;
    for j = 1:numel(tags)
        idx = find(strcmp({b1.tag}, tags{j}));
        r1  = b1(idx(1)).ratio;   r10 = b1(idx(end)).ratio;
        nl  = abs(r10/r1 - 1);
        nl_max = max(nl_max, nl);
        fprintf('%-10s  ratio(1)=%.6e  ratio(10)=%.6e  非线性度=%.3e\n', ...
            tags{j}, r1, r10, nl);
    end

    rel_naive_max  = max([b1.rel_naive]);
    rel_strict_max = max([b1.rel_strict]);
    b1_pass = (rel_strict_max < 1e-2) && (nl_max < 1e-2) && all([b1.ok_root]);

    % ── B3：多通道叠加 ──
    multi_tags = {'n-clad', 'n-wave', 'QW', 'p-wave', 'EBL', 'p-clad'};
    multi_a    = [3, 5, 10, 4, 2, 6];      % cm^-1，各 ≤10
    STAGE2_ALPHA = containers.Map(multi_tags, num2cell(multi_a));
    L   = stack_7layer(lambda);                       % 同上：改 α 后重建层结构
    fh  = @(x) dispersion_tmm(x, L, lambda, pol);
    [ne3, ~, Pl3, ck3] = stage2_s5_solve_te0(fh, L, lambda, pol, opts, yb, n_lo, n_hi);
    STAGE2_ALPHA = [];

    am3      = 2*k0*imag(ne3)*1e4;
    pred_sum = 0;  pred_strict = 0;
    for m = 1:numel(multi_tags)
        jj = find(strcmp(tags, multi_tags{m}));
        pred_sum    = pred_sum    + P0(jj)*multi_a(m);
        pred_strict = pred_strict + P0(jj)*multi_a(m)*real(n450_of_tag(tags{jj}))/nr0;
    end
    rel3_sum    = abs(am3/pred_sum - 1);
    rel3_strict = abs(am3/pred_strict - 1);

    fprintf('\n===== B3 多通道叠加（6 层同时赋 α ≤10）=====\n');
    for m = 1:numel(multi_tags)
        jj = find(strcmp(tags, multi_tags{m}));
        fprintf('%-10s α=%5.1f  P=%10.6f  贡献=%10.4e\n', ...
            multi_tags{m}, multi_a(m), P0(jj), P0(jj)*multi_a(m));
    end
    fprintf('α_mode(实算) = %.6e cm^-1\n', am3);
    fprintf('ΣP_j α_j       = %.6e  相对偏差 = %.3e\n', pred_sum,    rel3_sum);
    fprintf('ΣP_j α_j n_rj/n_eff = %.6e  相对偏差 = %.3e\n', pred_strict, rel3_strict);
    fprintf('根后验 max ratio = %.3e\n', max([ck3.ratio]));
    b3_pass = (rel3_sum < 2e-2) && all([ck3.ratio] < opts.chk_ratio);

    % ── 汇总 ──
    out = struct();
    out.ne0 = ne0; out.Gamma0 = G0; out.P0 = P0; out.alpha0 = a0;
    out.b1 = b1;
    out.b1_rel_naive_max  = rel_naive_max;
    out.b1_rel_strict_max = rel_strict_max;
    out.b1_nonlin_max     = nl_max;
    out.b1_pass = b1_pass;
    out.b3_alpha_mode = am3;
    out.b3_pred_sum = pred_sum;
    out.b3_pred_strict = pred_strict;
    out.b3_rel_sum = rel3_sum;
    out.b3_rel_strict = rel3_strict;
    out.b3_pass = b3_pass;
    out.pass = b1_pass && b3_pass;

    fprintf('\n===== S5 汇总 =====\n');
    fprintf('B1  max rel(朴素判据) = %.3e   max rel(严格判据) = %.3e   非线性度 = %.3e\n', ...
        rel_naive_max, rel_strict_max, nl_max);
    fprintf('B3  rel(ΣP_jα_j) = %.3e   rel(严格式) = %.3e\n', rel3_sum, rel3_strict);
    fprintf('S5  pass=%d  (B1=%d  B3=%d)\n', out.pass, b1_pass, b3_pass);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_s5_solve_te0（解 TE0：复本征值 + 模场 + Γ + 逐层功率占比）
% 输入：fh：色散函数; L：层结构; lambda：波长; pol：极化; opts：求解参数;
%       yb：层边界; n_lo/n_hi：搜索区间
% 输出：ne：TE0 的 n_eff（可复）; G：Γ; Pl：逐层功率占比; ck：根后验记录
% 说明：
%   1.取实部最大的根为 TE0；按实部降序（A1：不能按模排序）。
%   2.复用主流程的 mode_profile / normalize_mode / confinement_factor，不另立口径。
% ───────────────────────────────────────────────────────────────────────────────────
function [ne, G, Pl, ck] = stage2_s5_solve_te0(fh, L, lambda, pol, opts, yb, n_lo, n_hi)
    [roots_n, info] = find_roots_complex(fh, n_lo, n_hi, 6, opts);
    assert(~isempty(roots_n), 'stage2_s5_solve_te0: 未找到根');
    [~, ord] = sort(real(roots_n(:).'), 'descend');
    ne = roots_n(ord(1));

    [y, E, n_y] = mode_profile(ne, L, lambda, pol, opts);
    [E, ~]      = normalize_mode(y, E, n_y, pol);
    [G, det]    = confinement_factor(y, E, n_y, yb, L, pol);
    Pl          = det.P_layer;
    if isfield(info, 'checks')
        ck = info.checks(ord(1));
    else
        ck = stage2_check_one(fh, ne, opts);
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_s6_tmm_stability（S6：TMM 数值稳定性）
% 输入：无
% 输出：out：结构体，含三段扫描表、失真起点、是否失败
% 调用函数：stage2_s6_point
% 说明：
%   1.第一段按施工单：对称三层（n_core=2.455、n_clad=2.404、d=0.12 μm，取自 n-wave），
%     三层同赋吸收，α 从 1 扫到 1e5 cm⁻¹，搜根窗口是全窗口；
%   2.第二段测条件数：d=0.3 μm（仍单模），把搜根窗口收到解析根 ±0.02，
%     排除"模序错配"，直接看同一个根的精度随 |Im κd| 的变化；
%   3.第三段是反面对照：d=5 μm（多模），全窗口与窄窗口各跑一次，
%     用来证明"全窗口下 dn 变大"是多模选错根，不是矩阵病态；
%   4.误差：dn=|n_tmm−n_解析|、dα=2k0·1e4·|Im(n_tmm)−Im(n_解析)| [cm⁻¹]；
%   5."找不到根"与"精度失真"分开记：前者是求解窗口问题，不是矩阵条件数问题。
% ───────────────────────────────────────────────────────────────────────────────────
function out = stage2_s6_tmm_stability()
    lambda = stage2_lambda();
    opts   = struct('num_scan', 2e4, 'tol', 1e-13, 'maxit', 300, 'merge', 1e-6, ...
                    'chk_ratio', 0.5, 'chk_h', 10.^-(3:8));
    n_core = 2.455;
    n_clad = 2.404;

    % ── 第一段：施工单规定，d = 0.12 μm，全窗口扫 α ──
    d1    = 0.12;
    alist = [0 1 3 10 30 100 300 1e3 3e3 1e4 3e4 6e4 1e5];
    fprintf('\n===== S6 第一段：d=%.2f μm，全窗口，扫 α =====\n', d1);
    fprintf('%9s %10s %11s %13s %12s %12s %10s %5s\n', ...
        'α[cm^-1]', 'α·d', '|Im κd|', 'Re(n_tmm)', 'dn', 'dα[cm^-1]', 'dα/α', '根数');
    t1 = repmat(stage2_s6_point(n_core, n_clad, d1, 0, lambda, opts), 1, numel(alist));
    seed = 2.43;
    for i = 1:numel(alist)
        t1(i) = stage2_s6_point(n_core, n_clad, d1, alist(i), lambda, opts, [], seed);
        seed  = t1(i).ne_a;
        fprintf('%9g %10.3e %11.3e %13.9f %12.3e %12.3e %10.2e %5d\n', ...
            t1(i).alpha, t1(i).ad, t1(i).kdi, real(t1(i).ne_t), ...
            t1(i).dn, t1(i).dalpha, t1(i).rel, t1(i).nroots);
    end

    % ── 第二段：d = 0.3 μm 单模，窄窗口 ±0.02，压 |Im κd| ──
    d2     = 0.30;
    win2   = 0.02;
    alist2 = [100 300 1e3 3e3 1e4 3e4 6e4 1e5 2e5];
    fprintf('\n===== S6 第二段：d=%.2f μm，窄窗口 ±%.2f，加压 =====\n', d2, win2);
    fprintf('%9s %10s %11s %13s %12s %12s %10s %9s %5s\n', ...
        'α[cm^-1]', 'α·d', '|Im κd|', 'Re(n_tmm)', 'dn', 'dα/α', 'cond(M)', 'amp', '根数');
    t2 = repmat(stage2_s6_point(n_core, n_clad, d2, 0, lambda, opts, win2), 1, numel(alist2));
    seed = 2.43;
    for i = 1:numel(alist2)
        t2(i) = stage2_s6_point(n_core, n_clad, d2, alist2(i), lambda, opts, win2, seed);
        seed  = t2(i).ne_a;
        fprintf('%9g %10.3e %11.3e %13.9f %12.3e %12.2e %10.2e %9.2e %5d\n', ...
            t2(i).alpha, t2(i).ad, t2(i).kdi, real(t2(i).ne_t), ...
            t2(i).dn, t2(i).rel, t2(i).condM, t2(i).amp, t2(i).nroots);
    end

    % ── 第三段：反面对照，d=5 μm 多模，全窗口 vs 窄窗口 ──
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

    % ── 第四段：真实 7 层 + 文献量级 α（S7 的预演）──
    global STAGE2_ALPHA
    k0 = 2*pi/lambda;
    STAGE2_ALPHA = [];
    L0  = stack_7layer(lambda);
    yb0 = layer_bounds(L0);
    nmx = max(real([L0.n]));
    w0  = 1e-9*max(1,nmx);
    lo0 = max(real(L0(1).n), real(L0(end).n)) + w0;
    hi0 = nmx - w0;
    fh0 = @(x) dispersion_tmm(x, L0, lambda, 'TE');
    [neA, GA, PA, ckA] = stage2_s5_solve_te0(fh0, L0, lambda, 'TE', opts, yb0, lo0, hi0);

    STAGE2_ALPHA = containers.Map({'n-clad','p-clad','EBL'}, {36, 18, 36});  % 文献量级，仅作条件数预演
    L1  = stack_7layer(lambda);
    fh1 = @(x) dispersion_tmm(x, L1, lambda, 'TE');
    [neB, GB, PB, ckB] = stage2_s5_solve_te0(fh1, L1, lambda, 'TE', opts, layer_bounds(L1), lo0, hi0);
    [~, aux] = dispersion_tmm(neB, L1, lambda, 'TE');
    STAGE2_ALPHA = [];

    fprintf('\n===== S6 第四段：真实 7 层 + 量级 α（n-clad 36 / p-clad 18 / EBL 36 cm⁻¹）=====\n');
    fprintf('α=0    : Re(n)=%.12f  Γ=%.9f  ΣP−1=%.1e\n', real(neA), GA, abs(sum(PA)-1));
    fprintf('α≠0    : Re(n)=%.12f  Γ=%.9f  ΣP−1=%.1e  Im(n)=%.9f  α_mode=%.4f cm⁻¹\n', ...
        real(neB), GB, abs(sum(PB)-1), imag(neA*0+neB), 2*k0*1e4*imag(neB));
    fprintf('ΔRe(n)=%.3e  ΔΓ=%.3e  cond(M)=%.2f  根后验 max ratio=%.2e\n', ...
        abs(neB-neA), abs(GB-GA), cond([aux.M(1) aux.M(2); aux.M(3) aux.M(4)]), ...
        max([ckA.ratio, ckB.ratio]));

    % ── 判定 ──
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
    % ## fprintf('第二段：解析参考根处 |f| 最大 %.3e（应≈0，否则参考解不是根）\n', max([t2.res_a]));
    fprintf('第二段：解析参考根处 |f| 最大 %.3e（未归一化，不能单独据此判根；判根看下两行 ne_t 与 ne_a 的逐位比较）\n', max([t2.res_a]));
    fprintf('第二段 max |Im κd| = %.3f，失真起点 %s\n', max([t2.kdi]), ...
        ternary(isinf(onset_kd), '未出现（全段 dα/α ≤ 1e-6）', sprintf('|Im κd| ≈ %.3f', onset_kd)));

    out = struct();
    out.table1 = t1;  out.table2 = t2;
    out.onset_alpha = onset1;
    out.onset_kd    = onset_kd;
    out.max_rel_lo  = max(rel_lo);
    out.max_rel_all = max([[t1.rel], [t2.rel]]);
    out.max_kdi     = max([t2.kdi]);
    out.n_root_fail = sum(~[t1.found]) + sum(~[t2.found]);
    out.pass_target = (max(rel_lo) < 1e-6);
    out.onset       = onset1;
    out.max_rel_upto_100 = max(rel_lo);
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

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_s7_poly（四阶多项式取值，Horner 形式）
% 输入：a：5×1 系数（a(1)=a0 … a(5)=a4）; x：自变量
% 输出：v：多项式值
% 说明：
%   1.Table I 的 f(x)=Σ_{n=0..4} a_n x^n，x 取 λ[nm]；
%   2.用 Horner 形式，避免 x^4 在 680 nm 处的舍入放大。
% ───────────────────────────────────────────────────────────────────────────────────
function v = stage2_s7_poly(a, x)
    v = a(5);
    for k = 4:-1:1
        v = v*x + a(k);
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_s7_table1（Kioupakis 2010 Table I 的四阶多项式系数）
% 输入：无
% 输出：T：结构体
%       T.ph/T.de/T.al 各含 e_perp/e_par/h_perp/h_par 五元列向量
% 说明：
%   1.来源：Kioupakis, Rinke, Van de Walle, Appl. Phys. Express 3, 082101 (2010),
%     DOI 10.1143/APEX.3.082101, Table I（380–680 nm，300 K）；
%   2.perp = 光偏振 ⊥ c 轴，par = ∥ c 轴；e/h = 电子/空穴；
%   3.本模型 TE 模的 E 沿面内 x（⊥ c）⇒ 用 perp 列；TM 模的 E 主要沿 y（∥ c）⇒ 用 par 列；
%   4.系数按 PDF 排版逐格抄录（out\kioupakis2010_table1*.png 是抄录时的对照图），
%     抄写正确性由 stage2_s7_selftest 反查 Table II 校验，不靠"看着像"。
% ───────────────────────────────────────────────────────────────────────────────────
function T = stage2_s7_table1()
    T.ph.e_perp = [-2.629e1;  1.999e-1; -5.406e-4;  6.294e-7; -2.673e-10];
    T.ph.e_par  = [-2.858e1;  2.126e-1; -5.670e-4;  6.546e-7; -2.764e-10];
    T.ph.h_perp = [ 5.548;   -5.631e-2;  2.035e-4; -3.060e-7;  1.794e-10];
    T.ph.h_par  = [ 6.052e-1; -6.675e-3; 2.317e-5; -2.841e-8;  1.652e-11];

    T.de.e_perp = [-2.202e-3; 1.503e-5; -3.572e-8;  3.506e-11; -7.916e-15];
    T.de.e_par  = [-2.354e-3; 1.660e-5; -4.202e-8;  4.528e-11; -1.308e-14];
    T.de.h_perp = [ 5.059e-3; -3.991e-5; 1.011e-7; -8.548e-11;  4.357e-14];
    T.de.h_par  = [-1.007e-3; 9.129e-6; -3.153e-8;  4.853e-11; -2.458e-14];

    T.al.e_perp = [-7.675e1;  6.598e-1; -2.129e-3;  3.034e-6; -1.584e-9];
    T.al.e_par  = [ 1.752e1; -1.179e-1;  2.511e-4; -1.614e-7;  0];
    T.al.h_perp = [-5.347e1;  4.498e-1; -1.391e-3;  1.875e-6; -9.161e-10];
    T.al.h_par  = [-4.851e1;  3.727e-1; -1.061e-3;  1.355e-6; -6.400e-10];
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_s7_sigma（自由载流子吸收截面 [cm²]）
% 输入：lambda：波长 [μm]; pol：'TE'（E⊥c）或 'TM'（E∥c）;
%       carrier：'e' 或 'h'; xln：该层 In 组分（0 = GaN）
% 输出：s：σ = 10^-18·[f_el-ph + f_alloy·x(1−x)/(0.25·0.75)] [cm²]
% 说明：
%   1.charged-defect 通道记 0：原文式 (2) 要乘 Z²·n_I/10^18，缺陷密度未知，
%     且原文结论是该通道比声子辅助弱两个数量级以上——没有 n_I 就不假装知道；
%   2.alloy 通道按原文式 (3) 用 x(1−x)/[0.25(1−0.25)] 缩放到本层组分。
% ───────────────────────────────────────────────────────────────────────────────────
function s = stage2_s7_sigma(lambda, pol, carrier, xln)
    T  = stage2_s7_table1();
    nm = lambda * 1e3;
    if strcmpi(pol, 'TE'), sfx = 'perp'; else, sfx = 'par'; end
    key = sprintf('%s_%s', lower(carrier), sfx);
    fph = stage2_s7_poly(T.ph.(key), nm);
    fal = stage2_s7_poly(T.al.(key), nm) * (xln*(1-xln)) / (0.25*0.75);
    s   = 1e-18 * (fph + fal);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_s7_alpha_table（各层吸收系数表）
% 输入：lambda：波长 [μm]; pol; scenario：'min' | 'kio' | 'zero'
% 输出：M：containers.Map（层名 → α [cm⁻¹]）; 
%      tab：struct 数组（tag/N/carrier/x/alpha）
% 说明：
%   1.载流子浓度不是阶段 2 能算的量（属阶段 4/5），这里取 Kioupakis 2010 Table II 的量级，
%     按本结构层名映射：p 侧 3e19 / 1e19 / 2e19，n 侧 2e18，QW 2e19，衬底 0；这是声明过的假设；
%   2.'min' 场景用施工单给的两个量级值（EBL 36、p 包层 18 cm⁻¹），其余 0；
%   3.'kio' 场景用 Table I 参数化：α = N·σ(λ,pol,x)。
% ───────────────────────────────────────────────────────────────────────────────────
function [M, tab] = stage2_s7_alpha_table(lambda, pol, scenario)
    tags = {'substrate','n-clad','n-wave','QW','QW1','barrier','QW2','p-wave','EBL','p-clad'};
    Nn   = [ 0,          2e18,    2e18,    2e19, 2e19,  2e19,    2e19,  1e19,    2e19,  3e19];
    car  = {'e','e','e','h','h','h','h','h','h','h'};
    xln  = [ 0,          0,       0.04,    0.15, 0.15,  0,       0.15,  0.04,    0,     0];
    al   = zeros(1, numel(tags));
    switch lower(scenario)
        case 'kio'
            for k = 1:numel(tags)
                al(k) = Nn(k) * stage2_s7_sigma(lambda, pol, car{k}, xln(k));
            end
        case 'min'
            al(strcmp(tags, 'EBL'))    = 36;
            al(strcmp(tags, 'p-clad')) = 18;
        case 'zero'
            % 全 0，仅用于取基准 P_layer
        otherwise
            error('stage2_s7_alpha_table: 未知场景 %s', scenario);
    end
    M   = containers.Map(tags, num2cell(al));
    tab = struct('tag', tags, 'N', num2cell(Nn), 'carrier', car, ...
                 'x', num2cell(xln), 'alpha', num2cell(al));
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_s7_selftest（用 Table II 反查 Table I 系数抄录是否正确）
% 输入：无
% 输出：out：结构体
%       out.含逐层预测/原文对照与最大相对偏差
% 说明：
%   1.这是"系数抄对了没有"的把手：Table II 给了 Γ_i、N_i、Γ_iα_i，反推 α_i = Γ_iα_i/Γ_i，
%     再与 N_i·σ 比较；GaN 层应落在百分之几以内；
%   2.QW 行单独标注：它是 InGaN 且有量子阱态，本文只按合金散射缩放，若偏差大就如实列出，
%     不当成抄录错误（避免用 QW 的偏差否定其他四行）。
% ───────────────────────────────────────────────────────────────────────────────────
function out = stage2_s7_selftest()
    lambda = 0.52;
    rows = { ...
        'p-cladding', 0.23, 3e19, 'h', 0.00, 9.5; ...
        'p-guiding',  0.16, 1e19, 'h', 0.10, 2.9; ...
        'EBL',        0.03, 2e19, 'h', 0.00, 0.8; ...
        'QW',         0.03, 2e19, 'h', 0.30, 2.7; ...
        'n-guiding',  0.19, 2e18, 'e', 0.10, 0.4; ...
        'n-cladding', 0.27, 2e18, 'e', 0.00, 0.2 };

    fprintf('\n===== S7 自校验：Table I 参数化 vs 原文 Table II（520 nm）=====\n');
    fprintf('%-12s %6s %10s %9s %9s %9s %9s %9s\n', ...
        '层', 'Γ', 'N[cm^-3]', 'α_原文', 'α_预测', 'Γα_原文', 'Γα_预测', '|ΔΓα|');
    % 判据说明：Table II 的 Γα 只印到 0.1 cm^-1，末位一位的舍入就有 ±0.05，
    % 所以"α 的相对偏差"在小 Γα 行会被放大（n-clad 的 Γα=0.2，0.01 的舍入就是 5%）。
    % 正确的把手是直接比 Γα，容差取半个末位 0.06；QW 行单独列出、不参与判据。
    dev = [];
    for k = 1:size(rows, 1)
        tag = rows{k,1}; G = rows{k,2}; Nc = rows{k,3}; cr = rows{k,4};
        xl  = rows{k,5}; GA = rows{k,6};
        a_ref  = GA / G;
        s      = stage2_s7_sigma(lambda, 'TE', cr, xl);
        a_pred = Nc * s;
        dGA    = abs(G*a_pred - GA);
        if ~strcmp(tag, 'QW'), dev(end+1) = dGA; end   %#ok<AGROW>
        fprintf('%-12s %6.2f %10.2e %9.2f %9.2f %9.2f %9.2f %9.3f\n', ...
            tag, G, Nc, a_ref, a_pred, GA, G*a_pred, dGA);
    end
    out.max_dga = max(dev);
    out.pass    = out.max_dga <= 0.06;
    fprintf('非 QW 行 max|ΔΓα| = %.3f cm^-1（判据 0.06 = 原文一个末位的舍入）⇒ pass=%d\n', ...
        out.max_dga, out.pass);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_s7_solve（给定层结构解 TE0，返回 n_eff/Γ/P_layer）
% 输入：L：层结构; lambda; pol; opts
% 输出：ne, G, Pl, ck, ok：ok=false 表示没找到根（不抛异常，交给调用方记 NaN）
% 调用函数：stage2_empty_check [用于数据记录 (→ 后验检查) 的空数组]
% 说明：
%   1.搜根窗口每次按当前 L 重算——扫厚度/扫衬底折射率时窗口会变；
%   2.n_lo 用 max(衬底, 覆盖层) + tiny，n_hi 用 max(n) − tiny（与主函数同口径）。
% ───────────────────────────────────────────────────────────────────────────────────
function [ne, G, Pl, ck, ok] = stage2_s7_solve(L, lambda, pol, opts)
    ne = NaN; G = NaN; Pl = []; ck = stage2_empty_check(); ok = false;
    N  = numel(L);
    n_max = max(real([L.n]));
    tiny  = 1e-9 * max(1, n_max);
    n_lo  = max(real(L(1).n), real(L(N).n)) + tiny;
    n_hi  = n_max - tiny;
    if ~(n_hi > n_lo), return; end
    fh = @(x) dispersion_tmm(x, L, lambda, pol);
    try
        [ne, G, Pl, ck] = stage2_s5_solve_te0(fh, L, lambda, pol, opts, layer_bounds(L), n_lo, n_hi);
        ok = true;
    catch
        ok = false;
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_s7_channels（S7：四通道分解 + 厚度扫描 + 衬底折射率扫描）
% 输入：无（波长走 stage2_lambda）
% 输出：out：结构体
%       out.含分解表、扫描表与判据值
% 调用函数：stack_7layer，layer_bounds，stage2_s7_alpha_table，stage2_s7_solve
% 说明：
%   1.四通道：n 侧（substrate+n-clad+n-wave）/ 有源区（QW）/ p 侧（p-wave+EBL+p-clad）/ 接触金属（本结构无）；
%   2.分解给两种口径：字面 ΣP_jα_j 与严格 ΣP_jα_j·n_rj/Re(n_eff)（S5 结论）；
%   3.厚度扫描：n-wave 与 p-wave 同步加厚，看 α_i(t) 与 Γ(t) 的折衷（B5）；
%   4.衬底扫描：n_sub 从低于包层扫到高于包层，判 n_eff 是否越过 n_sub（B6）。
% ───────────────────────────────────────────────────────────────────────────────────
function out = stage2_s7_channels(lambda)
    % 函数：stage2_s7_channels（S7 四通道损耗分解 + 厚度扫描 + 衬底扫描）
    % 输入：lambda：波长 [μm]；省略时取 stage2_lambda()
    % 输出：out：decomp（各场景 α_mode/通道/P/α）、thick（厚度扫描）、
    %       substrate（衬底扫描）、nsub_th（束缚模截止的 n_sub）
    % 说明：
    %   场景 min = EBL 36 + p-clad 18 cm^-1（其余 0）；kio = Kioupakis 2010 Table I 参数化。
    %   注：n_of_lambda 目前不随 λ 色散（405/520 无可靠来源），随 λ 变的只有 k0 与 α。
    global STAGE2_ALPHA
    if nargin < 1 || isempty(lambda), lambda = stage2_lambda(); end
    k0     = 2*pi/lambda;
    pol    = 'TE';
    opts   = struct('num_scan', 2e4, 'tol', 1e-13, 'maxit', 200, 'merge', 1e-6, ...
                    'npts', 800, 'tail_dec', 20, 'tail_max', 50, 'num_ff', 2001, ...
                    'chk_ratio', 0.5, 'chk_h', 10.^-(3:8));

    % ── 基准：零吸收 ──
    STAGE2_ALPHA = [];
    L0 = stack_7layer(lambda);
    [ne0, G0, P0, ck0, ok0] = stage2_s7_solve(L0, lambda, pol, opts);
    assert(ok0, 'stage2_s7_channels: 基准未找到根');
    nr0  = real(ne0);
    tags = {L0.tag};

    fprintf('\n===== S7 基准（无吸收，%.0f nm，TE0）=====\n', lambda*1000);
    fprintf('Re(n_eff)=%.12f  Γ=%.9f  ΣP−1=%.2e  根后验 max ratio=%.2e\n', ...
        nr0, G0, sum(P0)-1, max([ck0.ratio]));

    scen = {'min', 'kio'};
    out.decomp = struct([]);
    for is = 1:numel(scen)
        [M, tab] = stage2_s7_alpha_table(lambda, pol, scen{is});
        STAGE2_ALPHA = M;
        L  = stack_7layer(lambda);                    % 改 α 后必须重建层结构（S5 教训）
        [ne, ~, Pl, ck, ok] = stage2_s7_solve(L, lambda, pol, opts);
        if ~ok, fprintf('场景 %s 未找到根\n', scen{is}); continue; end
        am    = 2*k0*imag(ne)*1e4;
        al    = cell2mat(values(M, tags));
        nr    = real(cellfun(@(t) n450_of_tag(t), tags));
        naive = sum(Pl .* al);
        stric = sum(Pl .* al .* nr) / real(ne);

        fprintf('\n===== S7 四通道分解：场景 %s =====\n', upper(scen{is}));
        fprintf('%-10s %8s %10s %10s %10s %8s\n', '层', 'P[%]', 'N[cm^-3]', 'α[cm^-1]', 'P·α', 'n_r/n_eff');
        % F7 修正：tab 是 10 层的全表，tags 是当前结构的子集（7 层没有 QW1/barrier/QW2），
        % 直接取 tab(j) 会错行（p-wave 会显示成 QW 的 N、p-clad 显示成 EBL 的 N）。
        % 必须按层名回查下标。
        for j = 1:numel(tags)
            kt = find(strcmp({tab.tag}, tags{j}), 1);
            fprintf('%-10s %8.4f %10.2e %10.3f %10.4f %8.4f\n', ...
                tags{j}, 100*Pl(j), tab(kt).N, al(j), Pl(j)*al(j), nr(j)/real(ne));
        end
        ch = struct('nside', 0, 'active', 0, 'pside', 0, 'contact', 0);
        ch.nside   = sum(Pl(ismember(tags, {'substrate','n-clad','n-wave'})) .* ...
                         al(ismember(tags, {'substrate','n-clad','n-wave'})));
        ch.active  = sum(Pl(strcmp(tags,'QW')) .* al(strcmp(tags,'QW')));
        ch.pside   = sum(Pl(ismember(tags, {'p-wave','EBL','p-clad'})) .* ...
                         al(ismember(tags, {'p-wave','EBL','p-clad'})));
        fprintf('通道：n 侧 %.4f / 有源区 %.4f / p 侧 %.4f / 接触 %.4f  [cm^-1]\n', ...
            ch.nside, ch.active, ch.pside, ch.contact);
        fprintf('α_mode(解) = %.6f   字面 ΣP·α = %.6f（偏差 %.2e）   严格 ΣP·α·n_r/n = %.6f（偏差 %.2e）\n', ...
            am, naive, abs(am/naive-1), stric, abs(am/stric-1));
        out.decomp(is).scenario = scen{is};
        out.decomp(is).alpha_mode = am;
        out.decomp(is).naive = naive;
        out.decomp(is).strict = stric;
        out.decomp(is).channels = ch;
        out.decomp(is).P = Pl;
        out.decomp(is).alpha = al;
    end

    % ── B5：厚度扫描（n-wave 与 p-wave 同步加厚），场景用 kio ──
    [M, ~] = stage2_s7_alpha_table(lambda, pol, 'kio');
    STAGE2_ALPHA = M;
    tlist = [0.04 0.06 0.08 0.10 0.12 0.15 0.18 0.22 0.26 0.30];
    sc = struct('t', {}, 'ne', {}, 'G', {}, 'am', {}, 'Gm500', {}, 'Gm1000', {});
    fprintf('\n===== B5 厚度扫描（n-wave = p-wave = t，场景 kio）=====\n');
    fprintf('%6s %14s %10s %12s %12s %12s\n', 't[μm]', 'Re(n_eff)', 'Γ', 'α_i[cm^-1]', 'Γ·500−α', 'Γ·1000−α');
    for it = 1:numel(tlist)
        STAGE2_ALPHA = M;
        L = stack_7layer(lambda);
        L(3).d = tlist(it);  L(5).d = tlist(it);
        [ne, G, ~, ~, ok] = stage2_s7_solve(L, lambda, pol, opts);
        if ~ok
            fprintf('%6.2f %14s %10s %12s\n', tlist(it), 'nofind', '-', '-');
            continue;
        end
        am = 2*k0*imag(ne)*1e4;
        sc(it).t = tlist(it); sc(it).ne = real(ne); sc(it).G = G; sc(it).am = am;
        sc(it).Gm500  = G*500  - am;
        sc(it).Gm1000 = G*1000 - am;
        fprintf('%6.2f %14.9f %10.6f %12.4f %12.2f %12.2f\n', ...
            tlist(it), real(ne), G, am, sc(it).Gm500, sc(it).Gm1000);
    end
    out.thick = sc;

    % ── B6：衬底折射率扫描 ──
    % 分两段：粗扫看趋势（0.02 步长），细扫定位"根消失"的那一格（0.002 步长）。
    % 注意：n_sub 一旦超过 Re(n_eff)，束缚模在实轴上就不存在了，本实现的搜索下界
    % n_lo=max(n_sub,n_clad) 会把根挤出窗口——这时的"无根"是窗口构造的结果，
    % 不能读成"模消失了"。真正的泄漏模要求复平面分支/PML，归阶段 3。
    nsub = [1.70:0.05:2.40, 2.400:0.002:2.440];
    ss = struct('nsub', {}, 'ne', {}, 'bound', {}, 'note', {});
    fprintf('\n===== B6 衬底折射率扫描（其余层不变，场景 kio）=====\n');
    fprintf('%7s %14s %8s %10s\n', 'n_sub', 'Re(n_eff)', '束缚?', '备注');
    for k = 1:numel(nsub)
        STAGE2_ALPHA = M;
        L = stack_7layer(lambda);
        L(1).n = complex(nsub(k), 0);
        [ne, ~, ~, ~, ok] = stage2_s7_solve(L, lambda, pol, opts);
        bound = false; note = '';
        if ok
            bound = real(ne) > nsub(k);
            if ~bound, note = 'n_eff<n_sub：衬底内为振荡解'; end
        else
            note = '窗口内无根';
        end
        ss(k).nsub = nsub(k); ss(k).ne = real(ne); ss(k).bound = bound; ss(k).note = note;
        fprintf('%7.3f %14.9f %8d %10s\n', nsub(k), real(ne), bound, note);
    end
    nb = find([ss.bound], 1, 'last');
    if isempty(nb)
        fprintf('阈值：扫描区间内未找到束缚根\n');
    else
        fprintf('阈值：n_sub 到 %.3f 仍是束缚模（Re(n_eff)=%.9f）', ...
            ss(nb).nsub, ss(nb).ne);
        if nb < numel(ss)
            fprintf('；%.3f 起窗口内无根 ⇒ 束缚模截止\n', ss(nb+1).nsub);
        else
            fprintf('；扫描区间内未截止\n');
        end
    end
    out.substrate = ss;
    out.nsub_th = ss(nb).nsub;
    STAGE2_ALPHA = [];
    STAGE2_ALPHA = [];
end



% ─────────────────────────────────────────────────────
% 函数：stage2_s8_three_lambda（同结构 405/450/520 nm 对比，只读趋势）
% 输入：无（波长表写死在函数体内）
% 输出：out(i)：lambda / ne / G / alpha(每层) / Pl / am / ch(通道)
% 调用函数：stage2_s7_alpha_table、stage2_s7_solve、stack_7layer、n450_of_tag
% 说明：
%   1) n(λ) 目前固定为 450 nm 的常数（405/520 无可靠色散来源），随 λ 变的只有
%      k0 与 Kioupakis Table I 给出的 α(λ)，所以 P_layer/Γ 的变化只来自 k0。
%   2) 因此 405/520 只能读趋势，不能当绝对 α_i 报出去（施工单 B8 纪律）。
% ─────────────────────────────────────────────────────
function out = stage2_s8_three_lambda()
    global STAGE2_ALPHA
    pol  = 'TE';
    opts = struct('num_scan', 2e4, 'tol', 1e-13, 'maxit', 200, 'merge', 1e-6, ...
                  'npts', 800, 'tail_dec', 20, 'tail_max', 50, 'num_ff', 2001, ...
                  'chk_ratio', 0.5, 'chk_h', 10.^-(3:8));
    lams = [0.405 0.450 0.520];
    tags = {'substrate','n-clad','n-wave','QW','p-wave','EBL','p-clad'};

    fprintf('\n===== S8 三波长对比（同结构，TE0）=====\n');
    fprintf('说明：n(λ) 固定为 450 nm 常数（405/520 无可靠色散来源）；α 用 Kioupakis Table I（380–680 nm）。\n');
    fprintf('      结论只作趋势，405/520 的绝对 α_i 不得外报。\n');

    out = struct('lambda', {}, 'ne', {}, 'G', {}, 'alpha', {}, 'Pl', {}, 'am', {}, 'ch', {});
    for i = 1:numel(lams)
        lam = lams(i);
        k0  = 2*pi/lam;
        [M, ~] = stage2_s7_alpha_table(lam, pol, 'kio');
        STAGE2_ALPHA = M;
        L = stack_7layer(lam);
        [ne, G, Pl, ~, ok] = stage2_s7_solve(L, lam, pol, opts);
        if ~ok
            fprintf('λ=%.0f nm：窗口内无根，跳过\n', lam*1000);
            continue
        end
        al    = cell2mat(values(M, tags));
        nr    = real(cellfun(@(t) n450_of_tag(t), tags));
        am    = 2*k0*imag(ne)*1e4;
        stric = sum(Pl .* al .* nr) / real(ne);

        ch = struct('nside', 0, 'active', 0, 'pside', 0);
        ch.nside  = sum(Pl(ismember(tags, {'substrate','n-clad','n-wave'})) .* ...
                        al(ismember(tags, {'substrate','n-clad','n-wave'})));
        ch.active = sum(Pl(strcmp(tags, 'QW')) .* al(strcmp(tags, 'QW')));
        ch.pside  = sum(Pl(ismember(tags, {'p-wave','EBL','p-clad'})) .* ...
                        al(ismember(tags, {'p-wave','EBL','p-clad'})));

        out(i).lambda = lam; out(i).ne = ne; out(i).G = G;
        out(i).alpha = al;  out(i).Pl = Pl; out(i).am = am; out(i).ch = ch;

        fprintf('\nλ=%3.0f nm  Re(n_eff)=%.9f  Γ=%.6f  α_mode=%.4f cm^-1（严格 %.4f，偏 %.1e）\n', ...
            lam*1000, real(ne), G, am, stric, abs(am/stric - 1));
        fprintf('  通道：n 侧 %.4f / 有源区 %.4f / p 侧 %.4f  [cm^-1]\n', ...
            ch.nside, ch.active, ch.pside);
        for j = 1:numel(tags)
            fprintf('    %-10s P=%8.4f%%   α=%9.3f cm^-1   P·α=%8.4f\n', ...
                tags{j}, 100*Pl(j), al(j), Pl(j)*al(j));
        end
    end
    STAGE2_ALPHA = [];
end
