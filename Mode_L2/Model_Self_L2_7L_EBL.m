function modes = Model_Self_L2_7L_EBL(action, opt)
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
%
%   —— 模型说明 ──────────────────────────────────────────────
%       1.加入损耗后模场几乎不变：α 只进折射率的虚部，量级 1e-5；
%         模场形状由实部台阶（Δn≈0.05）决定。
%       2.KIO是指文献'Determination of Internal Loss in Nitride 
%         Lasers from First Principles'的作者 Kioupakis，通过文献
%         中的四阶多项式，由层材料计算层吸收系数; MIN 是DS的一个最小
%         假设，仅设置了EBL和p包层
%       3.本阶段模型支持解算参数显式传入，兼容本地预设调用。检验脚本只通
%         过 Model_Self_L2_7L_EBL('solve', opt) 取数，不再直接改全局表。
%       4.本模型共求解两遍：先求解实折射率(即第一阶段)，再求解复折射率，以
%         一阶段结果为锚点相比较。因为损耗对模场的影响几乎不计，本模型产出
%         中，仅损耗-吸收系数 \alpha_i的数据源于第二遍求解，其余数据源自第
%         一遍求解。
%
%  ── 符号说明────────────────────────────────────────────────
%       !! 仍有疑惑的地方，或者需要后续注意的内容（严重）
%       &&& 需要再想想，当前无法解决或者难以弄清。
%       % ## 上一版代码内容，在更新后被注释
%       % !! 代码更新后，追加的注释。

    %% 初始化和动作判定
    global STAGE2_ALPHA

    clc;                                                        %清空命令行窗口
    close all;                                                  %关闭所有已打开的图窗
    % clear;                                                      %清除变量
        
    if nargin < 2 || isempty(opt)
        opt = struct();                                         % 无第二入参时给空结构体，旧调用全部兼容
    end

    % 模型选择（已有，或者自己再建；下述代码可自行添加或修改）
    % 没有输入或者输入为空，默认action为"show"。nargin为函数实际收到的输入参数数量。
    % 检查字符串中是否含有"save"和"9"，有则为真。转字符串转小写。
    if nargin < 1 || isempty(action)
        action = 'show'; 
    end
    do_save  = contains(lower(string(action)), 'save');
    use9     = contains(lower(string(action)), '9');
    usetest  = contains(lower(string(action)), 'test');

    % 根后验检查自测
    if contains(lower(string(action)), 'rootcheck')
        modes = stage2_rootcheck_selftest();
        return
    end

    % 吸收截面 σ [cm^2] 输出
    % strtrim(…), 去掉首尾空格
    if strcmpi(strtrim(string(action)), 'sigma')
        modes = stage2_sigma_api(opt);
        return
    end

    % !! 阶段2 B 步：参数化解算接口。
    % !! 用法：S = Model_Self_L2_7L_EBL('solve', struct('lambda',0.45,'pol','TE', ...
    % !!        'stack','7','scenario','kio'))  或 'alpha' 显式给层表。
    % !! 返回：每模一个结构体，含 ne / Gamma / P_layer / alpha / contrib / ch / ok。
    if contains(lower(string(action)), 'solve')
        modes = stage2_solve_api(opt);
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

    % 层结构记录，用于损耗重建
    if use9
        stack_name = '9';
    elseif usetest
        stack_name = 'test';
    else
        stack_name = '7';
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
        modes(j).yb       = yb;              % 各层的上下限
        modes(j).ff       = ff;              % 快轴远场（结构体—角轴/分布/FWHM/峰值角）
    end

    %% 复折射率求解 + 层损耗重建
    % 正虚部 = 损耗(重申)。
    % α_i = 2·k0·Im(n_eff)·1e4 = 1e4/Re(n_eff) * Σ_j P_j·α_j·n_rj
    % 总损耗求解以及层分解用零吸收解的 P_j（各层功率占比；视为微扰）
    % loss(j).复折射率架构X; j是模式
    %   loss(j).X.ai        有效折射率的虚部
    %   loss(j).X.lit       各层损耗(alay)*功率占比 → 求和总损耗
    %   loss(j).X.alay      各层损耗
    %   loss(j).X.aj        层分解，求各层对应模场分布的损耗系数
    %   loss(j).X.ch        层损耗总结为四通道，sub/n/act/p/other(contact,=0)
    % 为什么用实部求解的模式数？ &&&

    scen = {'kio', 'min'};
    sno  = {'KIO', 'MIN'};
    loss = struct([]);
    % 结构体数组的字段必须预先对齐：先用占位值把 kio/min 同时建出来，
    %&&& 需要同步修改下方的 loss.kio、loss.min
    % 否则第二轮给 loss(1) 加新字段时与 loss(2..n) 字段不一致会报错（D1 同类）
    one0 = struct('ai', NaN, 'lit', NaN, 'alay', [], 'aj', [], ...
                  'ch', struct('sub',NaN,'n',NaN,'act',NaN,'p',NaN,'contact',NaN));
    for j = 1:numel(modes)
        loss(j).kio = one0;
        loss(j).min = one0;
    end

    for is = 1:numel(scen)
        [Mx, tabx] = stage2_alpha_table(lambda, pol, scen{is});
        STAGE2_ALPHA = Mx;                                          % 获取吸收系数
        Lx  = stage2_build_stack(stack_name, lambda);               % 层结构重建(将吸收系数注入层结构)
        Sx  = stage2_solve_modes(Lx, lambda, pol, n_modes, opts);   % 复折射率模式求解
        alx = stage2_alpha_by_tag(Lx, tabx);                        % 根据层结构获取对齐的损耗数值，用于后续计算和赋值
                                                                
        nrx = real([Lx.n]);
        for j = 1:numel(modes)
            one = struct('ai', NaN, 'lit', NaN, 'alay', alx, 'aj', nan(1, numel(Lx)), ...
                         'ch', struct('sub',NaN,'n',NaN,'act',NaN,'p',NaN,'contact',NaN));
            if j <= numel(Sx) && Sx(j).ok
                Pl = modes(j).P_layer(:).';                         % 实部求解得到的各层功率占比
                w  = Pl .* nrx / real(modes(j).n_eff);              % 严格微扰权重
                one.aj  = w .* alx;                                 % 层分解，求各层对应模场分布的损耗系数
                one.ai  = 2*k0*imag(Sx(j).ne)*1e4;                  % 等效折射率虚部
                one.lit = sum(Pl .* alx);                           % 忽略折射率修正，各层损耗及功率占比求和总损耗
                                                                    % 与ai 比较，作为自检
                one.ch  = stage2_channels(one.aj, {Lx.tag});        % 层损耗总结为四通道
            end
            loss(j).(scen{is}) = one;
        end


        fprintf('\n── 损耗（场景 %s，%.0f nm，%s）──\n', sno{is}, lambda*1000, pol);
        for j = 1:numel(modes)
            fprintf('%s%d : α_i = %8.3f cm^-1   （字面 ΣP·α = %8.3f）   [n 侧 %6.3f / 有源 %6.3f / p 侧 %6.3f / 衬底 %6.3f]\n', ...
                pol, j-1, loss(j).(scen{is}).ai, loss(j).(scen{is}).lit, ...
                loss(j).(scen{is}).ch.n, loss(j).(scen{is}).ch.act, ...
                loss(j).(scen{is}).ch.p, loss(j).(scen{is}).ch.sub);
        end
    end

    for j = 1:numel(modes)
        modes(j).loss = loss(j);                                    % 将损耗挂到模型输出上
    end
    STAGE2_ALPHA = [];                                              % 清表，后面自检/图只走几何量
    
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

    %% 绘图1 —— 波导模场、远场
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
    rowH = 28;                                         % 单行高 [px]（FontSize 16 实测约 26，留 2 px 余量）
    gap  = 26;                                         % 两块表的间隙 [px]
    mar  = 14;                                         % 窗口内边距 [px]
    hT1  = (numel(modes) + 1) * rowH + 8;              % 表一高：数据行 + 表头
    hT2  = (N + 1)           * rowH + 8;               % 表二高：数据行 + 表头

    fig3 = figure('Visible', vis, 'Name', 'GaN LD - data tables', 'Color', 'w');
    fig3.InnerPosition = [30 60 w3, hT1 + hT2 + gap + 2*mar];    % 直接定客户区，免估标题栏

    % 表一：纵轴：模式;  横轴：n_eff, beta, Γ, theta_FWHM, 峰值角, 损耗,归一化自检数
    mCols = {'mode', 'n_eff', 'beta [1/um]', 'Gamma', ...
             'FWHM [deg]', 'peak [deg]', 'alpha_i KIO', 'alpha_i MIN', 'norm check'};
    mData = cell(numel(modes), numel(mCols));
    for k = 1:numel(modes)
        mData(k, :) = { sprintf('%s%d', pol, k-1), ...
                        sprintf('%.8f', real(modes(k).n_eff)), ...
                        sprintf('%.6f', real(modes(k).beta)), ...
                        sprintf('%.6f', modes(k).Gamma), ...
                        sprintf('%.2f',   modes(k).ff.fwhm*180/pi), ...
                        sprintf('%.2f',   modes(k).ff.th_pk*180/pi), ...
                        sprintf('%.3f',   modes(k).loss.kio.ai), ...
                        sprintf('%.3f',   modes(k).loss.min.ai), ...
                        sprintf('%.3e',   modes(k).norm_val) };
    end
    t1 = uitable(fig3, 'Data', mData, 'ColumnName', mCols, 'FontSize', 16, ...
                 'Units', 'pixels');


    % 表二： 纵轴：层;  横轴：厚度，折射率，各个模式的功率占比，损耗贡献
    lCols = [{'#', 'tag', 'd [um]', 'n'}, ...
             arrayfun(@(j) sprintf('P(%s%d) [%%]', pol, j-1), ...
                      1:numel(modes), 'UniformOutput', false), ...
             {'alpha KIO', 'alpha MIN'}];
    lData = cell(N, numel(lCols));
    for k = 1:N
        dstr = 'Inf';
        if isfinite(L(k).d), dstr = sprintf('%.4f', L(k).d); end
        lData(k, 1:4) = { sprintf('%d', k), L(k).tag, dstr, sprintf('%.4f', real(L(k).n)) };
        for j = 1:numel(modes)
            lData{k, 4+j} = sprintf('%.4f', 100*modes(j).P_layer(k));
        end
        % !! 阶段2 A 步：逐层 α_j
        lData{k, 4+numel(modes)+1} = sprintf('%.3f', modes(1).loss.kio.alay(k));
        lData{k, 4+numel(modes)+2} = sprintf('%.3f', modes(1).loss.min.alay(k));
    end
    t2 = uitable(fig3, 'Data', lData, 'ColumnName', lCols, 'FontSize', 16, ...
                 'Units', 'pixels');

    % 首次布局：列宽按窗口宽等分铺满、两块表上下排好；拉伸窗口时自动重排
    relayout_tables(fig3, t1, t2, mCols, lCols);
    fig3.SizeChangedFcn = @(s, ~) relayout_tables(s, t1, t2, mCols, lCols);


    %% 绘图2 —— 损耗分解（三维方块 + 堆叠条）
    % 仅看纵向数值，层深(即x-y平面的数值)无意义
    fig4 = figure('Visible', vis, 'Name', 'GaN LD - loss decomposition', 'Color', 'w', ...
                  'Position', [60 60 980 940]);
    S4 = struct();
    S4.modes = modes;
    S4.pol   = pol;
    S4.tag   = {L.tag};
    S4.d     = [L.d];
    S4.scen  = {'kio', 'min'};
    S4.sn    = {'KIO', 'MIN'};
    uicontrol(fig4, 'Style', 'text', 'String', 'mode', 'Units', 'pixels', ...
              'Position', [20 896 60 24], 'FontSize', 10, 'BackgroundColor', 'w');
    S4.pmMode = uicontrol(fig4, 'Style', 'popupmenu', 'Units', 'pixels', ...
              'Position', [78 896 130 26], 'FontSize', 10, ...
              'String', arrayfun(@(j) sprintf('%s%d', pol, j-1), 1:numel(modes), ...
                                 'UniformOutput', false), ...
              'Value', 1, 'Callback', @(s, ~) draw_loss3d(fig4));
    uicontrol(fig4, 'Style', 'text', 'String', 'scenario', 'Units', 'pixels', ...
              'Position', [240 896 70 24], 'FontSize', 10, 'BackgroundColor', 'w');
    S4.pmScen = uicontrol(fig4, 'Style', 'popupmenu', 'Units', 'pixels', ...
              'Position', [318 896 110 26], 'FontSize', 10, ...
              'String', {'KIO', 'MIN'}, 'Value', 1, ...
              'Callback', @(s, ~) draw_loss3d(fig4));
    S4.ax  = axes(fig4, 'Units', 'pixels', 'Position', [70 470 560 400]);
    S4.axb = axes(fig4, 'Units', 'pixels', 'Position', [85 70 540 330]);
    fig4.UserData = S4;
    draw_loss3d(fig4);

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
        info.cand = cand_all;                                       % 逐候选记录（被拒的也在）
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
%           detail.P_layer  各层功率(向量)占比
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
    rowH = 28;                                         % 单行高 [px]，与制表段同值
    gap  = 26;                                         % 两表间隙 [px]
    mar  = 14;                                         % 窗口内边距 [px]

    wp  = get(fig, 'InnerPosition');                   % 客户区当前几何
    w   = wp(3);
    h1  = (size(get(t1, 'Data'), 1) + 1) * rowH + 8;   % 表一高：数据行 + 表头
    h2  = (size(get(t2, 'Data'), 1) + 1) * rowH + 8;

    gut = 46;
    cw1 = col_widths(mCols, get(t1, 'Data'), w - 2*mar - gut, 16);
    cw2 = col_widths(lCols, get(t2, 'Data'), w - 2*mar - gut, 16);

    set(t1, 'ColumnWidth', num2cell(cw1), ...
            'Position',    [mar, 2*mar + h2 + gap, w - 2*mar, h1]);
    set(t2, 'ColumnWidth', num2cell(cw2), ...
            'Position',    [mar, mar,               w - 2*mar, h2]);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：col_widths（按内容算列宽，再等比缩放到可用宽度）
% 输入：cols：列名 cell; data：uitable 的 Data; avail：可用宽度 [px]; fs：字号 [pt]
% 输出：cw：各列宽 [px]
% 说明：
%   1.每列宽 = 最长内容字符数 × 0.72×字号 + 18 px 余量（FontSize 16 实测约 11 px/字符）；
%   2.始终按 avail/sum(need) 整体缩放——放大也缩，铺满整行；
%     若只放大不缩、再把余量塞给最后一列，最后一列会变超宽（第一版就是这样）；
%   3.舍入差给第一列，不塞给最后一列。
% ───────────────────────────────────────────────────────────────────────────────────
function cw = col_widths(cols, data, avail, fs)
    n    = numel(cols);
    need = zeros(1, n);
    for k = 1:n
        Lk = numel(cols{k});
        for r = 1:size(data, 1)
            v = data{r, k};
            if ischar(v), Lk = max(Lk, numel(v)); end
        end
        need(k) = Lk * 0.72*fs + 18;
    end
    need = max(need, 50);
    cw   = floor(need * (avail/sum(need)));
    cw(1) = cw(1) + (avail - sum(cw));
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：波长赋值
% 输入：无
% 输出：lambda [μm]
% 说明：
%   1.阶段2 规定 λ 只在这里改；改这一行即换波长，别处不许再出现 0.45。
%   2.!! 405/520 nm 的 n(λ) 与 α(λ) 在《阶段2 施工单》§九 里明确"无可引用来源"，
%     换波长时实部仍取 450 nm 表（n450_of_tag），结论只能写"仅限几何/趋势"。
%   3.本函数是为了减少后续波长的连续调用和赋值（到损耗那段再看看 &&&）
% ───────────────────────────────────────────────────────────────────────────────────
function lambda = stage2_lambda()
    lambda = 0.45;             
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：复折射率赋值
% 输入：lambda [μm]; tags：层名 cell（长度 = 层数）
% 输出：nv：与 tags 等长的复折射率行向量
% 调用函数：n450_of_tag(实折射率)，alpha_of_layer（折射率虚部）
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
% 输入：tag：层名(字符串); lambda [μm]
% 输出：a：该层吸收系数 [cm⁻¹]
% 说明：
%   1.由全局 containers.Map 进行赋值。
%   2.使用 Map，是因为 Struct 不允许使用'-'作为字段名。
%       STAGE2_ALPHA = containers.Map({'QW'}, {10});
%   3.containers.Map(keys,value)：key=层名，value=吸收系数 α [cm^-1]; 
%     空表示默认 0。键区分大小写
%   4.未使用 lambda，作为预留接口。
%   5.&&& 加层时，需要对 known 同步修改。
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
        a = STAGE2_ALPHA(tag);          % 扫描，注入损耗
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
    if ~isfield(opts, 'chk_res'),   opts.chk_res   = 1e-6; end   % 绝对残差门限（相对 fscale）

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

    % 不能只比 ratio，会放过"陡峭函数的非零局部极小"，必须再加一条绝对残差门限
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
% 函数：系数f的四阶多项式取值
% 输入：a：系数向量，a(1)=a0 … a(5)=a4; x：λ[nm]
% 输出：v：多项式的值
% 说明：
%   1.采取 Horner 形式: P(x)=a0+a1x+a2x^2+a3x^3+a4x^4
%                          =a0+x(a1+x(a2+x(a3+xa4)))
%   下面这个，避免放大误差——x^4 在 680 nm 处的舍入放大。
% ───────────────────────────────────────────────────────────────────────────────────
function v = stage2_poly(a, x)
    v = a(5);
    for k = 4:-1:1
        v = v*x + a(k);
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：Kioupakis 2010 Table I 的四阶多项式系数
% 输入：无
% 输出：T：结构体
%       T.ph/T.de/T.al 各含 e_perp/e_par/h_perp/h_par 五元列向量
% 说明：
%   1.DOI 10.1143/APEX.3.082101, Table I(380–680 nm，300 K) (2010)。
%   2.perp = 光偏振 ⊥ c 轴，par = ∥ c 轴；e/h = 电子/空穴；
%   3.本模型 TE 模 ⇒ perp; TM 模 ⇒ par；
%   4.抄写正确性由 test_stage2\test_s7_selftest.m 用原文 Table II 反查校验
% ───────────────────────────────────────────────────────────────────────────────────
function T = stage2_KIO_table1()
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
% 函数：吸收截面计算 [cm²]
% 输入：lambda：波长 [μm]; pol：'TE'（E⊥c轴）或 'TM'（E∥c轴）;
%       carrier：多子类型，'e' 或 'h'; xln：该层 In 组分（0 = GaN）
% 输出：s：吸收截面(σ [cm²])
% 说明：
%   1.!! 本模型主要考虑声子辅助自由载流子吸收和合金辅助吸收; 未考虑带电缺陷辅助
%       吸收，记为0，因为缺陷密度未知，且原文结论是该通道比声子辅助弱两个数量级以上。
%   2.alloy 通道按原文式 (3) 用 x(1−x)/[0.25(1−0.25)]。
% ───────────────────────────────────────────────────────────────────────────────────
function s = stage2_sigma(lambda, pol, carrier, xln)
    T  = stage2_KIO_table1();
    lambda_nm = lambda * 1e3;
    if strcmpi(pol, 'TE'), sfx = 'perp'; else, sfx = 'par'; end
    key = sprintf('%s_%s', lower(carrier), sfx);
    fph = stage2_poly(T.ph.(key), lambda_nm);
    fal = stage2_poly(T.al.(key), lambda_nm) * (xln*(1-xln)) / (0.25*0.75);
    s   = 1e-18 * (fph + fal);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：各层吸收系数表
% 输入：lambda：波长 [μm]; pol; scenario：采用的吸收系数定义'min' | 'kio' | 'zero'
% 输出：M：containers.Map（层名 → α [cm⁻¹]）; 
%      tab：struct 数组
%           tab.tag         层标签
%           tab.N           载流子浓度
%           tab.carrier     多子类型
%           tab.x           合金摩尔分数
%           tab.alpha       吸收系数
% 说明：
%   1.载流子浓度不是阶段 2 能算的量（属阶段 4/5），这里取 Kioupakis 2010 Table II 的量级，
%     按本结构层名映射：p 侧 3e19 / 1e19 / 2e19，n 侧 2e18，QW 2e19，衬底 0。
%   2.&&& 加入新的层名，需额外添加。
%   3.'min' 场景用施工单给的两个量级值（EBL 36、p 包层 18 cm⁻¹），其余 0；
%     'kio' 场景用 Table I 参数化：α = N·σ(λ,pol,x)。
% ───────────────────────────────────────────────────────────────────────────────────
function [M, tab] = stage2_alpha_table(lambda, pol, scenario)
    tags = {'substrate','n-clad','n-wave','QW','QW1','barrier','QW2','p-wave','EBL','p-clad'};
    Nn   = [ 0,          2e18,    2e18,    2e19, 2e19,  2e19,    2e19,  1e19,    2e19,  3e19];  % 载流子浓度
    car  = {'e','e','e','h','h','h','h','h','h','h'};                                           % 多子类型
    xln  = [ 0,          0,       0.04,    0.15, 0.15,  0,       0.15,  0.04,    0,     0];     % In的组分
    al   = zeros(1, numel(tags));
    switch lower(scenario)
        case 'kio'
            for k = 1:numel(tags)
                al(k) = Nn(k) * stage2_sigma(lambda, pol, car{k}, xln(k));
            end
        case 'min'
            al(strcmp(tags, 'EBL'))    = 36;
            al(strcmp(tags, 'p-clad')) = 18;
        case 'zero'
        otherwise
            error('stage2_alpha_table: 未知场景 %s', scenario);
    end
    M   = containers.Map(tags, num2cell(al));
    tab = struct('tag', tags, 'N', num2cell(Nn), 'carrier', car, ...
                 'x', num2cell(xln), 'alpha', num2cell(al));
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_solve_api（B 步：参数化解算接口，检验脚本的唯一入口）
% 输入：opt：结构体，字段全部可选
%         .lambda  波长 [μm]，默认 0.45
%         .pol     'TE' / 'TM'，默认 'TE'
%         .n_modes 求几个模，默认 3
%         .stack   '7' / '9' / 'test'，默认 '7'
%         .scenario 'kio' / 'min' / 'zero'，按 Kioupakis 参数化建 α 表
%         .alpha   显式 α 表：struct（层名→α [cm^-1]）或 containers.Map；给定时优先于 scenario
%         .thick   层厚覆盖：struct / containers.Map（层名→厚度 [μm]），用于 B5 的厚度扫描
%         .nsub    衬底折射率实部覆盖 [无量纲]，用于 B6 的衬底扫描；虚部保持原值
% 说明：thick / nsub 在层结构建好之后改写，再走一遍 build_layers 的断言（层厚必须有限正值）
% 输出：S：结构体数组，每模一个，字段
%         ne      复有效折射率
%         Gamma   限制因子（口径不变，同主流程）
%         P_layer 逐层功率占比
%         alpha   逐层吸收系数 [cm^-1]
%         contrib 逐层贡献 P·α·n_r/Re(ne) [cm^-1]，其和 = α_i
%         ch      四通道（sub/n/act/p/contact）
%         ok      该模是否解出：false 时数值字段为 NaN，但 S 仍有 n_modes 个元素
%      另附 S 之外的标量信息放在每个元素的 lambda/pol/tag 里
% 调用函数：stage2_alpha_table，stage2_build_stack，stage2_solve_modes，
%           stage2_alpha_by_tag，stage2_channels
% 说明：
%   1.全局表 STAGE2_ALPHA 在这里 set / clear，对调用方不再可见（退化为实现细节）；
%   2.每次调用都重建层结构（坑 H5），所以连续调用不同 α 表不会互相污染；
%   3.本函数只做"解 + 收集"，不做判定、不打印——判定留给 test_stage2 下的脚本。
% ───────────────────────────────────────────────────────────────────────────────────
function S = stage2_solve_api(opt)
    global STAGE2_ALPHA

    % 参数缺省
    lambda   = 0.45;
    pol      = 'TE';
    n_modes  = 3;
    stack    = '7';
    scenario = '';
    alpha    = [];
    if isfield(opt, 'lambda')   && ~isempty(opt.lambda),   lambda   = opt.lambda;   end
    if isfield(opt, 'pol')      && ~isempty(opt.pol),      pol      = opt.pol;      end
    if isfield(opt, 'n_modes')  && ~isempty(opt.n_modes),  n_modes  = opt.n_modes;  end
    if isfield(opt, 'stack')    && ~isempty(opt.stack),    stack    = opt.stack;    end
    if isfield(opt, 'scenario') && ~isempty(opt.scenario), scenario = opt.scenario; end
    if isfield(opt, 'alpha')    && ~isempty(opt.alpha),    alpha    = opt.alpha;    end
    thick = []; nsub = [];
    if isfield(opt, 'thick') && ~isempty(opt.thick), thick = opt.thick; end
    if isfield(opt, 'nsub')  && ~isempty(opt.nsub),  nsub  = opt.nsub;  end

    % α 表：显式给定时优先，否则按场景建
    if ~isempty(alpha)
        Mx = alpha;
        if isstruct(Mx)
            fn = fieldnames(Mx);
            Mm = containers.Map();
            for k = 1:numel(fn)
                Mm(fn{k}) = Mx.(fn{k});
            end
            Mx = Mm;
        end
        % 由 α 表反推逐层表（tag/alpha），供 stage2_alpha_by_tag 使用
        tg = keys(Mx); av = values(Mx);
        tabx = struct('tag', {}, 'alpha', {});
        for k = 1:numel(tg)
            tabx(k) = struct('tag', tg{k}, 'alpha', av{k});
        end
    elseif ~isempty(scenario)
        [Mx, tabx] = stage2_alpha_table(lambda, pol, scenario);
    else
        Mx = []; tabx = struct('tag', {}, 'alpha', {});
    end

    STAGE2_ALPHA = Mx;
    Lx  = stage2_build_stack(stack, lambda);

    % 层厚覆盖（B5）：按层名查下标，不能按下标给——7 层与 9 层的层序不同
    if ~isempty(thick)
        if isstruct(thick)
            tk0 = fieldnames(thick);
            tv0 = cellfun(@(k) thick.(k), tk0, 'UniformOutput', false);
            thick = containers.Map(tk0, tv0);
        end
        tk = keys(thick);
        for k = 1:numel(tk)
            idx = find(strcmp({Lx.tag}, tk{k}), 1);
            assert(~isempty(idx), 'stage2_solve_api: thick 里的层名 %s 不在当前结构里', tk{k});
            Lx(idx).d = thick(tk{k});
        end
        Lx = build_layers(Lx);      % 重跑断言：中间层厚度必须有限正值
    end

    % 衬底折射率覆盖（B6）：只换实部，虚部（吸收）保持不动
    if ~isempty(nsub)
        Lx(1).n = complex(nsub, imag(Lx(1).n));
        Lx = build_layers(Lx);
    end
    Sx  = stage2_solve_modes(Lx, lambda, pol, n_modes, opt);
    alx = stage2_alpha_by_tag(Lx, tabx);
    nrx = real([Lx.n]);
    STAGE2_ALPHA = [];

    S = struct('ne', {}, 'Gamma', {}, 'P_layer', {}, 'alpha', {}, ...
               'contrib', {}, 'ch', {}, 'ok', {}, 'lambda', {}, 'pol', {}, 'tag', {});

    % 求根失败（例如衬底折射率越过 Re(n_eff) 后窗口内无根）：仍返回 n_modes 个元素，
    % 只把 ok 置 false、数值置 NaN——调用方不必先判断"数组是不是空的"。
    if isempty(Sx)
        for j = 1:n_modes
            S(j) = struct('ne', complex(NaN, NaN), 'Gamma', NaN, ...
                          'P_layer', nan(1, numel(Lx)), 'alpha', alx, ...
                          'contrib', nan(1, numel(Lx)), ...
                          'ch', stage2_channels([], {Lx.tag}), 'ok', false, ...
                          'lambda', lambda, 'pol', pol, 'tag', {{Lx.tag}});
        end
        return
    end
    % 注意：'tag' 存的是 cell（层名表），必须用双层花括号包住，
    % 否则 struct(...) 会按 cell 的元素个数展开成 1×N 的结构体数组（D1 同类）。
    for j = 1:numel(Sx)
        one = struct('ne', Sx(j).ne, 'Gamma', Sx(j).G, 'P_layer', Sx(j).Pl, ...
                     'alpha', alx, 'contrib', nan(1, numel(Lx)), ...
                     'ch', stage2_channels([], {Lx.tag}), 'ok', Sx(j).ok, ...
                     'lambda', lambda, 'pol', pol, 'tag', {{Lx.tag}});
        if Sx(j).ok
            Pl = Sx(j).Pl(:).';
            one.contrib = Pl .* alx .* nrx / real(Sx(j).ne);
            one.ch      = stage2_channels(one.contrib, {Lx.tag});
        end
        S(j) = one;
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_sigma_api（'sigma' 动作：返回某一层吸收截面 σ [cm^2]）
% 输入：opt：结构体，字段全部可选
%         .lambda  波长 [μm]，默认 0.45
%         .pol     'TE'（E⊥c）/'TM'（E∥c），默认 'TE'
%         .carrier 'e' / 'h'，默认 'h'
%         .x       该层 In 组分（0 = GaN），默认 0
% 输出：c：结构体，字段 lambda/pol/carrier/x/elph/alloy/total，单位 cm^2
% 调用函数：stage2_sigma
% 说明：
%   1.这是模型自己的吸收模型出口，不是判定逻辑：判定（例如用原文 Table II 反查系数
%     抄录是否正确）写在 test_stage2 的脚本里；
%   2.defect 通道恒为 0（缺陷密度未知，见 stage2_sigma 的说明），故不单列。
% ───────────────────────────────────────────────────────────────────────────────────
function c = stage2_sigma_api(opt)
    lambda = 0.45; pol = 'TE'; carrier = 'h'; x = 0;
    if isfield(opt, 'lambda')  && ~isempty(opt.lambda),  lambda  = opt.lambda;  end
    if isfield(opt, 'pol')     && ~isempty(opt.pol),     pol     = opt.pol;     end
    if isfield(opt, 'carrier') && ~isempty(opt.carrier), carrier = opt.carrier; end
    if isfield(opt, 'x')       && ~isempty(opt.x),       x       = opt.x;       end

    T  = stage2_KIO_table1();
    nm = lambda * 1e3;
    if strcmpi(pol, 'TE'), sfx = 'perp'; else, sfx = 'par'; end
    key = sprintf('%s_%s', lower(carrier), sfx);
    elph  = 1e-18 * stage2_poly(T.ph.(key), nm);
    alloy = 1e-18 * stage2_poly(T.al.(key), nm) * (x*(1-x)) / (0.25*0.75);
    c = struct('lambda', lambda, 'pol', pol, 'carrier', carrier, 'x', x, ...
               'elph', elph, 'alloy', alloy, 'total', elph + alloy);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：根据已选预设结构(按名字)重建层结构
% 输入：name：'7' / '9' / 'test'; lambda：波长 [μm]
% 输出：L：build_layers 校验后的层结构（折射率已按 alpha_of_layer 带虚部）
% 调用函数：stack_7layer，stack_9layer_mqw，test_layer，
%          build_layers(校验层结构)
% 说明：
%   1.损耗重解时换的是全局表 STAGE2_ALPHA，必须重建层结构(即在求解中引入STAGE2_ALPHA);
%   2.&&& 堆栈选择与主函数 use9 / usetest 同源，避免两处判断各写一遍。(还没看明白)
% ───────────────────────────────────────────────────────────────────────────────────
function L = stage2_build_stack(name, lambda)
    switch name
        case '9'
            lay = stack_9layer_mqw(lambda);
        case 'test'
            lay = test_layer(lambda);
        otherwise
            lay = stack_7layer(lambda);
    end
    L = build_layers(lay);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：stage2_alpha_by_tag（把 α 表按层名对齐到层结构）
% 输入：L：层结构; tab：stage2_alpha_table 的逐层表（字段 tag/alpha）
% 输出：a：1×N 的逐层吸收系数 [cm^-1]
% 说明：
%   1.α 表按 10 个可能层名建，层结构只含其中 7 或 9 个 ⇒ 必须按名字查，不能用下标。
%   2.查不到（层不在表里）置 0，并在调用处由 ΣP·α 与 α_i 的比对暴露出来。
% ───────────────────────────────────────────────────────────────────────────────────
function a = stage2_alpha_by_tag(L, tab)
    a = zeros(1, numel(L));
    for k = 1:numel(L)
        idx = find(strcmp({tab.tag}, L(k).tag), 1);
        if ~isempty(idx)
            a(k) = tab(idx).alpha;
        end
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：将层损耗总结为四通道(衬底，n侧，有源区，p侧，其他(contact))
% 输入：aj：逐层贡献 P_j·α_j·n_rj/Re(n_eff) [cm^-1]; tags：层名 cell
% 输出：ch：结构体，字段 sub/n/act/p/contact，单位 cm^-1
% 说明：
%   1.通道口径与主输出的损耗分解一致：n 侧 = n-clad + n-wave；
%     p 侧 = p-wave + EBL + p-clad；有源区 = 全部 is_active 层并集；
%     衬底单列（它既不是包层也不是有源区）；其余归 contact（本结构为空）。
%   2.有源区把 barrier 也算进来，是为了让四通道之和 = α_i；
%     Γ 的口径不受影响（Γ 仍不含势垒，见 confinement_factor）。
% ───────────────────────────────────────────────────────────────────────────────────
function ch = stage2_channels(aj, tags)
    ch = struct('sub', 0, 'n', 0, 'act', 0, 'p', 0, 'contact', 0);
    if isempty(aj)
        return
    end
    is_sub = strcmp(tags, 'substrate');
    is_n   = ismember(tags, {'n-clad', 'n-wave'});
    is_act = ismember(tags, {'QW', 'QW1', 'barrier', 'QW2'});
    is_p   = ismember(tags, {'p-wave', 'EBL', 'p-clad'});
    is_oth = ~(is_sub | is_n | is_act | is_p);
    ch.sub     = sum(aj(is_sub));
    ch.n       = sum(aj(is_n));
    ch.act     = sum(aj(is_act));
    ch.p       = sum(aj(is_p));
    ch.contact = sum(aj(is_oth));
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：复折射率层结构模式求解
% 输入：L：层结构（折射率含虚部）; lambda; pol; n_modes：待求解模式数; opts
% 输出：S：结构体数组
%       S.ne    按实部降序的根(模式)    
%       S.G     限制因子
%       S.Pl    各层的功率
%       S.ck    用于复折射率求解时记录数据的空数组
%       S.ok    ok=false 表示该模没解出来
% 调用函数：find_roots_complex(求解色散零点)，mode_profile(模场重建)，
%          normalize_mode(模场归一化)，confinement_factor(求解限制因子)
% 说明：
%   1.与主流程同一套口径：按 有效折射率实部 降序排；
%   2.&&& Γ 与 P_layer 复用主流程函数，不另立算法；
%   3.窗口与主流程一致（n_lo 取两侧包层较高者），因此"窗口内无根"的含义也一致；
%   4.单个模失败不抛异常，只把 ok 置 false，由调用方决定怎么报。
%   5.&&& ck是干什么用的？
% ───────────────────────────────────────────────────────────────────────────────────
function S = stage2_solve_modes(L, lambda, pol, n_modes, opts)
    S = struct('ne', {}, 'G', {}, 'Pl', {}, 'ck', {}, 'ok', {});
    N     = numel(L);
    n_max = max(real([L.n]));
    tiny  = 1e-9 * max(1, n_max);
    n_lo  = max(real(L(1).n), real(L(N).n)) + tiny;
    n_hi  = n_max - tiny;
    if ~(n_hi > n_lo)
        return
    end
    yb = layer_bounds(L);
    fh = @(x) dispersion_tmm(x, L, lambda, pol);
    roots_n = find_roots_complex(fh, n_lo, n_hi, n_modes + 2, opts);
    if isempty(roots_n)
        return
    end
    % 按实部降序、截断
    [~, ord] = sort(real(roots_n(:).'), 'descend');
    roots_n  = roots_n(ord);
    if numel(roots_n) > n_modes
        roots_n = roots_n(1:n_modes);
    end

    for j = 1:numel(roots_n)
        s = struct('ne', roots_n(j), 'G', NaN, 'Pl', nan(1, N), ...
                   'ck', stage2_empty_check(), 'ok', false);
        try
            [y, E, n_y] = mode_profile(roots_n(j), L, lambda, pol, opts);
            [E, ~]      = normalize_mode(y, E, n_y, pol);
            [G, det]    = confinement_factor(y, E, n_y, yb, L, pol);
            s.G  = G;
            s.Pl = det.P_layer;
            s.ok = true;
        catch
            s.ok = false;
        end
        S(j) = s;
    end
end


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：draw_loss3d（按两个下拉框的状态重画三维方块与堆叠条）
% 输入：fig：图窗句柄（其 UserData 存 modes / tag / d / 下拉框与 axes 句柄）
% 输出：无（就地重画）
% 调用函数：box3
% 说明：
%   1.上格：方块高度 = 该层贡献 P_j·α_j·n_rj/Re(n_eff)，沿 z 累加，总高即 α_i；
%   2.下格：同一场景下所有模式的堆叠条，便于横向比 α_i；条顶标数值；
%   3.所有方块深度相同（只给立体感），不做厚度映射——三维图里按厚度画会让
%     QW（9 nm）与包层（1.2 μm）差 130 倍，QW 会压成一条线；
%   4.层名与数值放在右侧清单（annotation）：用数据坐标写会在 axes 边界被裁掉，
%     零贡献的层也照样列出，否则 MIN 场景看起来"少了几层"；
%   5.注意 sprintf 里要写 \alpha：单个 a 会被转成响铃字符（坑 H11）。
% ───────────────────────────────────────────────────────────────────────────────────
function draw_loss3d(fig)
    S  = fig.UserData;
    im = get(S.pmMode, 'Value');
    is = get(S.pmScen, 'Value');
    ax = S.ax;   axb = S.axb;
    Nl = numel(S.tag);
    cmap = parula(Nl);
    sc   = S.scen{is};

    % 上格：三维方块
    cla(ax); hold(ax, 'on'); box(ax, 'on');
    M  = S.modes(im).loss.(sc);
    aj = M.aj;
    dep = 0.35;  z0 = 0;
    for k = 1:Nl
        z1 = z0 + aj(k);
        if aj(k) > 0
            box3(ax, 0, 1, 0, dep, z0, z1, cmap(k, :));
        end
        z0 = z1;
    end
    xlim(ax, [0 1.35]); ylim(ax, [0 0.45]); zlim(ax, [0 max(z0, 1e-6)*1.12]);
    view(ax, [-40 24]);
    set(ax, 'XTick', [], 'YTick', []);
    zlabel(ax, 'P_j\alpha_j n_{r,j}/n_{eff}  [cm^{-1}]');
    xlabel(ax, 'schematic depth (not to scale)');
    title(ax, sprintf('%s%d   %s   \\alpha_i = %.3f cm^{-1}', ...
                      S.pol, im-1, S.sn{is}, M.ai));

    % 下格：同一场景下所有模式的堆叠条
    cla(axb); hold(axb, 'on'); grid(axb, 'on'); box(axb, 'on');
    nm = numel(S.modes);
    Y  = zeros(nm, Nl);
    for j = 1:nm
        Y(j, :) = S.modes(j).loss.(sc).aj;
    end
    hb = bar(axb, Y, 'stacked');
    for k = 1:Nl
        hb(k).FaceColor = cmap(k, :); hb(k).EdgeColor = 'none';
    end
    set(axb, 'XTick', 1:nm, 'XTickLabel', ...
        arrayfun(@(j) sprintf('%s%d', S.pol, j-1), 1:nm, 'UniformOutput', false));
    ylabel(axb, 'P_j\alpha_j n_{r,j}/n_{eff}  [cm^{-1}]');
    title(axb, sprintf('Loss decomposition   (all modes, %s;  bar height = \\alpha_i)', ...
                       S.sn{is}));
    ylim(axb, [0, 1.18*max(sum(Y, 2))]);
    for j = 1:nm
        text(axb, j, sum(Y(j, :)) + 0.03*max(Y(:)), ...
             sprintf('%.2f', S.modes(j).loss.(sc).ai), ...
             'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    legend(axb, S.tag, 'Location', 'eastoutside', 'FontSize', 8);
    hold(axb, 'off');

    % 右侧清单：逐层数值（零贡献层也列出）
    txt = [{'layer           cm^-1'}, ...
           arrayfun(@(k) sprintf('%-11s %7.3f', S.tag{k}, aj(k)), Nl:-1:1, ...
                    'UniformOutput', false), ...
           {'', sprintf('alpha_i      %7.3f', M.ai)}];
    delete(findall(fig, 'Tag', 'loss3d_list'));
    annotation(fig, 'textbox', [0.665 0.715 0.27 0.21], 'String', txt, ...
               'FontName', 'FixedWidth', 'FontSize', 9, 'EdgeColor', [0.75 0.75 0.75], ...
               'Tag', 'loss3d_list', 'Interpreter', 'none', 'FitBoxToText', 'off');
    hold(ax, 'off');
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：box3（画一个三维长方体，五个面 + 十二条棱）
% 输入：ax：坐标轴; x0/x1/y0/y1/z0/z1：范围; c：面色 [r g b]
% 输出：无
% 说明：
%   1.顶点按"下四点上四点"排列，面取 +x、+y、顶面、-x、-y 五个，不画底面；
%     底面朝下，视角在上方永远看不到，画了反而与下一块方块的顶面共面（z-fighting）；
%   2.十二条棱逐条用 plot3 画出：patch 只画它自己的五个面框，缺的棱会让方块
%     看起来像"没关上的抽屉"（用户 2026-10-22 指出）；远处的棱由深度缓冲遮挡。
% ───────────────────────────────────────────────────────────────────────────────────
function box3(ax, x0, x1, y0, y1, z0, z1, c)
    V = [x0 y0 z0; x1 y0 z0; x1 y1 z0; x0 y1 z0; ...
         x0 y0 z1; x1 y0 z1; x1 y1 z1; x0 y1 z1];
    F = [2 6 7 3; 3 7 8 4; 5 6 7 8; 1 5 8 4; 1 2 6 5];
    patch(ax, 'Vertices', V, 'Faces', F, 'FaceColor', c, ...
          'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.5);
    E = [1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5; 1 5; 2 6; 3 7; 4 8];
    for ie = 1:size(E, 1)
        plot3(ax, V(E(ie, :), 1), V(E(ie, :), 2), V(E(ie, :), 3), '-', ...
              'Color', [0.25 0.25 0.25], 'LineWidth', 0.5);
    end
end
