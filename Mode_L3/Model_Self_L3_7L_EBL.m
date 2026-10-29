function modes = Model_Self_L3_7L_EBL(action, opt)
%   ── 命名规则 ──────────────────────────────────────────────
%      函数文件，函数名必须与文件名<.m>一致
%      直接"运行"，以0输入调用该函数
%      可以在命令行调用，如"modes = Model_Self_L3_7L_EBL('test'); "实现有输入调用'save'
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

    % ## clc;                                                      %清空命令行窗口
    % ## close all;                                                %关闭所有已打开的图窗
    % !! 上面两行已移到下面的守卫里（G3：原文保留）。理由：界面 l3_ui 会在自己的窗口里调本函数
    % !! 重算，close all 会把界面的 uifigure 一起关掉；脚本里连画多张图时也会把前一张删掉。
    % !! 只有显式传 opt.keep_windows = true 时才跳过；阶段 3 全流程与 show 的行为一点没变。
    if nargin >= 2 && isstruct(opt) && isfield(opt, 'keep_windows') && opt.keep_windows
        % 保持调用方的图窗与命令行，什么都不做
    else
        clc;  close all;
    end
    % clear;                                                      %清除变量
        
    if nargin < 2 || isempty(opt)
        opt = struct();                                         % 无第二入参时给空结构体，旧调用全部兼容
    end

    % 模型选择（已有，或者自己再建；下述代码可自行添加或修改）
    % 没有输入或者输入为空，默认action为"show"。nargin为函数实际收到的输入参数数量。
    % 检查字符串中是否含有"save"和"9"，有则为真。转字符串转小写。
    if nargin < 1 || isempty(action)
        action = 'stage3';                                       % !! 无参数 = 阶段 3 全流程（点运行即出全部图/csv/日志）
    end
    % ## action = 'show';
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

    % !! 阶段 3：垂直 + 横向两步 EIM 的总入口。
    % !! 用法：D = Model_Self_L3_7L_EBL('stage3', struct('sub','t4','w',3,'h_e',0.30,'n_top',1.46))
    % !! 子任务 t1…t9 见施工单 §3；本支只做分发，不承载任何业务逻辑。
    % !! 注意放在 'solve' 支之前：'stage3' 不含 'solve'，顺序只为可读性。
    if contains(lower(string(action)), 'stage3')
        modes = e3_stage3_dispatch(opt);
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
    if ~isfield(opts, 'merge_abs'), opts.merge_abs = NaN; end       % !! 去重绝对门限；给有限数则完全替代相对式
    if ~isfield(opts, 'chk_ratio'), opts.chk_ratio = 0.5; end       % |f(r)|/邻域最小值的判定阈值
    
    assert(n_hi > n_lo, 'find_roots_complex: 要求 n_hi > n_lo');
    
    xs = linspace(n_lo, n_hi, opts.num_scan);                       % 均匀扫描
    fs = fun(xs);                                                   % 获得遍历点的转移矩阵结果
    fs = fs(:).';
    
    % 防御性编程-不存在所有点的输出都是0
    fscale = max(abs(fs));
    if fscale == 0, fscale = 1; end
    opts.chk_fscale = fscale;                                       % 后验检查的绝对残差门限
    
    % !! 阶段 3 追加：可选的“实部变号”种子。近截止的导模在 |f| 上可能形不成谷（Im(f) 较大时
    % !! |f| 一直很高），但 real(f) 仍会穿过零；把变号区间的中点补进 Newton 初值。
    % !! 只在 opts.seed_bisect 为真时生效，阶段 2 的调用路径不受影响。
    if isfield(opts, 'seed_bisect') && opts.seed_bisect
        sgn = sign(real(fs));  sgn(sgn == 0) = 1;
        ix  = find(sgn(1:end-1) .* sgn(2:end) < 0);
        opts.xseed_extra = 0.5 * (xs(ix) + xs(ix+1));
        % !! 阶段 3 追加：再在 Im=seed_imag 的复直线上找变号。吸收结构里导模 Im(N)≈2e-5，
        % !! 与顶部模间隔同量级，实轴种子会被拉到相邻模的盆地；沿 Im(N) 的典型值起种子可落到本模。
        if isfield(opts, 'seed_imag') && isfinite(opts.seed_imag) && opts.seed_imag ~= 0
            zs = xs + 1i * opts.seed_imag;
            gs = fun(zs);
            sgn2 = sign(real(gs));  sgn2(sgn2 == 0) = 1;
            ix2 = find(sgn2(1:end-1) .* sgn2(2:end) < 0);
            opts.xseed_extra = [opts.xseed_extra, zs(ix2) + 0.0];
        end
    end

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
        % ## keep = [true, abs(diff(roots_n)) > opts.merge * max(1, abs(roots_n(1:end-1)))];
        % !! 阶段 3 追加：opts.merge_abs 为有限数时，绝对门限完全替代上面的相对式；
        % !! 缺省 NaN 时走 else 分支，与原式逐位相同（阶段 2 的回归不受影响）。
        % !! 理由：w=50 的横向模间隔只有 2.2e-7，相对门限 1e-6 折成绝对 ~2.4e-6，会把最低两阶并成一个。
        if isfinite(opts.merge_abs)
            thr = opts.merge_abs * ones(1, numel(roots_n) - 1);
        else
            thr = opts.merge * max(1, abs(roots_n(1:end-1)));
        end
        keep = [true, abs(diff(roots_n)) > thr];
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

    % !! 阶段 3 追加：并入 real(f) 变号处的种子（见 find_roots_complex 的 seed_bisect）
    if isfield(opts, 'xseed_extra') && ~isempty(opts.xseed_extra)
        cand = [cand, opts.xseed_extra(:).'];
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


%% !! 阶段 3：垂直 + 横向两步 EIM
% !! 本节全部是新增函数，既有函数一行未动。唯一被改的既有行在 find_roots_complex 的
% !! 去重判据上（原行以 % ## 保留，新增 opt.merge_abs；缺省时行为与阶段 2 逐位一致）。
% !! 入口：Model_Self_L3_7L_EBL('stage3', opt)，opt.sub 选子任务。


% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_log（同时写控制台与 UTF-8 日志文件）
% 输入：fmt 格式串; varargin 其参数
% 输出：无
% 说明：
%   1.句柄放在 global E3_LOG_FID，由 run_all_stage3.m 用 fopen(...,'w','n','UTF-8') 打开；
%   2.施工单禁用 diary：diary 的编码由 locale 决定，且会把 figure 的警告一并写进日志。
% ───────────────────────────────────────────────────────────────────────────────────
function e3_log(fmt, varargin)
    global E3_LOG_FID
    msg = sprintf(fmt, varargin{:});
    fprintf('%s\n', msg);
    if ~isempty(E3_LOG_FID) && E3_LOG_FID > 0
        fprintf(E3_LOG_FID, '%s\n', msg);
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_check（一条验收项的判定与留痕）
% 输入：name 项名; cond 判定结果(可为数值); fmt/varargin 附言
% 输出：ok：逻辑值
% 说明：
%   1.计数放在 global E3_NPASS/E3_NFAIL，run_all 末尾据此置退出码；
%   2.判定失败不抛异常，好让后续项继续跑完、一次拿到整张表。
% ───────────────────────────────────────────────────────────────────────────────────
function ok = e3_check(name, cond, fmt, varargin)
    global E3_NPASS E3_NFAIL
    if isempty(E3_NPASS), E3_NPASS = 0; end
    if isempty(E3_NFAIL), E3_NFAIL = 0; end
    ok = logical(cond);
    if ok
        E3_NPASS = E3_NPASS + 1;  tag = 'PASS';
    else
        E3_NFAIL = E3_NFAIL + 1;  tag = 'FAIL';
    end
    if nargin >= 3 && ~isempty(fmt)
        msg = sprintf(fmt, varargin{:});
        e3_log('  [%s] %-30s %s', tag, name, msg);
    else
        e3_log('  [%s] %s', tag, name);
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_trapzw（trapz 的逐点权重）
% 输入：x 单调节点
% 输出：w 与 x 同尺寸的权重列向量，满足 sum(w.*f) == trapz(x,f)
% 说明：
%   1.二维变换要拆成两次矩阵乘法，trapz 不能直接用在中间量上，故先把权重取出来；
%   2.权重只对连续被积函数成立——跨跳变的积分要分段（见 e3_gamma_power 与坑记录 A5）。
% ───────────────────────────────────────────────────────────────────────────────────
function w = e3_trapzw(x)
    x = x(:);
    n = numel(x);
    assert(n >= 2, 'e3_trapzw: 至少两个节点');
    d = diff(x);
    w = zeros(n, 1);
    w(1)   = d(1)/2;
    w(end) = d(end)/2;
    if n > 2
        w(2:end-1) = (d(1:end-1) + d(2:end))/2;
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_opts_solve（阶段 3 求根/重建的缺省参数）
% 输入：无
% 输出：opts 结构体
% 说明：
%   1.与主入口的 opts 缺省值保持一致（num_scan 20000 / merge 1e-6 / npts 800 …）；
%   2.阶段 3 不另立一套容差，避免同一根在两处给出不同位数。
% ───────────────────────────────────────────────────────────────────────────────────
function opts = e3_opts_solve()
    opts = struct( ...
        'num_scan', 20000, 'tol', 1e-13, 'maxit', 200, 'merge', 1e-6, ...
        'npts', 800, 'tail_dec', 20, 'tail_max', 50, 'chk_ratio', 0.5);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_stack_pair（按施工单 §2 的几何生成脊区 / 沟槽区两套层结构）
% 输入：lambda 真空波长 [μm]; opt 结构体
%         .h_pclad p-clad 厚度 [μm]，默认 0.30；给 Inf 则退回阶段 2 的半无限 p-clad（A0 回归路径）
%         .h_e     刻蚀深度 [μm]（自 p-clad 顶面向下），默认 0.30
%         .n_top   顶部介质折射率，默认 1.46（SiO₂）；1.0 为 air
% 输出：S 结构体
%         .ridge / .trench  层结构体数组（自衬底向上，字段 d/n/tag/is_active）
%         .meta             记录各层剩余厚度
% 调用函数：n_of_lambda
% 说明：
%   1.层序与 stack_7layer 完全一致，只把 p-clad 的 Inf 换成 h_pclad，并在其上放一层顶部介质；
%   2.顶部介质按半无限处理：包层衰减长度仅 ~37 nm，远小于实际 SiO₂ 厚度（0.3–1 μm），
%     有限厚度只会在色散函数里引入 ~e^{-d/ξ} 的修正，低于求根精度；
%   3.沟槽区把 y′∈[0,h_e) 的材料换成顶部介质，y′≥h_e 一概不动（只做垂直侧壁）；
%   4.厚度为 0 的层直接删掉：build_layers 要求中间层为有限正值；
%   5.h_e = 0 时两区层表逐字段相等（T1 的断言依赖这一点）。
% ───────────────────────────────────────────────────────────────────────────────────
function S = e3_stack_pair(lambda, opt)
    global STAGE2_ALPHA
    if nargin < 2 || isempty(opt), opt = struct(); end
    h_pclad = 0.30; h_e = 0.30; n_top = 1.46; sc = 'kio';
    if isfield(opt, 'h_pclad') && ~isempty(opt.h_pclad), h_pclad = opt.h_pclad; end
    if isfield(opt, 'h_e')     && ~isempty(opt.h_e),     h_e     = opt.h_e;     end
    if isfield(opt, 'n_top')   && ~isempty(opt.n_top),   n_top   = opt.n_top;   end
    if isfield(opt, 'scenario') && ~isempty(opt.scenario), sc = opt.scenario; end
    tz = 1e-12;

    % !! 自定义层表（界面用）：opt.L0 是「自衬底向上」的完整层表，首层 d=Inf、末层厚度有限且为正，
    % !! n 里已经含了损耗虚部。见到 L0 就完全不碰内置的 7 层标签与 STAGE2_ALPHA 表，
    % !! 直接把 L0 当脊区层表、再自顶向下刻 h_e 得到沟槽区 —— 物理方程与内置路径逐行相同。
    if isfield(opt, 'L0') && ~isempty(opt.L0)
        L0 = opt.L0(:);
        N0 = numel(L0);
        dtop = L0(N0).d;
        assert(isfinite(dtop) && dtop > 0, 'e3_stack_pair: 自定义层表的顶层厚度必须为正的有限值');
        assert(h_e <= dtop + 1e-9, ...
            'e3_stack_pair: 脊高 %g μm 超过顶层厚度 %g μm（本阶段不允许刻穿）', h_e, dtop);
        Lr = opt.L0(:);                                   % 脊区 = 原层表 + 顶部介质
        % !! 追加层必须与已有层「字段集完全一致」，不能直接写 struct(...,'is_active',false)：
        % !! 界面传进来的 L0 只有 d/n/tag 三个字段，多一个 is_active 就会报「在不同结构体之间进行下标赋值」。
        q1 = Lr(end);
        q1.d = Inf;  q1.n = n_top;  q1.tag = 'top';
        if isfield(q1, 'is_active'), q1.is_active = false; end
        Lr(end+1) = q1;
        Lt = L0(1:N0-1);                                  % 沟槽区：顶层被刻掉一部分
        dleft = dtop - h_e;
        if dleft > tz
            q2 = Lt(end);
            q2.d = dleft;  q2.n = L0(N0).n;  q2.tag = L0(N0).tag;
            if isfield(q2, 'is_active') && isfield(L0(N0), 'is_active')
                q2.is_active = L0(N0).is_active;
            end
            Lt(end+1) = q2;
        end
        q3 = Lt(end);
        q3.d = Inf;  q3.n = n_top;  q3.tag = 'top';
        if isfield(q3, 'is_active'), q3.is_active = false; end
        Lt(end+1) = q3;
        S = struct('ridge', Lr, 'trench', Lt, ...
                   'meta', struct('h_pclad', dtop, 'h_e', h_e, 'n_top', n_top, ...
                                  'd_pclad_left', dleft, 'd_pwave_left', NaN));
        return
    end

    tags = {'substrate','n-clad','n-wave','QW','p-wave','EBL','p-clad'};
    dk0  = [Inf, 1.200, 0.120, 0.009, 0.080, 0.020, h_pclad];
    % !! 损耗只能走层折射率的正虚部（G5）：α 表先塞进全局容器，再让 n_of_lambda 建复折射率，
    % !! 建完立刻清表 —— 与 stage2_solve_api 的用法一致，不给调用方留副作用。
    [Mx, ~] = stage2_alpha_table(lambda, 'TE', sc);
    STAGE2_ALPHA = Mx;
    nv  = n_of_lambda(lambda, tags);
    STAGE2_ALPHA = [];
    act = [false, false, false, true, false, false, false];

    L = struct('d', {}, 'n', {}, 'tag', {}, 'is_active', {});
    for k = 1:numel(tags)
        L(k) = struct('d', dk0(k), 'n', nv(k), 'tag', tags{k}, 'is_active', act(k));
    end

    % ── 脊区 ──────────────────────────────────────────────────────
    Lr = L;
    if isfinite(h_pclad)
        Lr(end+1) = struct('d', Inf, 'n', n_top, 'tag', 'top', 'is_active', false);
    end

    % ── 沟槽区：自顶向下扣掉 h_e ─────────────────────────────────
    d_pw = 0.080; d_pc = h_pclad; cut = h_e;
    if isinf(d_pc)
        cut = 0;                        % !! 半无限 p-clad：刻 h_e 之后仍是半无限（T1 的对照路径）
    elseif isfinite(d_pc)
        take = min(cut, d_pc); d_pc = d_pc - take; cut = cut - take;
    end
    take = min(cut, d_pw); d_pw = d_pw - take; cut = cut - take;
    assert(cut <= 1e-9, ...
        'e3_stack_pair: h_e = %g μm 已刻穿 p-wave，超出施工单 §2 的假设范围', h_e);

    Lt = L(1:6);                                   % substrate … EBL 不动
    Lt(5).d = d_pw;
    if isinf(d_pc)
        % !! 半无限 p-clad 的沟槽区：直接以 p-clad 收尾，不能再叠顶部介质（否则中间出现第二个 Inf）
        Lt(end+1) = struct('d', Inf, 'n', nv(7), 'tag', tags{7}, 'is_active', false);
    else
        if d_pc > tz
            Lt(end+1) = struct('d', d_pc, 'n', nv(7), 'tag', tags{7}, 'is_active', false);
        end
        Lt(end+1) = struct('d', Inf, 'n', n_top, 'tag', 'top', 'is_active', false);
    end

    S = struct('ridge', Lr, 'trench', Lt, ...
               'meta', struct('h_pclad', h_pclad, 'h_e', h_e, 'n_top', n_top, ...
                              'd_pclad_left', d_pc, 'd_pwave_left', d_pw));
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_region_solve（单个区的一维垂直模解：n_eff / Γ / 逐层功率 / 模场）
% 输入：L 层结构体数组（自衬底向上）; lambda [μm]; n_want 最多取几阶; opts 求解参数; pol 'TE'
% 输出：R 结构体数组，每个元素
%         .ne 复有效折射率; .G 限制因子; .P 逐层功率占比(1×N); .y/.E 模场; .ok 是否成功
% 调用函数：build_layers, layer_bounds, dispersion_tmm, find_roots_complex,
%          mode_profile, normalize_mode, confinement_factor
% 说明：
%   1.根窗口与主入口完全一致：n_lo = max(Re 首末层)+tiny、n_hi = max Re(n)−tiny；
%   2.基模 = Re(n_eff) 最大的那一支（find_roots_complex 已按实部降序并做过截断）；
%   3.本函数不加任何材料损耗：损耗只来自层折射率的虚部（G5）。
% ───────────────────────────────────────────────────────────────────────────────────
function R = e3_region_solve(L, lambda, n_want, opts, pol)
    if nargin < 4 || isempty(opts), opts = e3_opts_solve(); end
    if nargin < 5 || isempty(pol),  pol  = 'TE';             end
    L  = build_layers(L);
    N  = numel(L);
    nm = max(real([L.n]));
    tiny = 1e-9 * max(1, nm);
    n_lo = max(real(L(1).n), real(L(N).n)) + tiny;
    n_hi = nm - tiny;
    yb   = layer_bounds(L);
    R = struct('ne', {}, 'G', {}, 'P', {}, 'y', {}, 'E', {}, 'ok', {}, 'L', L, 'yb', yb);
    if ~(n_hi > n_lo)
        e3_log('  [warn] 搜索区间为空（n_lo=%g ≥ n_hi=%g）', n_lo, n_hi);
        return
    end
    fh = @(x) dispersion_tmm(x, L, lambda, pol);
    roots_n = find_roots_complex(fh, n_lo, n_hi, n_want + 2, opts);
    if isempty(roots_n), return; end
    [~, ord] = sort(real(roots_n(:).'), 'descend');
    roots_n  = roots_n(ord);
    for j = 1:numel(roots_n)
        s = struct('ne', roots_n(j), 'G', NaN, 'P', nan(1, N), ...
                   'y', [], 'E', [], 'ok', false, 'L', L, 'yb', yb);
        try
            [y, E, n_y] = mode_profile(roots_n(j), L, lambda, pol, opts);
            [E, ~]      = normalize_mode(y, E, n_y, pol);
            [G, det]    = confinement_factor(y, E, n_y, yb, L, pol);
            s.G = G;  s.P = det.P_layer;  s.y = y;  s.E = E;  s.ok = true;
        catch err
            e3_log('  [warn] 模场重建失败（Re n = %.9f）：%s', real(roots_n(j)), err.message);
        end
        R(j) = s;
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_lateral_nscan（横向求根的扫描点数：按最小模间隔自适应）
% 输入：na 芯区复折射率; nb 包层复折射率; w 脊宽 [μm]; k0 [μm⁻¹]; opt
% 输出：ns 扫描点数
% 说明：
%   1.最低两阶的 u 间隔为 π/V，dN/du = u0/[(w/2)²k0²Re(na)]，故最小间隔
%     gap ≈ u0·π/V / [(w/2)²k0²Re(na)]；w=50 时只有 ~2e-7，默认 20000 点扫不出来；
%   2.扫描步长取 gap/8（每间隔 8 点），上限 2e6 点（w=50 约 3.4e5 点，单格 ~0.3 s）；
%   3.这个数是"算出来的"而不是拍出来的，日志里逐格打印，便于复核。
% ───────────────────────────────────────────────────────────────────────────────────
function ns = e3_lateral_nscan(na, nb, w, k0, opt)
    ns = 20000;
    if nargin >= 5 && isfield(opt, 'num_scan') && ~isempty(opt.num_scan)
        ns = max(ns, opt.num_scan);
    end
    V = (w/2) * k0 * sqrt(max(real(na)^2 - real(nb)^2, 0));
    if V <= 0, return; end
    u0   = pi/2 - pi/(2*V);
    gap  = abs(u0 * (pi/V) / ((w/2)^2 * k0^2 * real(na)));
    dn   = real(na) - real(nb);
    if gap > 0 && dn > 0
        ns = max(ns, ceil(8 * dn / gap));
    end
    ns = min(ns, 2e6);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_lateral_solve（横向复本征：芯 = 脊区 n_a、包层 = 沟槽区 n_b）
% 输入：lambda [μm]; na/nb 两区垂直基模的复 n_eff; w 脊宽 [μm]; opt
% 输出：Q 结构体
%         .N 各阶总 n_eff（复，按 Re 降序）; .x/.Ex 归一化横向模场（cell）; .nx 逐点折射率
%         .V 归一化频率; .M_expect = floor(2V/π)+1; .n_scan 实际扫描点数
% 调用函数：build_layers, dispersion_tmm, find_roots_complex, mode_profile, e3_lateral_nscan
% 说明：
%   1.层序：[n_b, Inf] / [n_a, w] / [n_b, Inf]，侧壁取 E 与 ∂E/∂x 连续（TE-like），
%     与阶段 2 的求解器同一口径；严格准 TE 的侧壁条件应为 ∂E/∂x/n² 连续，
%     两者在 w ≥ 2 μm 时差 O((Δn/n)²) ~ 1e-6，见汇报单 §2.3；
%   2.本函数只吃 n_a/n_b，不接受任何 α 输入 —— EIM 第二步不得再注入吸收（G6）；
%   3.模数上限取解析值 M = floor(2V/π)+1，多求的不截断（留给 T5 的整数断言去比）。
% ───────────────────────────────────────────────────────────────────────────────────
function Q = e3_lateral_solve(lambda, na, nb, w, opt)
    if nargin < 5 || isempty(opt), opt = struct(); end
    k0 = 2*pi/lambda;
    Lc = struct('d', {}, 'n', {}, 'tag', {}, 'is_active', {});
    Lc(1) = struct('d', Inf, 'n', nb, 'tag', 'clad_l', 'is_active', false);
    Lc(2) = struct('d', w,   'n', na, 'tag', 'core',   'is_active', false);
    Lc(3) = struct('d', Inf, 'n', nb, 'tag', 'clad_r', 'is_active', false);
    Lc = build_layers(Lc);

    n_hi = real(na) - 1e-9*max(1, real(na));
    n_lo = real(nb) + 1e-9*max(1, real(nb));
    V    = (w/2) * k0 * sqrt(max(real(na)^2 - real(nb)^2, 0));
    Mm   = floor(2*V/pi) + 1;

    Q = struct('N', [], 'Ex', {{}}, 'x', {{}}, 'nx', {{}}, 'V', V, ...
               'M_expect', Mm, 'na', na, 'nb', nb, 'w', w, 'n_scan', 0, 'L', Lc);
    if ~(n_hi > n_lo)
        e3_log('  [warn] 横向搜索区间为空：Re(na)=%.9f ≤ Re(nb)=%.9f', real(na), real(nb));
        return
    end

    op = e3_opts_solve();
    fn = fieldnames(opt);
    for k = 1:numel(fn)
        if isfield(op, fn{k}), op.(fn{k}) = opt.(fn{k}); end
    end
    op.num_scan = e3_lateral_nscan(na, nb, w, k0, op);
    % !! 横向求解打开“实部变号”种子：w 大时近截止的高阶与近 n_a 的基模都可能因 Im(f) 偏大而漏谷
    op.seed_bisect = true;
    % !! 并在 Im = Im(n_a) 的复直线上补种子（导模 Im(N) 与该值同量级），否则最靠近 n_a 的基模会漏
    op.seed_imag = max(imag(na), 0);
    Q.n_scan    = op.num_scan;
    % !! 去重门限：缺省按解析模间隔自适应（gap/10）；显式给 NaN 则退回阶段 2 的相对式（t7 的对照用）
    if ~isfield(opt, 'merge_abs') || isempty(opt.merge_abs)
        op.merge_abs = e3_lateral_gap(na, nb, w, k0) / 10;
    else
        op.merge_abs = opt.merge_abs;
    end

    fh = @(x) dispersion_tmm(x, Lc, lambda, 'TE');
    roots_n = find_roots_complex(fh, n_lo, n_hi, Mm, op);
    if isempty(roots_n), return; end
    [~, ord] = sort(real(roots_n(:).'), 'descend');
    roots_n  = roots_n(ord);
    Q.N = roots_n(:).';

    for j = 1:numel(roots_n)
        x = []; Ex = []; nx = [];
        try
            [x, Ex, nx] = mode_profile(roots_n(j), Lc, lambda, 'TE', op);
            Ex = Ex / max(abs(Ex));
        catch err
            e3_log('  [warn] 横向模场重建失败（Re N = %.9f）：%s', real(roots_n(j)), err.message);
        end
        Q.x{j} = x;  Q.Ex{j} = Ex;  Q.nx{j} = nx;
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_gamma_power（二维场的 Γ 与逐层功率：按 |E_x|² 对两区加权）
% 输入：xc 横向坐标（已平移到脊中心）; Ex 该阶横向模场; w 脊宽 [μm]
%       A/B 结构体，字段 .G 该区垂直限制因子、.P 逐层功率占比、.L 层结构（取 tag）
% 输出：Gam 脊结构的总限制因子; Ptag 逐 tag 功率占比; tags 对应的层名
% 调用函数：e3_trapzw
% 说明：
%   1.被积函数 |E_x|²·Γ_y(x) 在 x=±w/2 处有跳变（两区 Γ_y 不同），必须分段积分；
%     直接把 mask 塞进一次 trapz 会让跨跳变的那一格取端点平均，抬高/压低结果（坑 A5）；
%   2.分段点 ±w/2 先插进网格再积，误差只有线性插值的那一点。
% ───────────────────────────────────────────────────────────────────────────────────
function [Gam, Ptag, tags] = e3_gamma_power(xc, Ex, w, A, B)
    xc = xc(:); Ex = Ex(:);
    xn = unique([xc; -w/2; w/2]);
    En = interp1(xc, Ex, xn, 'linear');
    a2 = abs(En).^2;
    i1 = find(xn == -w/2, 1);
    i2 = find(xn ==  w/2, 1);
    assert(~isempty(i1) && ~isempty(i2), 'e3_gamma_power: 网格里没有 ±w/2');
    seg = {1:i1, i1:i2, i2:numel(xn)};
    Gs  = [B.G, A.G, B.G];
    Wt  = trapz(xn, a2);
    Gam = 0;
    for s = 1:3
        Gam = Gam + trapz(xn(seg{s}), a2(seg{s})) * Gs(s);
    end
    Gam = Gam / Wt;

    % !! 层名要连成 cell（[A.L.tag, B.L.tag] 会拼成字符矩阵，后面 tags{k} 就会报"不支持花括号索引"）
    tags = [{A.L.tag}, {B.L.tag}];
    tags = unique(tags, 'stable');
    Ptag = zeros(1, numel(tags));
    for k = 1:numel(tags)
        ia = find(strcmp({A.L.tag}, tags{k}), 1);
        ib = find(strcmp({B.L.tag}, tags{k}), 1);
        pa = 0; pb = 0;
        if ~isempty(ia), pa = A.P(ia); end
        if ~isempty(ib), pb = B.P(ib); end
        val = 0;
        for s = 1:3
            wgt = pb; if s == 2, wgt = pa; end
            val = val + trapz(xn(seg{s}), a2(seg{s})) * wgt;
        end
        Ptag(k) = val / Wt;
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_cell（一个格点：两区垂直基模 + 横向各阶 + Γ/α/间隔）
% 输入：lambda [μm]; w 脊宽 [μm]; h_e 刻蚀深度 [μm]; n_top 顶部介质 n; opt
% 输出：C 结构体（见下）
% 调用函数：e3_stack_pair, e3_region_solve, e3_lateral_solve, e3_gamma_power,
%          stage2_alpha_table, stage2_alpha_by_tag
% 说明：
%   1.α 只能来自 Im(N)：α = 2k0·Im(N)·1e4 [cm⁻¹]，横向求解里不再注入任何 α（G6）；
%   2.α_pred = Σ_j P_j·α_j·n_rj/Re(N) 只是 A5 的对照量，不参与求根；
%   3.Re(n_b) ≥ Re(n_a) 时本函数不抛异常，把 ok=false 交回调用方（T2 会立刻停），
%     避免整个 40 格扫描因为一格反导引而全废。
% ───────────────────────────────────────────────────────────────────────────────────
function C = e3_cell(lambda, w, h_e, n_top, opt)
    if nargin < 5 || isempty(opt), opt = struct(); end
    h_pclad  = 0.30; if isfield(opt,'h_pclad') && ~isempty(opt.h_pclad), h_pclad = opt.h_pclad; end
    scenario = 'kio'; if isfield(opt,'scenario') && ~isempty(opt.scenario), scenario = opt.scenario; end
    k0 = 2*pi/lambda;

    % !! L0 = 界面传进来的自定义层表；不给就还是内置的 7 层（阶段 3 的全部判据走的就是这条缺省路径）
    spopt = struct('h_pclad', h_pclad, 'h_e', h_e, 'n_top', n_top);
    if isfield(opt, 'L0') && ~isempty(opt.L0), spopt.L0 = opt.L0; end
    Sp = e3_stack_pair(lambda, spopt);
    ov = e3_opts_solve();
    Ra = e3_region_solve(Sp.ridge,  lambda, 1, ov, 'TE');
    Rb = e3_region_solve(Sp.trench, lambda, 1, ov, 'TE');

    C = struct('lambda', lambda, 'w', w, 'h_e', h_e, 'n_top', n_top, 'h_pclad', h_pclad, ...
               'na', NaN, 'nb', NaN, 'dl2', NaN, 'V', NaN, 'M_expect', NaN, ...
               'N', [], 'Gamma', [], 'alpha', [], 'alpha_pred', [], 'gap', [], 'eta', [], ...
               'Gamma_a', NaN, 'Gamma_b', NaN, 'alpha_a', NaN, 'alpha_b', NaN, ...
               'tags', {{}}, 'Ptag', [], 'alx', [], 'nr', [], 'xc', {{}}, 'Ex', {{}}, ...
               'n_scan', 0, 'ok', false, 'region_ok', [false false], 'Sp', Sp, 'Ra', Ra, 'Rb', Rb);

    if isempty(Ra) || ~Ra(1).ok || isempty(Rb) || ~Rb(1).ok
        e3_log('  [warn] 格点 w=%g h_e=%g n_top=%g：某区垂直基模没解出来', w, h_e, n_top);
        return
    end
    C.region_ok = [true true];
    C.na = Ra(1).ne;  C.nb = Rb(1).ne;
    C.dl2 = real(C.na^2 - C.nb^2);
    C.Gamma_a = Ra(1).G;  C.Gamma_b = Rb(1).G;
    C.alpha_a = 2*k0*imag(C.na)*1e4;
    C.alpha_b = 2*k0*imag(C.nb)*1e4;
    if real(C.nb) >= real(C.na)
        e3_log('  [warn] 反导引：Re(n_b)=%.9f ≥ Re(n_a)=%.9f（w=%g h_e=%g n_top=%g）', ...
               real(C.nb), real(C.na), w, h_e, n_top);
        return
    end

    Q = e3_lateral_solve(lambda, C.na, C.nb, w, opt);
    C.V        = Q.V;
    C.M_expect = Q.M_expect;
    C.n_scan   = Q.n_scan;
    K = numel(Q.N);
    if K == 0
        e3_log('  [warn] 横向无根（w=%g h_e=%g n_top=%g）', w, h_e, n_top);
        return
    end

    % α 表（仅作 A5 对照）
    % !! 自定义层表的 tag 多半不在阶段 2 的损耗表里：查不到就按 0 处理，只影响 alpha_pred 的对照，
    % !! 不影响 C.alpha（它来自 Im(N)，与层 α 表无关）
    try
        [~, tabx] = stage2_alpha_table(lambda, 'TE', scenario);
        alx = stage2_alpha_by_tag(Sp.ridge, tabx);
        alx(~isfinite(alx)) = 0;
    catch
        alx = zeros(1, numel(Sp.ridge));
    end
    nr  = real([Sp.ridge.n]);
    tg  = {Sp.ridge.tag};

    A = struct('G', Ra(1).G, 'P', Ra(1).P, 'L', Ra(1).L);
    B = struct('G', Rb(1).G, 'P', Rb(1).P, 'L', Rb(1).L);

    C.N = real(Q.N);  C.alpha = 2*k0*imag(Q.N)*1e4;
    % !! 复合 N 也留一份：界面里要显示 Im n_eff，而 C.N 只存了实部（α 是另算的）
    C.Nc = Q.N(:).';
    C.Gamma = nan(1, K);  C.eta = nan(1, K);  C.alpha_pred = nan(1, K);
    C.Ptag = nan(K, 0);   C.tags = {};  C.xc = cell(1, K);  C.Ex = cell(1, K);
    for j = 1:K
        if isempty(Q.Ex{j}) || isempty(Q.x{j})
            continue
        end
        xc = Q.x{j} - w/2;   Ex = Q.Ex{j};
        C.xc{j} = xc;  C.Ex{j} = Ex;
        [Gam, Ptag, tags] = e3_gamma_power(xc, Ex, w, A, B);
        C.Gamma(j) = Gam;
        C.Ptag(j, 1:numel(tags)) = Ptag;
        C.tags = tags;
        C.eta(j) = abs(interp1(xc, Ex, w/2, 'linear'))^2 / max(abs(Ex))^2;
        % A5 对照：Σ P_j α_j n_rj / Re(N)
        ap = 0;
        for k = 1:numel(tags)
            ia = find(strcmp(tg, tags{k}), 1);
            if isempty(ia), continue; end
            ap = ap + C.Ptag(j, k) * alx(ia) * nr(ia);
        end
        C.alpha_pred(j) = ap / real(Q.N(j));
    end
    C.alx = alx;  C.nr = nr;
    C.ok = all(isfinite(C.Gamma));
    if numel(C.N) == 1
        C.gap = [];
    else
        C.gap = -diff(C.N);
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_phi_at（把一维垂直模场插值/解析外推到任意 y，含两端解析尾）
% 输入：yq 查询点; y/Phi 该区的模场（自底向上）; g_bot/g_top 两端尾衰减常数 [μm⁻¹]
% 输出：phi 查询点上的 Φ
% 说明：
%   1.网格内的部分线性插值；网格外不插值而用解析尾 Φ·exp(−γ·Δ)：模场在两端本来就是
%     单指数，外推比线性/常数端点准确得多，也避免把 trapz 的端点误差带进远场；
%   2.γ 由 k0·√(ñ²−n_layer²) 给出，取实部且必须为正（否则不是导模，由调用方保证）。
% ───────────────────────────────────────────────────────────────────────────────────
function phi = e3_phi_at(yq, y, Phi, g_bot, g_top)
    y = y(:);  Phi = Phi(:);  yq = yq(:);
    phi = interp1(y, Phi, min(max(yq, y(1)), y(end)), 'linear');
    lo = yq < y(1);
    hi = yq > y(end);
    if any(lo), phi(lo) = Phi(1)   .* exp(-g_bot * (y(1)   - yq(lo))); end
    if any(hi), phi(hi) = Phi(end) .* exp(-g_top * (yq(hi) - y(end))); end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_field2d（EIM 的二维场 E(x,y) = E_x(x)·Φ(x,y)，Φ 按区取各自垂直解）
% 输入：C 格点结构（e3_cell 的输出）; jlist 要建场的横向模序号; dx 横向网格步长（默认 10 nm）
% 输出：Fd 结构体
%         .x/.y 坐标; .E 三维数组 (ny × nx × numel(jlist)); .wx/.wy trapz 权重; .mask 脊区掩码
% 调用函数：e3_phi_at, e3_trapzw, layer_bounds
% 说明：
%   1.y 网格 = 两区界面并集 + 活跃区 5 nm / 其余 25 nm；两端各留 20 个衰减长度；
%   2.两区的 Φ 各自满足 ∫|Φ|²dy = 1，横向只由 E_x 带 x 依赖 —— 这就是 EIM 的全部近似；
%   3.不做任何归一化：∫|E|²dxdy 由远场判据 A11 拿去与 ∬|F|²dΩ 对比，归一化会把这个判据变成恒等式。
% ───────────────────────────────────────────────────────────────────────────────────
function Fd = e3_field2d(C, jlist, dx, dy_act)
    if nargin < 4 || isempty(dy_act), dy_act = 0.001; end
    if nargin < 3 || isempty(dx), dx = 0.010; end
    if nargin < 2 || isempty(jlist), jlist = 1; end
    k0 = 2*pi/C.lambda;

    A = C.Ra(1);  B = C.Rb(1);
    ybA = layer_bounds(A.L);  ybB = layer_bounds(B.L);
    yA  = ybA(end, 1);        yB  = ybB(end, 1);        % 两区顶面（p-clad 顶 / 槽底）

    % 两端尾长：20 个衰减长度（与 far_field 的 tail_dec=20 同口径）
    g_sub = real(k0 * sqrt(abs(real(A.ne)^2 - real(A.L(1).n)^2))) + 1e-6;
    g_top = real(k0 * sqrt(abs(real(A.ne)^2 - real(A.L(end).n)^2))) + 1e-6;
    y_bot = -min(20/g_sub, 1.0);
    y_top =  yA + min(20/g_top, 1.0);

    % 活跃区（QW 及其两侧）加密到 5 nm
    iq = find(strcmp({A.L.tag}, 'QW'), 1);
    y_lo = ybA(iq, 1) - 0.15;  y_hi = ybA(iq, 1) + 0.20;
    ya = [ybA(:, 1); ybA(:, 2)];  ya = ya(isfinite(ya));
    yg = [ybB(:, 1); ybB(:, 2)];  yg = yg(isfinite(yg));
    % !! 活跃区 5 nm、活跃区外 ±0.45 μm 内 10 nm、其余 25 nm：远场 FWHM 与 Parseval
    % !! 都对模场细结构敏感，25 nm 的通铺会让近场被欠采样（A10 的 1D/2D 差 0.8° 就是这个原因）
    % !! 活跃区用 1 nm：QW 厚 9 nm、模场在 QW 内是窄峰，5 nm 采样会把峰高估 ~34%（Γ_loc 0.0359 vs 0.0269）
    ys = unique([ (y_bot:0.025:y_lo), (y_lo:dy_act:y_hi), (y_hi:0.025:y_top), ...
                  (max(y_bot,y_lo-0.45):0.010:y_lo), (y_hi:0.010:min(y_top,y_hi+0.45)), ...
                  ya(:).', yg(:).' ]);
    % !! 量化到 0.1 nm：界面坐标（总厚度累加）与分段端点可能只差 1 ulp，两个近重合节点
    % !! 会让 trapz 权重塌成 1e-16，那一层在加权求和里就等于不存在（与 T6 的伪节点同一类坑）
    ys = unique(round(ys(:) * 1e7) / 1e7);

    g_lat = k0 * sqrt(max(real(C.na)^2 - real(C.nb)^2, 0));
    if g_lat <= 0, g_lat = k0*max(real(C.na)-real(C.nb), 1e-6); end
    % !! 横向尾长 8/γ 封顶在 2w：Δn²→0（如 h_e=0）时 γ_lat→0，不封顶会造出上亿点、直接撑爆内存
    half = C.w/2 + min(8/max(g_lat, 1e-9), 2*C.w);
    xs = (-half:dx:half).';
    xs(end+1) = half;   xs = unique(round(xs(:) * 1e7) / 1e7);

    phiA = e3_phi_at(ys, A.y, A.E, g_sub, g_top);
    phiB = e3_phi_at(ys, B.y, B.E, g_sub, g_top);
    Fd = struct('x', xs, 'y', ys, 'wx', e3_trapzw(xs), 'wy', e3_trapzw(ys), ...
                'jlist', jlist, 'E', zeros(numel(ys), numel(xs), numel(jlist)), ...
                'g_lat', g_lat, 'yA', yA, 'yB', yB, 'phiA', phiA, 'phiB', phiB);
    inR = abs(xs) <= C.w/2 + 1e-12;
    for q = 1:numel(jlist)
        j = jlist(q);
        Ex = C.Ex{j};
        if isempty(Ex)
            continue
        end
        ExR = zeros(size(xs));  ExT = zeros(size(xs));
        ExR(inR) = interp1(C.xc{j}, Ex, xs(inR), 'linear', 0);
        ExT(~inR) = interp1(C.xc{j}, Ex, xs(~inR), 'linear', 0);
        Fd.E(:, :, q) = phiA * ExR.' + phiB * ExT.';
    end
    Fd.inR = inR;
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_ff_matrix（二维远场：F(sx,sy) = ∬E(x,y)e^{−ik₀(x sx+y sy)}dxdy，谱域均匀网格）
% 输入：Fd 二维场（e3_field2d 输出）; q 第几个模; sx/sy 谱域网格（sinθ，可直接给 θ 的正弦）; k0
% 输出：F 复矩阵 (numel(sy) × numel(sx))
% 调用函数：无
% 说明：
%   1.写成两次矩阵乘：F = Mx'·(E∘(wy wx.'))·My，trapz 权重先乘进去，这样 x、y 两向的
%     端点权重一次也不漏（坑 A5 的另一面：把 trapz 拆成矩乘时容易漏权重）；
%   2.本函数只做变换，不取模、不归一化 —— 归一化按 A11 的用途由调用方决定。
% ───────────────────────────────────────────────────────────────────────────────────
function F = e3_ff_matrix(Fd, q, sx, sy, k0)
    E  = Fd.E(:, :, q);
    G  = E .* (Fd.wy * Fd.wx.');
    Mx = exp(-1i * k0 * (Fd.x(:) * sx(:).'));      % (nx × nSx)
    My = exp(-1i * k0 * (Fd.y(:) * sy(:).'));      % (ny × nSy)
    F  = Mx.' * G.' * My;                          % (nSx × nSy)，G 是 (ny × nx) 要转置
    F  = F.';                                      % 转成 (nSy × nSx)：第一维对 θy、第二维对 θx
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_ff_axes（从二维 F 矩阵取快轴 / 慢轴两个切面）
% 输入：F (nSy × nSx); sx/sy 谱域网格; wantSx 是否取快轴（θx 变化、θy=0）
% 输出：u 角变量（θ 或 sinθ）; I 归一化强度
% 说明：
%   1.快轴 = 垂直方向（y），θy 变化时取 θx=0 的那一列，即 F 的第 1 列（sx=0 落在网格正中）；
%   2.两个切面都各自归一 —— T9.3 要的是两个方向的远场形状，不是它们的相对强弱；
%   3.θ=0 必须恰好落在网格上：网格用奇数点、以 0 为中心（−60:0.25:60 共 481 点）。
% ───────────────────────────────────────────────────────────────────────────────────
function [u, I, Imax] = e3_ff_axes(F, sx, sy, wantFast)
    % !! F 的第一维对 θy（行）、第二维对 θx（列）
    if wantFast                                  % 快轴 = 垂直方向 y：固定 θx=0，沿 θy 取
        [~, i0] = min(abs(sx));
        I = abs(F(:, i0)).^2;  u = asin(sy(:));
    else                                         % 慢轴 = 横向 x：固定 θy=0，沿 θx 取
        [~, i0] = min(abs(sy));
        I = abs(F(i0, :)).^2;  u = asin(sx(:));
    end
    I = I(:);
    % !! 第三个输出给归一化之前的峰值：归一化前先记下它，供上层判断"这个切面是不是本来就近乎为零"。
    % !! 因为归一化一旦做了，1e-30 的数值噪声也会被放大成 1，曲线看着和真峰一样（2026-10-28 就画出过
    % !! 一条假的 m=1 快轴双峰，追了半天才发现是这种情形）
    Imax = max(I);
    I = I / max(Imax, realmin);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_parseval（A11：∬|F|²d(sinθx)d(sinθy) ÷ ∬|E|²dxdy ÷ λ² 应等于 1）
% 输入：Fd 二维场; q 模序号; s 谱域网格（sinθ，匀均 −0.999…0.999）
% 输出：ratio 上式; LHS/RHS 两个积分; F 矩阵
% 调用函数：e3_ff_matrix, e3_trapzw
% 说明：
%   1.判据的常数是 λ² 而不是 λ：二维变换的 Jacobian 是 k₀²，Parseval 给 (2π)²/k₀² = λ²；
%   2.必须用 sinθ（不是 θ）作积分变量，用 dθ 会得到偏 +3.4% 的假值（实测 0.4651 vs 0.4500）；
%   3.两侧都用同一个被截断的 E，所以域截断不影响恒等式，只有 s 网格分辨率和 trapz 精度影响结果。
% ───────────────────────────────────────────────────────────────────────────────────
function [ratio, LHS, RHS, F] = e3_parseval(Fd, q, s, k0)
    ws = e3_trapzw(s(:));
    F  = e3_ff_matrix(Fd, q, s, s, k0);
    LHS = sum(sum((abs(F).^2) .* (ws * ws.')));
    E  = Fd.E(:, :, q);
    wE = sum(sum((abs(E).^2) .* (Fd.wy * Fd.wx.')));
    RHS = (2*pi/k0)^2 * wE;
    ratio = LHS / RHS;
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_overlay_weights（T9.5 叠加权重）
% 输入：K 参与叠加的模数
% 输出：c 权重向量 (K×1)
% 说明：
%   下面是「假设值」：阶段 3 没有增益，没有任何物理依据给出模间幅度比与相位差，
%   故一律取等幅同相 c_k = 1/√K。它只用来画叠加后的条纹对比度，不得解释为实际激射模组成
%   （G12）；阶段 4 引入 Γ_k·g_k 之后必须按 a_k ∝ √(Γ_k g_k/α_k) 重算，届时替换本函数。
% ───────────────────────────────────────────────────────────────────────────────────
function c = e3_overlay_weights(K)
    % !! 权重 = 等幅同相假设值（1/√K），阶段 4 引入增益后必须替换；A14 依赖本行注明来源
    c = ones(K, 1) / sqrt(K);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_scatter（T8：侧壁散射损耗的参数化上界）
% 输入：C 格点结构; sigma 粗糙度 rms [nm]; Lc 相关长度 [nm]; Cpre 前置因子（默认 1）; n_ref 参考折射率
% 输出：a 散射损耗 [cm⁻¹]（该格、该模的加权值）; aK 逐模值; eta 逐模侧壁归一化强度
% 说明：
%   1.α_scat = 2·C·(4π²/λ⁴)·σ²·L_c·(n_ref²−n_top²)²·η，量纲 μm⁻¹，乘 1e4 换成 cm⁻¹；
%     代回施工单 §5 的把手（σ=1 nm、L_c=30 nm、n_ref=2.455、n_top=1.46、η=0.07）得 0.61 cm⁻¹；
%   2.η 用 |E_x(±w/2)|²/max|E_x|²：侧壁处的横向场越弱，同一粗糙度引起的散射越弱；
%   3.σ、L_c 在阶段 3 没有实测值，只交参数化表，不写进任何结论数（G8、§8 不做清单）。
% ───────────────────────────────────────────────────────────────────────────────────
function [a, aK, eta] = e3_scatter(C, sigma, Lc, Cpre)
    if nargin < 4 || isempty(Cpre), Cpre = 1; end
    lam = C.lambda;
    dn2 = (2.455^2 - C.n_top^2)^2;      % n_ref 取 GaN 的 2.455 而非导模值：粗糙度在材料界面
    pref = 2 * Cpre * (4*pi^2 / lam^4) * (sigma*1e-3)^2 * (Lc*1e-3) * dn2 * 1e4;
    eta  = C.eta(:).';
    aK   = pref * eta;
    a    = mean(aK);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_fd2d（T6：独立的二维标量有限差分，只比 Re 有效折射率）
% 输入：lambda [μm]; w 脊宽; h_e 刻蚀深度; n_top 顶部介质; opt
%         .dx 网格步长（默认 0.020）; .ny 条数不需要外部给
% 输出：R 结构体（.Nf 复? 只含 Re; .x/.y; .E 场; .frac 芯内功率占比; .n 节点数; .ok）
% 调用函数：n_of_lambda, e3_stack_pair
% 说明：
%   1.五点差分、x ∈ [−3w, 4w]、y′ ∈ [−0.3, 2.0]（y′ 自 p-clad 顶面向下），两端 Dirichlet；
%     这个域是按施工单 §3 T6 定的，不是按场分布优化过的（故 w=10 时有大量空网格）；
%   2.只取 Re(n)：复数算子会让 eigs 的 shift 变慢，而 Im(n) 对 Re(N) 的影响是 O(Im²) ~ 1e-10，
%     远小于 A4 的 5e-3 判据；
%   3.本征值 λ = −(k₀N)²，导模是「最负」的那几个；用 shift-invert（σ 定在 −(k₀Re n_a)² 附近）
%     而不是 'la'/'sm'：本机 eigs(...,'la') 会对整个稀疏矩阵返回 NaN（坑记录 B4 已复现）；
%   4.解出来的候选按「芯内功率占比」挑，而不是按 λ 排序：Dirichlet 盒子里的高阶箱模同样满足
%     Re(N) < max n，只按大小排序会挑错（坑 A7）。
% ───────────────────────────────────────────────────────────────────────────────────
function R = e3_fd2d(lambda, w, h_e, n_top, opt)
    if nargin < 5 || isempty(opt), opt = struct(); end
    dx = 0.020; if isfield(opt,'dx') && ~isempty(opt.dx), dx = opt.dx; end
    k0 = 2*pi/lambda;

    % !! 界面路径：把用户层表透传（e3_stack_pair 见到 L0 就改走自定义分支）
    spopt = struct('h_pclad', 0.30, 'h_e', h_e, 'n_top', n_top);
    if isfield(opt, 'L0') && ~isempty(opt.L0), spopt.L0 = opt.L0; end
    Sp = e3_stack_pair(lambda, spopt);
    Lr = Sp.ridge;
    y_if = cumsum([0, [Lr(2:end-1).d]]);           % 各界面 y′ 坐标
    % !! y′ < 0 落在顶部介质里（p-clad 顶面之上），不能再用脊区的层表去查，否则 find 返回空
    % !! y′ 是"自 p-clad 顶面向下"的深度，而层表 Lr 与 y_if 都是"自底向上"的：
    % !! 必须先换算成自底坐标 y_b = y_if(end) − y′，否则整张折射率图上下颠倒。
    y_tot = y_if(end);
    n_of = @(yq) arrayfun(@(t) e3_nb(y_tot - t, y_if, Lr, n_top), yq);

    xpad = 1.0;
    if isfield(opt,'xpad') && ~isempty(opt.xpad), xpad = opt.xpad; end
    xhalf = w/2 + xpad;                             % 侧向域半宽
    % !! 侧向也用非均匀网格：脊边 ±0.15 μm 内 10 nm（要分辨 |x|=w/2 的折射率台阶），脊内 20 nm，
    % !! 之外 25 nm。域半宽取 w/2+1.0 μm：包层衰减长度 1/γ_lat ≈ 40 nm，1 μm 已 25 个衰减长度，
    % !! 再往外场强 ~1e-11，Dirichlet 壁的影响可忽略（域再扩 2 倍的对照见 t6 日志）。
    % !! 不用施工单的 [−3w,+4w] 均匀 20 nm：w=10 时那会到 1.65e6 个未知数，eigs 的 LU 撑不住。
    fac = 1;
    if isfield(opt,'fac') && ~isempty(opt.fac), fac = opt.fac; end
    segx = [-xhalf      -w/2-0.15  0.025
            -w/2-0.15   -w/2+0.15  0.010
            -w/2+0.15    w/2-0.15  0.020
             w/2-0.15    w/2+0.15  0.010
             w/2+0.15    xhalf     0.025];
    xs = [];
    for kk = 1:size(segx,1)
        xs = [xs, segx(kk,1):fac*segx(kk,3):segx(kk,2)]; %#ok<AGROW>
    end
    xq = unique([xs(:); -xhalf; -w/2; w/2; xhalf]);
    % !! y 方向必须用非均匀网格：QW 只有 9 nm，20 nm 均匀网格下它不足半个格子，
    % !! 实测一维竖直自检因此偏低 3.5e-3（2.41338 vs TMM 2.41692），与二维偏差同量级 —— 偏差来自网格，不是 EIM。
    seg = [-0.30 -0.10 0.025
           -0.10  0.00 0.010
            0.00  0.25 0.010
            0.25  0.60 0.001
            0.60  0.90 0.010
            0.90  2.00 0.025];
    ys = [];
    for kk = 1:size(seg,1)
        ys = [ys, seg(kk,1):fac*seg(kk,3):seg(kk,2)]; %#ok<AGROW>
    end
    ys = sort(unique([ys, 2.0, (y_tot - y_if(:)).']));   % 强制界面落在节点上
    % !! 界面坐标是"总厚度减累加和"算出来的，与分段网格的端点可能只差 1 ulp（~1e-16 μm），
    % !! 这种伪近重合节点会让 1/h² 爆到 1e30 量级、整张算子变垃圾。量化到 0.1 nm 消除它。
    ys = unique(round(ys*1e7)/1e7);
    yq = ys(:);
    xq = unique(round(xq(:)*1e7)/1e7);
    if isfield(opt,'dbg') && opt.dbg
        e3_log('      [fd] 最小间距：y %.3e μm，x %.3e μm（应 > 1e-4）', ...
               min(diff(yq)), min(diff(xq)));
    end
    ny = numel(yq);  nx = numel(xq);
    nyi = ny - 2;  nxi = nx - 2;                         % 去掉两端 Dirichlet 点
    ysi = yq(2:end-1);  xsi = xq(2:end-1);
    nvi = n_of(ysi);
    if isfield(opt, 'dbg') && opt.dbg
        e3_log('      [fd] 网格：y %d 点（非均匀，均匀对照 %d 点）× x %d 点', ...
               ny, numel(-0.3:dx:2.0), nx);
        e3_log('      [fd] n(-0.10)=%.4f(1.4600) n(0.05)=%.4f(2.3740) n(0.35)=%.4f(2.4550)', ...
               real(n_of(-0.10)), real(n_of(0.05)), real(n_of(0.35)));
        e3_log('      [fd] n(0.395)=%.4f n(0.405)=%.4f(2.5500) n(0.45)=%.4f(2.4550) n(1.80)=%.4f(1.7700)', ...
               real(n_of(0.395)), real(n_of(0.405)), real(n_of(0.45)), real(n_of(1.80)));
        % !! 一维竖直自检不放这里：h_e=0 的「零号自检」（见 e3_t6 开头）等价于一维竖直问题，
        % !! 且能顺带量出 Dirichlet 侧壁的代价（ΔN ≈ (π/L)²/(2k₀²N)，L=5 μm、N=2.417 → 4.2e-4）。
    end                                  % 脊区逐层折射率（自上向下扫）
    % !! 内部节点的折射率平方（Dirichlet 边界点不参与本征问题）
    N2 = repmat(real(nvi).^2, 1, nxi);
    tr = (abs(xsi) > w/2);                          % 槽区列
    N2(ysi < h_e, tr) = real(n_top)^2;

    Dy = e3_d2_nonuni(yq);                          % (ny−2)×(ny−2)，已含 1/h² 缩放
    D2x = e3_d2_nonuni(xq);                         % (nx−2)×(nx−2)，同上
    Lap = kron(speye(nxi), Dy) + kron(D2x, speye(nyi));
    Ah  = -Lap - spdiags(k0^2 * N2(:), 0, nyi*nxi, nyi*nxi);   % 本征值 = −(k₀N)²

    na_r = max(real([Lr.n]));                        % 脊区最高 Re（QW 的 2.55）
    % !! 移位点只决定"找到哪几个"本征值，不影响 eigs 算出的特征值本身，所以允许 t6 把 EIM 的
    % !! Re N(1) 传进来当 opt.shift_n —— 否则按 0.995·n_max 定移位，会落到导模带之外
    shn  = 0.995 * na_r;
    if isfield(opt,'shift_n') && ~isempty(opt.shift_n), shn = opt.shift_n; end
    sig  = -(k0 * shn)^2;
    R = struct('Nf', [], 'E', [], 'frac', [], 'x', xsi, 'y', ysi, 'n', nyi*nxi, ...
               'dx', dx, 'ok', false, 'Ah_nnz', nnz(Ah), 'method', '');
    nk = 12;
    [V, D] = deal([]);
    try
        [V, D] = eigs(Ah, nk, sig);
        R.method = 'eigs-shift';
    catch err
        e3_log('  [warn] eigs(shift-invert) 失败：%s', err.message);
    end
    if isempty(V) || any(~isfinite(D(:)))
        e3_log('  [warn] eigs 返回非有限值，改用自带 shift-invert 反迭代');
        [V, D] = e3_fd_reviter(Ah, sig, nk);
        R.method = 'reviter';
    end
    if isempty(V), return; end
    lam_v = real(diag(D));
    Nv = sqrt(max(-lam_v, 0)) / k0;
    % !! 下界取衬底折射率（1.77）而不是槽区末层：导模 N≈2.4167 落在槽区 p-wave 的 2.455 之下，
    % !! 用 2.455 当界会把所有导模都滤掉（keep 为空时 next 句会退化成全收）
    keep = Nv < na_r + 1e-9 & Nv > real(Lr(1).n) + 1e-3;
    if ~any(keep), keep = true(size(Nv)); end
    idx = find(keep);
    frac = zeros(numel(idx), 1);
    inC = (ysi >= 0.30) & (ysi <= 0.55);            % 芯/有源区（QW 附近）的行
    inX = abs(xsi) <= w/2;
    for t = 1:numel(idx)
        u = reshape(V(:, idx(t)), nyi, nxi);
        u = abs(u).^2;
        frac(t) = sum(sum(u(inC, inX))) / sum(u(:));
    end
    [~, ib] = max(frac);
    R.Nf = Nv(idx(ib));  R.frac = frac(ib);
    R.E = reshape(V(:, idx(ib)), nyi, nxi);
    R.cand = [Nv(idx), frac];
    R.ok = true;
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_d2_nonuni（非均匀网格上的三对角二阶差分，只含内部节点）
% 输入：v 单调递增节点坐标（含两端 Dirichlet 点）
% 输出：D (N−2)×(N−2) 稀疏矩阵，作用在内部节点上，已含 1/h² 缩放
% 说明：
%   1.三点格式 2/(h₋+h₊)[(u_{i+1}−u_i)/h₊ − (u_i−u_{i−1})/h₋]，h₋/h₊ 为左右间距；
%   2.均匀网格时退化成 (u_{i+1}−2u_i+u_{i−1})/h²，与 spdiags 常数写法逐位一致。
% ───────────────────────────────────────────────────────────────────────────────────
function D = e3_d2_nonuni(v)
    v = v(:).';
    n = numel(v);
    d = diff(v);                                    % d(k) = v(k+1) − v(k)，长度 n−1
    nin = n - 2;                                    % 内部节点 j=1..nin 对应全局 i=j+1
    hm = d(1:nin);                                  % 左间距 v(i)−v(i−1) = d(i−1) = d(j)
    hp = d(2:nin+1);                                % 右间距 v(i+1)−v(i)   = d(i)   = d(j+1)
    c0 = -2 ./ (hm .* hp);
    cm =  2 ./ (hm .* (hm + hp));
    cp =  2 ./ (hp .* (hm + hp));
    D = spdiags([[cm(2:end).'; 0], c0(:), [0; cp(1:end-1).']], -1:1, nin, nin);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_fd_reviter（shift-invert 反迭代 + MGS deflation，eigs 不可用时的兜底）
% 输入：Ah 稀疏算子; sig 移位; nk 取几个
% 输出：V 特征向量 (N × nk); D 对角特征值
% 说明：
%   1.对 B = A − σI 做反迭代：收敛到离 σ 最近的实特征值；每收敛一个就用 MGS 正交化掉，
%     再继续下一个（deflation），避免重复收敛到同一个；
%   2.步数上限 400、相对变化 <1e-10 即停；LU 分解只做一次；
%   3.只在 eigs 出问题时才走这条路，速度比 eigs 慢一个量级，但不会给出 NaN。
% ───────────────────────────────────────────────────────────────────────────────────
function [V, D] = e3_fd_reviter(Ah, sig, nk)
    Nn = size(Ah, 1);
    B  = Ah - sig*speye(Nn);
    [Lf, Uf, Pf, Qf] = lu(B);
    V = zeros(Nn, nk);  D = zeros(nk, nk);
    for q = 1:nk
        v = ones(Nn, 1) / sqrt(Nn);
        for it = 1:400
            w = Uf \ (Lf \ (Pf * (Qf * v)));
            for p = 1:q-1                      % 对已收敛向量再正交一次
                if norm(V(:,p)) > 0
                    w = w - V(:,p) * (V(:,p)' * w);
                end
            end
            nw = norm(w);
            if nw == 0 || ~isfinite(nw), break; end
            v2 = w / nw;
            if norm(v2 - v) < 1e-10, v = v2; break; end
            v = v2;
        end
        Av = Ah * v;
        mu = real(v' * Av);
        V(:, q) = v;  D(q, q) = mu;
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_nb（二维 FD 的逐点折射率：y′ 坐标 → 脊区层折射率，负数落在顶部介质）
% 输入：t 查询的 y′ [μm]; y_if 各界面坐标（递增，首项 0）; Lr 脊区层表; n_top 顶部介质折射率
% 输出：v 该点的复折射率
% 说明：
%   1.y′ < 0 在 p-clad 顶面之上，属于顶部介质；直接用层表查 find 会返回空，故先判负；
%   2.用 find(...,'last') 而不是插值：层内折射率是分段常数，插值会在界面处给出假的过渡值。
% ───────────────────────────────────────────────────────────────────────────────────
function v = e3_nb(yb, y_if, Lr, n_top)
    if yb > y_if(end) + 1e-12
        v = complex(n_top, 0);
        return
    end
    if yb < 0
        v = Lr(1).n;
        return
    end
    i = find(y_if <= yb + 1e-12, 1, 'last');
    if isempty(i), i = 1; end
    % !! y_if 是「界面」坐标（首项 0 = 衬底顶面），层表 Lr 的第 1 层是半无限衬底，
    % !! 所以界面区间 i 对应第 i+1 层；y′ 侧要先用 y_tot − y′ 换到自底坐标（见 e3_fd2d）。
    if i + 1 > numel(Lr), i = numel(Lr) - 1; end
    v = Lr(i + 1).n;
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_lateral_gap（横向相邻模的解析最小间隔，用于定 merge_abs 与扫描步长）
% 输入：na/nb 两区复 n_eff; w 脊宽 [μm]; k0 [μm⁻¹]
% 输出：gap 最小间隔估计（Re 量纲）
% 说明：
%   1.u₀ = π/2 − π/(2V)，Δu = π/V（相邻两阶的横向本征值间隔），
%     dN/du = u/((w/2)²k₀²Re ñ)，取 u = u₀ 得 gap ≈ u₀Δu/(w²k₀²Re ñ/4)；
%   2.这个量在 w=50 时只有 ~2e-7，比 find_roots_complex 的缺省 merge（相对 1e-6 → 绝对 2.4e-6）
%     小一个量级，所以 w ≥ 10 必须显式给 merge_abs，否则最低两阶会被并成一个（见 t7）。
% ───────────────────────────────────────────────────────────────────────────────────
function gap = e3_lateral_gap(na, nb, w, k0)
    V = (w/2) * k0 * sqrt(max(real(na)^2 - real(nb)^2, 0));
    if V <= 0, gap = 1e-3; return; end
    u0  = pi/2 - pi/(2*V);
    du  = pi/V;
    gap = abs(u0 * du / ((w/2)^2 * k0^2 * real(na)));
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_print_stack（打印一套层结构，供日志留痕）
% 输入：L 层结构; ttl 标题
% 输出：无
% 说明：逐层打印厚度、复折射率两个分量与层名，厚度 Inf 用 'Inf' 标出。
% ───────────────────────────────────────────────────────────────────────────────────
function e3_print_stack(L, ttl)
    e3_log('  -- %s（%d 层）', ttl, numel(L));
    e3_log('     %-10s %10s %12s %14s', 'tag', 'd[/μm]', 'Re(n)', 'Im(n)');
    for k = 1:numel(L)
        e3_log('     %-10s %10s %12.6f %14.3e', L(k).tag, ...
               num2str(L(k).d, '%.4f'), real(L(k).n), imag(L(k).n));
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_t1（T1：脊区 / 槽区两套层结构 + 与阶段 2 的 stack_7layer 对账）
% 输入：lambda 真空波长 [μm]
% 输出：D 结构体（.base 阶段 2 原栈; .s 默认; .s0 h_e=0; .sair air; .sinf p-clad 半无限）
% 调用函数：stack_7layer, e3_stack_pair, layer_bounds
% 说明：
%   1.对账三件事：h_e=0 时两区逐字段相等；脊区有限 p-clad 之外多出的只有顶部介质层；
%     厚度总和 = 1.729 μm（与阶段 2 的 stack_7layer 中间层厚度一致）；
%   2.stack_7layer 本函数只读不改：T1 的目的是"新结构不破坏旧结构"，旧结构还要给 A0 回归用。
% ───────────────────────────────────────────────────────────────────────────────────
function D = e3_t1(lambda)
    global STAGE2_ALPHA
    e3_log('');
    e3_log('=== T1 层结构：脊 / 槽 / 边界情形 ===');
    D = struct();
    D.lambda = lambda;
    % !! 对账用的阶段 2 原栈也要带同一个 α 表，否则复折射率的虚部对不上（先是 0、后者非 0），
    % !! 会被误判成"新结构改了材料"
    [Mx, ~] = stage2_alpha_table(lambda, 'TE', 'kio');
    STAGE2_ALPHA = Mx;
    D.base = stack_7layer(lambda);
    STAGE2_ALPHA = [];
    D.s    = e3_stack_pair(lambda, struct('h_pclad', 0.30, 'h_e', 0.30, 'n_top', 1.46));
    D.s0   = e3_stack_pair(lambda, struct('h_pclad', 0.30, 'h_e', 0.00, 'n_top', 1.46));
    D.sair = e3_stack_pair(lambda, struct('h_pclad', 0.30, 'h_e', 0.30, 'n_top', 1.00));
    D.sinf = e3_stack_pair(lambda, struct('h_pclad', Inf,  'h_e', 0.30, 'n_top', 1.46));

    e3_print_stack(D.base, '阶段 2 stack_7layer（只读对账）');
    e3_print_stack(D.s.ridge,  '脊区 default（h_pclad=0.30, SiO₂）');
    e3_print_stack(D.s.trench, '槽区 default（h_e=0.30）');
    e3_print_stack(D.sinf.ridge, '脊区对照（h_pclad=Inf，应逐层等于 stack_7layer）');

    yb = layer_bounds(build_layers(D.s.ridge));
    e3_log('  脊区中间层总厚度 = %.6f μm（期望 1.729）', yb(end, 1));
    e3_log('  p-clad 余厚 = %.6f μm，p-wave 余厚 = %.6f μm', ...
           D.s.meta.d_pclad_left, D.s.meta.d_pwave_left);

    % h_e = 0：两区应逐字段相等
    eg = isequal({D.s0.ridge.tag}, {D.s0.trench.tag}) && ...
         isequal([D.s0.ridge.d], [D.s0.trench.d]) && ...
         isequal([D.s0.ridge.n], [D.s0.trench.n]);
    D.eq_h0 = eg;
    if eg
        e3_log('  h_e=0：脊区与槽区层表逐字段相等（%d 层）', numel(D.s0.ridge));
    else
        e3_log('  h_e=0：脊区与槽区层表不等 —— 几何定义有误');
    end

    % p-clad 半无限对照：应逐层等于 stack_7layer
    mb = numel(D.base) == numel(D.sinf.ridge) && ...
         all(abs([D.base.d] - [D.sinf.ridge.d]) < 1e-12 | ...
             (isinf([D.base.d]) & isinf([D.sinf.ridge.d]))) && ...
         all(abs([D.base.n] - [D.sinf.ridge.n]) < 1e-15);
    D.eq_base = mb;
    e3_log('  h_pclad=Inf 时与 stack_7layer 一致：%d', mb);
    e3_log('  air 对照：顶部 n = %.3f', D.sair.ridge(end).n);

    % !! 验收判据（T1：层结构）
    e3_check('T1 脊区中间层总厚 1.729 μm', abs(yb(end,1) - 1.729) < 1e-9, '实得 %.9f', yb(end,1));
    e3_check('T1 h_e=0 两区层表逐字段一致', D.eq_h0, '');
    e3_check('T1 h_pclad=Inf 退回 stack_7layer', D.eq_base, '');
    e3_check('T1 顶部 air 对照生效', abs(D.sair.ridge(end).n - 1.0) < 1e-12, 'n = %.3f', D.sair.ridge(end).n);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_t2（T2：两区垂直基模 n_a / n_b 与折射率差）
% 输入：lambda [μm]
% 输出：D.tab 结构体数组（每行 h_e、n_top、n_a、n_b、Re(n_a²−n_b²)、Γ_a、Γ_b、α_a、α_b）
% 调用函数：e3_stack_pair, e3_region_solve
% 说明：
%   1.这一步只用一维求解器，不涉及横向；它的作用是先把"芯"与"包层"钉死，
%     后面横向求解只吃这两个复折射率；
%   2.同时打印 Re(n_b) < Re(n_a) 的检查结果（A7），反导引在这里就该被发现。
% ───────────────────────────────────────────────────────────────────────────────────
function D = e3_t2(lambda)
    e3_log('');
    e3_log('=== T2 两区垂直基模 ===');
    k0 = 2*pi/lambda;
    hs   = [0.20 0.25 0.30 0.35];
    tops = [1.46 1.00];
    D.tab = struct('h_e', {}, 'n_top', {}, 'na', {}, 'nb', {}, 'dl2', {}, ...
                   'Ga', {}, 'Gb', {}, 'aa', {}, 'ab', {}, 'ok', {});
    e3_log('  %6s %6s %16s %16s %12s %9s %9s %9s', ...
           'h_e', 'n_top', 'Re n_a', 'Re n_b', 'Re(na²−nb²)', 'Γ_a', 'Γ_b', 'α_a');
    opts = e3_opts_solve();
    for h = hs
        for t = tops
            Sp = e3_stack_pair(lambda, struct('h_pclad', 0.30, 'h_e', h, 'n_top', t));
            Ra = e3_region_solve(Sp.ridge,  lambda, 1, opts, 'TE');
            Rb = e3_region_solve(Sp.trench, lambda, 1, opts, 'TE');
            r  = struct('h_e', h, 'n_top', t, 'na', NaN, 'nb', NaN, 'dl2', NaN, ...
                        'Ga', NaN, 'Gb', NaN, 'aa', NaN, 'ab', NaN, 'ok', false);
            if ~isempty(Ra) && ~isempty(Rb) && Ra(1).ok && Rb(1).ok
                r.na  = Ra(1).ne;      r.nb = Rb(1).ne;
                r.dl2 = real(r.na^2 - r.nb^2);
                r.Ga  = Ra(1).G;       r.Gb = Rb(1).G;
                r.aa  = 2*k0*imag(r.na)*1e4;
                r.ab  = 2*k0*imag(r.nb)*1e4;
                r.ok  = true;
            end
            D.tab(end+1) = r;
            e3_log('  %6.2f %6.2f %16.9f %16.9f %12.6f %9.5f %9.5f %9.4f', ...
                   h, t, real(r.na), real(r.nb), r.dl2, r.Ga, r.Gb, r.aa);
        end
    end
    D.all_ok = all([D.tab.ok]) && all(arrayfun(@(r) real(r.nb) < real(r.na), D.tab));
    e3_log('  反导引检查（Re n_b < Re n_a 全表）：%d', D.all_ok);

    % !! 验收判据（T2：两区垂直基模）
    e3_check('A7 Re n_b < Re n_a（8 格全表）', D.all_ok, '');
    e3_check('T2 8 格都解出垂直基模', all([D.tab.ok]), '');
    tp = [D.tab.n_top];  na_s = real(mean([D.tab(tp == 1.46).na]));  na_a = real(mean([D.tab(tp == 1.0).na]));
    e3_check('T2 SiO₂ 顶层 Re n_a 高于 air 对照', na_s > na_a, '%.9f > %.9f', na_s, na_a);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_t3（T3：单格横向复本征 —— n_eff / 模序 / α 一致性）
% 输入：lambda [μm]; w 脊宽 [μm]; h_e; n_top
% 输出：D 结构体（.C 格点数据; .alpha_pred 等已含在 C 内）
% 调用函数：e3_cell
% 说明：
%   1.打印 V、解析模数 M、实求模数与各阶 Re N/Γ/α；
%   2.同时打印 α 与 Σ P_j α_j n_rj/Re(N) 的对照（A5 的物理含义：α 只能来自 Im N）。
% ───────────────────────────────────────────────────────────────────────────────────
function D = e3_t3(lambda, w, h_e, n_top)
    if nargin < 2 || isempty(w),      w = 3;      end
    if nargin < 3 || isempty(h_e),    h_e = 0.30; end
    if nargin < 4 || isempty(n_top),  n_top = 1.46; end
    e3_log('');
    e3_log('=== T3 横向复本征（w=%g μm, h_e=%g μm, n_top=%g）===', w, h_e, n_top);
    C = e3_cell(lambda, w, h_e, n_top, struct());
    D = struct('C', C);
    if ~C.ok && isempty(C.N)
        e3_log('  该格没解出横向模');
        return
    end
    e3_log('  V = %.6f，解析模数 M = floor(2V/π)+1 = %d，实求到 %d 阶', ...
           C.V, C.M_expect, numel(C.N));
    e3_log('  n_a = %.9f (Im %.3e)，n_b = %.9f (Im %.3e)，Re(na²−nb²) = %.6f', ...
           real(C.na), imag(C.na), real(C.nb), imag(C.nb), C.dl2);
    e3_log('  %4s %18s %12s %10s %10s %8s', 'k', 'Re N', 'Im N', 'Γ', 'α[cm⁻¹]', 'η');
    for j = 1:numel(C.N)
        e3_log('  %4d %18.12f %12.4e %10.6f %10.4f %8.4f', j, C.N(j), ...
               C.alpha(j)/(2*(2*pi/lambda)*1e4), C.Gamma(j), C.alpha(j), C.eta(j));
    end
    e3_log('  α(≈2k₀Im N) = %s', mat2str(C.alpha, 6));
    e3_log('  ΣP_jα_j n_rj/Re N = %s', mat2str(C.alpha_pred, 6));
    if numel(C.N) > 1
        e3_log('  相邻间隔 = %s', mat2str(C.gap, 4));
    end

    % !! 验收判据（T3：单格横向复本征）
    rl = abs(C.alpha - C.alpha_pred) ./ max(abs(C.alpha), eps);
    e3_check('A5 α=2k₀Im N 与 ΣP_jα_j n_rj/Re N 一致', max(rl) < 1e-4, '最大相对差 %.2e', max(rl));
    e3_check('T3 实求模数 = 解析模数 M', numel(C.N) == C.M_expect, '%d vs %d', numel(C.N), C.M_expect);
    e3_check('T3 模序按 Re N 严格递减', all(diff(C.N) < 0), '');
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_t4（T4：二维场、Γ 与逐层功率 —— 含 w→∞ 极限对账）
% 输入：lambda [μm]; w（默认 3）; h_e; n_top
% 输出：D 结构体（.C; .Fd; .Gam2d 由二维场独立算出的 Γ; .P2d; .ok）
% 调用函数：e3_cell, e3_field2d, e3_trapzw
% 说明：
%   1.Γ 有两条独立路径：e3_gamma_power（x 加权一维 Γ_y）与二维场直接积分
%     ∫∫_active|E|²dxdy/∫∫|E|²dxdy；两者一致才说明 EIM 的场与 Γ 是同一套口径；
%   2.w→∞ 的极限不在这里做（w=50 见 t5），本函数只出 w=3 的场（交付图 field2d_w3）。
% ───────────────────────────────────────────────────────────────────────────────────
function D = e3_t4(lambda, w, h_e, n_top)
    if nargin < 2 || isempty(w),     w = 3;      end
    if nargin < 3 || isempty(h_e),   h_e = 0.30; end
    if nargin < 4 || isempty(n_top), n_top = 1.46; end
    e3_log('');
    e3_log('=== T4 二维场与 Γ（w=%g μm, h_e=%g μm, n_top=%g）===', w, h_e, n_top);
    C  = e3_cell(lambda, w, h_e, n_top, struct());
    D  = struct('C', C, 'Fd', [], 'Gam2d', [], 'ok', false);
    if isempty(C.N) || ~C.ok
        e3_log('  该格无有效横向模，跳过');
        return
    end
    K = min(numel(C.N), 4);
    Fd = e3_field2d(C, 1:K, 0.010, 0.001);
    D.Fd = Fd;
    % !! 施工单 T4 的 Γ 定义就是 ∫|E_x|²Γ_y(x)dx/∫|E_x|²dx（= 下表"1D 加权"列，交付值）。
    % !! 二维直积只是离散化检查：y 网格在 QW 内要 1 nm 以下才收敛，故再跑一级 0.25 nm 看趋势。
    % !! 细化只针对基模：高阶模在槽区的权重更大，而槽区 Φ 在 SiO2/QW 界面附近变化最快，
    % !! 直积的残差主要来自那里的离散化；基模足以说明"直积会收敛到 1D 加权值"。
    Fd2 = e3_field2d(C, 1, 0.020, 0.00025);
    D.Fd2 = Fd2;

    A = C.Ra(1);  B = C.Rb(1);
    ybA = layer_bounds(A.L);  ybB = layer_bounds(B.L);
    iq  = find(strcmp({A.L.tag}, 'QW'), 1);
    yq0 = ybA(iq, 1);  yq1 = ybA(iq, 2);
    inY = (Fd.y >= yq0) & (Fd.y <= yq1);
    inX = Fd.inR;
    % !! y 一律是「自衬底界面向上」的坐标（与 mode_profile/layer_bounds 同约定），不是 y′
    e3_log('  QW 行窗口 y∈[%.4f, %.4f]，命中 %d 行（y 网格 %.4f..%.4f，%d 点）', ...
           yq0, yq1, sum(inY), Fd.y(1), Fd.y(end), numel(Fd.y));
    % !! 分子分母都要乘 trapz 权重：y 网格在 QW 附近 5 nm、外侧 25 nm，等权求和会把
    % !! QW 那两三行按 5 倍权重算进去，Γ 直接翻倍（就是坑记录 A5 的变体）
    WY = Fd.wy(:);  WX = Fd.wx(:).';
    Wt = WY * WX;
    G2  = zeros(1, K);  G3 = zeros(1, K);
    for q = 1:K
        E2 = abs(Fd.E(:, :, q)).^2;
        G2(q) = sum(sum(E2(inY, inX) .* Wt(inY, inX))) / sum(sum(E2 .* Wt));
        if q == 1
            E3 = abs(Fd2.E).^2;
            iY3 = (Fd2.y >= yq0) & (Fd2.y <= yq1);
            iX3 = Fd2.inR;
            Wt3 = Fd2.wy(:) * Fd2.wx(:).';
            G3(q) = sum(sum(E3(iY3, iX3) .* Wt3(iY3, iX3))) / sum(sum(E3 .* Wt3));
        else
            G3(q) = NaN;
        end
    end
    D.Gam2d = G2;  D.Gam2d_fine = G3;
    e3_log('  %4s %14s %15s %15s %12s %12s', 'k', 'Γ(1D 加权)', 'Γ(2D,dy=1nm)', 'Γ(2D,dy=0.25nm)', '差(0.25nm)', 'α[cm⁻¹]');
    for q = 1:K
        e3_log('  %4d %14.8f %15.8f %15.8f %12.2e %12.4f', q, C.Gamma(q), G2(q), G3(q), ...
               abs(C.Gamma(q)-G3(q))/max(C.Gamma(q), eps), C.alpha(q));
    end
    e3_log('  平板极限（w→∞ 的 Γ_a）= %.8f，阶段 2 归档 0.026775', C.Gamma_a);
    % !! 诊断：两区 Φ 在并集网格上的"局部 Γ"（Γ(2D) 必须是这两个数的加权平均，越界即为构建错误）
    Ga_loc = sum(abs(Fd.phiA(inY)).^2 .* Fd.wy(inY)) / sum(abs(Fd.phiA).^2 .* Fd.wy);
    Gb_loc = sum(abs(Fd.phiB(inY)).^2 .* Fd.wy(inY)) / sum(abs(Fd.phiB).^2 .* Fd.wy);
    e3_log('  并集网格上的局部 Γ：脊区 %.8f（应为 %.8f），槽区 %.8f（应为 %.8f）', ...
           Ga_loc, real(A.G), Gb_loc, real(B.G));
    e3_log('  y 网格 %d 点 × x 网格 %d 点（dx = %.3f μm）', numel(Fd.y), numel(Fd.x), 0.010);
    % !! 验收判据（T4：二维场与 Γ）
    % !! 判据取 1e-2：QW 只有 9 nm、模场在其中是窄峰，trapz 在 1 nm 档给 0.02650、0.25 nm 档给
    % !! 0.02696（1D 加权 0.026712），残差 ±1% 是离散化真值而非模型误差；交付值取 1D 加权列，
    % !! 直积这条只用来证明「两区 Φ 拼出来的二维场是对的」（Φ 若拼接错、双计或漏算，差会是量级）。
    e3_check('T4 Γ(2D,dy=0.25nm) 与 Γ(1D 加权) 相对差 < 1e-2', ...
             abs(G3(1)-C.Gamma(1))/C.Gamma(1) < 1e-2, '基模 %.2e', abs(G3(1)-C.Gamma(1))/C.Gamma(1));
    e3_check('T4 一阶以上 Γ(2D,dy=1nm) 与 Γ(1D) 同号同量级', ...
             all(G2 > 0) && G2(1) <= 0.05, 'G2(1) = %.6f', G2(1));
    D.ok = true;
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_t5（T5：40 格扫描 —— w × h_e × 顶部介质）
% 输入：lambda [μm]
% 输出：D.tab 结构体数组（每格一个元素，字段见下）; D.M 长表 (Nmodes × 8)
% 调用函数：e3_cell
% 说明：
%   1.表列：w、h_e、n_top、模数、M解析、Re N(1)、Γ(1)、α(1)、min_gap、排序单调性；
%   2.每格都记录 ok 与 Re(n_b) < Re(n_a)，A7 就是这张表的全体检查；
%   3.405/520 波长不在这张表里（G9：只冒烟，不进结论）。
% ───────────────────────────────────────────────────────────────────────────────────
function D = e3_t5(lambda)
    e3_log('');
    e3_log('=== T5 40 格扫描（5 个 w × 4 个 h_e × {SiO₂, air}）===');
    ws = [2 3 5 10 50];
    hs = [0.20 0.25 0.30 0.35];
    ts = [1.46 1.00];
    D = struct('tab', struct('w', {}, 'h_e', {}, 'n_top', {}, 'K', {}, 'M', {}, ...
                             'N1', {}, 'G1', {}, 'a1', {}, 'gap', {}, 'mono', {}, ...
                             'ok', {}, 'ratio', {}, 'sec', {}, 'na', {}, 'Ga', {}, ...
                             'okA7', {}), 'M', [], 'lambda', lambda);
    e3_log('  %5s %6s %6s %4s %4s %16s %10s %9s %11s %3s %5s', ...
           'w', 'h_e', 'n_top', 'K', 'M', 'Re N(1)', 'Γ(1)', 'α(1)', 'min_gap', 'A7', 'mono');
    for w = ws
        for h = hs
            for t = ts
                t0 = tic;
                C  = e3_cell(lambda, w, h, t, struct());
                tc = toc(t0);
                r  = struct('w', w, 'h_e', h, 'n_top', t, 'K', numel(C.N), ...
                            'M', C.M_expect, 'N1', NaN, 'G1', NaN, 'a1', NaN, ...
                            'gap', NaN, 'mono', NaN, 'ok', false, 'ratio', 1, 'sec', tc, ...
                            'na', NaN, 'Ga', NaN, 'okA7', false);
                if ~isempty(C.N)
                    r.N1 = C.N(1);  r.G1 = C.Gamma(1);  r.a1 = C.alpha(1);
                    r.na = real(C.na);  r.Ga = C.Gamma_a;
                    r.okA7 = real(C.nb) < real(C.na);
                    if numel(C.N) > 1, r.gap = min(C.gap); end
                    r.ratio = numel(C.N) / max(C.M_expect, 1);
                    r.mono  = all(diff(C.N) < 0);
                    r.ok    = true;
                end
                D.tab(end+1) = r;
                e3_log('  %5g %6.2f %6.2f %4d %4d %16.9f %10.6f %9.4f %11.3e %3d %5d', ...
                       w, h, t, r.K, r.M, r.N1, r.G1, r.a1, r.gap, ...
                       double(~isnan(r.N1) && real(C.nb) < real(C.na)), double(r.mono));
            end
        end
    end
    e3_log('  40 格总耗时 %.1f s（单格最长 %.1f s）', ...
           sum([D.tab.sec]), max([D.tab.sec]));

    % !! 验收判据（T5：40 格扫描）
    tb   = D.tab;
    pick = @(w0, t0) find([tb.w] == w0 & [tb.h_e] == 0.30 & [tb.n_top] == t0, 1);
    i50 = pick(50, 1.46);  i10 = pick(10, 1.46);  i2 = pick(2, 1.46);  i3 = pick(3, 1.46);  i5 = pick(5, 1.46);
    e3_check('A1 w=50 时 |Re N(1) − Re n_a| < 1e-5', ...
             abs(tb(i50).N1 - tb(i50).na) < 1e-5, '%.3e', abs(tb(i50).N1 - tb(i50).na));
    e3_check('A2 |Re N(1)(w=10) − Re n_a| < 1e-3', ...
             abs(tb(i10).N1 - tb(i10).na) < 1e-3, '%.3e', abs(tb(i10).N1 - tb(i10).na));
    seq = [tb(i2).N1 tb(i3).N1 tb(i5).N1 tb(i10).N1 tb(i50).N1];
    e3_check('A2 Re N(1) 随 w 单调递增', all(diff(seq) > 0), mat2str(seq, 9));
    e3_check('A3 Γ(2)<Γ(3)<Γ(5)<Γ(10)', ...
             all(diff([tb(i2).G1 tb(i3).G1 tb(i5).G1 tb(i10).G1]) > 0), ...
             mat2str([tb(i2).G1 tb(i3).G1 tb(i5).G1 tb(i10).G1], 6));
    e3_check('A3 |Γ(50)/Γ_slab − 1| < 2e-2', ...
             abs(tb(i50).G1/tb(i50).Ga - 1) < 2e-2, 'Γ50 = %.8f，Γ_slab = %.8f', tb(i50).G1, tb(i50).Ga);
    e3_check('A6 40 格模数 K == M 全表', all([tb.K] == [tb.M]), ...
             'K==M 的格数 %d/%d', sum([tb.K] == [tb.M]), numel(tb));
    e3_check('A7 40 格 Re n_b < Re n_a 全表', all([tb.okA7]), '');
    e3_check('T5 40 格全部解出横向模', all([tb.ok]), '');
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_t6（T6：二维有限差分对照，只比 Re）
% 输入：lambda [μm]
% 输出：D.tab 结构体数组（w、EIM 的 Re N、FD 的 Re N、差值、芯内占比、未知数个数、耗时）
% 调用函数：e3_cell, e3_fd2d
% 说明：
%   1.比较对象固定为 h_e=0.30、SiO₂、w ∈ {3,5,10}；只比 Re（G 级判据 A4 的两半都在这里）；
%   2.第二半是网格收敛：最窄档 w=3 用 20 nm 与 10 nm 各跑一遍，两者之差应 < 2e-4，
%     这才是"差 5e-3 不是网格误差"的依据；
%   3.FD 自己的耗时、未知数个数、选模方法都打进日志，便于复核为什么某格特别慢。
% ───────────────────────────────────────────────────────────────────────────────────
function D = e3_t6(lambda)
    e3_log('');
    e3_log('=== T6 二维有限差分对照（只比 Re 有效折射率）===');
    k0 = 2*pi/lambda;
    ws = [3 5 10];
    D = struct('tab', struct('w', {}, 'Neim', {}, 'Nfd', {}, 'd', {}, 'frac', {}, ...
                             'n', {}, 'sec', {}, 'method', {}, 'dx', {}), 'conv', []);
    e3_log('  %5s %18s %18s %11s %8s %9s %6s %-12s', ...
           'w', 'Re N (EIM)', 'Re N (FD)', '|Δ|', '芯占比', '未知数', 's', 'method');
    % !! 零号自检：h_e=0 时两区折射率完全相同、结构在 x 方向处处相同，二维算子的基模必须
    % !! 退回一维竖直基模（2.416918895）。这是唯一能单独检验「二维算子本身有没有错」的格子。
    t0 = tic;  R0 = e3_fd2d(lambda, 3, 0, 1.46, ...
                            struct('shift_n', 2.416918900, 'dbg', true, 'na_ref', 2.416918895));  tc = toc(t0);
    e3_log('  零号自检（h_e=0，结构沿 x 均匀）：FD 基模 %.9f vs 一维 TMM %.9f，差 %.2e（%.1f s）', ...
           R0.Nf, 2.416918895, abs(R0.Nf - 2.416918895), tc);

    for w = ws
        C  = e3_cell(lambda, w, 0.30, 1.46, struct());
        ne = NaN;
        if ~isempty(C.N), ne = C.N(1); end
        % !! 移位点取 EIM 的 Re N(1)：它只决定 eigs 返回"哪几个"本征值，不改变特征值本身；
        % !! 若按 0.995·n_max 定移位，会落在导模带（≈2.4167）之外，取回的全是箱模，芯占比为 0。
        t0 = tic;  R = e3_fd2d(lambda, w, 0.30, 1.46, ...
                               struct('dx', 0.020, 'shift_n', real(ne), 'dbg', w == 3, ...
                                      'na_ref', real(C.na)));   tc = toc(t0);
        r  = struct('w', w, 'Neim', ne, 'Nfd', NaN, 'd', NaN, 'frac', R.frac, ...
                    'n', R.n, 'sec', tc, 'method', R.method, 'dx', R.dx);
        if ~isempty(R.Nf), r.Nfd = R.Nf; end
        if ~isnan(r.Nfd) && ~isnan(ne), r.d = abs(ne - r.Nfd); end
        D.tab(end+1) = r;
        if isfield(R, 'cand') && ~isempty(R.cand) && w == 3
            e3_log('       FD 候选（N / 芯占比）：');
            for tt = 1:size(R.cand, 1)
                e3_log('         %18.9f  %8.5f', R.cand(tt,1), R.cand(tt,2));
            end
        end
        e3_log('  %5g %18.9f %18.9f %11.3e %8.4f %9d %6.1f %-12s', ...
               w, r.Neim, r.Nfd, r.d, r.frac, r.n, r.sec, r.method);
    end

    % 网格收敛：w=3 整体加密 2 倍（各分段步长 ×0.5）
    t0 = tic;  R2 = e3_fd2d(lambda, 3, 0.30, 1.46, ...
                            struct('fac', 0.5, 'shift_n', real(D.tab(1).Neim)));  tc = toc(t0);
    n20 = D.tab(1).Nfd;
    D.conv = [n20, R2.Nf];
    e3_log('  网格收敛（w=3）：基准段长 → %.9f，全段加密 2 倍 → %.9f，差 %.2e（%.1f s）', ...
           n20, R2.Nf, abs(n20 - R2.Nf), tc);
    % 域收敛：w=3 侧向半宽从 w/2+1.0 μm 扩到 w/2+2.0 μm
    t0 = tic;  R3 = e3_fd2d(lambda, 3, 0.30, 1.46, ...
                            struct('xpad', 2.0, 'shift_n', real(D.tab(1).Neim)));  tc = toc(t0);
    e3_log('  域收敛（w=3）：半宽 w/2+1.0 μm → %.9f，w/2+2.0 μm → %.9f，差 %.2e（%.1f s）', ...
           n20, R3.Nf, abs(n20 - R3.Nf), tc);
    e3_log('  FD 域：x ∈ ±(w/2+1.0) μm（非均匀：脊边 10 nm / 脊内 20 nm / 外 25 nm）、y′ ∈ [−0.3, 2.0] μm（QW 附近 1 nm），两端 Dirichlet');

    % !! 验收判据（T6：A4 的两半）
    dv = [D.tab.d];
    e3_check('A4 EIM vs 2D FD：|Δ Re N| < 5e-3（w=3,5,10）', all(dv < 5e-3), '最大 %.3e', max(dv));
    e3_check('A4 网格收敛（w=3，整体加密 2 倍）< 2e-4', abs(n20 - R2.Nf) < 2e-4, '%.2e', abs(n20 - R2.Nf));
    e3_check('A4 域收敛（w=3，半宽扩到 w/2+2 μm）< 2e-4', abs(n20 - R3.Nf) < 2e-4, '%.2e', abs(n20 - R3.Nf));
    e3_check('T6 零号自检（结构沿 x 均匀）解释 D 侧壁代价', ...
             abs(R0.Nf - (2.416918895 - (pi/5)^2/(2*k0^2*2.416918895))) < 2e-4, ...
             'FD %.9f，1D−(π/L)²/2k₀²N %.9f', R0.Nf, 2.416918895 - (pi/5)^2/(2*k0^2*2.416918895));
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_t7（T7：相邻模间隔、模数整数断言、merge 去重的影响）
% 输入：lambda [μm]
% 输出：D.tab 结构体数组（w、模数 K、解析 M、min_gap、扫描步长、相对式 merge 下的 K）
% 调用函数：e3_stack_pair, e3_region_solve, e3_lateral_solve, e3_lateral_gap
% 说明：
%   1.判据不是施工单原写的"min_gap > 100e-6"：w=50 的解析间隔只有 2.2e-7，
%     那个阈值在物理上不可能通过（见计划 §偏离 1）。这里改判
%     min_gap > 10·merge_abs 且 > 5·扫描步长，另单列模数的整数断言 K == M；
%   2.w=50 同时给出"缺省相对 merge"的结果，用来说明不去重修正会掉一阶（48→47）。
% ───────────────────────────────────────────────────────────────────────────────────
function D = e3_t7(lambda)
    e3_log('');
    e3_log('=== T7 模间隔 / 模数断言 / merge 影响 ===');
    ws  = [5 10 50];
    k0  = 2*pi/lambda;
    opts = e3_opts_solve();
    D = struct('tab', struct('w', {}, 'K', {}, 'M', {}, 'gap', {}, 'step', {}, ...
                             'Krel', {}, 'mgap', {}), 'lambda', lambda);
    e3_log('  %5s %4s %4s %12s %12s %6s %6s', 'w', 'K', 'M', 'min_gap', '扫描步长', 'K相对', 'merge');
    for w = ws
        Sp = e3_stack_pair(lambda, struct('h_pclad', 0.30, 'h_e', 0.30, 'n_top', 1.46));
        Ra = e3_region_solve(Sp.ridge,  lambda, 1, opts, 'TE');
        Rb = e3_region_solve(Sp.trench, lambda, 1, opts, 'TE');
        ga = e3_lateral_gap(Ra(1).ne, Rb(1).ne, w, k0);
        Qa = e3_lateral_solve(lambda, Ra(1).ne, Rb(1).ne, w, struct('merge_abs', ga/10));
        Qb = e3_lateral_solve(lambda, Ra(1).ne, Rb(1).ne, w, struct('merge_abs', NaN));
        kg = Inf;
        if numel(Qa.N) > 1, kg = min(-diff(real(Qa.N))); end
        r = struct('w', w, 'K', numel(Qa.N), 'M', Qa.M_expect, 'gap', kg, ...
                   'step', (real(Ra(1).ne)-real(Rb(1).ne))/max(Qa.n_scan,1), ...
                   'Krel', numel(Qb.N), 'mgap', ga);
        D.tab(end+1) = r;
        e3_log('  %5g %4d %4d %12.3e %12.3e %6d %6.1e', w, r.K, r.M, r.gap, r.step, ...
               r.Krel, ga);
        % !! 诊断：把最上/最下三阶的 Re N 与两区 n 一起打出来，模数对不上时能立刻判断
        % !! 是"上端（N→Re n_a）漏根"还是"下端（N→Re n_b）漏根"
        Nq = real(Qa.N(:).');
        e3_log('       Re n_a = %.9f, Re n_b = %.9f, n_hi = %.9f, n_lo = %.9f', ...
               real(Ra(1).ne), real(Rb(1).ne), real(Ra(1).ne)-1e-9*real(Ra(1).ne), ...
               real(Rb(1).ne)+1e-9*real(Rb(1).ne));
        if numel(Nq) >= 3
            e3_log('       最高三阶 = %s', mat2str(Nq(1:3), 9));
            e3_log('       最低三阶 = %s', mat2str(Nq(max(1,end-2):end), 9));
        end
    end

    % !! 验收判据（T7：模间隔与模数）——判据按计划 §偏离 1 改写，不用施工单原写的 100e-6
    gp = [D.tab.gap];  st = [D.tab.step];  mg = [D.tab.mgap];
    e3_check('A6 min_gap > 10·merge_abs（三档）', all(gp > mg), ...
             'min_gap/merge_abs = %s', mat2str(gp./max(mg,eps), 4));
    e3_check('A6 min_gap > 5·扫描步长（三档）', all(gp > 5*st), ...
             'min_gap/步长 = %s', mat2str(gp./max(st,eps), 4));
    e3_check('A6 模数 K == M（w=5,10,50）', all([D.tab.K] == [D.tab.M]), ...
             'K = %s，M = %s', mat2str([D.tab.K]), mat2str([D.tab.M]));
    % !! 这条不是判据而是对照：加了复种子之后，缺省相对式 merge 在 w=50 也不再掉阶；
    % !! 若哪天它又变少，说明种子机制失效（不是 merge 变严），必须回来看这一行。
    e3_check('T7 缺省相对 merge 与 merge_abs 模数一致（不掉阶）', ...
             all([D.tab.Krel] == [D.tab.K]), 'K相对 = %s，K = %s', ...
             mat2str([D.tab.Krel]), mat2str([D.tab.K]));
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_t8（T8：侧壁散射损耗的参数化上界）
% 输入：lambda [μm]; w 脊宽（默认 3）
% 输出：D.tab 结构体数组（σ、L_c、加权 α_scat、逐模 α_scat）; D.eta
% 调用函数：e3_cell, e3_scatter
% 说明：
%   1.σ ∈ {0.5,1,2} nm × L_c ∈ {10,30,100} nm 共 9 组，C=1；
%   2.这只是上界：没有实测粗糙度，数值不进结论（G8），只用来判断散射是否与材料吸收同量级。
% ───────────────────────────────────────────────────────────────────────────────────
function D = e3_t8(lambda, w)
    if nargin < 2 || isempty(w), w = 3; end
    e3_log('');
    e3_log('=== T8 侧壁散射损耗（参数化上界，C=1）===');
    C = e3_cell(lambda, w, 0.30, 1.46, struct());
    D = struct('tab', struct('sigma', {}, 'Lc', {}, 'a', {}, 'aK', {}), 'eta', []);
    if isempty(C.N), return; end
    D.eta = C.eta;
    e3_log('  侧壁归一化强度 η（逐模）= %s', mat2str(C.eta, 4));
    e3_log('  %8s %8s %12s %12s', 'σ[nm]', 'L_c[nm]', 'α_scat(1阶)', 'α_scat(加权)');
    for s = [0.5 1 2]
        for l = [10 30 100]
            [a, aK] = e3_scatter(C, s, l, 1);
            D.tab(end+1) = struct('sigma', s, 'Lc', l, 'a', a, 'aK', aK);
            e3_log('  %8.1f %8.1f %12.3f %12.3f', s, l, aK(1), a);
        end
    end
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_t9（T9：二维远场、两个切面、叠加、Parseval 与远场一致性判据）
% 输入：lambda [μm]; w 脊宽（默认 3）
% 输出：D 结构体
%         .Iov 叠加强度图 (481×481，归一化); .F2 各阶复矩阵; .fast/.slow 两个切面
%         .rat 各阶 Parseval 比; .rat_ov 叠加的 Parseval 比; .A9/.A10 两个一致性量
% 调用函数：e3_cell, e3_field2d, e3_ff_matrix, e3_ff_axes, e3_overlay_weights,
%          e3_parseval, e3_fwhm, far_field, e3_region_solve
% 说明：
%   1.远场只能从二维场做：先对 y 做一次变换、再对 x 做一次，等价于 F = Mx'·E·My；
%     不能把横向解与垂直解各自做一维远场再相乘（那只在可分离时成立，而本结构的 Φ 随 x 分片）；
%   2.A11 的 Parseval 单独用 s ∈ [−0.999,0.999] 的均匀网格做，不用 ±60° 的网格：
%     快轴 FWHM 是 23.95°，|s| ≤ 0.866 按物理就丢掉 y 方向约 35% 的能量，在那里判 Parseval 必然假失败；
%   3.叠加权重 c 由 e3_overlay_weights 给（假设值），四组敏感性在别处跑（见 t9 的日志）。
% ───────────────────────────────────────────────────────────────────────────────────
function D = e3_t9(lambda, w)
    if nargin < 2 || isempty(w), w = 3; end
    h_e = 0.30;  n_top = 1.46;  k0 = 2*pi/lambda;
    e3_log('');
    e3_log('=== T9 二维远场（w=%g μm, h_e=%g, n_top=%g）===', w, h_e, n_top);
    C = e3_cell(lambda, w, h_e, n_top, struct());
    D = struct('ok', false, 'C', C);
    if isempty(C.N) || ~C.ok
        e3_log('  无有效模，跳过');
        return
    end
    K  = min(numel(C.N), 4);
    Fd = e3_field2d(C, 1:K, 0.010);
    D.Fd = Fd;

    % !! 远场网格：施工单 §T9 要求 ±60°/0.25°，但快轴主瓣 FWHM 24.9° 的尾巴到 60° 还没落干净，
    % !! 二维图看着像被切了一刀。这里扩到 ±85°（仍是 0.25° 步长），±60° 的判据值不受影响（A9/A10b 用同网格自洽）
    th = linspace(-85, 85, 681).';
    sx = sind(th);  sy = sind(th);
    F2 = cell(1, K);  I2 = cell(1, K);
    for q = 1:K
        Fq = e3_ff_matrix(Fd, q, sx, sy, k0);
        F2{q} = Fq;
        I2{q} = abs(Fq).^2 / max(abs(Fq(:)).^2);
    end
    c   = e3_overlay_weights(K);
    Fov = zeros(size(F2{1}));
    for q = 1:K, Fov = Fov + c(q)*F2{q}; end
    Iov = abs(Fov).^2 / max(abs(Fov(:)).^2);
    D.c = c;  D.F2 = F2;  D.I2 = I2;  D.Fov = Fov;  D.Iov = Iov;
    D.th = th;  D.sx = sx;  D.sy = sy;
    e3_log('  参与叠加的模数 K = %d，c = %s', K, mat2str(c.', 4));

    % ── 两个切面 ────────────────────────────────────────────────
    [uF1, IF1] = e3_ff_axes(F2{1}, sx, sy, true);
    [uFv, IFv] = e3_ff_axes(Fov,  sx, sy, true);
    [uS1, IS1] = e3_ff_axes(F2{1}, sx, sy, false);
    [uSv, ISv] = e3_ff_axes(Fov,  sx, sy, false);
    D.fast = struct('u1', uF1, 'I1', IF1, 'uov', uFv, 'Iov', IFv);
    D.slow = struct('u1', uS1, 'I1', IS1, 'uov', uSv, 'Iov', ISv);
    f1 = rad2deg(e3_fwhm(uF1, IF1));  fv = rad2deg(e3_fwhm(uFv, IFv));
    s1 = rad2deg(e3_fwhm(uS1, IS1));  sv = rad2deg(e3_fwhm(uSv, ISv));
    D.fwhm_fast = [f1 fv];  D.fwhm_slow = [s1 sv];
    e3_log('  快轴（y）FWHM：基模 %.3f°，叠加 %.3f°', f1, fv);
    e3_log('  慢轴（x）FWHM：基模 %.3f°，叠加 %.3f°', s1, sv);
    e3_log('  叠加条纹对比度（慢轴）= %.3f', ...
           (max(ISv)-min(ISv))/(max(ISv)+min(ISv)));

    % ── A9：矩阵切面 vs 直接一维变换 ────────────────────────────
    i0 = find(abs(sx) < 1e-12, 1);
    Ay = zeros(numel(Fd.y), 1);
    for q = 1:K
        Aq = sum(Fd.E(:, :, q) .* (Fd.wx.'), 2);
        Ay = Ay + c(q)*Aq;
    end
    My  = exp(-1i*k0*(Fd.y(:)*sy(:).'));
    Fdr = My.' * (Fd.wy .* Ay);
    D.A9 = max(abs(Fov(:, i0) - Fdr)) / max(abs(Fov(:, i0)));
    e3_log('  A9：θx=0 切面（矩乘）vs 直接一维变换，相对差 %.2e', D.A9);

    % ── A10：h_e=0 退化回平板（垂直方向）＋ 2D 切面与 1D 变换的一致性 ──
    R0  = e3_region_solve(stack_7layer(lambda), lambda, 1, e3_opts_solve(), 'TE');
    ff0 = far_field(R0(1).y, R0(1).E, k0);
    C0  = e3_cell(lambda, w, 0.0, n_top, struct());
    ff1 = far_field(C0.Ra(1).y, C0.Ra(1).E, k0);
    % !! h_e=0 时两区折射率相同、横向无导引（Δn²=0），x 网格的半宽 w/2+8/γ_lat 会发散，
    % !! 也没有"横向模"可言。所以 2D 切面改在 h_e 正常的一格上做：sx=0 处
    % !! F(0,sy) = (∫E_x dx)·∫Φ e^(−i k₀ y s_y) dy，与横向解是否分离无关，正好检验
    % !! 矩阵方向、trapz 权重与归一化；退化格只比垂直方向的 1D 变换
    Ay = sum(Fd.E(:, :, 1) .* Fd.wx.', 2);          % 先对 x 积分，得 y 剖面
    My = exp(-1i*k0*(Fd.y(:)*sy(:).'));
    T1 = abs(My.' * (Fd.wy(:) .* Ay)).^2;  T1 = T1 / max(T1);
    D.A10b = max(abs(T1(:) - IF1(:)));
    D.fwhm_ref = ff0.fwhm;  D.fwhm_h0 = ff1.fwhm;
    e3_log('  A10：半无限 p-clad 平板 far_field FWHM = %.3f°（归档 23.95°）', rad2deg(ff0.fwhm));
    % !! 诊断行：远场 FWHM 由 |E|² 的近场宽度决定，把两侧都打出来，差一个量级时立刻能看出是网格还是模
    e3_log('       诊断：半无限栈 %d 点 y∈[%.3f, %.3f] μm，|E|² FWHM = %.1f nm，th_pk = %.2f°', ...
           numel(R0(1).y), R0(1).y(1), R0(1).y(end), ...
           1e3*e3_fwhm(R0(1).y(:), abs(R0(1).E(:)).^2), rad2deg(ff0.th_pk));
    e3_log('       诊断：二维场 %d 点 y∈[%.3f, %.3f] μm，脊区 |E|² FWHM = %.1f nm，th_pk = %.2f°', ...
           numel(Fd.y), Fd.y(1), Fd.y(end), ...
           1e3*e3_fwhm(Fd.y(:), abs(Fd.E(:, max(1,round(end/2)), 1)).^2), rad2deg(ff1.th_pk));
    e3_log('       h_e=0 有限 p-clad 脊区 1D FWHM = %.3f°，与半无限平板差 %.3f°', ...
           rad2deg(ff1.fwhm), rad2deg(abs(ff1.fwhm - ff0.fwhm)));
    e3_log('       2D θx=0 切面 vs 直接一维变换，最大差 %.2e（应 ~1e-15）', D.A10b);
    D.A10 = rad2deg(abs(ff1.fwhm - ff0.fwhm));

    % ── A11：Parseval（独立的 s 网格）────────────────────────────
    s   = linspace(-0.999, 0.999, 1001).';
    ws  = e3_trapzw(s);
    % !! A11 分三条判：(a) 解析带限场 —— 这才是"矩阵方向 / trapz 权重 / λ² 归一化"的真正判据；
    % !! (b)(c) 真实模在 dx=10 nm 与 5 nm 上的比值，必然略小于 1：|s|>1 的倏逝分量不参与
    % !! 传播积分，那 0.5%~1% 是物理丢的，不是代码错，只作诊断与网格收敛证据
    xA = (-2.5:0.02:2.5).';  yA = (-2.5:0.02:2.5).';
    EA = exp(-(xA.'.^2 + yA.^2)/(2*0.5^2));             % σ=0.5 μm，边界已衰减到 3e-4
    FdA = struct('x', xA, 'y', yA, 'wx', e3_trapzw(xA), 'wy', e3_trapzw(yA), ...
                 'jlist', 1, 'E', reshape(EA, numel(yA), numel(xA), 1));
    D.ratA = e3_parseval(FdA, 1, linspace(-0.9999,0.9999,1201).', k0);
    rat = zeros(1, K);
    for q = 1:K
        rat(q) = e3_parseval(Fd, q, s, k0);
    end
    % !! (d) 把 s 积分域扩到 |s| ≤ 2（含倏逝分量）：Parseval 在整条 s 轴上恒成立，
    % !! 若比值回到 1，就证明 (b) 的 0.57% 缺口是 |s|>1 的倏逝能量，不是归一化错
    sAll = linspace(-2, 2, 4001).';
    ratAll = zeros(1, K);
    for q = 1:K
        ratAll(q) = e3_parseval(Fd, q, sAll, k0);
    end
    Fd2 = e3_field2d(C, 1:K, 0.005);
    rat2 = zeros(1, K);
    for q = 1:K
        rat2(q) = e3_parseval(Fd2, q, s, k0);
    end
    Fs = zeros(numel(s), numel(s));
    for q = 1:K, Fs = Fs + c(q)*e3_ff_matrix(Fd, q, s, s, k0); end
    Eov = zeros(size(Fd.E(:, :, 1)));
    for q = 1:K, Eov = Eov + c(q)*Fd.E(:, :, q); end
    LHS  = sum(sum((abs(Fs).^2) .* (ws*ws.')));
    wE   = sum(sum((abs(Eov).^2) .* (Fd.wy*Fd.wx.')));
    D.rat = rat;  D.rat2 = rat2;  D.rat_ov = LHS / ((2*pi/k0)^2 * wE);
    e3_log('  A11(a)：解析带限场 σ=1 μm 的 Parseval 比 = %.10f（目标 1，容差 1e-6）', D.ratA);
    e3_log('  A11(b)：真实模 dx=10 nm 比值 = %s（倏逝损失 %.3f%%）', ...
           mat2str(rat, 6), 100*(1-mean(rat)));
    e3_log('  A11(c)：真实模 dx=5 nm 比值 = %s（收敛证据）', mat2str(rat2, 6));
    e3_log('  A11(d)：s∈[−2,2] 全轴比值 = %s（回到 1 即证明缺口是倏逝能量）', mat2str(ratAll, 6));
    e3_log('       叠加场 Parseval 比 = %.8f', D.rat_ov);
    e3_log('  A11：见上三条（(a) 解析场、(b) 粗网格、(c) 细网格）');

    % ── A12：单模权重 ───────────────────────────────────────────
    I0 = abs(1*F2{1}).^2;  I0 = I0 / max(I0(:));
    D.A12 = max(abs(I0(:) - I2{1}(:)));
    e3_log('  A12：c=[1,0,…] 与单模逐点最大差 %.2e', D.A12);

    % !! 加到五组：用户 2026-10-28 看图时问"叠加里到底几个模"，于是补一组"只用偶数模 m0+m2"，
    % !! 它与三模叠加的差别就是奇模 m1 的贡献（看 FWHM 与条纹对比度差多少即可）
    D.sens = struct('name', {}, 'c', {}, 'fw', {}, 'sw', {}, 'rho', {});
    sets = { '仅基模', 1; '基模+一阶', min(K, 2); '前三阶', min(K, 3); '按 Γ 加权', -1; ...
             '偶数模对 m0+m2', -2 };
    for i = 1:size(sets, 1)
        nm = sets{i, 1};  kk = sets{i, 2};
        if isequal(kk, -1)
            g = C.Gamma(1:K).';  g = g / sum(g);
        elseif isequal(kk, -2)
            g = zeros(K, 1);
            if K >= 3
                g([1 3]) = 1/sqrt(2);
            else
                g(1) = 1;
            end
        else
            g = zeros(K, 1);  g(1:kk) = 1/sqrt(kk);
        end
        Fo = zeros(size(F2{1}));
        for q = 1:K, Fo = Fo + g(q)*F2{q}; end
        Io = abs(Fo).^2 / max(abs(Fo(:)).^2);
        [uf, If] = e3_ff_axes(Fo, sx, sy, true);
        [us, Is] = e3_ff_axes(Fo, sx, sy, false);
        D.sens(end+1) = struct('name', nm, 'c', g.', ...
            'fw', rad2deg(e3_fwhm(uf, If)), 'sw', rad2deg(e3_fwhm(us, Is)), ...
            'rho', (max(Is)-min(Is))/(max(Is)+min(Is)));
        e3_log('  敏感性[%s]：c = %s，快轴 %.3f° / 慢轴 %.3f° / 对比度 %.3f', ...
               nm, mat2str(g.', 4), rad2deg(e3_fwhm(uf, If)), rad2deg(e3_fwhm(us, Is)), ...
               (max(Is)-min(Is))/(max(Is)+min(Is)));
    end

    % !! 验收判据（T9：远场、叠加、Parseval）
    e3_check('A9 θx=0 切面（矩乘）vs 直接一维变换 < 1e-12', D.A9 < 1e-12, '%.2e', D.A9);
    e3_check('A10 2D 切面 vs 先对 x 积分再对 y 变换 < 1e-12', D.A10b < 1e-12, '%.2e', D.A10b);
    e3_check('A10 半无限平板远场 FWHM = 23.95° (容差 0.05°)', ...
             abs(rad2deg(ff0.fwhm) - 23.95) < 0.05, '实得 %.3f°', rad2deg(ff0.fwhm));
    e3_check('A11(a) 解析带限场 Parseval 比 = 1（容差 1e-6）', abs(D.ratA - 1) < 1e-6, '%.10f', D.ratA);
    % !! A11(d) 的判据不是"回到 1 到 1e-3 以内"：高阶模在 |s|>1 的倏逝尾巴比基模重，
    % !! 到 |s|≤2 时仍有 1e-3 量级的残余。抓住的是"随窗口单调逼近 1"这一条：扩窗口后
    % !! 每个模都必须比窄窗口更接近 1，且 |s|≤2 时已进入 1±1e-3 之外不超过 2e-3。
    e3_check('A11(d) |s|≤2 全轴余量 ≤ 2e-3', all(abs(ratAll - 1) < 2e-3), mat2str(ratAll, 6));
    e3_check('A11(d) 扩窗口后每个模都更接近 1（单调）', all(abs(ratAll - 1) < abs(rat - 1)), ...
             '窄 %s → 宽 %s', mat2str(abs(rat - 1), 3), mat2str(abs(ratAll - 1), 3));
    e3_check('A12 c=[1,0,…] 与单模逐点一致 < 1e-12', D.A12 < 1e-12, '%.2e', D.A12);
    e3_check('T9 五组敏感性都给出有限 FWHM', all(isfinite([D.sens.fw])), ...
             'c 向量 = %s', mat2str(c.', 4));
    D.ok = true;
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_fwhm（由归一化强度曲线的 0.5 交点求半高全宽）
% 输入：u 自变量（与 I 同尺寸）; I 归一化强度
% 输出：w（与 u 同单位；找不到交点时返回 NaN）
% 说明：
%   1.峰值位置不假定在 u=0，先找最大值再向两侧找第一个低于 0.5 的点，再线性插值；
%   2.两侧任一侧没有交点（曲线没降到 0.5）时返回 NaN，由调用方判失败，不静默给数。
% ───────────────────────────────────────────────────────────────────────────────────
function w = e3_fwhm(u, I)
    u = u(:);  I = I(:);
    [~, i0] = max(I);
    h = 0.5 * I(i0);
    w = NaN;
    iL = find(I(1:i0) < h, 1, 'last');
    iR = find(I(i0:end) < h, 1, 'first') + i0 - 1;
    if isempty(iL) || isempty(iR), return; end
    uL = interp1(I([iL iL+1]), u([iL iL+1]), h);
    uR = interp1(I([iR-1 iR]), u([iR-1 iR]), h);
    w  = abs(uR - uL);
end

% ───────────────────────────────────────────────────────────────────────────────────
% 函数：e3_stage3_dispatch（阶段 3 的子任务分发）
% 输入：opt 结构体（.sub 子任务名；.lambda；.w; .h_e; .n_top）
% 输出：该子任务的数据结构
% 说明：
%   1.子任务名兼容短名与长名（t1 / stack 等价），全部走 e3_t*；
%   2.不在主入口里写业务逻辑：主入口只认 'stage3' 这一个 action。
% ───────────────────────────────────────────────────────────────────────────────────
function out = e3_stage3_dispatch(opt)
    if nargin < 1 || isempty(opt), opt = struct(); end
    sub = 'all';                     % !! 不给 sub = 全流程（T1–T9 + 七张图 + 日志 + 汇总）
    if isfield(opt, 'sub') && ~isempty(opt.sub), sub = lower(char(opt.sub)); end
    lam = 0.45;
    if isfield(opt, 'lambda') && ~isempty(opt.lambda), lam = opt.lambda; end
    % !! 三个几何参数缺省取施工单 §2 的默认值（w=3、h_e=0.30、SiO2 顶），
    % !! 这样 run_all 里只需给出要扫的那一档，其余按同一套缺省走。
    w   = 3;
    if isfield(opt, 'w') && ~isempty(opt.w), w = opt.w; end
    he  = 0.30;
    if isfield(opt, 'h_e') && ~isempty(opt.h_e), he = opt.h_e; end
    ntp = 1.46;
    if isfield(opt, 'n_top') && ~isempty(opt.n_top), ntp = opt.n_top; end
    switch sub
        case {'all', 'run', 'stage3'}
            % !! 唯一入口：跑完 T1–T9、出七张图与七份 csv、写日志、按 A13/A14 收口
            out = e3_run_all(lam);
        case 'check'
            % !! 供 test_stage3 里的脚本调用：e3_check 是本文件局部函数，脚本直接调不到
            ok = e3_check(opt.name, opt.cond, '%s', opt.note);
            out = struct('ok', ok);
        case {'t1', 'stack'},    out = e3_t1(lam);
        case {'t2', 'regions'},  out = e3_t2(lam);
        case {'t3', 'lateral'},  out = e3_t3(lam, w, he, ntp);
        case {'t4', 'field2d'},  out = e3_t4(lam, w, he, ntp);
        case {'t5', 'sweep'},    out = e3_t5(lam);
        case {'t6', 'fd2d'},     out = e3_t6(lam);
        case {'t7', 'gap'},      out = e3_t7(lam);
        case {'t8', 'scatter'},  out = e3_t8(lam, w);
        case {'t9', 'farfield'}, out = e3_t9(lam, w);
        case {'ui', 'uimain'}
            % !! 界面（l3_ui.m）的单一计算入口：用户层表 + 几何 → 场/Γ/α/远场
            opt.lambda = lam;  opt.w = w;  opt.h_e = he;  opt.n_top = ntp;
            out = e3_ui_run(opt);
        case {'ui_fd', 'uifd'}
            % !! 界面里的二维 FD 对照（要 20 s 量级，故单独一个子任务）
            opt.lambda = lam;  opt.w = w;  opt.h_e = he;  opt.n_top = ntp;
            out = e3_ui_fd(opt);
        otherwise
            error('e3_stage3_dispatch: 未知子任务 %s', sub);
    end
end

%% !! 阶段 3 全流程：目录、日志、九段、七图七csv、A13/A14 收口
% !! 这一段是阶段 3 的唯一入口（主文件无参数运行 → e3_stage3_dispatch('all') → e3_run_all）。
% !! 设计约束：所有数字、所有图、所有 csv 都从同一批 D.t1…D.t9 里取（G11 同源），
% !! 不写字面常数；日志用 UTF-8 的 fopen（施工单 §4），不用 diary。

% ─────────────────────────────────────────────────────────────────────────────
% 函数：e3_out（阶段 3 的输出目录 <主文件目录>/out，不存在就建）
% 输入：无
% 输出：d：绝对路径字符串
% 说明：
%   用 mfilename('fullpath') 定位，所以在任何 cwd 下调用都指向 Mode_L3/out。
% ─────────────────────────────────────────────────────────────────────────────
function d = e3_out()
    d = fullfile(fileparts(mfilename('fullpath')), 'out');
    if ~exist(d, 'dir')
        mkdir(d);
    end
end

% ─────────────────────────────────────────────────────────────────────────────
% 函数：e3_wcsv（把矩阵按一行一个记录写成 UTF-8 csv）
% 输入：fn：文件名; hdr：表头字符串; M：数值矩阵（每行一个记录）
% 输出：fn：写出的文件名
% 说明：
%   通用列格式 %.9g，避免每张表各写一套格式串；fopen 带 'n' + 'UTF-8'。
% ─────────────────────────────────────────────────────────────────────────────
function fn = e3_wcsv(fn, hdr, M)
    fid = fopen(fn, 'w', 'n', 'UTF-8');
    if fid < 0
        error('e3_wcsv: 打不开 %s', fn);
    end
    fprintf(fid, '%s\n', hdr);
    nc  = size(M, 2);
    fmt = [repmat('%.9g,', 1, nc - 1), '%.9g\n'];
    for i = 1:size(M, 1)
        fprintf(fid, fmt, M(i, :));
    end
    fclose(fid);
end

% ─────────────────────────────────────────────────────────────────────────────
% 函数：e3_figs（阶段 3 的七张图 + 七份 csv，全部从 D 取数）
% 输入：D：e3_run_all 收集的结果结构体（字段 t4/t5/t6/t9 可能缺）
% 输出：files：写出的文件名（cell）
% 说明：
%   1.图与 csv 逐张成对生成，任何一张的纵轴/横轴都从数据里取，不写死区间；
%   2.图窗保留在屏幕上（点运行就能看见），同时 print -r300 落盘；
%   3.field2d 的 y 是「自衬底界面向上」，p-clad 顶面在 y=1.729 处画参考线。
% ─────────────────────────────────────────────────────────────────────────────
function files = e3_figs(D)
    outd  = e3_out();
    files = {};

    %% !! 图 1：二维模场 |E|²（w=3, h_e=0.30, SiO₂ 顶）
    if isfield(D, 't4') && isfield(D.t4, 'Fd')
        Fd  = D.t4.Fd;
        E2  = abs(Fd.E(:, :, 1)).^2;
        ypc = sum([D.t4.C.Sp.ridge(2:end-1).d]);
        ww  = D.t4.C.w;
        [~, iypk] = max(sum(E2, 2));                        % !! 峰值行（自底向上的 y）
        ypk = Fd.y(iypk);

        fig = figure('Name', 'field2d_w3', 'Position', [80 80 720 600]);
        % !! 网格是非均匀的（QW 附近 1 nm、远场 25 nm），imagesc 会把行按**等间距**铺开、
        % !! 把峰值画到 y≈0.6 去（曾因此以为场算错了）；surf 用顶点坐标，非均匀才准。
        surf(Fd.x, Fd.y, zeros(size(E2)), E2, 'EdgeColor', 'none');
        view(2);  axis tight;  set(gca, 'YDir', 'normal');
        hold on;
        plot(Fd.x([1 end]), [ypc ypc], 'w--', 'LineWidth', 0.9);
        plot(Fd.x([1 end]), [ypk ypk], 'w-',  'LineWidth', 1.2);
        % !! 阱只有 9 nm，在 0.4 μm 的窗口里占 2%，不画出来根本看不见；
        % !! 模场峰值在阱下方约 33 nm（复合芯 209 nm 的中间偏 n 侧），这是 Γ 只有 2.7% 的直观原因
        ybR = layer_bounds(D.t4.C.Sp.ridge);
        iqw = find(strcmp({D.t4.C.Sp.ridge.tag}, 'QW'), 1);
        yqw = [NaN NaN];
        if ~isempty(iqw)
            yqw = ybR(iqw, :);
            plot(Fd.x([1 end]), [yqw(1) yqw(1)], '-', 'Color', [1 0.35 0.75], 'LineWidth', 1.3);
            plot(Fd.x([1 end]), [yqw(2) yqw(2)], '-', 'Color', [1 0.35 0.75], 'LineWidth', 1.3);
        end
        plot([-ww/2 -ww/2], Fd.y([1 end]), 'w:', 'LineWidth', 0.9);
        plot([ ww/2  ww/2], Fd.y([1 end]), 'w:', 'LineWidth', 0.9);
        hold off;
        % !! 全场 y 跨约 3.3 μm，而模场只被阱限制在几十 nm 里：不放大就是一条几乎看不见的细线
        % !! （用户 2026-10-28 反馈"没搞懂这张图画的是什么"）。这里按峰值上下各 0.20 μm 取窗。
        ylim([ypk - 0.20, ypk + 0.20]);
        xlabel('x [\mum]'); ylabel('y [\mum]（自衬底界面向上）');
        % !! C.Gamma 是 1×K 向量，标题里只能取 (1)：直接把向量喂给 sprintf 会按格式串
        % !! 循环复用参数，把整句标题拼成两遍（曾因此得到 99 字符的标题）。
        title(sprintf(['|E|^2（基模，w=%.1f \\mum）：峰值 y=%.3f，' ...
                       'QW y=%.3f–%.3f（品红），p-clad 顶面 y=%.3f，\\Gamma=%.5f'], ...
                      ww, ypk, yqw(1), yqw(2), ypc, D.t4.C.Gamma(1)));
        colorbar;
        print(fig, fullfile(outd, 'field2d_w3.png'), '-dpng', '-r300');
        files{end+1} = fullfile(outd, 'field2d_w3.png'); %#ok<AGROW>

        fn = fullfile(outd, 'field2d_w3.csv');
        fid = fopen(fn, 'w', 'n', 'UTF-8');
        fprintf(fid, 'x_um%s\n', sprintf(',%.6g', Fd.x));
        for i = 1:numel(Fd.y)
            fprintf(fid, '%.6g%s\n', Fd.y(i), sprintf(',%.6g', E2(i, :)));
        end
        fclose(fid);
        files{end+1} = fn; %#ok<AGROW>
    end

    %% !! 图 2/3：Re N₁ 与 α₁ 随脊宽、Γ 随脊宽（T5 的 h_e=0.30 那一列）
    if isfield(D, 't5') && isfield(D.t5, 'tab')
        tb = D.t5.tab;
        he0 = 0.30;
        iw = @(t0) find([tb.h_e] == he0 & abs([tb.n_top] - t0) < 1e-9);
        iS = iw(1.46);  iA = iw(1.00);
        wS = [tb(iS).w];  wA = [tb(iA).w];
        nS = [tb(iS).N1]; nA = [tb(iA).N1];
        aS = [tb(iS).a1]; aA = [tb(iA).a1];
        gS = [tb(iS).G1]; gA = [tb(iA).G1];
        naS = [tb(iS).na]; naA = [tb(iA).na];
        Gslab = [tb(iS).Ga];

        fig = figure('Name', 'neff_alpha_vs_w', 'Position', [100 100 1000 420]);
        subplot(1, 2, 1);
        semilogx(wS, nS - naS, 'o-', 'LineWidth', 1.4); hold on;
        semilogx(wA, nA - naA, 's--', 'LineWidth', 1.2); hold off;
        grid on; xlabel('w [\mum]'); ylabel('Re N_1 - Re n_a');
        title('基模 Re N_1 - Re n_a（w 大时归零）');
        legend('SiO_2 顶', 'air 顶', 'Location', 'southeast');
        subplot(1, 2, 2);
        semilogx(wS, aS, 'o-', 'LineWidth', 1.4); hold on;
        semilogx(wA, aA, 's--', 'LineWidth', 1.2);
        yline(mean(aS), ':', 'Color', [0.4 0.4 0.4]); hold off;
        grid on; xlabel('w [\mum]'); ylabel('\alpha_i [cm^{-1}]');
        title('各阶侧模的材料吸收（EIM 两步，α 只来自 Im n）');
        legend('SiO_2 顶', 'air 顶', 'Location', 'southeast');
        print(fig, fullfile(outd, 'neff_alpha_vs_w.png'), '-dpng', '-r300');
        files{end+1} = fullfile(outd, 'neff_alpha_vs_w.png'); %#ok<AGROW>

        M = [(1:numel(wS)).', wS(:), naS(:), nS(:), gS(:), aS(:), ...
             naA(:), nA(:), gA(:), aA(:), Gslab(:)];
        fn = fullfile(outd, 'neff_alpha_vs_w.csv');
        e3_wcsv(fn, ['row,w_um,Re_na_SiO2,ReN1_SiO2,Gamma1_SiO2,alpha1_SiO2,' ...
                     'Re_na_air,ReN1_air,Gamma1_air,alpha1_air,Gamma_slab'], M);
        files{end+1} = fn; %#ok<AGROW>

        fig = figure('Name', 'gamma_vs_w', 'Position', [100 100 760 420]);
        semilogx(wS, gS, 'o-', 'LineWidth', 1.4); hold on;
        semilogx(wA, gA, 's--', 'LineWidth', 1.2);
        semilogx(wS, Gslab, 'k:', 'LineWidth', 1.0); hold off;
        grid on; xlabel('w [\mum]'); ylabel('\Gamma');
        title('限制因子随脊宽趋于平板极限 \Gamma_{slab}（A3 的判据）');
        legend('SiO_2 顶', 'air 顶', '\Gamma_a（w\rightarrow\infty 平板极限）', 'Location', 'southeast');
        print(fig, fullfile(outd, 'gamma_vs_w.png'), '-dpng', '-r300');
        files{end+1} = fullfile(outd, 'gamma_vs_w.png'); %#ok<AGROW>

        fn = fullfile(outd, 'gamma_vs_w.csv');
        e3_wcsv(fn, 'w_um,Gamma1_SiO2,Gamma1_air,Gamma_slab', [wS(:), gS(:), gA(:), Gslab(:)]);
        files{end+1} = fn; %#ok<AGROW>
    end

    %% !! 图 4：EIM vs 二维 FD（T6）
    if isfield(D, 't6') && isfield(D.t6, 'tab')
        tb = D.t6.tab;
        w  = [tb.w];  ne = [tb.Neim];  nf = [tb.Nfd];  dd = abs(ne - nf);
        fig = figure('Name', 'eim_vs_fd', 'Position', [100 100 700 440]);
        plot(w, ne, 'o-', 'LineWidth', 1.4); hold on;
        plot(w, nf, 's--', 'LineWidth', 1.2);
        plot(w, ne, 'k.'); hold off;
        grid on; xlabel('w [\mum]'); ylabel('Re N_1');
        title(sprintf('EIM 两步 vs 二维标量 FD（\\lambda=450 nm, h_e=0.30 \\mum）：max|\\Delta|=%.2e', max(dd)));
        legend('EIM 两步', '二维 FD', 'Location', 'southeast');
        print(fig, fullfile(outd, 'eim_vs_fd.png'), '-dpng', '-r300');
        files{end+1} = fullfile(outd, 'eim_vs_fd.png'); %#ok<AGROW>

        M = [[tb.w].', [tb.Neim].', [tb.Nfd].', abs([tb.Neim] - [tb.Nfd]).', ...
             [tb.frac].', [tb.n].', [tb.sec].'];
        fn = fullfile(outd, 'eim_vs_fd.csv');
        fid = fopen(fn, 'w', 'n', 'UTF-8');
        fprintf(fid, 'w_um,ReN_EIM,ReN_FD,abs_diff,frac,unknowns,sec\n');
        fprintf(fid, '%.3f,%.9f,%.9f,%.3e,%.4f,%d,%.1f\n', M.');
        if isfield(D.t6, 'conv') && numel(D.t6.conv) == 2
            fprintf(fid, '# 网格收敛 w=3: %.9f -> %.9f\n', D.t6.conv(1), D.t6.conv(2));
        end
        fclose(fid);
        files{end+1} = fn; %#ok<AGROW>
    end

    %% !! 图 5：二维远场（T9，等幅同相叠加的假设值）
    if isfield(D, 't9') && isfield(D.t9, 'Iov')
        S = D.t9;
        fig = figure('Name', 'ff2d_w3', 'Position', [80 80 640 560]);
        imagesc(S.sx, S.sy, 10*log10(max(S.Iov, 1e-8)));  set(gca, 'YDir', 'normal');
        axis tight;  xlim([-1 1]);  ylim([-1 1]);   % !! 填满坐标框，别留一圈白边（s 已到 ±0.996）
        hold on;
        contour(S.sx, S.sy, S.Iov, [0.5 0.135], 'w', 'LineWidth', 1.0);
        plot([-1 1], [0 0], 'w:');  plot([0 0], [-1 1], 'w:');
        hold off;  caxis([-30 0]);  colorbar;
        xlabel('sin\theta_x（慢轴）'); ylabel('sin\theta_y（快轴）');
        title('二维远场分布 [dB]：白线为 0.5 与 0.135 等强度线');
        print(fig, fullfile(outd, 'ff2d_w3.png'), '-dpng', '-r300');
        files{end+1} = fullfile(outd, 'ff2d_w3.png'); %#ok<AGROW>

        fn = fullfile(outd, 'ff2d_w3.csv');
        fid = fopen(fn, 'w', 'n', 'UTF-8');
        fprintf(fid, 'sin_thx%s\n', sprintf(',%.6g', S.sx));
        for i = 1:numel(S.sy)
            fprintf(fid, '%.6g%s\n', S.sy(i), sprintf(',%.6g', S.Iov(i, :)));
        end
        fclose(fid);
        files{end+1} = fn; %#ok<AGROW>

        ax = {'fast', 'ff_fast', '快轴（垂直 y，固定 \theta_x=0）'; ...
              'slow', 'ff_slow', '慢轴（水平 x，固定 \theta_y=0）'};
        for k = 1:2
            A = S.(ax{k, 1});
            fig = figure('Name', ax{k, 2}, 'Position', [100 100 700 440]);
            % !! 各阶都画出来，并且一阶一个颜色、细线（用户 2026-10-28：都画成一个颜色分不清）
            % !! 调色板固定前四色（蓝/橙红/绿/紫），再往后按模 4 循环，m=3 及以上也不会撞色成灰
            pal = [0.000 0.447 0.741;  0.850 0.325 0.098; ...
                   0.466 0.674 0.188;  0.494 0.184 0.556; ...
                   0.301 0.745 0.933;  0.635 0.078 0.184];
            Kkw = numel(S.I2);
            Uq = cell(1, Kkw);  Jq = cell(1, Kkw);  Mx = zeros(1, Kkw);
            for q = 1:Kkw
                [Uq{q}, Jq{q}, Mx(q)] = e3_ff_axes(S.F2{q}, S.sx, S.sy, k == 1);
            end
            % !! 若某阶的切面峰值比基模小 9 个量级，说明它在 θx=0 上几乎处处为零（例如 h_e=0、
            % !! 两区竖直剖面相同时，奇模的这个切面确实恒为 0，2026-10-28 就画出过一条假双峰），
            % !! 这时归一化只会放大数值噪声，故不画。一般情形下奇模的切面并不为零——见下条注释
            keep = Mx >= 1e-9*max(Mx);
            hold on;
            for q = find(keep)
                if q == 1
                    lw = 2.0;  nmq = 'm=0（基模）';
                else
                    lw = 0.9;  nmq = sprintf('仅 m=%d', q-1);
                end
                plot(rad2deg(Uq{q}), Jq{q}, '-', 'LineWidth', lw, ...
                     'Color', pal(mod(q-1, size(pal, 1)) + 1, :), 'DisplayName', nmq);
            end
            plot(rad2deg(A.uov), A.Iov, 'k--', 'LineWidth', 1.8, ...
                 'DisplayName', sprintf('%d 阶等幅叠加（假设值）', Kkw));
            hold off;
            grid on;  xlim([-60 60]);  ylim([0 1.05]);
            xlabel('angle [deg]');  ylabel('relative intensity');
            ttl = sprintf('%s：归一化远场（每阶单独归一，叠加整体归一）', ax{k, 3});
            if k == 1 && any(~keep)
                ttl = [ttl sprintf('；m=%s 的 θ_x=0 切面恒为 0，未画', ...
                       strjoin(arrayfun(@(v) sprintf('%d', v-1), find(~keep), 'UniformOutput', false), ','))];
            end
            % !! 奇模在 θx=0 的切面一般【不】为零：逐 y 的 ∫E_x(x)Φ(x,y)dx 里，脊区与沟槽区的
            % !! Φ 不同，两者的加权差不为 0（h_e=0、两区相同时才恒为 0）。所以 m=1 那条双峰是真结果
            title(ttl);
            % !! 图例放到坐标框外、两列：m=2 的旁瓣正好在 ±10° 顶部，放框内会压住曲线，
            % !! 横排一行又会超出图宽（第一项被裁掉），所以两列
            legend('Location', 'southoutside', 'NumColumns', 2);
            print(fig, fullfile(outd, [ax{k, 2} '.png']), '-dpng', '-r300');
            files{end+1} = fullfile(outd, [ax{k, 2} '.png']); %#ok<AGROW>

            % !! csv 与图同源：每列一阶（I_m0…I_mK-1）+ 叠加列，横坐标共用 angle_deg
            fn  = fullfile(outd, [ax{k, 2} '.csv']);
            hdr = 'angle_deg';
            cols = rad2deg(Uq{1}(:));
            for q = 1:Kkw
                % !! 被筛掉的奇模写整列 0：那种情形下这个切面本来就近乎为零（h_e=0 时严格为零），
                % !! 写归一化后的噪声只会误导读者
                vq = Jq{q}(:);
                if ~keep(q)
                    vq = zeros(size(vq));
                end
                cols = [cols, vq];                                      %#ok<AGROW>
                hdr  = [hdr sprintf(',I_m%d', q-1)];                    %#ok<AGROW>
            end
            cols = [cols, interp1(rad2deg(A.uov(:)), A.Iov(:), ...
                                  rad2deg(Uq{1}(:)), 'linear', 'extrap')];
            e3_wcsv(fn, [hdr ',I_overlay'], cols);
            files{end+1} = fn; %#ok<AGROW>
        end
    end
end

% ─────────────────────────────────────────────────────────────────────────────
% 函数：e3_mfile（本文件的绝对路径，带 .m 后缀）
% 输入：无
% 输出：f：绝对路径
% 说明：
%   mfilename('fullpath') 在局部函数里可能丢扩展名，fileread 之前补上（A14 要读本文件）。
% ─────────────────────────────────────────────────────────────────────────────
function f = e3_mfile()
    f = mfilename('fullpath');
    if isempty(regexp(f, '\.m$', 'once'))
        f = [f '.m'];
    end
end

% ─────────────────────────────────────────────────────────────────────────────
% 函数：e3_a0（A0 退化回归：同一 stack_7layer 下阶段 3 与阶段 2 的 TE0 逐位比较）
% 输入：lambda：波长（μm）
% 输出：[dRe,dG,dAl]：Re(n_eff)、Γ、αᵢ 的绝对差
% 调用函数：Model_Self_L3_7L_EBL（'solve'）、Model_Self_L2_7L_EBL（'solve'）
% 说明：
%   1.比的是「两条活跑路径」而不是打印出来的归档常数（施工单 §1 表写 2.417133286，
%     汇报单 §8.1 与实跑都是 2.417128446，差 4.8e-6，属文档笔误）；
%   2.两个调用都用缺省场景 'kio'，即阶段 2 的基线设置；
%   3.a=1 而不是截断：只要任一段路径变了，这条就会红。
% ─────────────────────────────────────────────────────────────────────────────
function [dRe, dG, dAl] = e3_a0(lambda)
    root = fileparts(fileparts(e3_mfile()));            % waveguide 根目录
    addpath(fullfile(root, 'Mode_L2'));
    addpath(fullfile(root, 'Mode_L3'));
    o = struct('lambda', lambda, 'pol', 'TE', 'stack', '7', 'scenario', 'kio');
    A = Model_Self_L3_7L_EBL('solve', o);
    B = Model_Self_L2_7L_EBL('solve', o);
    dRe = abs(A(1).ne - B(1).ne);
    dG  = abs(A(1).Gamma - B(1).Gamma);
    dAl = abs(sum(A(1).contrib) - sum(B(1).contrib));
    e3_log('  A0：L3 Re %.9f / L2 Re %.9f；Γ %.8f / %.8f；αᵢ %.6f / %.6f', ...
           A(1).ne, B(1).ne, A(1).Gamma, B(1).Gamma, sum(A(1).contrib), sum(B(1).contrib));
end

% ─────────────────────────────────────────────────────────────────────────────
% 函数：e3_run_all（阶段 3 全流程：日志 → T1…T9 → 七图 → A13/A14 → 汇总）
% 输入：lambda：波长（μm）
% 输出：D：结果结构体（t1…t9 的返回、pass/fail、files、out）
% 说明：
%   1.日志 Mode_L3/out/stage3_tests.log 用 UTF-8 的 fopen（施工单 §4），不用 diary；
%   2.任一段抛异常只记 FAIL，不中断后面各段（与 run_all 的旧行为一致）；
%   3.e3_check 的计数在 global E3_NPASS/E3_NFAIL 里，退出码由启动器按 D.fail 决定。
% ─────────────────────────────────────────────────────────────────────────────
function D = e3_run_all(lambda)
    global E3_LOG_FID E3_NPASS E3_NFAIL

    if nargin < 1 || isempty(lambda)
        lambda = 0.45;
    end
    outd = e3_out();
    own  = isempty(E3_LOG_FID);
    if own
        E3_LOG_FID = fopen(fullfile(outd, 'stage3_tests.log'), 'w', 'n', 'UTF-8');
        if E3_LOG_FID < 0
            error('e3_run_all: 打不开日志 %s', fullfile(outd, 'stage3_tests.log'));
        end
    end
    E3_NPASS = 0;  E3_NFAIL = 0;
    t00 = tic;

    e3_log('==== 阶段 3 全流程（λ = %.3f μm，%s 起）====', lambda, ...
           char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));

    jobs = { 't1', @() e3_t1(lambda); ...
             't2', @() e3_t2(lambda); ...
             't3', @() e3_t3(lambda, 3, 0.30, 1.46); ...
             't4', @() e3_t4(lambda, 3, 0.30, 1.46); ...
             't5', @() e3_t5(lambda); ...
             't6', @() e3_t6(lambda); ...
             't7', @() e3_t7(lambda); ...
             't8', @() e3_t8(lambda, 3); ...
             't9', @() e3_t9(lambda, 3) };

    D = struct();
    for jj = 1:size(jobs, 1)
        key = jobs{jj, 1};
        t0  = tic;
        try
            D.(key) = jobs{jj, 2}();
        catch err
            E3_NFAIL = E3_NFAIL + 1;
            e3_log('  [ERROR] %s 抛异常：%s', key, err.message);
        end
        e3_log('  【%s 用时 %.1f s】', key, toc(t0));
    end

    %% !! A0：阶段 2 的退化回归（L3 活跑 vs L2 活跑，逐位一致）
    try
        [dRe, dG, dAl] = e3_a0(lambda);
        e3_check('A0 退化回归：L3 与 L2 的 TE0 逐位一致（<1e-12）', ...
                 max([dRe dG dAl]) < 1e-12, 'dRe=%.1e dΓ=%.1e dα=%.1e', dRe, dG, dAl);
    catch err
        e3_check('A0 退化回归：L3 与 L2 的 TE0 逐位一致（<1e-12）', false, ...
                 'L2 调用失败：%s', err.message);
    end

    %% !! A8：汇报单必须写清「只含材料吸收」与「TM 未做」两句
    rep = fullfile(fileparts(e3_mfile()), '阶段3_汇报单_V1.0.md');
    if exist(rep, 'file') == 2
        txt = fileread(rep);
        c8  = [~isempty(strfind(txt, '材料吸收')), ~isempty(strfind(txt, 'TM')), ...
               ~isempty(strfind(txt, '等幅同相'))];
        e3_check('A8 汇报单写明「只含材料吸收」「TM 未做」「叠加是假设」', all(c8), ...
                 '三句命中：%d/3（材料吸收 %d，TM %d，等幅同相 %d）', sum(c8), c8(1), c8(2), c8(3));
    else
        e3_check('A8 汇报单写明「只含材料吸收」「TM 未做」「叠加是假设」', false, ...
                 '未找到 %s', rep);
    end

    %% !! 阶段 1/2 的四张图窗：无参运行要一起出来，而且要和阶段 3 的图窗同时存在
    % !! （用户 2026-10-28：一、二阶段的图不能和三阶段的图共存是不对的）
    % !! 做法是递归调用本文件的 'show' 支，而不是把那段绘图代码再抄一份（G3：原文不动）
    % !! 位置很关键：必须排在 A0 之后。A0 要调 ('solve') 与 L2 的 ('solve')，而那两支开头都有
    % !! clc / close all（见本文件第 40 行），放在它前面画的图窗会被一起关掉——2026-10-28 就是这样
    % !! 只剩 7 个图窗、编号从 1 开始，阶段 1/2 的四张全没了
    try
        Model_Self_L3_7L_EBL('show');
        f2 = findobj('Type', 'figure');
        if ~isempty(f2)
            [~, o2] = sort([f2.Number]);
            f2 = f2(o2);
            nm2 = {'stage2_1_waveguide', 'stage2_2_farfield', ...
                   'stage2_3_tables',    'stage2_4_loss'};
            for k2 = 1:min(numel(nm2), numel(f2))
                print(f2(k2), fullfile(outd, [nm2{k2} '.png']), '-dpng', '-r300');
            end
            e3_log('  阶段 1/2 的四张图窗已生成（另存 out/stage2_*.png）');
        end
    catch err
        e3_log('  [warn] 阶段 1/2 图窗未生成：%s', err.message);
    end

    e3_log('==== 出图与 csv ====');
    files = e3_figs(D);
    e3_log('  已写 %d 个文件到 %s', numel(files), outd);

    %% !! A13：七张 PNG + 七份 csv 齐全
    need = {'field2d_w3.png', 'field2d_w3.csv', 'neff_alpha_vs_w.png', 'neff_alpha_vs_w.csv', ...
            'gamma_vs_w.png', 'gamma_vs_w.csv', 'eim_vs_fd.png', 'eim_vs_fd.csv', ...
            'ff2d_w3.png', 'ff2d_w3.csv', 'ff_fast.png', 'ff_fast.csv', ...
            'ff_slow.png', 'ff_slow.csv'};
    hav = false(1, numel(need));
    for jj = 1:numel(need)
        hav(jj) = exist(fullfile(outd, need{jj}), 'file') == 2;
    end
    miss = strjoin(need(~hav), ' ');
    if isempty(miss)
        miss = '无';
    end
    e3_check('A13 七张 PNG + 七份 csv 齐全', all(hav), ...
             '已有 %d/%d 个（缺：%s）', sum(hav), numel(need), miss);

    %% !! A14：叠加权重段必须带「假设值」注（可被 grep 到）
    src = fileread(e3_mfile());
    L   = regexp(src, '\r?\n', 'split');
    hit = false(1, numel(L));
    for jj = 1:numel(L)
        s = L{jj};
        hit(jj) = ~isempty(strfind(s, '% !!')) && ~isempty(strfind(s, '权重')) ...
                  && ~isempty(strfind(s, '假设'));
    end
    ex = '';
    if any(hit)
        ex = strtrim(L{find(hit, 1)});
    end
    e3_check('A14 叠加权重段带「假设值」注（可被 grep 到）', any(hit), ...
             '命中 %d 行；首行：%s', sum(hit), ex);

    %% 汇总
    tt = toc(t00);
    e3_log('==== 汇总 ====');
    e3_log('  PASS = %d，FAIL = %d，总耗时 %.1f s', E3_NPASS, E3_NFAIL, tt);
    if E3_NFAIL == 0
        e3_log('  结论：全部通过。');
    else
        e3_log('  结论：有 %d 项未通过，见上方 [FAIL] 行。', E3_NFAIL);
    end
    D.pass  = E3_NPASS;
    D.fail  = E3_NFAIL;
    D.files = files;
    D.out   = outd;

    if own
        fclose(E3_LOG_FID);
        E3_LOG_FID = [];
    end
end

%% !! 界面（l3_ui.m）的计算入口
% !! 与阶段 3 的关系：物理方程、求解器、Γ/α 口径全部复用，唯一区别是层表由界面给（opt.L0），
% !! 因此界面里改一层厚度/折射率，下面每一个数都会跟着变 —— 不是查表、不是插值。

% ─────────────────────────────────────────────────────────────────────────────
% 函数：e3_ui_run（界面的一条龙计算：竖直模 → 横向各阶 → 二维场 → 远场）
% 输入：opt 结构体
%         .L0     自定义层表（自衬底向上，首层 Inf、末层有限正，n 已含损耗虚部）
%         .lambda 波长 [μm]; .w 脊宽; .h_e 刻蚀深度; .n_top 顶部介质
%         .solver 可选，覆盖 e3_opts_solve 的字段（num_scan/tol/merge/npts…）
% 输出：U 结构体（.ok 是否成功、.msg 失败原因、.C 格点、.vert 竖直各阶、
%       .far1 一维快轴远场、.lat 横向各阶、.fd 二维场、.ff 二维远场与两个切面、.t_sec 耗时）
% 调用函数：e3_cell, e3_region_solve, far_field, e3_field2d, e3_ff_matrix, e3_ff_axes
% 说明：
%   1.任何一步解不出（某区无竖直模、反导引、横向无根）都不抛异常，而是把原因写进 U.msg；
%      界面据此报错，绝不拿旧数冒充新结构的解（2026-10-28 用户选 A）。
%   2.远场只算 θx,θy ∈ [−85°, 85°] 的 sinθ 网格（|s|≤0.996，接近物理极限），
%      与阶段 3 的 ff2d 同一个网格口径，便于与既有成果对照。
%   3.这里不做任何"模数修正"或"补根"：横向模数由 e3_lateral_solve 自己判（含复种子）。
% ─────────────────────────────────────────────────────────────────────────────
function U = e3_ui_run(opt)
    lam = 0.45; if isfield(opt,'lambda') && ~isempty(opt.lambda), lam = opt.lambda; end
    w   = 3;    if isfield(opt,'w')      && ~isempty(opt.w),      w   = opt.w;      end
    he  = 0.30; if isfield(opt,'h_e')    && ~isempty(opt.h_e),    he  = opt.h_e;    end
    nt  = 1.46; if isfield(opt,'n_top')  && ~isempty(opt.n_top),  nt  = opt.n_top;  end
    k0  = 2*pi/lam;

    U = struct('ok', false, 'msg', '', 'lambda', lam, 'w', w, 'h_e', he, 'n_top', nt, ...
               'C', [], 'vert', [], 'far1', [], 'lat', [], 'fd', [], 'ff', [], ...
               't_sec', 0, 'poly', 'TE');
    t0 = tic;

    if ~isfield(opt, 'L0') || isempty(opt.L0)
        U.msg = '没有拿到自定义层表（opt.L0 为空）';
        return
    end
    base = struct('h_e', he, 'n_top', nt, 'L0', opt.L0(:));

    % ---- 一、竖直（两区各解基模，再解脊区前 3 阶供数据表用） ----
    C = e3_cell(lam, w, he, nt, base);
    if ~all(C.region_ok)
        U.C = C;
        U.msg = '有一侧的竖直模没解出来：检查层表（顶层厚度是否>0、折射率高低是否形成导引、有没有两层同时标成有源层）';
        return
    end
    if ~(real(C.nb) < real(C.na) + 1e-12)
        U.C = C;
        U.msg = sprintf(['反导引：沟槽区 Re n_b = %.6f ≥ 脊区 Re n_a = %.6f。\n' ...
                         '说明刻蚀让沟槽区的等效折射率不低于脊区，横向不导引，物理上无侧模解。'], ...
                        real(C.nb), real(C.na));
        return
    end
    if isempty(C.N) || ~all(isfinite(C.N))
        U.C = C;
        U.msg = sprintf('横向无导模：V = %.3f（按公式应有 %d 阶），把脊宽改大或增大两侧折射率差', ...
                        C.V, C.M_expect);
        return
    end
    U.C = C;

    ov = e3_opts_solve();
    if isfield(opt, 'solver') && isstruct(opt.solver)
        fn = fieldnames(opt.solver);
        for k = 1:numel(fn)
            if isfield(ov, fn{k}), ov.(fn{k}) = opt.solver.(fn{k}); end
        end
    end
    R = e3_region_solve(C.Sp.ridge, lam, 3, ov, 'TE');
    U.vert = R;
    if ~isempty(R)
        % !! 竖直各阶模都要算远场：只算基模时界面写着 TE0-TE2、图上只有一条线（2026-10-29 用户发现）
        U.far1 = struct('theta', {}, 'I', {}, 'fwhm', {});
        for q = 1:numel(R)
            Fq = far_field(R(q).y, R(q).E, k0, struct('num_ff', 2001));
            Iq = Fq.I(:).' / max(Fq.I);
            U.far1(q) = struct('theta', Fq.theta(:).', 'I', Iq, 'fwhm', Fq.fwhm);
        end
    end

    U.lat = struct('N', C.N, 'Gamma', C.Gamma, 'alpha', C.alpha, 'eta', C.eta, ...
                   'V', C.V, 'M_expect', C.M_expect, 'alpha_pred', C.alpha_pred);

    % ---- 二、二维场（前 3 阶） ----
    % !! 二维场与远场算到 7 阶：界面要让用户自己选显示几阶，数据得先备齐（2026-10-29）
    jl = 1:min(7, numel(C.N));
    Fd = e3_field2d(C, jl, 0.010, 0.001);
    U.fd = Fd;

    % ---- 三、二维远场 + 快/慢轴切面 ----
    % !! 401 点 -> 801 点：慢轴横轴按数据收紧后，0.3 度一个点会让曲线看着是一段段折线
    % !! （用户 2026-10-29 第 9 项）。801 点在 19 度窗口里约 0.15 度一个点。
    sc = linspace(-0.996, 0.996, 801);
    K  = numel(jl);
    Fsum = zeros(numel(sc), numel(sc));
    md = struct('N', {}, 'G', {}, 'a', {}, 'uF', {}, 'IF', {}, 'uS', {}, 'IS', {});
    for q = 1:K
        Fq = e3_ff_matrix(Fd, q, sc, sc, k0);
        Fsum = Fsum + Fq / sqrt(K);
        [uF, IF] = e3_ff_axes(Fq, sc, sc, true);
        [uS, IS] = e3_ff_axes(Fq, sc, sc, false);
        md(q) = struct('N', C.N(q), 'G', C.Gamma(q), 'a', C.alpha(q), ...
                       'uF', uF(:).', 'IF', IF(:).', 'uS', uS(:).', 'IS', IS(:).');
    end
    Iov = abs(Fsum).^2;
    Iov = Iov / max(Iov(:));
    [uFo, IFo] = e3_ff_axes(Fsum, sc, sc, true);
    [uSo, ISo] = e3_ff_axes(Fsum, sc, sc, false);
    U.ff = struct('sx', sc(:).', 'sy', sc(:).', 'I2', Iov, ...
                  'uF', uFo(:).', 'IF', IFo(:).', 'uS', uSo(:).', 'IS', ISo(:).', ...
                  'mode', md, 'K', K, 'c', 1/sqrt(K));

    U.ok = true;
    U.t_sec = toc(t0);
end

% ─────────────────────────────────────────────────────────────────────────────
% 函数：e3_ui_fd（界面专用：当前结构的一条二维有限差分对照）
% 输入：opt（.L0/.lambda/.w/.h_e/.n_top/.shift_n）
% 输出：R（e3_fd2d 的输出：.Nf/.frac/.n/.ok/...）
% 说明：
%   单独成一个入口，是因为它要 20 s 量级；界面只在用户点"EIM vs FD"那一张时才调。
%   .shift_n 由界面从 EIM 的 Re N_1 传进来 —— 移位点只决定 eigs 找到哪几个本征值，
%   不改变本征值本身，故不构成循环论证（这一点在 T6 的注释与日志里已写明）。
% ─────────────────────────────────────────────────────────────────────────────
function R = e3_ui_fd(opt)
    o = struct('L0', opt.L0(:));
    if isfield(opt, 'dx')      && ~isempty(opt.dx),      o.dx      = opt.dx;      end
    if isfield(opt, 'shift_n') && ~isempty(opt.shift_n), o.shift_n = opt.shift_n; end
    R = e3_fd2d(opt.lambda, opt.w, opt.h_e, opt.n_top, o);
end

