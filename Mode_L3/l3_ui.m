function h = l3_ui(mode)
% ───────────────────────────────────────────────────────────────────────────────────
% 函数：l3_ui（阶段 3 参数界面：三列结构编辑 + 三维预览 + 按当前结构重算的成果画布）
% 输入：mode：可选。不给 = 打开界面；'selftest' = 无界面地把 11 张图各画一遍并落盘（自检用）
% 输出：h：结构体。h.fig 为 uifigure 句柄；h.handles 为控件句柄；h.refresh/@h.state 供自动化测试
% 调用函数：Model_Self_L3_7L_EBL（'stage3' 的 sub='ui' / 'ui_fd'）
% 说明：
%   1. 第一列管结构：波导类型、层表（可增删，首末层不可删）、慢/快/纵轴长、脊形特有的脊宽与脊高。
%   2. 第二列管选中层的参数：tag、折射率、损耗、厚度、是否 QW；左上角「高级设置」改 opt 与顶部介质。
%   3. 第三列三维预览：每层不同颜色，可拖动旋转、滚轮缩放。
%   4. 任何输入变化立刻校验，非法时 Sure 置灰并在第一列底部写原因，不弹窗、不阻塞。
%   5. Sure 之后成果画布的每一张图都是用当前结构现算的（用户 2026-10-28 选 A：结构真的进计算，
%      无解就报错），不读任何旧图；矩形波导没有横向限制，横向相关的几张会写明不适用。
% ───────────────────────────────────────────────────────────────────────────────────

    if nargin >= 1 && (ischar(mode) || isstring(mode)) && strcmpi(char(mode), 'selftest')
        h = l3_selftest();
        return
    end

    S = lui_live(lui_default());            % !! 补上拖动排序与三维相机这些"运行时"字段
    Ucur = [];                              % 当前结构算出来的成果（结构一改就作废）
    rowBtn = gobjects(0);                   % 层表每行的"名称按钮"，拖动排序要拿它的屏幕坐标
    hfig = uifigure('Name', 'GaN 蓝光激光器 · 波导结构设置与阶段成果', ...
                    'Position', [50 30 1620 950], 'Color', [0.96 0.96 0.97]);
    % !! 层表拖动排序要靠图窗级的鼠标事件（uibutton 只有"松开"回调，拿不到"按下"）
    hfig.WindowButtonDownFcn = @(~, ~) on_down();
    hfig.WindowButtonMotionFcn = @(~, ~) on_motion();
    hfig.WindowButtonUpFcn = @(~, ~) on_up();
    % !! 滚轮缩放也自己接（不要 zoomInteraction：它会把相机推进盒子里）
    hfig.WindowScrollWheelFcn = @(~, e) on_wheel(e);

    %% ===== 根布局 =====
    root = uigridlayout(hfig, [2 1]);
    root.RowHeight = {'1x', 54};
    root.ColumnWidth = {'1x'};
    root.Padding = [10 8 10 8];

    mainGrid = uigridlayout(root, [1 3]);
    mainGrid.Layout.Row = 1;
    mainGrid.ColumnWidth = {360, 360, '1x'};
    mainGrid.ColumnSpacing = 10;

    %% ===== 第一列：结构 =====
    p1 = uipanel(mainGrid, 'Title', '① 波导结构', 'FontWeight', 'bold', 'FontSize', 12);
    g1 = uigridlayout(p1, [6 1]);
    g1.RowHeight = {30, '1x', 30, 152, 24, 60};
    g1.Padding = [8 8 8 8];
    g1.RowSpacing = 6;

    ddKind = uidropdown(g1, 'Items', {'脊形波导', '矩形波导'}, 'Value', '脊形波导', ...
                        'ValueChangedFcn', @(s, ~) on_kind(s.Value));

    pLay = uipanel(g1, 'Title', '层表（自上而下：顶层 … 衬底）', 'FontSize', 11);
    pLay.Scrollable = 'on';
    glRows = uigridlayout(pLay, [3 4]);
    glRows.RowHeight = repmat({28}, 1, 3);
    % !! 第 1 列是拖动手柄；拖动排序只对"中间层"开放（顶层与衬底固定）
    glRows.ColumnWidth = {20, 24, '1x', 26};
    glRows.Padding = [4 4 4 4];
    glRows.RowSpacing = 3;

    btnAdd = uibutton(g1, 'Text', '＋ 添加一层（插在顶层下方）', ...
                      'ButtonPushedFcn', @(s, ~) on_add());

    gp = uigridlayout(g1, [5 2]);
    gp.RowHeight = repmat({26}, 1, 5);
    gp.ColumnWidth = {'fit', '1x'};
    gp.ColumnSpacing = 6;
    gp.RowSpacing = 4;
    uilabel(gp, 'Text', '慢轴长 [μm]');
    efLx = uieditfield(gp, 'numeric', 'Value', S.Lx, 'Limits', [1 1e5], ...
                       'ValueChangedFcn', @(s, ~) on_geo('Lx', s.Value));
    uilabel(gp, 'Text', '快轴长 [μm]');
    efLy = uieditfield(gp, 'numeric', 'Value', S.Ly, 'Limits', [0.5 1e4], ...
                       'ValueChangedFcn', @(s, ~) on_geo('Ly', s.Value));
    uilabel(gp, 'Text', '纵轴长 [μm]');
    efLz = uieditfield(gp, 'numeric', 'Value', S.Lz, 'Limits', [10 1e6], ...
                       'ValueChangedFcn', @(s, ~) on_geo('Lz', s.Value));
    labRW = uilabel(gp, 'Text', '脊宽 [μm]');
    efRW = uieditfield(gp, 'numeric', 'Value', S.wr, 'Limits', [0.2 1e4], ...
                       'ValueChangedFcn', @(s, ~) on_geo('wr', s.Value));
    labRH = uilabel(gp, 'Text', '脊高 [μm]');
    efRH = uieditfield(gp, 'numeric', 'Value', S.he, 'Limits', [0.001 1e3], ...
                       'ValueChangedFcn', @(s, ~) on_geo('he', s.Value));

    lblSummary = uilabel(g1, 'Text', '', 'FontSize', 10, 'FontColor', [0.35 0.35 0.40], ...
                         'VerticalAlignment', 'center');

    lblMsg = uilabel(g1, 'Text', '', 'FontSize', 11, 'FontWeight', 'bold', ...
                     'FontColor', [0.75 0.10 0.10], 'WordWrap', 'on', ...
                     'VerticalAlignment', 'top');

    %% ===== 第二列：层参数 =====
    p2 = uipanel(mainGrid, 'Title', '② 层参数（在①里点一层来选中）', 'FontWeight', 'bold', 'FontSize', 12);
    g2 = uigridlayout(p2, [3 1]);
    g2.RowHeight = {42, '1x', 56};
    g2.Padding = [8 6 8 6];

    g2top = uigridlayout(g2, [1 2]);
    g2top.ColumnWidth = {'1x', 104};
    g2top.Padding = [0 0 0 0];
    lblSel = uilabel(g2top, 'Text', '', 'FontSize', 12, 'FontWeight', 'bold');
    btnAdv = uibutton(g2top, 'Text', '高级设置…', 'FontSize', 10, ...
                      'ButtonPushedFcn', @(s, ~) open_advanced());

    gp2 = uigridlayout(g2, [7 2]);
    % !! ②列原来下半截全空：网格 6 行固定 28 px 撑不满面板，看着像"错位"。
    % !! 补第 7 行（1x 自适应）放「本层上下文」，空白就被有用的信息填掉了。
    gp2.RowHeight = [repmat({28}, 1, 6), {'1x'}];
    gp2.ColumnWidth = {'fit', '1x'};
    gp2.ColumnSpacing = 6;
    gp2.RowSpacing = 8;
    uilabel(gp2, 'Text', 'tag（层名）');
    efTag = uieditfield(gp2, 'text', 'Value', '', 'ValueChangedFcn', @(s, ~) on_lay('tag', s.Value));
    uilabel(gp2, 'Text', '折射率 n');
    efN = uieditfield(gp2, 'numeric', 'Value', 1, 'Limits', [0.5 6], ...
                      'ValueChangedFcn', @(s, ~) on_lay('n', s.Value));
    uilabel(gp2, 'Text', '损耗 α [cm⁻¹]');
    efA = uieditfield(gp2, 'numeric', 'Value', 0, 'Limits', [0 1e5], ...
                      'ValueChangedFcn', @(s, ~) on_lay('alpha', s.Value));
    uilabel(gp2, 'Text', '厚度 d [μm]');
    efD = uieditfield(gp2, 'numeric', 'Value', 0.1, 'Limits', [1e-4 1e4], ...
                      'ValueChangedFcn', @(s, ~) on_lay('d', s.Value));
    uilabel(gp2, 'Text', '有源层（QW）');
    cbAct = uicheckbox(gp2, 'Text', '这一层是发光层（Γ 只对有源层积分）', 'FontSize', 11, ...
                       'ValueChangedFcn', @(s, ~) on_lay('active', s.Value));
    uilabel(gp2, 'Text', '');
    lblDnote = uilabel(gp2, 'Text', '', 'FontSize', 10, 'FontColor', [0.40 0.40 0.45], 'WordWrap', 'on');
    lblCtx = uilabel(gp2, 'Text', '', 'FontSize', 10, 'FontColor', [0.22 0.22 0.28], ...
                     'WordWrap', 'on', 'VerticalAlignment', 'top');
    lblCtx.Layout.Row = 7;
    lblCtx.Layout.Column = [1 2];

    lblHint2 = uilabel(g2, 'Text', ['提示：衬底与顶层的厚度固定为无限大（三维图里衬底按 10% 高度示意）。' ...
                                    '顶层厚度就是本阶段允许刻的最大深度，不允许刻穿。'], ...
                       'FontSize', 10, 'FontColor', [0.40 0.40 0.45], 'WordWrap', 'on');

    %% ===== 第三列：三维预览 =====
    p3 = uipanel(mainGrid, 'Title', '③ 波导截面示意（拖动旋转，＋/－ 缩放截面）', ...
                 'FontWeight', 'bold', 'FontSize', 12);
    g3 = uigridlayout(p3, [3 1]);
    g3.RowHeight = {'1x', 30, 34};
    g3.Padding = [6 6 6 6];
    g3.RowSpacing = 6;
    ax3d = uiaxes(g3);
    ax3d.Box = 'on';
    ax3d.XGrid = 'on';  ax3d.YGrid = 'on';  ax3d.ZGrid = 'on';
    ax3d.FontSize = 10;
    ax3d.Color = [1 1 1];
    try
        % !! 不用内置的 rotate/zoomInteraction（在用户机器上转不动、滚轮还会把相机推进盒子）；
        % !! 旋转与缩放全部由图窗级的 on_down/on_motion/on_wheel 自己实现。
        ax3d.Interactions = [];
    catch
    end

    zrow = uigridlayout(g3, [1 5]);
    zrow.ColumnWidth = {'1x', 52, 84, 52, '1x'};
    zrow.Padding = [0 0 0 0];
    uibutton(zrow, 'Text', '', 'Enable', 'off');
    uibutton(zrow, 'Text', '－', 'FontSize', 13, 'Tooltip', '缩小（放大看细节用 ＋）', ...
             'ButtonPushedFcn', @(s, ~) on_zoom(1.25));
    lblZoom = uilabel(zrow, 'Text', '100%', 'HorizontalAlignment', 'center', 'FontSize', 11);
    uibutton(zrow, 'Text', '＋', 'FontSize', 13, 'Tooltip', '放大截面', ...
             'ButtonPushedFcn', @(s, ~) on_zoom(1/1.25));
    uibutton(zrow, 'Text', '复位', 'FontSize', 10, 'Tooltip', '回到 100% 与初始朝向', ...
             'ButtonPushedFcn', @(s, ~) on_view_reset());

    lbl3d = uilabel(g3, 'Text', '', 'FontSize', 10, 'FontColor', [0.40 0.40 0.45], ...
                    'WordWrap', 'on');

    %% ===== 底部按钮 =====
    bot = uigridlayout(root, [1 5]);
    bot.Layout.Row = 2;
    bot.ColumnWidth = {'1x', 190, 170, 120, '1x'};
    bot.Padding = [0 0 0 0];
    btnSure = uibutton(bot, 'Text', 'Sure（确认结构，进入成果画布）', 'FontSize', 12, ...
                       'FontWeight', 'bold', 'BackgroundColor', [0.80 0.89 0.98], ...
                       'ButtonPushedFcn', @(s, ~) on_sure());
    btnSure.Layout.Column = 2;
    btnClear = uibutton(bot, 'Text', 'Clear（清空并复位）', 'FontSize', 12, ...
                        'ButtonPushedFcn', @(s, ~) on_clear());
    btnClear.Layout.Column = 3;
    btnExit = uibutton(bot, 'Text', 'Exit（退出）', 'FontSize', 12, ...
                       'ButtonPushedFcn', @(s, ~) delete(hfig));
    btnExit.Layout.Column = 4;

    %% ===== 成果画布 =====
    repGrid = uigridlayout(root, [2 1]);
    repGrid.Layout.Row = 1;
    repGrid.RowHeight = {46, '1x'};
    repGrid.Padding = [0 0 0 0];
    repGrid.Visible = 'off';

    gHead = uigridlayout(repGrid, [1 7]);
    gHead.ColumnWidth = {'1x', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
    lblRepHead = uilabel(gHead, 'Text', '', 'FontSize', 12, 'FontWeight', 'bold');
    lblRepHead.Layout.Column = 1;
    % !! 用户 2026-10-29 的要求：告诉他一共有几个横向模，显示几阶由他定（上限 7）；
    % !! 损耗分解那张要能选模式，所以再给一个「模式」下拉。
    lblNS = uilabel(gHead, 'Text', '显示阶数', 'FontSize', 11, ...
                    'HorizontalAlignment', 'right', 'Tooltip', '横向模最多显示到第几阶（1-7）');
    lblNS.Layout.Column = 2;
    spnN = uispinner(gHead, 'Limits', [1 7], 'Value', 3, 'Step', 1, 'FontSize', 11, ...
                     'Tooltip', '横向模最多显示到第几阶（1-7）', ...
                     'ValueChangedFcn', @(s, ~) on_nshow(s.Value));
    spnN.Layout.Column = 3;
    lblMS = uilabel(gHead, 'Text', '模式', 'FontSize', 11, ...
                    'HorizontalAlignment', 'right', 'Tooltip', '第 5、11 项按这一阶出图');
    lblMS.Layout.Column = 4;
    ddMode = uidropdown(gHead, 'Items', {'m=0'}, 'Value', 'm=0', 'FontSize', 11, ...
                        'Tooltip', '第 5、11 项按这一阶出图', ...
                        'ValueChangedFcn', @(s, ~) on_msel(s.Value));
    ddMode.Layout.Column = 5;
    btnRecomp.Layout.Column = 6;
    btnBack.Layout.Column = 7;
    btnRecomp = uibutton(gHead, 'Text', '按当前结构重算', 'FontSize', 11, ...
                         'ButtonPushedFcn', @(s, ~) force_recompute());
    btnBack = uibutton(gHead, 'Text', '← 返回结构设置', 'FontSize', 11, ...
                       'ButtonPushedFcn', @(s, ~) on_back());

    gRep = uigridlayout(repGrid, [1 2]);
    gRep.ColumnWidth = {350, '1x'};

    pList = uipanel(gRep, 'Title', '成果按钮（一次只出一张，全部按当前结构现算）', 'FontSize', 11);
    pList.Scrollable = 'on';
    glRep = uigridlayout(pList, [1 1]);
    glRep.RowHeight = {26};
    glRep.Padding = [6 6 6 6];

    pView = uipanel(gRep, 'Title', '图像', 'FontSize', 11);
    gView = uigridlayout(pView, [2 1]);
    gView.RowHeight = {'1x', 30};
    axRep = uiaxes(gView);
    axRep.XTick = [];  axRep.YTick = [];
    axRep.Color = [1 1 1];
    lblRepMsg = uilabel(gView, 'Text', '', 'FontSize', 11, 'FontColor', [0.35 0.35 0.40], ...
                        'WordWrap', 'on');

    hRep = gobjects(0, 1);
    kCur = 1;   % !! 当前在看的第几项（切换显示阶数/模式时用它重画）
    build_report_list();
    refresh();

    h = struct('fig', hfig, 'handles', struct('ddKind', ddKind, 'efLx', efLx, 'efLy', efLy, ...
             'efLz', efLz, 'efRW', efRW, 'efRH', efRH, 'efTag', efTag, 'efN', efN, ...
             'efA', efA, 'efD', efD, 'cbAct', cbAct, 'btnSure', btnSure, 'btnAdd', btnAdd, ...
             'btnAdv', btnAdv, 'lblMsg', lblMsg, 'ax3d', ax3d, 'axRep', axRep, ...
             'lblRepMsg', lblRepMsg), ...
             'refresh', @refresh, 'state', @() S, 'report', @(k) show_item(k), ...
             'compute', @() compute_now());

    %% ================= 嵌套：界面逻辑 =================

    function on_kind(v)
        if strcmp(v, '脊形波导'), S.kind = 'ridge'; else, S.kind = 'rect'; end
        invalidate();
        refresh();
    end

    function on_geo(f, v)
        S.(f) = v;
        invalidate();
        refresh();
    end

    function on_add()
        k = numel(S.layers);
        newL = struct('tag', sprintf('layer%d', k - 1), 'n', 2.400, 'd', 0.050, ...
                      'alpha', 0, 'active', false);
        S.layers = [S.layers(1:k-1), newL, S.layers(k)];
        S.sel = k;
        invalidate();
        refresh();
    end

    function on_del(i)
        if i == 1 || i == numel(S.layers) || numel(S.layers) <= 3
            return
        end
        S.layers(i) = [];
        S.sel = min(S.sel, numel(S.layers));
        invalidate();
        refresh();
    end

    function on_select(i)
        S.sel = i;
        refresh();
    end

    function on_lay(f, v)
        i = S.sel;
        if i < 1 || i > numel(S.layers), return, end
        switch f
            case 'tag',    S.layers(i).tag = char(v);
            case 'n',      S.layers(i).n = v;
            case 'alpha',  S.layers(i).alpha = v;
            case 'd',      S.layers(i).d = v;
            case 'active'
                S.layers(i).active = logical(v);
                if v
                    for j = 1:numel(S.layers)
                        if j ~= i, S.layers(j).active = false; end
                    end
                end
        end
        invalidate();
        refresh();
    end

    function invalidate()
        Ucur = [];
    end

    function on_sure()
        [ok, ~] = check_state();
        if ~ok, return, end
        S.page = 'report';
        mainGrid.Visible = 'off';
        repGrid.Visible = 'on';
        lblRepHead.Text = sprintf('成果画布 · %s / %d 层 / 脊宽 %.2f μm / 脊高 %.3f μm / 顶部 n=%.2f', ...
                                  kind_name(), numel(S.layers), S.wr, S.he, S.n_top);
        lblRepMsg.Text = '正在按当前结构计算…';
        drawnow;
        if compute_now()
            show_item(1);
        end
    end

    function ok = compute_now()
        lblRepMsg.Text = '正在按当前结构计算（竖直 → 横向 → 二维场 → 远场）…';
        drawnow;
        try
            U = lui_compute(S);
        catch err
            U = struct('ok', false, 'msg', ['计算异常：' err.message]);
        end
        Ucur = U;
        ok = U.ok;
        if ~ok
            cla(axRep);
            text(axRep, 0.5, 0.5, {'无法求解', '', U.msg}, 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', [0.85 0.35 0.30], ...
                 'Interpreter', 'none');
            axRep.XTick = [];  axRep.YTick = [];
            lblRepMsg.Text = ['按当前结构无解：' strrep(U.msg, newline, ' ')];
        else
            % !! 算完立刻把「最多几阶」告诉用户，并把控件范围对齐到实际阶数
            nmax = max(1, min(7, numel(U.lat.N)));
            ddMode.Items = arrayfun(@(q) sprintf('m=%d', q - 1), 1:nmax, 'UniformOutput', false);
            if ~any(strcmp(ddMode.Items, ddMode.Value))
                ddMode.Value = ddMode.Items{1};
            end
            S.opt.msel = max(0, min(S.opt.msel, nmax - 1));
            spnN.Limits = [1 max(2, nmax)];
            S.opt.nshow = max(1, min(S.opt.nshow, nmax));
            spnN.Value = S.opt.nshow;
            lblRepMsg.Text = sprintf(['已按当前结构算完（%.1f s）。横向共 %d 阶（理论上限 %d 阶），' ...
                                      '当前显示前 %d 阶。左侧点按钮出图。'], ...
                                     U.t_sec, nmax, U.lat.M_expect, S.opt.nshow);
        end
    end

    % !! 切换「显示阶数」：只重画当前这一张，不重算结构
    function on_nshow(v)
        S.opt.nshow = max(1, min(7, round(v)));
        if ~isempty(Ucur) && Ucur.ok, show_item(kCur); end
    end

    % !! 切换「模式」：第 5 项（二维模场）与第 11 项（损耗分解）按这一阶出图
    function on_msel(txt)
        v = sscanf(txt, 'm=%d');
        if isempty(v), v = 0; end
        S.opt.msel = max(0, v);
        if ~isempty(Ucur) && Ucur.ok, show_item(kCur); end
    end

    function force_recompute()
        invalidate();
        if compute_now(), show_item(1); end
    end

    function on_back()
        S.page = 'main';
        repGrid.Visible = 'off';
        mainGrid.Visible = 'on';
    end

    function on_clear()
        % !! 复位结构时，运行时字段（拖动排序、三维相机）必须一起补回来：少了它们，
        % !! 图窗的鼠标回调每次移动都会报「无法识别的字段名称 "camDrag"」（用户 2026-10-29 报过）。
        S = lui_live(lui_default());
        Ucur = [];
        ddKind.Value = '脊形波导';
        efLx.Value = S.Lx;  efLy.Value = S.Ly;  efLz.Value = S.Lz;
        efRW.Value = S.wr;  efRH.Value = S.he;
        on_back();
        refresh();
    end

    function nm = kind_name()
        if strcmp(S.kind, 'ridge'), nm = '脊形波导'; else, nm = '矩形波导'; end
    end

    %% ---------- 校验 ----------
    function [ok, msg] = check_state()
        bad = {};
        if ~(isfinite(S.Lx) && S.Lx > 0), bad{end+1} = '慢轴长必须为正的有限值'; end
        if ~(isfinite(S.Ly) && S.Ly > 0), bad{end+1} = '快轴长必须为正的有限值'; end
        if ~(isfinite(S.Lz) && S.Lz > 0), bad{end+1} = '纵轴长必须为正的有限值'; end
        if numel(S.layers) < 3, bad{end+1} = '至少要 3 层（衬底 + 芯/包层 + 顶层）'; end
        if ~isinf(S.layers(1).d), bad{end+1} = '第 1 层（衬底）的厚度必须是无限大'; end
        if ~(isfinite(S.layers(end).d) && S.layers(end).d > 0)
            bad{end+1} = '顶层厚度必须是正的有限值（脊高不能超过它）';
        end
        Tm = 0;
        for j = 2:numel(S.layers)-1
            d = S.layers(j).d;
            if ~(isfinite(d) && d > 0)
                bad{end+1} = sprintf('第 %d 层（%s）厚度必须为正的有限值', j, S.layers(j).tag);
            end
            if ~(isfinite(S.layers(j).n) && S.layers(j).n >= 1 && S.layers(j).n <= 5)
                bad{end+1} = sprintf('第 %d 层（%s）折射率 n=%.3f 不合理（应在 1–5）', j, S.layers(j).tag, S.layers(j).n);
            end
            if ~(S.layers(j).alpha >= 0)
                bad{end+1} = sprintf('第 %d 层（%s）损耗必须 ≥ 0', j, S.layers(j).tag);
            end
            Tm = Tm + d;
        end
        Tg = Tm + S.layers(end).d;
        if Tm > S.Ly
            bad{end+1} = sprintf('中间各层厚度之和 %.4f μm 超过快轴长 %.4f μm', Tm, S.Ly);
        end
        if Tg < 0.05*S.Ly
            bad{end+1} = sprintf('结构总厚 %.4f μm 相对快轴长 %.4f μm 太小，三维图看不清', Tg, S.Ly);
        end
        tg = {S.layers.tag};
        if numel(unique(tg)) ~= numel(tg)
            bad{end+1} = '层名（tag）不能重复，否则损耗表会串层';
        end
        if sum([S.layers.active]) > 1
            bad{end+1} = '有源层最多一层';
        end
        if strcmp(S.kind, 'ridge')
            if ~(S.wr > 0 && S.wr <= S.Lx)
                bad{end+1} = sprintf('脊宽 %.4f μm 必须在 (0, 慢轴长 %.4f μm] 内', S.wr, S.Lx);
            end
            dtop = S.layers(end).d;
            if ~(S.he > 0)
                bad{end+1} = '脊高必须为正';
            elseif S.he > dtop + 1e-12
                bad{end+1} = sprintf('脊高 %.4f μm 超过顶层厚度 %.4f μm（本阶段不允许刻穿）', S.he, dtop);
            end
            if S.he > S.Ly
                bad{end+1} = sprintf('脊高 %.4f μm 超过快轴长 %.4f μm', S.he, S.Ly);
            end
        end
        if ~(isfinite(S.n_top) && S.n_top >= 1 && S.n_top <= 4)
            bad{end+1} = '顶部介质折射率应在 1–4';
        end
        if ~(S.opt.lambda > 0.2 && S.opt.lambda < 1.5)
            bad{end+1} = '波长应在 0.2–1.5 μm';
        end
        ok = isempty(bad);
        if ok
            msg = '';
        else
            msg = strjoin(bad, '；');
        end
    end

    %% ---------- 刷新 ----------
    function refresh()
        [ok, msg] = check_state();
        vis = 'on';
        if strcmp(S.kind, 'rect'), vis = 'off'; end
        labRW.Visible = vis;  efRW.Visible = vis;
        labRH.Visible = vis;  efRH.Visible = vis;

        if ok
            btnSure.Enable = 'on';
            lblMsg.Text = '结构合法，可以点 Sure。';
            lblMsg.FontColor = [0.10 0.50 0.15];
        else
            btnSure.Enable = 'off';
            lblMsg.Text = msg;
            lblMsg.FontColor = [0.75 0.10 0.10];
        end

        rebuild_rows();

        i = min(max(1, S.sel), numel(S.layers));
        S.sel = i;
        L = S.layers(i);
        lblSel.Text = sprintf('第 %d 层：%s', i, L.tag);
        efTag.Value = L.tag;  efN.Value = L.n;  efA.Value = L.alpha;
        if isinf(L.d)
            efD.Value = 0.001;  efD.Enable = 'off';
            lblDnote.Text = '这一层是半无限层，厚度不可改。';
        else
            efD.Value = L.d;    efD.Enable = 'on';
            lblDnote.Text = sprintf('中间层厚度之和 %.4f / 快轴长 %.4f μm。', sum_mid(), S.Ly);
        end
        cbAct.Value = L.active;
        if i == 1
            cbAct.Enable = 'off';
        else
            cbAct.Enable = 'on';
        end

        % !! 第 7 行「本层上下文」：层序 + 上下折射率 + 这层在竖直方向上的角色。
        % !! 用户看三维图时最常问的就是"哪层在限制光"，这里直接给判断依据（n 与两侧比较）。
        N = numel(S.layers);
        if i > 1, sUp = sprintf('上方 n = %.4f', S.layers(i-1).n); else, sUp = '上方：无（最外层）'; end
        if i < N, sDn = sprintf('下方 n = %.4f', S.layers(i+1).n); else, sDn = '下方：无（衬底）'; end
        sRole = '';
        if i > 1 && i < N
            nb = [S.layers(i-1).n, S.layers(i+1).n];
            if L.n > max(nb) + 1e-9
                sRole = '本层 n 高于上下两侧 → 竖直方向形成折射率峰，光会被限制在这一层附近。';
            elseif L.n < min(nb) - 1e-9
                sRole = '本层 n 低于上下两侧 → 包层角色，竖直方向不限制光。';
            else
                sRole = '本层 n 夹在上下两侧之间 → 不构成折射率峰（是台阶的一侧）。';
            end
        elseif i == 1
            sRole = '最外层为覆盖介质（顶部 n_top），位于整个结构之上。';
        else
            sRole = '最外层为衬底，位于整个结构之下（半无限）。';
        end
        lblCtx.Text = sprintf(['第 %d 层 / 共 %d 层（自上而下）\n%s　　本层 n = %.4f　　%s\n%s'], ...
                              i, N, sUp, L.n, sDn, sRole);

        lblSummary.Text = sprintf('%d 层 · 中间层 %.3f / 总厚 %.3f / 顶层可刻 %.3f μm', ...
                                  numel(S.layers), sum_mid(), sum_mid() + S.layers(end).d, S.layers(end).d);
        draw3d();
    end

    function t = sum_mid()
        t = 0;
        for j = 2:numel(S.layers)-1, t = t + S.layers(j).d; end
    end

    function rebuild_rows()
        delete(glRows.Children);
        N = numel(S.layers);
        glRows.RowHeight = repmat({28}, 1, N);
        rowBtn = gobjects(1, N);
        for r = 1:N
            i = N - r + 1;
            L = S.layers(i);
            % !! 拖动手柄：只有中间层显示 ≡（顶层与衬底不能拖，留空占位对齐）
            if i > 1 && i < N
                gp = uilabel(glRows, 'Text', '≡', 'FontSize', 13, ...
                             'HorizontalAlignment', 'center', ...
                             'Tooltip', '按住这一行上下拖动，可调整该层的位置');
            else
                gp = uilabel(glRows, 'Text', '');
            end
            gp.Layout.Row = r;  gp.Layout.Column = 1;
            sw = uilabel(glRows, 'Text', '', 'BackgroundColor', lui_color(i), ...
                         'VerticalAlignment', 'center');
            sw.Layout.Row = r;  sw.Layout.Column = 2;
            if isinf(L.d)
                dtxt = '∞';
            else
                dtxt = sprintf('%.3f', L.d);
            end
            mk = '';
            fw = 'normal';
            if i == S.sel, mk = '▸ '; fw = 'bold'; end
            nm = sprintf('%s%s    d=%s μm   n=%.3f', mk, L.tag, dtxt, L.n);
            b = uibutton(glRows, 'Text', nm, 'HorizontalAlignment', 'left', ...
                         'FontSize', 11, 'FontWeight', fw, ...
                         'ButtonPushedFcn', @(s, ~) on_select(i));
            b.Layout.Row = r;  b.Layout.Column = 3;
            rowBtn(r) = b;
            if i == 1 || i == N || N <= 3
                % !! 不可删的层（表层 / 衬底 / 只剩 3 层时）原来画一个 Enable='off' 的空按钮，
                % !! 看着像坏控件；直接不显示。
                db = uibutton(glRows, 'Text', '');
                db.Visible = 'off';
            else
                db = uibutton(glRows, 'Text', '✕', 'FontSize', 10, ...
                              'Tooltip', sprintf('删除第 %d 层（%s）', i, L.tag), ...
                              'ButtonPushedFcn', @(s, ~) on_del(i));
            end
            db.Layout.Row = r;  db.Layout.Column = 4;
        end
    end

    %% ---------- 层表拖动排序（用户 2026-10-29 要求） ----------
    function p = fig_xy()
        % !! uifigure 的 CurrentPoint 给的是图内像素坐标；个别版本取不到（返回空），
        % !! 那种情况下拖动只能停用，并明确告诉用户改用「上移 / 下移」按钮，不做静默失败。
        p = hfig.CurrentPoint;
        if isempty(p) || numel(p) < 2 || any(~isfinite(p(1:2)))
            p = [NaN NaN];
            if S.dragOK
                S.dragOK = false;
                lblMsg.Text = '这个 MATLAB 版本取不到鼠标坐标，拖动排序不可用：请用「上移 / 下移」按钮。';
            end
        else
            p = p(1:2);
        end
    end

    function q = pos_in_fig(h)
        % 组件的 Position 是相对父容器的，逐级加上去才得到图窗坐标
        q = h.Position;
        p = h.Parent;
        while ~isempty(p) && p ~= hfig
            pp = get(p, 'Position');
            q(1) = q(1) + pp(1);  q(2) = q(2) + pp(2);
            p = p.Parent;
        end
    end

    function r = row_at(p)
        r = 0;
        for rr = 1:numel(rowBtn)
            if ~isvalid(rowBtn(rr)), continue, end
            q = pos_in_fig(rowBtn(rr));
            if p(1) >= q(1) && p(1) <= q(1) + q(3) && p(2) >= q(2) && p(2) <= q(2) + q(4)
                r = rr;  return
            end
        end
    end

    function tf = inrect(p, q)
        tf = p(1) >= q(1) && p(1) <= q(1) + q(3) && p(2) >= q(2) && p(2) <= q(2) + q(4);
    end

    function on_wheel(e)
        % !! 只在鼠标位于三维 axes 内时缩放（免得在别处滚一下整个图就变了）
        c = 0;
        try
            c = e.VerticalScrollCount;
        catch
        end
        if c == 0, return, end
        p = fig_xy();
        if any(isnan(p)), return, end
        if ~inrect(p, pos_in_fig(ax3d)), return, end
        on_zoom(1.18 ^ c);      % 向上滚（c<0）→ 缩小视窗范围 = 放大
    end

    function on_down()
        S.dragSrc = 0;  S.dragMoved = false;  S.camDrag = false;
        p = fig_xy();
        if any(isnan(p)), return, end
        r = row_at(p);
        if r > 0
            N = numel(S.layers);  i = N - r + 1;
            if i == 1 || i == N, return, end      % 顶层与衬底不能拖
            S.dragSrc = i;  S.dragP0 = p;
            rowBtn(r).BackgroundColor = [0.86 0.93 1.00];
            return
        end
        % !! 旋转自己实现：uifigure 里 ax3d.Interactions 的 rotate 在用户机器上不生效
        % !! （2026-10-29 反馈"不可以旋转"），所以在图窗级事件里按像素位移改相机姿态。
        if inrect(p, pos_in_fig(ax3d))
            S.camDrag = true;  S.camP0 = p;
            S.camR0 = S.R;
        end
    end

    function on_motion()
        if S.camDrag
            p = fig_xy();
            if any(isnan(p)), return, end
            d = p - S.camP0;
            % !! 球面轨道（trackball）：把鼠标位移拆成两个绕【屏幕轴】的转动，合成到 S.R 上。
            % !! 不用 view(az,el)——正对截面时 el=90 是万向锁，改方位角只会绕屏幕法线打转，
            % !! 用户反馈的"左右拉像是在原地转、看不到侧面"就是这个原因。
            a =  d(1) * 0.28 * pi/180;     % 横向：绕屏幕竖直轴（右拉→转出左侧面）
            b =  d(2) * 0.28 * pi/180;     % 纵向：上拉→模型跟着往上转（即从下方看）
            Ry = [cos(a) 0 -sin(a); 0 1 0; sin(a) 0 cos(a)];
            Rx = [1 0 0; 0 cos(b) -sin(b); 0 sin(b) cos(b)];
            S.R = S.camR0 * Rx * Ry;
            apply_cam();
            S.needview = false;
            drawnow limitrate
            return
        end
        if S.dragSrc == 0, return, end
        p = fig_xy();
        if any(isnan(p)), return, end
        if norm(p - S.dragP0) > 12, S.dragMoved = true; end
    end

    function on_up()
        if S.camDrag
            S.camDrag = false;
            return
        end
        if S.dragSrc == 0, return, end
        iSrc = S.dragSrc;  moved = S.dragMoved;
        S.dragSrc = 0;  S.dragMoved = false;
        p = fig_xy();
        r = 0;
        if ~any(isnan(p)), r = row_at(p); end
        if ~moved || r == 0
            refresh();                            % 只是点了一下：恢复高亮
            return
        end
        N = numel(S.layers);  iT = N - r + 1;
        iT = min(max(iT, 2), N - 1);
        if iT == iSrc
            refresh();  return
        end
        perm = 1:N;  x = perm(iSrc);
        perm(iSrc) = [];
        iT = min(iT, numel(perm) + 1);
        perm = [perm(1:iT-1), x, perm(iT:end)];
        S.layers = S.layers(perm);
        S.sel = iT;
        invalidate();
        refresh();
    end

    %% ---------- 三维 ----------
    function draw3d()
        cla(ax3d);
        N = numel(S.layers);
        Lx = S.Lx;  Lz = S.Lz;
        % !! 三维纵深取 0.25·Lx（方案 A「立方体示意」，用户 2026-10-29 选）：深度与慢轴同量级，
        % !! 才有立体感；它与真实 Lz = 1500 μm 无关，真实值写在 z 轴标签里。
        Zd = 0.25 * max(S.Lx, 1);
        y0 = 0.10;  y1 = 0.88;
        dtop = S.layers(end).d;
        % !! 竖直方向按 √(厚度) 压缩。之前按线性厚度画：9 nm 的 QW 只有 4 px 宽的一条线，
        % !! 在三维图上等于"没有这一层"（用户看不出来 QW 在哪）。√ 压缩后 QW≈24 px、EBL≈36 px，
        % !! 层序与"谁薄"仍看得出来；精确厚度以①列的层表为准（三位小数）。
        wmid = sqrt(max([S.layers(2:N-1).d], 1e-9));
        wtop = sqrt(max(dtop, 1e-9));
        wsum = sum(wmid) + wtop;
        band = y1 - y0;
        hmid = band * wmid / wsum;
        htop = band * wtop / wsum;


        lui_box(ax3d, [-Lx/2 Lx/2], [0 y0], [0 Zd], [0.80 0.80 0.80], 1);
        text(ax3d, -Lx/2*0.92, y0/2, Zd*1.02, S.layers(1).tag, 'FontSize', 9, ...
             'Color', [0.25 0.25 0.25], 'Interpreter', 'none');

        y = y0;
        % !! 先只画盒子、把层名记下来，等所有层的位置都定了再统一排标签（见下面的"标尺"）
        lab  = struct('tag', {}, 'act', {}, 'yc', {});
        for j = 2:N-1
            hh = hmid(j-1);
            lui_box(ax3d, [-Lx/2 Lx/2], [y y+hh], [0 Zd], lui_color(j), 1);
            lab(end+1) = struct('tag', S.layers(j).tag, 'act', S.layers(j).active, 'yc', y + hh/2);
            y = y + hh;
        end

        % !! 视角抬高到 −28° 后，薄层在屏幕上的高度被压得很小，层名会互相压住（p-wave/QW/n-wave
        % !! 挤成三行）。改成"标尺"排法：标签按最小间距 sep 自上而下排开，各自用一条细引线
        % !! 指回本层的中心 —— 引线短、且不会指错层，视觉上比就地写更清楚。
        if ~isempty(lab)
            [~, ord] = sort([lab.yc]);
            sep = 0.050;
            ylq = [lab(ord).yc];
            for q = 2:numel(ylq)
                if ylq(q) - ylq(q-1) < sep
                    ylq(q) = ylq(q-1) + sep;
                end
            end
            ov = ylq(end) - (y1 - 0.012);
            if ov > 0, ylq = ylq - ov; end
            un = (y0 + 0.012) - ylq(1);
            if un > 0, ylq = ylq + un; end
            for q = 1:numel(ord)
                jq = ord(q);
                if lab(jq).act
                    tc = [1 1 1];  tw = 'bold';
                else
                    tc = [0.15 0.15 0.15];  tw = 'normal';
                end
                xl0 = -Lx/2*0.92;
                % !! 必须用 line() 而不是 plot3()：plot3 是高层绘图命令，hold 为 off 时会把轴上的
                % !! 既有对象全部清掉 —— 2026-10-29 那次「下半截不见了」就是它：7 个层盒子画完后，
                % !! 第一次 plot3 引线把衬底与 5 个中间层全擦了，只剩后画的顶层 4 个 patch。
                line(ax3d, [xl0 + 0.035*Lx, xl0 + 0.20*Lx], [ylq(q), ylq(q)], ...
                     [Zd*1.02, Zd*1.02], 'Color', [0.35 0.35 0.35], 'LineWidth', 0.4);
                text(ax3d, xl0, ylq(q), Zd*1.02, lab(jq).tag, ...
                     'FontSize', 9, 'Color', tc, 'FontWeight', tw, 'Interpreter', 'none');
            end
        end

        % !! 顶层（p-clad）的脊：先把这一层按 sqrt 权重拆成"脊下整层"与"脊身"两段，
        % !! 脊身的视觉高度与沟槽里被刻成顶部介质的那段同高 —— 这样 h_e = 顶层厚度时脊面
        % !! 与槽面严格齐平（物理上正是"刻到 EBL"），h_e 调小才凸出来。
        hef = min(S.he, dtop);
        wu = sqrt(max(hef, 1e-9));
        wd = sqrt(max(dtop - hef, 1e-9));
        hup = htop * wu / (wu + wd);
        if strcmp(S.kind, 'ridge')
            hbase = max(htop - hup, 1e-4);
            lui_box(ax3d, [-Lx/2 Lx/2], [y y+hbase], [0 Zd], lui_color(N), 1);
            lui_box(ax3d, [-S.wr/2 S.wr/2], [y+hbase y+htop], [0 Zd], lui_color(N), 1);
            if Lx > S.wr
                lui_box(ax3d, [-Lx/2 -S.wr/2], [y+hbase y+htop], [0 Zd], [0.70 0.85 0.70], 0.28);
                lui_box(ax3d, [S.wr/2 Lx/2], [y+hbase y+htop], [0 Zd], [0.70 0.85 0.70], 0.28);
            end
            text(ax3d, 0, y + hbase + hup/2, Zd*1.02, sprintf('脊 %s', S.layers(N).tag), ...
                 'FontSize', 9, 'Color', [0.15 0.15 0.15], 'Interpreter', 'none');
            if hef >= dtop - 1e-9
                % !! 默认 h_e = 0.300 μm 正好等于顶层 p-clad 的厚度，沟槽把 p-clad 刻光，
                % !! 于是脊面与槽面同高、看不见"凸起"。这不是画错，必须写在图注里。
                lbl3d.Text = sprintf(['脊宽 %.2f μm（慢轴 %.1f μm 的 %.1f%%）· 脊高 %.3f μm · 顶层可刻 %.3f μm · ' ...
                                      '竖轴按 √(厚度) 压缩；当前脊高 = 顶层厚度 → 脊面与槽面齐平（槽底已到 EBL）'], ...
                                     S.wr, Lx, 100*S.wr/Lx, S.he, dtop);
            else
                lbl3d.Text = sprintf(['脊宽 %.2f μm（慢轴 %.1f μm 的 %.1f%%）· 脊高 %.3f μm · 顶层可刻 %.3f μm · ' ...
                                      '竖轴按 √(厚度) 压缩'], S.wr, Lx, 100*S.wr/Lx, S.he, dtop);
            end
        else
            lui_box(ax3d, [-Lx/2 Lx/2], [y y+htop], [0 Zd], lui_color(N), 1);
            text(ax3d, -Lx/2*0.92, y + htop/2, Zd*1.02, S.layers(N).tag, 'FontSize', 9, ...
                 'Color', [0.15 0.15 0.15], 'Interpreter', 'none');
            lbl3d.Text = sprintf('矩形波导：%d 层平铺，无侧向限制（慢轴 %.1f μm 仅供展示）· 竖轴按 √(厚度) 压缩', N, Lx);
        end

        % !! 初始等轴测（-37.5°, 22°）：层序与脊块一眼能看出来；想看正对截面自己拖过去。
        xl = Lx / 2 * S.zoom;
        yc = 0.49;
        yh = 0.5 * S.zoom;
        % !! 盒子四周留 9% 余量：等轴测下几何体的右上角会贴到框边，看着被切掉（用户 2026-10-29）
        padv = 0.09;
        xlim(ax3d, [-xl xl] + [-1 1]*padv*2*xl);
        ylim(ax3d, [yc - yh, yc + yh] + [-1 1]*padv*2*yh);
        zlim(ax3d, [0 Zd] + [-1 1]*padv*Zd);
        % !! 竖直方向的"厚度"是归一化过的（真实 1.7 μm 对 10 μm 宽会压成一张纸），
        % !! 所以用 pbaspect 把 y 撑到与 x 同高，看起来才是块而不是片（不改数据，只改显示比例）
        try
            pbaspect(ax3d, [1 0.95 0.5]);
        catch
        end
        % !! x 轴的数字默认画在远端顶棱上（三维轴总把刻度放在离相机最远的那条棱），
        % !! 用户 2026-10-29 要求放到下面。做法是只藏数字、保留刻度线与网格，
        % !! 再把数字自己画到底面近棱外侧：
        % !!   · 不能用 `XTick=[]`，那会把 XLimMode 悄悄改回 auto，紧接着的 plot3 按数据
        % !!     重算范围，整张图缩成 [-1,1]（实测过）；
        % !!   · `XTickLabel={}` 只影响标签，不动范围，是安全的。
        ax3d.XTickLabel = {};
        delete(findall(ax3d, 'Tag', 'uiaxlab'));
        yloT = yc - yh - padv*2*yh;
        yhiT = yc + yh + padv*2*yh;
        dyT = yhiT - yloT;
        for v = linspace(-xl, xl, 5)
            text(ax3d, v, yloT + 0.012*dyT, Zd*1.02, sprintf('%g', v), 'FontSize', 9, ...
                 'Color', [0.15 0.15 0.15], 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'top', 'Tag', 'uiaxlab');
        end
        % !! xlabel 也跟着数字走：轴标题本身留在顶棱上会跟底下的数字脱节
        ax3d.XLabel.String = '';
        text(ax3d, 0, yloT - 0.052*dyT, Zd*1.02, '慢轴 x [\mum]', 'FontSize', 10, ...
             'Color', [0.15 0.15 0.15], 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'top', 'Tag', 'uiaxlab');
        ylabel(ax3d, '快轴 y（归一化）');
        zlabel(ax3d, sprintf('纵轴 z（示意；实际 L_z = %.0f μm）', Lz));
        % !! 只在首次绘制和按「复位」时设一次视角：draw3d 每次重画都动相机会把用户刚转到
        % !! 的角度顶回去（上一版"根本转不动"就是这个）；之后交给交互，重画不动相机。
        if ~isfield(S, 'needview') || S.needview
            % !! 默认视角必须是「快轴朝上、层序竖着摞」（用户 2026-10-29：view(az,el) 无论怎么给，
            % !! z 都会被摆成竖直的，层序看着是横的——因为 view 的 up 永远是 +z）。
            % !! 这里改成 up = +y（快轴）：先 view(0,90) 把 up 与距离定下来，再用 R 把相机绕到
            % !! 右前上方，得到"能看到上表面"的立体感。R 也是拖动的起点。
            view(ax3d, 0, 90);
            S.R = lui_home_rot();  S.dc = [];
            S.needview = false;
            doCam = true;
        elseif max(abs(S.R(:) - reshape(eye(3), [], 1))) > 0
            doCam = true;     % 用户转过角度：重画后按 S.R 把相机摆回去
        else
            doCam = false;
        end
        lblZoom.Text = sprintf('%.0f%%', 100/S.zoom);
        % !! 不打灯光：默认 SpecularStrength 0.9 会把每层照成纯白（上一版「第三列彻底变白」就是它）
        try
            lighting(ax3d, 'none');
        catch
        end
        % !! 相机必须在所有 xlim/ylim/zlim 与 pbaspect 定稿之后再摆（先摆会被重算覆盖）。
        % !! 这里原来漏了这一步，只有拖动时才 apply_cam，所以默认视角一直停在 view(0,90)
        % !! = 整个面正对用户（用户 2026-10-29 的反馈）。
        if doCam
            apply_cam();
        end


    end

    %% ---------- 三维缩放 / 复位 ----------
    function on_zoom(f)
        % !! 缩放改的是 x–y 的坐标范围，不是相机距离：相机永远待在盒外，
        % !! 所以不会像滚轮那样"往里钻"进结构内部。
        S.zoom = min(max(S.zoom * f, 0.08), 1.0);
        draw3d();
    end

    function on_view_reset()
        S.zoom = 1.0;
        S.R = lui_home_rot();  S.dc = [];
        S.needview = true;   % !! 复位按钮显式要求回到默认视角（快轴朝上 + 轻微倾斜）
        draw3d();
    end

    % ─────────────────────────────────────────────────────
    % 函数：lui_home_rot（默认视角的旋转矩阵）
    % 输入：无
    % 输出：R（3×3，把相机方向 [0;0;1] 与 up [0;1;0] 一起转）
    % 说明：
    %   !! 绕 y 转 −25°、绕 x 倾 −28°：从右前上方看，能看到上表面（−12° 太平，用户要求再抬高），
    %   !! 但 y（快轴）仍然基本竖直 —— 层序是竖着摞的，不会看着横过来。
    % ─────────────────────────────────────────────────────
    function R = lui_home_rot()
        % !! 参数按「仰角 / 偏航」定义，不用旋转矩阵的符号：仰角为正 = 相机在水平面之上、
        % !! 从上往下看。当前值由用户 2026-10-29 亲手设定：仰角 3°、偏航 -20°
        % !! （3° 近似正对前脸，层带是完整长条、顶面只剩窄边）。
        % !! 内部 Rx 用负号：b 为负表示把模型朝观察者方向扳倒，等价于相机升高。
        elev = 3;   yaw = -20;
        a = yaw*pi/180;   b = -elev*pi/180;
        Ry = [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)];
        Rx = [1 0 0; 0 cos(b) -sin(b); 0 sin(b) cos(b)];
        R = Ry*Rx;
    end

    % ─────────────────────────────────────────────────────
    % 函数：apply_cam（按旋转矩阵 S.R 摆放三维相机，球面轨道）
    % 输入：无（读 S.R / S.dc / ax3d 当前状态）
    % 输出：无（直接 set campos / camup）
    % 说明：
    %   !! 为什么不用 view(az,el)：正对截面时 el = 90 是万向锁，那里改方位角
    %   !! 只会绕屏幕法线转，画面上就是"原地打转、看不到侧面"（用户 2026-10-29 反馈）。
    %   !! 这里把相机方向取成 R·[0;0;1]，up 取成 R·[0;1;0]——up 跟着一起转，
    %   !! 所以永远不会出现滚转（画面不会斜），横向拖动是真正的"绕竖直轴转出侧面"。
    %   !! 距离只在第一次（view(0,90) 之后）读一次存进 S.dc，缩放走坐标范围不走相机。
    % ─────────────────────────────────────────────────────
    function apply_cam()
        % !! Camera*Mode 一律留给 auto：钉成 manual 后 pbaspect 改盒形时相机不再重算，
        % !! 结构会被推到画面外（2026-10-29 实测整张空白）；CameraViewAngle 钉死更糟，
        % !! 会沿用上一次的小角度、缩成看不见的点。
        if isempty(S.dc) || ~isfinite(S.dc) || S.dc <= 0
            S.dc = norm(campos(ax3d) - camtarget(ax3d));
        end
        if ~isfinite(S.dc) || S.dc <= 0
            S.dc = 4 * max(S.Lx, 1);
        end
        t = camtarget(ax3d);
        v = (S.R * [0; 0; 1]).';
        campos(ax3d, t + v * S.dc);
        camup(ax3d, (S.R * [0; 1; 0]).');
        % !! 把相机钉住：uifigure 布局完成（或用户拉窗口）时，auto 模式会按 az/el 重算 campos，
        % !! 把刚设的方向与 up 覆盖回 view(0,90)。这几个 manual 必须在坐标范围与 pbaspect
        % !! 全部定稿之后才设（否则会重现坑 G6：内容被推出画面）。
        try
            ax3d.CameraPositionMode = 'manual';
            ax3d.CameraTargetMode   = 'manual';
            ax3d.CameraUpVectorMode = 'manual';
        catch
        end
    end

    %% ---------- 高级设置 ----------
    function open_advanced()
        f2 = uifigure('Name', '高级设置（opt 与顶部介质）', 'Position', [400 300 470 470]);
        try
            f2.WindowStyle = 'modal';
        catch
        end
        g = uigridlayout(f2, [11 2]);
        g.RowHeight = [repmat({30}, 1, 10), {44}];
        g.ColumnWidth = {'fit', '1x'};
        g.Padding = [12 12 12 12];
        g.RowSpacing = 8;

        uilabel(g, 'Text', '波长 λ [μm]');
        e1 = uieditfield(g, 'numeric', 'Value', S.opt.lambda, 'Limits', [0.2 1.5]);
        uilabel(g, 'Text', '偏振');
        e2 = uidropdown(g, 'Items', {'TE（本阶段只做准 TE）'}, 'Value', 'TE（本阶段只做准 TE）', 'Enable', 'off');
        uilabel(g, 'Text', '顶部介质 n（1.46=SiO₂）');
        e3 = uieditfield(g, 'numeric', 'Value', S.n_top, 'Limits', [1 4]);
        uilabel(g, 'Text', '根扫描点数 num_scan');
        e4 = uieditfield(g, 'numeric', 'Value', S.opt.num_scan, 'Limits', [1000 1e6], 'RoundFractionalValues', 'on');
        uilabel(g, 'Text', '收敛容差 tol');
        e5 = uieditfield(g, 'numeric', 'Value', S.opt.tol, 'Limits', [1e-16 1e-6]);
        uilabel(g, 'Text', '去重门限 merge');
        e6 = uieditfield(g, 'numeric', 'Value', S.opt.merge, 'Limits', [1e-12 1e-2]);
        uilabel(g, 'Text', '模场网格点数 npts');
        e7 = uieditfield(g, 'numeric', 'Value', S.opt.npts, 'Limits', [200 20000], 'RoundFractionalValues', 'on');
        uilabel(g, 'Text', '尾部衰减长度 tail_dec');
        e8 = uieditfield(g, 'numeric', 'Value', S.opt.tail_dec, 'Limits', [5 200], 'RoundFractionalValues', 'on');
        uilabel(g, 'Text', '远场采样数 num_ff');
        e9 = uieditfield(g, 'numeric', 'Value', S.opt.num_ff, 'Limits', [201 20001], 'RoundFractionalValues', 'on');

        gb = uigridlayout(f2, [1 2]);
        gb.Layout.Row = 11;
        gb.Layout.Column = [1 2];
        gb.ColumnWidth = {'1x', '1x'};
        uibutton(gb, 'Text', '应用', 'FontWeight', 'bold', ...
                 'ButtonPushedFcn', @(~, ~) apply_adv(e1, e3, e4, e5, e6, e7, e8, e9, f2));
        uibutton(gb, 'Text', '取消', 'ButtonPushedFcn', @(~, ~) delete(f2));
    end

    function apply_adv(e1, e3, e4, e5, e6, e7, e8, e9, f2)
        S.opt.lambda   = e1.Value;
        S.n_top        = e3.Value;
        S.opt.num_scan = round(e4.Value);
        S.opt.tol      = e5.Value;
        S.opt.merge    = e6.Value;
        S.opt.npts     = round(e7.Value);
        S.opt.tail_dec = round(e8.Value);
        S.opt.num_ff   = round(e9.Value);
        delete(f2);
        invalidate();
        refresh();
    end

    %% ---------- 成果画布 ----------
    function build_report_list()
        % !! 一项一行按钮（26 px + 4 px 间隔），说明改挂 Tooltip；原来"按钮 + 下面两行说明"
        % !! 的排法遇到换行就错位（用户 2026-10-28 反馈第 6 条）。选中项的说明显示在画布下方。
        items = lui_items();
        n = numel(items);
        glRep.RowHeight = repmat({28}, 1, n);
        glRep.RowSpacing = 5;
        hRep = gobjects(n, 1);
        for k = 1:n
            hRep(k) = uibutton(glRep, 'Text', sprintf('%d. %s', k, items(k).name), ...
                               'HorizontalAlignment', 'left', 'FontSize', 11, ...
                               'Tooltip', items(k).desc, ...
                               'ButtonPushedFcn', @(s, ~) show_item(k));
            hRep(k).Layout.Row = k;
            hRep(k).Layout.Column = 1;
        end
    end

    function show_item(k)
        items = lui_items();
        if k < 1 || k > numel(items), return, end
        kCur = k;
        if isempty(Ucur)
            if ~compute_now(), return, end
        end
        if ~Ucur.ok
            return
        end
        % !! 第 6、10 项要算好几秒，先清屏写"正在计算"再挂一个不确定进度条，
        % !! 否则用户以为程序卡死（用户 2026-10-29 反馈）。
        dlg = [];
        if any(k == [6 10])
            cla(axRep);
            text(axRep, 0.5, 0.55, items(k).name, 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', [0.25 0.30 0.45]);
            text(axRep, 0.5, 0.42, '正在计算，请稍候…', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', [0.45 0.45 0.50]);
            axRep.XTick = [];  axRep.YTick = [];
            drawnow;
            try
                dlg = uiprogressdlg(hfig, 'Title', '正在计算', ...
                                    'Message', [items(k).name ' …'], 'Indeterminate', 'on');
                drawnow;
            catch
                dlg = [];
            end
        end
        try
            lblRepMsg.Text = sprintf('第 %d 项：%s · %s', k, items(k).name, items(k).desc);
            drawnow;
            l3_draw(axRep, k, Ucur, S);
            lblRepMsg.Text = sprintf('第 %d 项：%s · %s', k, items(k).name, items(k).desc);
        catch err
            lblRepMsg.Text = sprintf('第 %d 项画图失败：%s', k, err.message);
        end
        if ~isempty(dlg)
            try, close(dlg); catch, end
        end
    end
end

%% ================= 以下是局部函数（不依赖界面状态，可在无界面下测试） =================

% ─────────────────────────────────────────────────────────────────────────────
% 函数：lui_default（缺省状态：本工程真实的 7 层，损耗取阶段 2 归档表）
% 输入：无
% 输出：S：状态结构体
% 说明：
%   损耗列直接抄 stage2_alpha_table(0.45,'TE','kio') 的逐层值（0 / 1.173209 / 1.401295 /
%   25.352696 / 10.034002 / 17.790425 / 26.685638 cm⁻¹），不是估的。
% ─────────────────────────────────────────────────────────────────────────────
function S = lui_default()
    tags  = {'substrate','n-clad','n-wave','QW','p-wave','EBL','p-clad'};
    dd    = [Inf,     1.200,   0.120,  0.009, 0.080,  0.020, 0.300];
    nn    = [1.770,   2.404,   2.455,  2.550, 2.455,  2.343, 2.374];
    al    = [0,       1.173209, 1.401295, 25.352696, 10.034002, 17.790425, 26.685638];
    ac    = [false, false, false, true, false, false, false];
    L = struct('tag', {}, 'n', {}, 'd', {}, 'alpha', {}, 'active', {});
    for j = 1:7
        L(j) = struct('tag', tags{j}, 'n', nn(j), 'd', dd(j), 'alpha', al(j), 'active', ac(j));
    end
    S = struct();
    S.kind = 'ridge';
    S.Lx = 10.000;  S.Ly = 2.500;  S.Lz = 1500.000;
    S.wr = 3.000;   S.he = 0.300;
    S.zoom = 1.0;   S.R = eye(3);   S.dc = [];   S.needview = true;   % 三维：缩放倍数 + 相机旋转矩阵（复位＝单位阵 = 正对截面）
    S.sel = 4;
    S.opt = struct('lambda', 0.45, 'pol', 'TE', 'num_scan', 20000, 'tol', 1e-13, ...
                   'merge', 1e-6, 'npts', 800, 'tail_dec', 20, 'num_ff', 2001, ...
                   'nshow', 3, 'msel', 0);
    S.n_top = 1.46;
    S.layers = L;
    S.page = 'main';
end

% ─────────────────────────────────────────────────────────────────────────────
% 函数：lui_live（给状态补上"运行时"字段：拖动排序与三维相机）
% 输入：S 结构状态（来自 lui_default 或任何重建的）
% 输出：S 结构状态（补全字段后）
% 说明：
%   !! 为什么单独一个函数：`S = lui_default()` 会抹掉这些字段，而图窗的
%   !! WindowButtonMotionFcn 一直在调用 on_motion，缺一个字段就每次移动都报错
%   !! （「无法识别的字段名称 camDrag」，用户 2026-10-29 报过）。凡是从默认值
%   !! 重建状态的地方，都要经过这里。
% ─────────────────────────────────────────────────────────────────────────────
function S = lui_live(S)
    S.dragSrc   = 0;        % 正在拖动的源层号（0 = 没在拖）
    S.dragMoved = false;    % 已移动超过阈值（区分"点一下"与"拖动"）
    S.dragOK    = true;     % 本版本能否取到鼠标坐标
    S.dragP0    = [0 0];
    S.camDrag   = false;    % 是否正在转三维相机
    S.camP0     = [0 0];
    S.camR0     = eye(3);
    if ~isfield(S, 'R') || isempty(S.R) || ~isequal(size(S.R), [3 3])
        S.R = eye(3);       % 相机旋转矩阵（单位阵 = 正对 x–y 截面）
    end
    if ~isfield(S, 'dc'), S.dc = []; end
    S.needview  = true;
    if ~isfield(S, 'zoom') || ~isfinite(S.zoom) || S.zoom <= 0
        S.zoom = 1.0;
    end
end

% ─────────────────────────────────────────────────────────────────────────────
% 函数：lui_stack（把界面层表变成模型要的 L0：自衬底向上、n 带损耗虚部）
% 输入：S 状态
% 输出：L0 结构体数组（.d/.n/.tag/.is_active）
% 说明：
%   损耗走 n 的正虚部（G5）：Im n = α / (2 k₀ · 1e4)，与 n_of_lambda 里同一式子。
% ─────────────────────────────────────────────────────────────────────────────
function L0 = lui_stack(S)
    k0 = 2*pi/S.opt.lambda;
    N = numel(S.layers);
    L0 = struct('d', {}, 'n', {}, 'tag', {}, 'is_active', {});
    for j = 1:N
        L = S.layers(j);
        L0(j) = struct('d', L.d, 'n', complex(L.n, L.alpha/(2*k0*1e4)), ...
                       'tag', L.tag, 'is_active', logical(L.active));
    end
end

% ─────────────────────────────────────────────────────────────────────────────
% 函数：lui_compute（按当前结构算一遍：竖直 → 横向 → 二维场 → 远场）
% 输入：S 状态
% 输出：U（模型 e3_ui_run 的返回）
% 说明：
%   矩形波导没有侧向限制，传 slab=true 让模型只做竖直部分，不去硬解一个不存在的横向本征问题。
% ─────────────────────────────────────────────────────────────────────────────
function U = lui_compute(S)
    isSlab = strcmp(S.kind, 'rect');
    he = S.he;
    if isSlab, he = 0; end
    U = Model_Self_L3_7L_EBL('stage3', struct( ...
            'sub', 'ui', 'L0', lui_stack(S), 'lambda', S.opt.lambda, ...
            'w', S.wr, 'h_e', he, 'n_top', S.n_top, 'slab', isSlab, ...
            'solver', S.opt, 'keep_windows', true));
end

% ─────────────────────────────────────────────────────────────────────────────
% 函数：lui_color（第 i 层的颜色）
% 输入：i 层号
% 输出：c 1×3 RGB
% ─────────────────────────────────────────────────────────────────────────────
function c = lui_color(i)
    pal = [0.20 0.45 0.75; 0.95 0.60 0.25; 0.35 0.70 0.40; 0.85 0.30 0.30; ...
           0.55 0.40 0.80; 0.90 0.78 0.20; 0.25 0.68 0.72; 0.80 0.45 0.70];
    c = pal(mod(i - 1, size(pal, 1)) + 1, :);
end

% ─────────────────────────────────────────────────────────────────────────────
% 函数：lui_box（在三维轴上画一个长方体）
% 输入：ax 轴; xr/yr/zr 三个方向的 [起 止]; col RGB; alpha 透明度
% 输出：无
% 说明：
%   六个面用一个 patch 画（8 顶点 + 6 面），比循环画 6 个 surf 快，光照也连续。
% ─────────────────────────────────────────────────────────────────────────────
function lui_box(ax, xr, yr, zr, col, alpha)
    x1 = xr(1); x2 = xr(2); y1 = yr(1); y2 = yr(2); z1 = zr(1); z2 = zr(2);
    V = [x1 y1 z1; x2 y1 z1; x2 y2 z1; x1 y2 z1; ...
         x1 y1 z2; x2 y1 z2; x2 y2 z2; x1 y2 z2];
    F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    patch(ax, 'Vertices', V, 'Faces', F, 'FaceColor', col, ...
          'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.5, ...
          'FaceAlpha', alpha, 'FaceLighting', 'none');
end

% ─────────────────────────────────────────────────────────────────────────────
% 函数：lui_items（成果画布的 11 个按钮：名称、说明、是否只在脊形下有意义）
% 输入：无
% 输出：items 结构体数组（.name/.desc/.slabOk）
% 说明：
%   slabOk = false 的几条与横向限制有关，矩形波导下会直接写明"不适用"，不画假图。
% ─────────────────────────────────────────────────────────────────────────────
function items = lui_items()
    items = struct('name', {}, 'desc', {}, 'slabOk', {});
    items(1) = mk('竖直结构剖面与模场', 'n(y) 台阶 + |Φ|²（有源层用色带标出）', true);
    items(2) = mk('整层竖直方向的远场（快轴）', '各阶竖直模的 1D 远场与半高全宽（阶数随结构）', true);
    items(3) = mk('竖直各阶数据表', 'Re n_eff、Γ、αᵢ（α 与层损耗的加权式对照）', true);
    items(4) = mk('横向模场 E_x(x)', '各阶 E_x 的强度分布，虚线是脊边 ±w/2', false);
    items(5) = mk('二维模场 |E(x,y)|²', 'E_x(x)·Φ(x,y)，品红线是 QW', false);
    items(6) = mk('n_eff / Γ / α vs 脊宽', '在当前脊宽的 0.5–2 倍内扫 5 个点', false);
    items(7) = mk('二维远场 θx–θy', '强度（dB），白线是 0.5 与 0.135 等高线', false);
    items(8) = mk('快轴远场（2D 场）', '各阶单画 + 等幅同相叠加（黑虚线）', false);
    items(9) = mk('慢轴远场（2D 场）', '各阶单画 + 等幅同相叠加（黑虚线）', false);
    items(10) = mk('EIM 与二维有限差分对比', '当前结构的 Re N：两步 EIM vs 五点差分（约 20 s）', false);
    items(11) = mk('损耗分解', '逐层 P_j·α_j·n_rj/Re N 的贡献（基模）', true);

    function s = mk(nm, ds, ok)
        s = struct('name', nm, 'desc', ds, 'slabOk', ok);
    end
end

% ─────────────────────────────────────────────────────────────────────────────
% 函数：l3_draw（总入口：按 k 选一张图画到 ax 上）
% 输入：ax 轴句柄（uiaxes 或普通 axes）; k 第几张; U 成果; S 状态
% 输出：无
% 说明：
%   画图一律用 ax 显式指定目标轴，所以同一套代码在界面（uiaxes）与自检（普通 axes）里通用。
% ─────────────────────────────────────────────────────────────────────────────
function l3_draw(ax, k, U, S)
    cla(ax, 'reset');
    % !! 每张图都从干净状态开始：上一张可能留下 axis tight / log 轴 / 手动坐标范围，
    % !! 不清就会串到下一张（用户 2026-10-28 反馈第 8 条"二维远场没法看"就是这么来的）。
    try
        axis(ax, 'normal');
    catch
    end
    ax.Color = [1 1 1];
    ax.XScale = 'linear';  ax.YScale = 'linear';
    ax.XDir = 'normal';    ax.YDir = 'normal';
    ax.CLimMode = 'auto';
    ax.XLimMode = 'auto';  ax.YLimMode = 'auto';  ax.ZLimMode = 'auto';
    ax.XTickMode = 'auto';  ax.YTickMode = 'auto';
    ax.XGrid = 'off';  ax.YGrid = 'off';
    ax.XTickLabelRotation = 0;
    if ~isempty(ax.Legend), delete(ax.Legend); end
    % !! colorbar 挂在图窗上、不是 axes 的 child，cla 清不掉；每画一张就会多叠一条
    % !! 旧色标（用户 2026-10-29 截图里色标同时出现 0..1 与 -30..0 两套刻度）。
    try
        fig0 = ancestor(ax, 'figure');
        delete(findall(fig0, 'Type', 'ColorBar'));
        tt = findall(fig0, 'Type', 'text');
        for q = 1:numel(tt)
            if isequal(tt(q).Parent, fig0), delete(tt(q)); end
        end
    catch
    end
    % !! 上一张若是 axis off 的文本表（第 3、10 项），坐标轴会一直藏着；显式打开
    ax.Visible = 'on';  ax.Box = 'on';
    ax.XColor = [0.15 0.15 0.15];  ax.YColor = [0.15 0.15 0.15];
    if strcmp(S.kind, 'rect')
        it = lui_items();
        if ~it(k).slabOk
            text(ax, 0.5, 0.5, {'矩形波导没有侧向限制', '', '这一张需要横向本征解，故不适用'}, ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                 'FontSize', 13, 'Color', [0.80 0.55 0.30]);
            return
        end
    end
    switch k
        case 1,  d_vert(ax, U, S);
        case 2,  d_far1(ax, U, S);
        case 3,  d_table(ax, U, S);
        case 4,  d_lat(ax, U, S);
        case 5,  d_field2d(ax, U, S);
        case 6,  d_sweep(ax, U, S);
        case 7,  d_ff2d(ax, U, S);
        case 8,  d_ffax(ax, U, S, true);
        case 9,  d_ffax(ax, U, S, false);
        case 10, d_fd(ax, U, S);
        case 11, d_loss(ax, U, S);
    end
end

function d_vert(ax, U, S)
% !! 用户 2026-10-29 的意见（第 1 张图）：层名不要压在折射率折线上，改到横轴下方排一行；
% !! QW、EBL 这类薄层在横轴上写不下，用箭头指到 n(y) 的平台上（QW 折射率高，从右上
% !! 2 点钟方向指下来；EBL 折射率低，从下方指上去），字号放大一档。
% !! 注意：yyaxis 有"当前作用侧"，text/plot 的 y 值按当前侧的量程解释。画完折射率必须
% !! 切回 left，否则层名会被丢到右轴量程（1.7–2.6）之外、整行消失（2026-10-29 踩过）。
    R = U.vert(1);
    y = R.y(:).';  I = abs(R.E(:).').^2;  I = I / max(I);
    [~, im] = max(I);
    tags = {R.L.tag};
    ylo = min(y);  yhi = max(y);
    nv = zeros(1, numel(y));
    for k = 1:numel(y)
        j = find(R.yb(:,1) <= y(k) + 1e-12, 1, 'last');
        if isempty(j), j = 1; end
        nv(k) = real(R.L(j).n);
    end
    nlo = min(nv);  nhi = max(nv);
    rng = [nlo - 0.05*(nhi-nlo), nhi + 0.12*(nhi-nlo)];
    ylL = [0 1.14];

    hold(ax, 'on');
    yyaxis(ax, 'left');
    plot(ax, y, I, '-', 'Color', [0.00 0.30 0.75], 'LineWidth', 2.0);
    ylim(ax, ylL);  yticks(ax, 0:0.25:1);
    ylabel(ax, '|Φ|^2（归一化）');
    ax.YColor = [0.00 0.30 0.75];

    yyaxis(ax, 'right');
    stairs(ax, y, nv, '-', 'Color', [0.85 0.45 0.10], 'LineWidth', 2.4);
    ylim(ax, rng);
    ylabel(ax, '折射率 n(y)（实部）');
    ax.YColor = [0.75 0.40 0.05];
    yyaxis(ax, 'left');

    xlim(ax, [ylo yhi]);
    dy = diff(ylL);  xr = [ylo yhi];  dxr = diff(xr);

    iq = find(strcmp(tags, 'QW'), 1);
    if ~isempty(iq)
        plot(ax, [R.yb(iq,1) R.yb(iq,1)], ylL, 'm:', 'LineWidth', 1.2);
        plot(ax, [R.yb(iq,2) R.yb(iq,2)], ylL, 'm:', 'LineWidth', 1.2);
    end

    % !! 宽层：层名放横轴下方，拥挤时两行交替
    thin = [];  nrow = 0;  lastC = -Inf;
    for j = 1:numel(R.L)
        y1 = R.yb(j,1);  y2 = R.yb(j,2);
        if ~isfinite(y1) || ~isfinite(y2), continue, end
        if (y2 - y1) < 0.055
            thin(end+1) = j; %#ok<AGROW>
            continue
        end
        yc = (y1 + y2) / 2;
        if (yc - lastC) < 0.18 && nrow == 0
            nrow = 1;
        else
            nrow = 0;
        end
        lastC = yc;
        yt = 0.075 - 0.043*nrow;
        text(ax, yc, yt, tags{j}, 'FontSize', 10.5, 'Clipping', 'off', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
             'Color', [0.45 0.22 0.02], 'Interpreter', 'none');
    end

    % !! 薄层：箭头指认（折射率换算到左轴坐标）
    for q = 1:numel(thin)
        j = thin(q);
        yc = (R.yb(j,1) + R.yb(j,2)) / 2;
        nq = real(R.L(j).n);
        yv = (nq - rng(1)) / diff(rng) * dy + ylL(1);
        jb = max(j-1, 1);
        if nq > real(R.L(jb).n)
            p1 = [yc, yv - 0.012*dy];
            p0 = [yc + 0.130*dxr, yv + 0.030*dy];
        else
            p1 = [yc, yv - 0.010*dy];
            p0 = [yc + 0.010*dxr, yv - 0.290*dy];
        end
        vert_arrow(ax, p0, p1, [0.55 0.10 0.45], dxr, dy);
        text(ax, p0(1) + 0.030*dxr, p0(2), ...
             {tags{j}, sprintf('n = %.3f · %g nm', nq, round(1000*(R.yb(j,2)-R.yb(j,1))))}, ...
             'FontSize', 10, 'Clipping', 'off', 'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'middle', 'Color', [0.55 0.10 0.45], 'Interpreter', 'none');
    end
    hold(ax, 'off');

    xlabel(ax, 'y [μm]（自衬底界面向上）', 'FontSize', 10.5);
    title(ax, {sprintf('竖直剖面 · 基模 Γ = %.5f · Re n_{eff} = %.6f · 场峰 y = %.3f μm', ...
                       R.G, real(R.ne), y(im)), ...
               sprintf('自上而下：%s', strjoin(tags(end:-1:1), ' / '))}, ...
          'Interpreter', 'none', 'FontSize', 10.5);
    grid(ax, 'on');
end

function vert_arrow(ax, p0, p1, col, dx, dy)
% ─────────────────────────────────────────────────────
% 函数：vert_arrow（在数据坐标里画一根带三角箭头的短线）
% 输入：ax 坐标轴；p0 起点 [x y]；p1 尖端 [x y]（数据坐标）
%       col 颜色；dx/dy 横纵坐标跨度（把像素尺寸折算成数据尺寸）
% 输出：无
% 说明：
%   箭头几何按像素给（头长 16 px、半宽 6 px）再折算回数据单位；直接用数据单位
%   画头会因为横纵比例悬殊（3.3 μm 对 1.14）而变形成一条色块。
%   uiaxes 的 Position 是归一化的，必须用 getpixelposition 才拿得到真实像素。
% ─────────────────────────────────────────────────────
    pp = getpixelposition(ax, true);
    wp = max(pp(3), 1);  hp = max(pp(4), 1);
    sx = dx / wp;  sy = dy / hp;
    vx = (p1(1) - p0(1)) / sx;  vy = (p1(2) - p0(2)) / sy;
    L = hypot(vx, vy);  if L < 1e-9, return, end
    ux = vx / L;  uy = vy / L;
    bx = p1(1) - 16*ux*sx;  by = p1(2) - 16*uy*sy;
    qx = -uy;  qy = ux;
    pa = [p1(1), p1(2); ...
          bx + 6*qx*sx, by + 6*qy*sy; ...
          bx - 6*qx*sx, by - 6*qy*sy];
    plot(ax, [p0(1) bx], [p0(2) by], '-', 'Color', col, 'LineWidth', 1.3, ...
         'Clipping', 'off');
    patch(ax, 'XData', pa(:,1), 'YData', pa(:,2), 'FaceColor', col, ...
          'EdgeColor', col, 'Clipping', 'off');
end
function d_far1(ax, U, S)
    hold(ax, 'on');
    K = numel(U.far1);
    pal = lui_pal();
    for j = 1:K
        c = pal(mod(j-1, size(pal,1)) + 1, :);
        plot(ax, rad2deg(U.far1(j).theta), U.far1(j).I, '-', 'Color', c, ...
             'LineWidth', 2.0 - 0.6*(j > 1));
    end
    hold(ax, 'off');
    xlim(ax, [-85 85]);  ylim(ax, [0 1.05]);
    xlabel(ax, 'θ_y [deg]（快轴方向出射角）', 'FontSize', 10.5);  ylabel(ax, '相对强度（各阶自归一）', 'FontSize', 10.5);
    lg = arrayfun(@(j) sprintf('TE%d  FWHM %.2f°', j-1, rad2deg(U.far1(j).fwhm)), ...
                  1:K, 'UniformOutput', false);
    legend(ax, lg, 'Location', 'northeast');
    % !! 用户 2026-10-29：原标题「竖直模场快轴远场」看不懂。改成先说清楚画的是什么方向的远场。
    title(ax, {sprintf('竖直方向（快轴）远场 · %s 各阶竖直模', S.opt.pol), '每条曲线各自归一到峰值 1；纯竖直结构、不含横向模叠加（对比第 8 项）'}, 'Interpreter', 'none', 'FontSize', 10.5);
    grid(ax, 'on');
end

function d_table(ax, U, S)
% !! 2026-10-29 改版：原来是把几行等宽文本丢在空白画布上，下半张全空。
% !! 现在做成"带表头的表格 + 浅灰底 + 下方注释块"，并把结构参数一并列出。
    axis(ax, 'off');
    C = U.C;
    K = numel(C.N);
    hold(ax, 'on');

    xc = [0.06 0.15 0.37 0.50 0.66 0.82];
    hdr = {'m', 'Re n_eff', 'Im n_eff', 'Gamma', 'alpha_i [1/cm]', 'alpha_p [1/cm]'};
    y0 = 0.86;  dy = 0.058;

    % 表头底色
    patch(ax, [0.02 0.98 0.98 0.02], [y0+0.030 y0+0.030 y0-0.030 y0-0.030], ...
          [0.90 0.92 0.96], 'EdgeColor', 'none');
    for c = 1:numel(hdr)
        text(ax, xc(c), y0, hdr{c}, 'FontName', 'Microsoft YaHei', 'FontSize', 10.5, ...
             'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
    end
    plot(ax, [0.02 0.98], [y0-dy+0.008 y0-dy+0.008], '-', 'Color', [0.45 0.45 0.50], 'LineWidth', 1.0);

    for j = 1:K
        yr = y0 - dy * j;
        if mod(j, 2) == 0
            patch(ax, [0.02 0.98 0.98 0.02], [yr+0.030 yr+0.030 yr-0.030 yr-0.030], ...
                  [0.965 0.970 0.980], 'EdgeColor', 'none');
        end
        Im = NaN;
        if isfield(C, 'Nc') && numel(C.Nc) >= j, Im = imag(C.Nc(j)); end
        Gm = C.Gamma(min(j, numel(C.Gamma)));
        a1 = C.alpha(min(j, numel(C.alpha)));
        ap = NaN;
        if isfield(C, 'alpha_pred') && numel(C.alpha_pred) >= j, ap = C.alpha_pred(j); end
        vals = {sprintf('m=%d', j-1), sprintf('%.9f', real(C.N(j))), sprintf('%.3e', Im), ...
                sprintf('%.6f', Gm), sprintf('%.4f', a1), sprintf('%.4f', ap)};
        for c = 1:numel(vals)
            text(ax, xc(c), yr, vals{c}, 'FontName', 'Consolas', 'FontSize', 12, ...
                 'VerticalAlignment', 'middle', 'Interpreter', 'none');
        end
    end
    ybot = y0 - dy * K - 0.030;
    plot(ax, [0.02 0.98], [ybot ybot], '-', 'Color', [0.45 0.45 0.50], 'LineWidth', 1.0);

    % 说明块
    note = {sprintf('结构：%d 层，脊宽 %.3f μm，脊高 %.3f μm，λ = %.3f μm，顶部介质 n = %.2f，横向网格 dx = 10 nm', ...
                    numel(S.layers), S.wr, S.he, S.opt.lambda, S.n_top), ...
            sprintf('侧向：Re n_a = %.9f，Re n_b = %.9f，V = %.4f，公式给 M = %d 阶，实解 %d 阶', ...
                    real(C.na), real(C.nb), C.V, C.M_expect, K), ...
            'alpha_p = 逐层 P_j·alpha_j·n_rj / Re N 的加权和（与直接算出的 alpha_i 对照）'};
    yt = ybot - 0.030;  yb2 = ybot - 0.185;
    patch(ax, [0.02 0.98 0.98 0.02], [yt yt yb2 yb2], ...
          [0.975 0.975 0.965], 'EdgeColor', [0.80 0.80 0.75]);
    for k = 1:numel(note)
        text(ax, 0.035, yt - 0.042*(k-1) - 0.018, note{k}, 'FontName', 'Microsoft YaHei', ...
             'FontSize', 10.5, 'VerticalAlignment', 'middle', 'Interpreter', 'none');
    end
    if strcmp(S.kind, 'ridge'), knm = '脊形波导'; else, knm = '矩形波导'; end
    title(ax, sprintf('竖直各阶数据（m = 0 … %d）· %s', K-1, knm), ...
          'FontName', 'Microsoft YaHei', 'FontSize', 11.5, 'Interpreter', 'none');
    hold(ax, 'off');
end

function d_lat(ax, U, S)
    hold(ax, 'on');
    pal = lui_pal();
    % !! 只画用户选的阶数，但标题必须写出「一共几阶」，否则看图的人会以为只有 K 阶
    K = min(max(1, S.opt.nshow), numel(U.C.xc));
    for j = 1:K
        x = U.C.xc{j};  E = U.C.Ex{j};
        if isempty(x), continue; end
        c = pal(mod(j-1, size(pal,1)) + 1, :);
        plot(ax, x, abs(E).^2 / max(abs(E).^2), '-', 'Color', c, 'LineWidth', 2.0 - 0.6*(j > 1));
    end
    xline(ax, -S.wr/2, 'k--');  xline(ax, S.wr/2, 'k--');
    hold(ax, 'off');
    xl = 2*S.wr;
    xlim(ax, [-xl xl]);
    xlabel(ax, 'x [μm]（芯中心为 0）');  ylabel(ax, '|E_x|² (归一化)');
    lg = arrayfun(@(j) sprintf('m=%d  Γ=%.5f', j-1, U.C.Gamma(min(j, numel(U.C.Gamma)))), ...
                  1:K, 'UniformOutput', false);
    legend(ax, lg, 'Location', 'southoutside', 'NumColumns', 2);
    title(ax, {sprintf('横向模场 E_x(x) · 脊宽 %.2f μm', S.wr), ...
               sprintf(['当前结构共 %d 阶（由 V = %.3f 定），这里显示前 %d 阶；' ...
                        '要看更多用表头「显示阶数」'], numel(U.C.xc), U.C.V, K)}, ...
          'Interpreter', 'none', 'FontSize', 10.5);
    grid(ax, 'on');
end

function d_field2d(ax, U, S)
    Fd = U.fd;
    qm = min(max(1, S.opt.msel + 1), size(Fd.E, 3));
    E2 = abs(Fd.E(:,:,qm)).^2;
    E2 = E2 / max(E2(:));
    surf(ax, Fd.x, Fd.y, zeros(size(E2)), E2, 'EdgeColor', 'none');
    view(ax, 0, 90);  axis(ax, 'tight');
    xlim(ax, [-3*S.wr/2, 3*S.wr/2]);
    yl = [min(Fd.y) max(Fd.y)];
    iq = find(strcmp(U.C.tags, 'QW'), 1);
    hold(ax, 'on');
    yp = peak_row(Fd.y, E2);
    ylim(ax, [yp-0.25, yp+0.25]);
    yq = U.C.Ra(1).yb;
    m = find(strcmp({U.C.Ra(1).L.tag}, 'QW'), 1);
    if ~isempty(m)
        yline(ax, [yq(m,1) yq(m,2)], 'm-', 'LineWidth', 1.4);
    end
    % !! 脊边位置要画出来：不然看不出场是被 |x| < w/2 的那一段夹住的（2026-10-29 看图后补）
    yy = ylim(ax);
    plot(ax, [S.wr/2 S.wr/2], yy, 'w--', 'LineWidth', 1.1);
    plot(ax, [-S.wr/2 -S.wr/2], yy, 'w--', 'LineWidth', 1.1);
    text(ax, 0, yy(2) - 0.03*diff(yy), sprintf('脊宽 w = %.2f μm（白虚线＝脊边）', S.wr), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'w', 'FontSize', 9);
    hold(ax, 'off');
    colormap(ax, parula(256));
    cb = colorbar(ax);
    cb.Label.String = '|E|^2（归一化）';
    cb.Label.FontSize = 10;
    xlabel(ax, 'x [μm]');  ylabel(ax, 'y [μm]（自衬底界面向上）');
    title(ax, sprintf('二维模场 |E|²（m=0）· 峰值 y=%.3f μm · 品红线是 QW 上下沿 · Γ=%.5f', ...
                      yp, U.C.Gamma(1)), 'Interpreter', 'none');
end

function yp = peak_row(y, E2)
    [~, r] = max(sum(E2, 2));
    yp = y(min(r, numel(y)));
end

function d_sweep(ax, U, S)
% !! 2026-10-29 改版：原来把 Γ、α 的数字打印在坐标区左下角，压在刻度和轴标签上；
% !! 现在 Γ 走右轴画曲线、α 写进标题，坐标区里只放数据，并在当前 w 处加红圈。
    w0 = S.wr;
    ws = w0 * [0.5 0.75 1 1.25 1.5];
    ws = min(ws, S.Lx);
    Nn = nan(size(ws));  Gg = nan(size(ws));  Aa = nan(size(ws));
    for j = 1:numel(ws)
        Sj = S;  Sj.wr = ws(j);
        try
            Uj = lui_compute(Sj);
            if Uj.ok
                Nn(j) = Uj.C.N(1);  Gg(j) = Uj.C.Gamma(1);  Aa(j) = Uj.C.alpha(1);
            end
        catch
        end
    end
    dN = real(Nn) - real(U.C.na);
    yyaxis(ax, 'left');
    h1 = semilogx(ax, ws, dN, 'o-', 'LineWidth', 1.9, 'MarkerFaceColor', 'w');
    hold(ax, 'on');
    i0 = find(abs(ws - w0) < 1e-9, 1);
    if ~isempty(i0)
        semilogx(ax, ws(i0), dN(i0), 'o', 'MarkerSize', 10, 'LineWidth', 2.2, ...
                 'Color', [0.85 0.10 0.10]);
    end
    hold(ax, 'off');
    ylabel(ax, 'Re N_1 - Re n_a');
    ax.YColor = [0.00 0.30 0.75];
    yyaxis(ax, 'right');
    h2 = plot(ax, ws, Gg, 's--', 'LineWidth', 1.5, 'MarkerSize', 6);
    ylabel(ax, 'Gamma_1 (基模限制因子)');
    ax.YColor = [0.75 0.40 0.05];
    yyaxis(ax, 'left');
    xlabel(ax, '脊宽 w [μm]');
    grid(ax, 'on');
    try
        legend(ax, [h1 h2], {'Re N_1 - Re n_a (左轴)', 'Gamma_1 (右轴)'}, ...
               'Location', 'northwest', 'FontSize', 9);
    catch
    end
    title(ax, sprintf(['脊宽扫描 %.2f ~ %.2f μm，红圈 = 当前 w = %.2f μm\n' ...
                       'alpha_i (基模) = %s cm^{-1}'], min(ws), max(ws), w0, ...
                      strjoin(arrayfun(@(v) sprintf('%.4f', v), Aa, 'UniformOutput', false), ' / ')), ...
          'Interpreter', 'tex', 'FontSize', 10.5);
end
function d_ff2d(ax, U, S)
    I = 10*log10(max(U.ff.I2, 1e-8));
    imagesc(ax, U.ff.sx, U.ff.sy, I);
    set(ax, 'YDir', 'normal');
    hold(ax, 'on');
    contour(ax, U.ff.sx, U.ff.sy, U.ff.I2, [0.5 0.135], 'k-', 'LineWidth', 1.1);
    plot(ax, [min(U.ff.sx) max(U.ff.sx)], [0 0], 'w:', 'LineWidth', 0.8);
    plot(ax, [0 0], [min(U.ff.sy) max(U.ff.sy)], 'w:', 'LineWidth', 0.8);
    hold(ax, 'off');
    % !! 显式把坐标范围钉在数据范围上：万一上一张留了手动 XLim，imagesc 会被挤成一条竖线
    xlim(ax, [min(U.ff.sx) max(U.ff.sx)]);
    ylim(ax, [min(U.ff.sy) max(U.ff.sy)]);
    caxis(ax, [-30 0]);  colormap(ax, parula(256));
    % !! 加了色标：没有色标时"归一化到 0 dB"这句话在图上无法兑现（2026-10-29 看图后补）
    cb = colorbar(ax);
    cb.Label.String = '相对强度 [dB]';
    cb.Label.FontSize = 10;
    xlabel(ax, 'sin\theta_x（慢轴）');  ylabel(ax, 'sin\theta_y（快轴）');
    title(ax, sprintf(['二维远场 · %d 阶等幅同相叠加（假设值）· 归一化到峰值 0 dB\n' ...
                       '黑线：0.5 与 0.135 等高线；θ_x = θ_y = 0 见白色点线'], U.ff.K), ...
          'Interpreter', 'tex');
end

function d_ffax(ax, U, S, isFast)
    hold(ax, 'on');
    pal = lui_pal();
    for j = 1:numel(U.ff.mode)
        if isFast
            u = U.ff.mode(j).uF;  I = U.ff.mode(j).IF;
        else
            u = U.ff.mode(j).uS;  I = U.ff.mode(j).IS;
        end
        c = pal(mod(j-1, size(pal,1)) + 1, :);
        plot(ax, rad2deg(u), I, '-', 'Color', c, 'LineWidth', 2.0 - 0.6*(j > 1));
    end
    if isFast
        uo = U.ff.uF;  Io = U.ff.IF;
        nm = '快轴（固定 θx=0）';
        xl = 'θy [deg]';
    else
        uo = U.ff.uS;  Io = U.ff.IS;
        nm = '慢轴（固定 θy=0）';
        xl = 'θx [deg]';
    end
    plot(ax, rad2deg(uo), Io, 'k--', 'LineWidth', 1.8);
    hold(ax, 'off');
    % !! 横轴按数据自己撑到哪就画到哪：一律 ±85° 时慢轴只有中间 ±20° 有内容，
    % !! 两边大片空白（用户 2026-10-29 第 9 项反馈）。判据用"到场值 0.5%"，再留 12% 余量。
    xd = rad2deg(uo(:));
    xw = 0;
    for j = 1:numel(U.ff.mode)
        if isFast
            Ij = U.ff.mode(j).IF(:);
        else
            Ij = U.ff.mode(j).IS(:);
        end
        keep = Ij >= 0.02;
        if any(keep), xw = max(xw, max(abs(xd(keep)))); end
    end
    keep = Io(:) >= 0.02;
    if any(keep), xw = max(xw, max(abs(xd(keep)))); end
    if xw <= 0, xw = 20; end
    xw = min(max(1.10*xw, 3), 85);
    xlim(ax, [-xw xw]);  ylim(ax, [0 1.05]);
    xlabel(ax, xl);  ylabel(ax, '相对强度');
    lg = [arrayfun(@(j) sprintf('m=%d', j-1), 1:numel(U.ff.mode), 'UniformOutput', false), ...
          {sprintf('%d 阶等幅叠加（假设值）', U.ff.K)}];
    legend(ax, lg, 'Location', 'northeast');
    title(ax, sprintf('%s远场 · 脊宽 %.2f μm', nm, S.wr), 'Interpreter', 'none');
    grid(ax, 'on');
end

function d_fd(ax, U, S)
% !! 2026-10-29 改版：原来只有一屏等宽文本、下半张全空。
% !! 现在做成"参数表 + 说明"两段式。注意右列必须是纯 ASCII——Consolas 没有中文字形，
% !! 中英混排会出方块（这张图第一版就踩过）。
    axis(ax, 'off');
    text(ax, 0.5, 0.55, '正在做二维有限差分（约 20 s）…', 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'FontSize', 13);
    drawnow;
    try
        R = Model_Self_L3_7L_EBL('stage3', struct('sub', 'ui_fd', 'L0', lui_stack(S), ...
                'lambda', S.opt.lambda, 'w', S.wr, 'h_e', S.he, 'n_top', S.n_top, ...
                'shift_n', real(U.C.N(1)), 'keep_windows', true));
    catch err
        cla(ax);
        text(ax, 0.5, 0.5, ['二维有限差分失败：' err.message], 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'FontSize', 11);
        return
    end
    cla(ax);
    axis(ax, 'off');
    if ~isfield(R, 'Nf') || ~isfinite(R.Nf)
        text(ax, 0.5, 0.5, '二维有限差分没有解出导模', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'FontSize', 12);
        return
    end
    Ne   = real(U.C.N(1));
    d    = abs(Ne - R.Nf);
    crit = 5e-3;
    if d <= crit, okstr = '通过'; else, okstr = '不通过'; end
    rl = d / max(Ne, eps);

    left = { 'EIM（两步有效折射率，本阶段交付值）'
             '二维有限差分（独立校核：非均匀网格 + Dirichlet 边界）'
             sprintf('两者之差（相对差 %.2e）', rl)
             sprintf('A4 判据 |Δ| ≤ 5e-3 —— %s，余量 %.0f 倍', okstr, crit/max(d, eps))
             '二维 FD 规模（未知数 / 芯内功率占比）' };
    right = { sprintf('Re N_1 = %.9f', Ne)
              sprintf('Re N_1 = %.9f', R.Nf)
              sprintf('|Delta| = %.3e', d)
              sprintf('%.2e / 5e-03 = %.3f', d, d/crit)
              sprintf('%d / %.3f', R.n, R.frac) };

    patch(ax, [0.04 0.96 0.96 0.04], [0.97 0.97 0.40 0.40], [0.975 0.978 0.985], ...
          'EdgeColor', [0.80 0.82 0.86]);
    y0 = 0.905;  dy = 0.105;
    for k = 1:numel(left)
        yk = y0 - dy * (k-1);
        text(ax, 0.07, yk, left{k}, 'FontName', 'Microsoft YaHei', 'FontSize', 11, ...
             'VerticalAlignment', 'middle', 'Interpreter', 'none');
        text(ax, 0.63, yk, right{k}, 'FontName', 'Consolas', 'FontSize', 11.5, ...
             'VerticalAlignment', 'middle', 'Interpreter', 'none');
    end

    text(ax, 0.07, 0.30, {'EIM：先解竖直方向的复模，得到脊区 n_a 与沟槽区 n_b，', ...
                          '    再以这两个折射率解一次横向模——本阶段交付的 n_eff/Γ/α 都出自它。', ...
                          '二维 FD：直接在 (x, y) 上解标量 Helmholtz 方程，不解两次一维问题，', ...
                          '    因此它与 EIM 的偏差是"两条独立路线之差"，不是自证。'}, ...
         'Units', 'normalized', 'FontName', 'Microsoft YaHei', 'FontSize', 10.5, ...
         'VerticalAlignment', 'top', 'Interpreter', 'none');
    title(ax, sprintf('EIM 与二维有限差分对比（脊宽 w = %.3f μm，脊高 %.3f μm）', S.wr, S.he), ...
          'FontName', 'Microsoft YaHei', 'FontSize', 11.5, 'Interpreter', 'none');
end
function d_loss(ax, U, S)
    qm = min(max(1, S.opt.msel + 1), size(U.C.Ptag, 1));
    P = U.C.Ptag(qm, :);
    tags = U.C.tags;
    alx = U.C.alx;
    nr = U.C.nr;
    n = min(numel(P), numel(alx));
    contrib = zeros(1, n);
    for j = 1:n
        contrib(j) = P(j) * alx(j) * nr(j) / real(U.C.N(qm));
    end
    b = bar(ax, contrib, 'FaceColor', 'flat');
    % !! 层名直接写在 x 轴上（用户：损耗图没有层的名字），有源层用暖色突出
    for j = 1:n
        if strcmpi(tags{j}, 'QW')
            b.CData(j, :) = [0.85 0.30 0.30];
        else
            b.CData(j, :) = [0.30 0.55 0.85];
        end
        text(ax, j, contrib(j) + max(contrib)*0.04, sprintf('%.2f', contrib(j)), ...
             'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.20 0.20 0.20]);
    end
    xlim(ax, [0.3 n+0.7]);
    ylim(ax, [0 min(1.25*max(contrib), 1.25*max(contrib))]);
    set(ax, 'XTick', 1:n, 'XTickLabel', tags(1:n), 'XTickLabelRotation', 30);
    ylabel(ax, 'α_i 贡献 [cm^{-1}]');
    title(ax, sprintf(['损耗分解（m=%d）：α_i = Σ_j P_j α_j n_{rj} / Re N ' ...
                       '· 合计 %.4f cm⁻¹（直接算得 %.4f）'], qm - 1, sum(contrib), ...
                      U.C.alpha(qm)), 'Interpreter', 'none', 'FontSize', 10);
    grid(ax, 'on');
end

function pal = lui_pal()
    pal = [0 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188; ...
           0.494 0.184 0.556; 0.301 0.745 0.933; 0.635 0.078 0.184];
end

% ─────────────────────────────────────────────────────────────────────────────
% 函数：l3_selftest（无界面自检：按缺省结构算一遍，把 11 张图各画到普通 axes 上并存盘）
% 输入：无
% 输出：U：模型算出来的成果（便于命令行检查数值）
% 说明：
%   本机 matlab -batch 环境创建 uifigure 会挂（无桌面/无 CEF），但普通 figure 正常，
%   所以把画图逻辑全放进接受 ax 的局部函数，再用普通 axes 跑一遍 —— 画的代码与界面逐行相同。
% ─────────────────────────────────────────────────────────────────────────────
function U = l3_selftest()
    S = lui_default();
    outd = fullfile(fileparts(mfilename('fullpath')), 'out');
    if ~exist(outd, 'dir'), mkdir(outd); end
    U = lui_compute(S);
    fprintf('[selftest] ok=%d  t=%.2f s  msg=%s\n', U.ok, U.t_sec, U.msg);
    if ~U.ok
        return
    end
    fprintf('[selftest] Re N = %s\n', mat2str(U.C.N, 9));
    fprintf('[selftest] Gamma = %s\n', mat2str(U.C.Gamma, 6));
    fprintf('[selftest] alpha = %s\n', mat2str(U.C.alpha, 5));
    items = lui_items();
    for k = 1:11
        f = figure('Visible', 'off', 'Position', [80 60 1000 620]);
        ax = axes(f);
        l3_draw(ax, k, U, S);
        print(f, fullfile(outd, sprintf('_ui_%02d.png', k)), '-dpng', '-r110');
        close(f);
        fprintf('[selftest] 第 %d 张：%s\n', k, items(k).name);
    end
end
