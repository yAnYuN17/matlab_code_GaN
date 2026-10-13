function startup()
%STARTUP  阶段1「非对称多层平板波导本征方程」环境初始化
%
%   用法：把 Mode_L1 设为当前目录后执行
%       run('startup.m')
%
%   ── 全项目单位约定（唯一定义处，不得在别处另立门户）─────────────
%       长度        μm              层厚 d、波长 λ、模场横坐标 y
%       折射率      无量纲           n = n' + 1i·n''
%       真空波数    k0 = 2π/λ        [μm^-1]
%       传播常数    β  = k0·n_eff    [μm^-1]     转 cm^-1 乘 1e4
%       吸收系数    α  [cm^-1] = 2·k0·Im(n_eff)·1e4
%
%   ── 符号与时间约定（阶段2 起至关重要，先在此写死）───────────────
%       时间因子取 exp(-iωt)，沿 +z 传播取 exp(+iβz)：
%           E(y,z,t) = φ(y)·exp(i(βz - ωt)),   β = k0·n_eff
%       于是光强 I(z) ∝ exp(-2·Im(β)·z)，得到
%           α = 2·Im(β) = 2·k0·Im(n_eff)
%       结论：Im(n_eff) > 0 表示【损耗】，Im(n_eff) < 0 表示【净增益】。
%       阶段1 全部取 Im(n_eff) = 0，上式不生效，但约定先行写死。
%   ────────────────────────────────────────────────────────────
%
%   为什么这条注释必须存在：阶段2 引入复折射率后，若符号约定与前序
%   不一致，"增益算成损耗"这类错误会静默通过全部数值检查。

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'layers'));
addpath(fullfile(root, 'tests'));

fprintf('[阶段1] 环境就绪：%s\n', root);
fprintf('[阶段1] 单位：长度 μm ｜ β μm^-1 ｜ α cm^-1\n');
fprintf('[阶段1] 约定：exp(-iωt)·exp(+iβz)，Im(n_eff) > 0 即损耗\n');
end
