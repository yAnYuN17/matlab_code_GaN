function layers = layers_three(n_sub, n_core, n_cov, d_core, i_active)
%LAYERS_THREE  构造严格三层非对称平板波导的层结构
%
%   layers = layers_three(n_sub, n_core, n_cov, d_core)
%   layers = layers_three(n_sub, n_core, n_cov, d_core, i_active)
%
%   层序（y 从小到大）：
%       L(1)  下包层 substrate   d = Inf   n = n_sub
%       L(2)  芯层   core        d = d_core n = n_core
%       L(3)  上包层 cover       d = Inf   n = n_cov
%
%   可选参数 i_active : 芯层内哪些部分算有源区。三种用法
%       （省略）            芯层不标记为有源区
%       'all'              整个芯层标记为有源区
%
%   说明：这是所有验证工作的基准结构。三层非对称平板波导有闭式色散方程
%   （见 dispersion_analytic3.m），是唯一能被独立核对的解析基准。
%
%   注意 n_sub 与 n_cov 谁大谁小都可以：当 n_cov > n_sub 时，模场会偏向
%   上侧，V 数与不对称度的定义要随之对调。代码内部按"较高的那个包层"
%   统一处理，不要在任何地方预设 n_sub > n_cov。

if nargin < 5, i_active = false; end

layers = struct( ...
    'd',         {Inf,      d_core,      Inf}, ...
    'n',         {n_sub,    n_core,      n_cov}, ...
    'tag',       {'substrate', 'core',   'cover'}, ...
    'is_active', {false,    false,       false});

if ischar(i_active) && strcmpi(i_active, 'all')
    layers(2).is_active = true;
elseif islogical(i_active)
    layers(2).is_active = i_active;
end
end
