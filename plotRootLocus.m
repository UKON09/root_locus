function [r, k, k_critical, asymptotes, angles] = plotRootLocus(sys)
% PLOTROOTLOCUS 绘制系统根轨迹图
%   [r, k, k_critical, asymptotes, angles] = plotRootLocus(sys)
%
% 输入参数：
%   sys          - 开环传递函数 (tf 或 zpk 对象)
%
% 输出参数：
%   r          - 根轨迹上的极点位置
%   k          - 对应的增益向量
%   k_critical - 临界增益值数组（自动计算的虚轴交点）
%   asymptotes - 渐近线信息结构体
%   angles     - 出射角和入射角信息结构体

% 获取零极点
[sys_zeros, sys_poles, ~] = zpkdata(sys, 'v');
n = length(sys_poles);  % 极点数
m = length(sys_zeros);  % 零点数

% 计算并输出零极点与实轴的夹角
fprintf('\n=== 零极点与实轴夹角 ===\n');

% 计算极点与实轴夹角
fprintf('极点与实轴夹角：\n');
for i = 1:length(sys_poles)
    pole = sys_poles(i);
    pole_angle = angle(pole) * 180 / pi;  % 转换为度数
    if imag(pole) == 0
        fprintf('极点 %d: %.4f (实数极点，夹角 = 0°)\n', i, real(pole));
    else
        fprintf('极点 %d: %.4f%+.4fj，夹角 = %.1f°\n', i, real(pole), imag(pole), pole_angle);
    end
end

% 计算零点与实轴夹角（如果存在零点）
if ~isempty(sys_zeros)
    fprintf('\n零点与实轴夹角：\n');
    for i = 1:length(sys_zeros)
        zero = sys_zeros(i);
        zero_angle = angle(zero) * 180 / pi;  % 转换为度数
        if imag(zero) == 0
            fprintf('零点 %d: %.4f (实数零点，夹角 = 0°)\n', i, real(zero));
        else
            fprintf('零点 %d: %.4f%+.4fj，夹角 = %.1f°\n', i, real(zero), imag(zero), zero_angle);
        end
    end
else
    fprintf('\n系统无零点\n');
end

%% 计算复数零极点与其他点的夹角
fprintf('\n=== 复数零极点与其他点的夹角分析 ===\n');

% 分析复数极点与其他点的夹角
complex_poles_idx = find(imag(sys_poles) ~= 0);
if ~isempty(complex_poles_idx)
    fprintf('复数极点与其他点的夹角：\n');
    for i = 1:length(complex_poles_idx)
        idx = complex_poles_idx(i);
        pole = sys_poles(idx);
        fprintf('\n极点 %d (%.4f%+.4fj) 的角度分析：\n', idx, real(pole), imag(pole));
        
        % 与其他极点的夹角
        fprintf('  与其他极点的夹角：\n');
        for j = 1:length(sys_poles)
            if j ~= idx
                other_pole = sys_poles(j);
                angle_diff = angle(pole - other_pole) * 180 / pi;
                distance = abs(pole - other_pole);
                if imag(other_pole) == 0
                    fprintf('    到极点 %d (%.4f): 角度 = %.1f°, 距离 = %.4f\n', ...
                            j, real(other_pole), angle_diff, distance);
                else
                    fprintf('    到极点 %d (%.4f%+.4fj): 角度 = %.1f°, 距离 = %.4f\n', ...
                            j, real(other_pole), imag(other_pole), angle_diff, distance);
                end
            end
        end
        
        % 与零点的夹角（如果存在）
        if ~isempty(sys_zeros)
            fprintf('  与零点的夹角：\n');
            for j = 1:length(sys_zeros)
                zero = sys_zeros(j);
                angle_diff = angle(pole - zero) * 180 / pi;
                distance = abs(pole - zero);
                if imag(zero) == 0
                    fprintf('    到零点 %d (%.4f): 角度 = %.1f°, 距离 = %.4f\n', ...
                            j, real(zero), angle_diff, distance);
                else
                    fprintf('    到零点 %d (%.4f%+.4fj): 角度 = %.1f°, 距离 = %.4f\n', ...
                            j, real(zero), imag(zero), angle_diff, distance);
                end
            end
        end
    end
end

% 分析复数零点与其他点的夹角
complex_zeros_idx = find(imag(sys_zeros) ~= 0);
if ~isempty(complex_zeros_idx)
    fprintf('\n复数零点与其他点的夹角：\n');
    for i = 1:length(complex_zeros_idx)
        idx = complex_zeros_idx(i);
        zero = sys_zeros(idx);
        fprintf('\n零点 %d (%.4f%+.4fj) 的角度分析：\n', idx, real(zero), imag(zero));
        
        % 与极点的夹角
        fprintf('  与极点的夹角：\n');
        for j = 1:length(sys_poles)
            pole = sys_poles(j);
            angle_diff = angle(zero - pole) * 180 / pi;
            distance = abs(zero - pole);
            if imag(pole) == 0
                fprintf('    到极点 %d (%.4f): 角度 = %.1f°, 距离 = %.4f\n', ...
                        j, real(pole), angle_diff, distance);
            else
                fprintf('    到极点 %d (%.4f%+.4fj): 角度 = %.1f°, 距离 = %.4f\n', ...
                        j, real(pole), imag(pole), angle_diff, distance);
            end
        end
        
        % 与其他零点的夹角
        fprintf('  与其他零点的夹角：\n');
        for j = 1:length(sys_zeros)
            if j ~= idx
                other_zero = sys_zeros(j);
                angle_diff = angle(zero - other_zero) * 180 / pi;
                distance = abs(zero - other_zero);
                if imag(other_zero) == 0
                    fprintf('    到零点 %d (%.4f): 角度 = %.1f°, 距离 = %.4f\n', ...
                            j, real(other_zero), angle_diff, distance);
                else
                    fprintf('    到零点 %d (%.4f%+.4fj): 角度 = %.1f°, 距离 = %.4f\n', ...
                            j, real(other_zero), imag(other_zero), angle_diff, distance);
                end
            end
        end
    end
end

% 创建新图形窗口
figure;

% 绘制根轨迹
rlocus(sys);
hold on;

% 添加网格和标签
grid on;
title('系统根轨迹图 (增强版)');
xlabel('实轴 (Re)');
ylabel('虚轴 (Im)');

% 添加阻尼比和自然频率网格
sgrid;

% 获取根轨迹数据
[r, k] = rlocus(sys);

%% 1. 自动计算临界增益（虚轴交点）
k_critical = [];
fprintf('\n=== 临界增益计算 ===\n');

% 构造特征方程 1 + K*G(s) = 0
% 即 den(s) + K*num(s) = 0
[num, den] = tfdata(sys, 'v');

% 寻找虚轴上的根
syms s K_sym
char_eq = poly2sym(den, s) + K_sym * poly2sym(num, s);

% 设 s = jw，分离实部和虚部
syms w real
s_jw = 1i * w;
char_eq_jw = subs(char_eq, s, s_jw);

% 分离实部和虚部
real_part = real(char_eq_jw);
imag_part = imag(char_eq_jw);

try
    % 求解虚部为0的频率
    w_solutions = solve(imag_part == 0, w, 'Real', true);
    w_solutions = double(w_solutions);
    w_solutions = w_solutions(w_solutions >= 0 & isreal(w_solutions)); % 只取非负实数解
    
    % 对每个频率求对应的K值
    for i = 1:length(w_solutions)
        w_val = w_solutions(i);
        K_val = -double(subs(real_part, [w, K_sym], [w_val, 0])) / double(subs(diff(real_part, K_sym), w, w_val));
        
        if K_val > 0 % 只考虑正增益
            k_critical = [k_critical, K_val];
            fprintf('虚轴交点: s = ±%.4fj, 临界增益: K = %.4f\n', w_val, K_val);
            
            % 在图上标记交点
            plot([0, 0], [w_val, -w_val], 'mo', 'MarkerSize', 10, 'LineWidth', 3);
            text(0.1, w_val, sprintf('K=%.2f', K_val), 'FontSize', 10, 'Color', 'magenta');
            if w_val > 0
                text(0.1, -w_val, sprintf('K=%.2f', K_val), 'FontSize', 10, 'Color', 'magenta');
            end
        end
    end
catch
    fprintf('无法自动计算临界增益，可能系统无虚轴交点\n');
end

if isempty(k_critical)
    fprintf('系统在正增益下无虚轴交点\n');
end

%% 2. 计算分离点
fprintf('\n=== 分离点计算 ===\n');
separation_points = [];
separation_gains = [];

% 方法1：使用特征方程的导数
% 分离点条件：dK/ds = 0，即 d/ds[1 + K*G(s)] = 0
% 这等价于：G'(s) + K*G'(s)*G(s) = 0
% 即：dG(s)/ds = 0 或者 G(s) = -1/K (这是根轨迹条件)

[num, den] = tfdata(sys, 'v');

% 计算 G(s) = num(s)/den(s) 的导数
% G'(s) = [num'(s)*den(s) - num(s)*den'(s)] / [den(s)]^2

% 计算多项式的导数
num_deriv = polyder(num);
den_deriv = polyder(den);

% 计算分子：num'(s)*den(s) - num(s)*den'(s)
% 使用conv进行多项式乘法
if ~isempty(num_deriv) && ~isempty(den)
    term1 = conv(num_deriv, den);
else
    term1 = 0;
end

if ~isempty(num) && ~isempty(den_deriv)
    term2 = conv(num, den_deriv);
else
    term2 = 0;
end

% 确保两个多项式长度相同以便相减
max_len = max(length(term1), length(term2));
if length(term1) < max_len
    term1 = [zeros(1, max_len - length(term1)), term1];
end
if length(term2) < max_len
    term2 = [zeros(1, max_len - length(term2)), term2];
end

numerator_deriv = term1 - term2;

% 分离点条件：numerator_deriv = 0
% 求解多项式方程的根
if length(numerator_deriv) > 1
    separation_candidates = roots(numerator_deriv);
    
    % 只保留实数解
    real_candidates = separation_candidates(abs(imag(separation_candidates)) < 1e-10);
    real_candidates = real(real_candidates);
    
    if ~isempty(real_candidates)
        fprintf('分离点候选点：\n');
        
        for i = 1:length(real_candidates)
            s_point = real_candidates(i);
            
            % 计算该点处的增益值
            G_value = polyval(num, s_point) / polyval(den, s_point);
            K_value = -1 / G_value;
            
            % 检查增益是否为正实数
            if isreal(K_value) && K_value > 0 && isfinite(K_value)
                separation_points = [separation_points, s_point];
                separation_gains = [separation_gains, K_value];
                
                fprintf('  s = %.4f, K = %.4f\n', s_point, K_value);
                
                % 在图上标记分离点
                plot(s_point, 0, 'bs', 'MarkerSize', 12, 'MarkerFaceColor', 'blue', 'LineWidth', 2);
                text(s_point, 0.3, sprintf('分离点\nK=%.2f', K_value), ...
                     'FontSize', 9, 'Color', 'blue', 'FontWeight', 'bold', ...
                     'HorizontalAlignment', 'center');
            end
        end
        
        if isempty(separation_points)
            fprintf('候选点中无有效分离点（要求K>0且为实数）\n');
        end
    else
        fprintf('无实数分离点\n');
    end
else
    fprintf('无法计算分离点（系统阶数过低）\n');
end

% 方法2：验证分离点（可选）
if ~isempty(separation_points)
    fprintf('\n分离点验证：\n');
    for i = 1:length(separation_points)
        s_point = separation_points(i);
        K_value = separation_gains(i);
        
        % 在分离点处，特征方程应该有重根
        % 计算特征方程：den(s) + K*num(s) = 0
        char_poly = den + K_value * num;
        char_roots = roots(char_poly);
        
        % 检查是否有重根（在分离点附近）
        distances = abs(char_roots - s_point);
        close_roots = char_roots(distances < 0.1);
        
        fprintf('  分离点 s=%.4f 处的特征根：', s_point);
        for j = 1:length(close_roots)
            if abs(imag(close_roots(j))) < 1e-10
                fprintf(' %.4f', real(close_roots(j)));
            else
                fprintf(' %.4f%+.4fj', real(close_roots(j)), imag(close_roots(j)));
            end
        end
        fprintf('\n');
    end
end

% 方法3：使用根轨迹实轴规则验证
fprintf('\n实轴上根轨迹区段分析：\n');
all_critical_points = [real(sys_poles)', real(sys_zeros)'];
all_critical_points = sort(all_critical_points);

if length(all_critical_points) > 1
    fprintf('实轴上的关键点：');
    fprintf(' %.4f', all_critical_points);
    fprintf('\n');
    
    % 检查每个区间
    for i = 1:length(all_critical_points)-1
        test_point = (all_critical_points(i) + all_critical_points(i+1)) / 2;
        
        % 计算该点右侧的极点和零点数量
        poles_right = sum(real(sys_poles) > test_point);
        zeros_right = sum(real(sys_zeros) > test_point);
        
        % 根轨迹存在条件：右侧极点数-零点数为奇数
        if mod(poles_right - zeros_right, 2) == 1
            fprintf('区间 [%.4f, %.4f] 上有根轨迹\n', all_critical_points(i), all_critical_points(i+1));
            
            % 检查该区间内是否有分离点
            interval_sep_points = separation_points(separation_points > all_critical_points(i) & ...
                                                   separation_points < all_critical_points(i+1));
            if ~isempty(interval_sep_points)
                fprintf('  该区间内的分离点：');
                fprintf(' %.4f', interval_sep_points);
                fprintf('\n');
            end
        end
    end
    
    % 检查无穷远区间
    if ~isempty(all_critical_points)
        % 左侧无穷远区间
        test_point = all_critical_points(1) - 1;
        poles_right = sum(real(sys_poles) > test_point);
        zeros_right = sum(real(sys_zeros) > test_point);
        if mod(poles_right - zeros_right, 2) == 1
            fprintf('区间 (-∞, %.4f] 上有根轨迹\n', all_critical_points(1));
        end
        
        % 右侧无穷远区间
        test_point = all_critical_points(end) + 1;
        poles_right = sum(real(sys_poles) > test_point);
        zeros_right = sum(real(sys_zeros) > test_point);
        if mod(poles_right - zeros_right, 2) == 1
            fprintf('区间 [%.4f, +∞) 上有根轨迹\n', all_critical_points(end));
        end
    end
end

%% 3. 绘制渐近线
fprintf('\n=== 渐近线计算 ===\n');

if n > m
    % 渐近线数量
    num_asymptotes = n - m;
    
    % 渐近线与实轴的交点（质心）
    sigma_a = (sum(sys_poles) - sum(sys_zeros)) / (n - m);
    
    % 渐近线角度
    asymptote_angles = zeros(1, num_asymptotes);
    for i = 1:num_asymptotes
        asymptote_angles(i) = (2*i - 1) * 180 / (n - m);
    end
    
    % 存储渐近线信息
    asymptotes.center = sigma_a;
    asymptotes.angles = asymptote_angles;
    asymptotes.count = num_asymptotes;
    
    fprintf('渐近线中心: σ = %.4f\n', sigma_a);
    fprintf('渐近线角度: ');
    fprintf('%.1f° ', asymptote_angles);
    fprintf('\n');
    
    % 绘制渐近线
    axis_limits = axis;
    line_length = max(abs(axis_limits)) * 2;
    
    for i = 1:num_asymptotes
        angle_rad = asymptote_angles(i) * pi / 180;
        x_end = sigma_a + line_length * cos(angle_rad);
        y_end = line_length * sin(angle_rad);
        
        % 绘制渐近线
        plot([sigma_a - line_length * cos(angle_rad), x_end], ...
             [-line_length * sin(angle_rad), y_end], ...
             'k--', 'LineWidth', 1.5);
    end
else
    asymptotes = struct('center', [], 'angles', [], 'count', 0);
    fprintf('n ≤ m，无渐近线\n');
end

%% 4. 计算出射角和入射角
fprintf('\n=== 出射角/入射角计算 ===\n');

angles = struct();
angles.departure = [];  % 出射角
angles.arrival = [];    % 入射角

% 计算复数极点的出射角
complex_poles = sys_poles(imag(sys_poles) ~= 0);
if ~isempty(complex_poles)
    fprintf('复数极点出射角:\n');
    angles.departure = zeros(1, length(complex_poles));
    
    for i = 1:length(complex_poles)
        pole = complex_poles(i);
        
        % 计算角度贡献
        angle_sum = 0;
        
        % 来自其他极点的角度贡献
        other_poles = sys_poles(sys_poles ~= pole);
        for j = 1:length(other_poles)
            angle_sum = angle_sum + angle(pole - other_poles(j));
        end
        
        % 来自零点的角度贡献
        for j = 1:length(sys_zeros)
            angle_sum = angle_sum - angle(pole - sys_zeros(j));
        end
        
        % 出射角 = 180° - (其他极点角度之和 - 零点角度之和)
        departure_angle = 180 - angle_sum * 180 / pi;
        
        % 确保角度在 [-180°, 180°] 范围内
        while departure_angle > 180
            departure_angle = departure_angle - 360;
        end
        while departure_angle <= -180
            departure_angle = departure_angle + 360;
        end
        
        angles.departure(i) = departure_angle;
        
        fprintf('极点 %.4f%+.4fj: 出射角 = %.1f°\n', ...
                real(pole), imag(pole), departure_angle);
        
        % 在图上绘制出射角箭头
        arrow_length = 0.5;
        arrow_angle = departure_angle * pi / 180;
        arrow_end_x = real(pole) + arrow_length * cos(arrow_angle);
        arrow_end_y = imag(pole) + arrow_length * sin(arrow_angle);
        
        quiver(real(pole), imag(pole), ...
               arrow_length * cos(arrow_angle), arrow_length * sin(arrow_angle), ...
               0, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.3);
        
        % 添加角度标注
        text(real(pole) + 0.1, imag(pole) + 0.1, ...
             sprintf('%.0f°', departure_angle), ...
             'FontSize', 9, 'Color', 'green', 'FontWeight', 'bold');
    end
end

% 计算复数零点的入射角
complex_zeros = sys_zeros(imag(sys_zeros) ~= 0);
if ~isempty(complex_zeros)
    fprintf('复数零点入射角:\n');
    angles.arrival = zeros(1, length(complex_zeros)); 
    
    for i = 1:length(complex_zeros)
        zero = complex_zeros(i);
        
        % 计算角度贡献
        angle_sum = 0;
        
        % 来自极点的角度贡献
        for j = 1:length(sys_poles)
            angle_sum = angle_sum + angle(zero - sys_poles(j));
        end
        
        % 来自其他零点的角度贡献
        other_zeros = sys_zeros(sys_zeros ~= zero);
        for j = 1:length(other_zeros)
            angle_sum = angle_sum - angle(zero - other_zeros(j));
        end
        
        % 入射角 = -(极点角度之和 - 其他零点角度之和)
        arrival_angle = -angle_sum * 180 / pi;
        
        % 确保角度在 [-180°, 180°] 范围内
        while arrival_angle > 180
            arrival_angle = arrival_angle - 360;
        end
        while arrival_angle <= -180
            arrival_angle = arrival_angle + 360;
        end
        
        angles.arrival(i) = arrival_angle;
        
        fprintf('零点 %.4f%+.4fj: 入射角 = %.1f°\n', ...
                real(zero), imag(zero), arrival_angle);
        
        % 在图上绘制入射角箭头
        arrow_length = 0.5;
        arrow_angle = arrival_angle * pi / 180;
        
        quiver(real(zero) - arrow_length * cos(arrow_angle), ...
               imag(zero) - arrow_length * sin(arrow_angle), ...
               arrow_length * cos(arrow_angle), arrow_length * sin(arrow_angle), ...
               0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.3);
        
        % 添加角度标注
        text(real(zero) + 0.1, imag(zero) + 0.1, ...
             sprintf('%.0f°', arrival_angle), ...
             'FontSize', 9, 'Color', 'red', 'FontWeight', 'bold');
    end
end

%% 添加极零点标记
h_poles = plot(real(sys_poles), imag(sys_poles), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
if ~isempty(sys_zeros)
    h_zeros = plot(real(sys_zeros), imag(sys_zeros), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
end

% 添加文本说明替代图例，避免与控制系统工具箱冲突
text_x = min(real(sys_poles)) - 2;
text_y = max(imag(sys_poles)) + 1;
text(text_x, text_y, '× 极点', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
if ~isempty(sys_zeros)
    text(text_x, text_y - 0.5, '○ 零点', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
end
if ~isempty(k_critical)
    text(text_x, text_y - 1, '● 虚轴交点', 'FontSize', 10, 'Color', 'magenta', 'FontWeight', 'bold');
end
if ~isempty(separation_points)
    text(text_x, text_y - 1.5, '■ 分离点', 'FontSize', 10, 'Color', 'blue', 'FontWeight', 'bold');
end
if n > m
    text(text_x, text_y - 2, '-- 渐近线', 'FontSize', 10, 'Color', 'black', 'FontWeight', 'bold');
end

% 更新标题显示临界增益和分离点信息
title_str = '根轨迹图';
if ~isempty(k_critical) && ~isempty(separation_points)
    if isscalar(k_critical) && isscalar(separation_points)
        title_str = sprintf('根轨迹图 [临界增益: K_c=%.3f, 分离点: s=%.3f]', k_critical(1), separation_points(1));
    else
        title_str = sprintf('根轨迹图 [临界增益: %d个, 分离点: %d个]', length(k_critical), length(separation_points));
    end
elseif ~isempty(k_critical)
    if isscalar(k_critical)
        title_str = sprintf('根轨迹图 [临界增益: K_c = %.4f]', k_critical(1));
    else
        title_str = sprintf('根轨迹图 [临界增益: K_c = %s]', mat2str(k_critical, 4));
    end
elseif ~isempty(separation_points)
    if isscalar(separation_points)
        title_str = sprintf('根轨迹图 [分离点: s = %.4f]', separation_points(1));
    else
        title_str = sprintf('根轨迹图 [分离点: %d个]', length(separation_points));
    end
end
title(title_str);

hold off;

%% 输出总结
fprintf('\n=== 分析总结 ===\n');
fprintf('系统阶数: %d (极点数: %d, 零点数: %d)\n', n, n, m);
if ~isempty(k_critical)
    fprintf('临界增益数量: %d\n', length(k_critical));
end
if ~isempty(separation_points)
    fprintf('分离点数量: %d\n', length(separation_points));
end
if n > m
    fprintf('渐近线数量: %d\n', num_asymptotes);
end
fprintf('复数极点数量: %d\n', length(complex_poles));
fprintf('复数零点数量: %d\n', length(complex_zeros));

end
