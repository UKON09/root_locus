% 清除工作空间
clc;
clear;
close all;


% 定义系统模型
%示例：开环传递函数
%      sys(s) =             1 
%                ------------------------
%                s(s + 20)(s^2  - 2s + 2)
% 方法一：已知零点、极点和增益
sys_zeros = [];
sys_poles = [0 -20 -2-4i -2+4i];
sys_gain = 1;
G = zpk(sys_zeros,sys_poles,sys_gain);
% 方法二：只知道开环传递函数
% s = tf('s');
% G = 1 / (s * (s + 20) * (s^2 + 4*s + 20));


% 调用根轨迹函数
[r, k, k_crit, asymp, ang] = plotRootLocus(G);
