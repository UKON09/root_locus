% 清除工作空间
clc;
clear;
close all;


%% 定义系统模型
%示例：开环传递函数
%      sys(s) =             10 
%                ------------------------
%                 s(s + 20)(s² - 2s + 2)


%% ❗❗定义系统增益，不然临界增益会少kin倍

kin = 10;

%% 方法一：已知零点、极点和增益

sys_zeros = [];
sys_poles = [0 -20 -2-4i -2+4i];
sys_gain = kin;
G = zpk(sys_zeros,sys_poles,sys_gain);

%% 方法二：只知道开环传递函数多项式

%分母多项式为s⁴ + 24s³ + 100s² + 400s
% num = kin;  % 分子多项式系数
% den = [1 24 100 400 0];  % 分母多项式系数
% G = tf(num, den);

%% 方法三：直接写
% s = tf('s');
% G = 1 / (s * (s + 20) * (s^2 + 4*s + 20));

%% ====== 绘制根轨迹（右半平面不能出现零极点） ======
[r, k, k_crit, asymp, ang] = plotRootLocus(G, kin);
