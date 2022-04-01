%%-------------------------------------------------------------------------
% 作者：       赵敏琨
% 日期：       2021年5月
% 说明：       辨识系统传递函数
% 版本：       MATLAB R2018a
% 方法：   1)采用Hankel矩阵法;  2)采用脉冲响应的差分法.
% 要求：   采样周期T0分别为0.2、0.5、0.8秒
% 传递函数形式为：
%              b0 * s + b1
%   ---------------------------------
%   a0 * s^3 + a1 * s^2 + a2 * s + a3
% 数据集如下：
%     b0   b1   a0   a1   a2   a3
% 1    1   10    1    3    2   10
% 2    1   15    1    5    4   15
% 3    1   20    1    3    5   20
% 4    1   25    1    6    8   25
% 5    1   30    1    7    9   30
% 6    1   40    1    5    3   40
%%-------------------------------------------------------------------------
%% 输入
clc 
close all
clear
% 写入数据集
DataSet = [ 1   10    1    3    2   10
            1   15    1    5    4   15
            1   20    1    3    5   20
            1   25    1    6    8   25
            1   30    1    7    9   30
            1   40    1    5    3   40  ];
b0 = DataSet(:, 1); b1 = DataSet(:, 2);
a0 = DataSet(:, 3); a1 = DataSet(:, 4); a2 = DataSet(:, 5); a3 = DataSet(:, 6);

while 1
    k = input('选择数据编号：k=');    %手动输入选择数据编号
    if k < 1 || k > 6
        disp('非法编号，请重新输入')
    else
        break;
    end
end
while 1
    T0 = input('采样周期：T0=');      %手动输入采样周期
    if T0 < 0 || T0 > 3
        disp('非法采样周期，请重新输入')
    else
        break;
    end
end
while 1
    TSim = input('仿真时间：TSim=');  %手动输入仿真时间
    if TSim < 10 || TSim > 100
        disp('非法仿真时间，请重新输入')
    else
        break;
    end
end
while 1                             %手动选择方法
    Method = input('脉冲响应辨识方法（HK/DE）：', 's');     %Hankel矩阵法/Difference Equation差分方程法
    switch Method
        case 'HK'
            break;
        case 'DE'
            while 1
                m = input('起始拍(第0拍,m=1)：m=');  %手动输入起始拍(第0拍，m=1)
                if m < 1 || mod(m,1) ~= 0
                    disp('非法起始拍，请重新输入')
                else
                    break;
                end
            end
            break;
        otherwise
            disp('非法，请重新输入')
    end
end
% k=5, T0=0.2效果好

%% Hankel矩阵算法 
if Method == 'HK'  
    num = [b0(k) b1(k)];
    den = [a0(k) a1(k) a2(k) a3(k)];
    sys = tf(num, den);             %sys为实际的传递函数
    sysd = c2d(sys, T0, 'zoh');            %传递函数离散化
    disp('实际传递函数为：')
    sys
%     disp('系统极点为：')
%     pole(sys)
%     disp('离散化实际传函为：')
%     sysd
    [g, gt] = impulse(sysd);
    H = [g(1+1) g(2+1) g(3+1)
        g(2+1) g(3+1) g(4+1)
        g(3+1) g(4+1) g(5+1)];
    if det(H) == 0
        error('Hankel矩阵奇异')
    else
        A = H^(-1) * [-g(4+1); -g(5+1); -g(6+1)];
        B = [1 0 0; A(3) 1 0; A(2) A(3) 1] * [g(1+1); g(2+1); g(3+1)];
        numd = B'*T0;   %乘以T0补偿求脉冲传递函数G(z^(-1))=c*(zI-A)^(-1)*b过程中，由采样时间引起的误差
        dend = [1 A(3) A(2) A(1)];
        sysd_identi = tf(numd, dend, T0);   %创建1个采样时间为T0的离散传递函数
    end
    sys_identi = d2c(sysd_identi, 'zoh');  %sys_identi为辨识出的传递函数
%     disp('辨识所得离散传函为：')
%     sysd_identi
    disp('辨识所得连续传函为：')
    sys_identi
    figure('Name','单位脉冲响应曲线')
    impulse(sysd, 'r-', TSim), hold on
    impulse(sysd_identi, 'b-.', TSim)
    legend({'实际脉冲传递函数', '辨识脉冲传递函数'}, 'Location', 'best')
    figure('Name','单位阶跃响应曲线')
    step(sys, 'r-', TSim), hold on
    step(sys_identi, 'b-.', TSim)
    legend({'实际传递函数', '辨识传递函数'}, 'Location', 'best')
    
    % 准则判断误差分析
    [h, ht] = step(sys, TSim+10);
    [h_identi, ht_identi] = step(sys_identi, TSim+10);
    figure('Name','误差曲线')
    plot(ht, h_identi - h, 'g-'), grid on
    axis tight
    xlabel('Time（seconds)'), ylabel('Amplitude')
    eval(['title(''采样时间 T0=' num2str(T0) 's 时的辨识误差曲线'')'])    
end

%% 差分方程算法
if Method == 'DE'
    num = [b0(k) b1(k)];
    den = [a0(k) a1(k) a2(k) a3(k)];
    sys = tf(num, den);             %sys为实际的传递函数
    disp('实际传递函数为：')
    sys
%     disp('系统极点为：')
%     pole(sys)
    [g, gt] = impulse(sys, 0:T0:TSim);
    % m为起始拍，第0拍，m=1
    A = [g(m+1) g(m+2) g(m+3); g(m+2) g(m+3) g(m+4); g(m+3) g(m+4) g(m+5)];
    B = [-g(m); -g(m+1); -g(m+2)];
    a = A^(-1)*B;                   %待定系数a
    p = [a(3) a(2) a(1) 1];
    x = roots(p);                   %求解x为特征方程的单根
    s = log(x)/T0;                  %求解s为传递函数的极点
    c = ([(x.^(m-1)).'; (x.^m).'; (x.^(m+1)).'])^(-1)*[g(m); g(m+1); g(m+2)]; %部分分式分子
    H = [];
    [num1, den1] = residue(c, s, H); %将部分分式展开式转换
    num2 = real(num1);              %返回复数阵列每个元素的实部
    den2 = real(den1);
    sys_identi = tf(num2, den2);    %sys_identi为辨识出的传递函数
    disp('辨识所得连续传函为：')
    sys_identi
    figure('Name','单位脉冲响应曲线')
    impulse(sys, 'r-', TSim), hold on     %脉冲响应
    plot(gt, g, 'bo')
    legend({'单位脉冲响应', '单位脉冲响应点'}, 'Location', 'best')
    figure('Name','单位阶跃响应曲线')
    step(sys, 'r-', TSim), hold on
    step(sys_identi, 'b-.', TSim)
    legend({'实际传递函数', '辨识传递函数'}, 'Location', 'best')
     
    % 准则判断误差分析
    [h, ht] = step(sys, TSim+10);
    [h_identi, ht_identi] = step(sys_identi, TSim+10);
    figure('Name','误差曲线')
    plot(ht, h_identi - h, 'g-'), grid on
    axis tight
    xlabel('Time（seconds)'), ylabel('Amplitude')
    eval(['title(''采样时间 T0=' num2str(T0) 's 时的辨识误差曲线'')'])    
end