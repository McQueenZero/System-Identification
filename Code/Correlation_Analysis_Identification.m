%%-------------------------------------------------------------------------
% 作者：       赵敏琨
% 日期：       2021年6月
% 说明：       相关分析法辨识脉冲响应，进而辨识系统传递函数
% 版本：       MATLAB R2018a
% 辨识方法：
% 脉冲响应：   1)白噪声;  2)M序列;  3)逆M序列  
% 传递函数：   1)采用Hankel矩阵法;  2)采用脉冲响应的差分法.
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
%% 数据输入、数据选择、方法选择
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
while 1                             %手动选择方法
    gMethod = input('脉冲响应辨识方法(WN/M/RM)：', 's');     %白噪声/M序列/逆M序列
    switch gMethod
        case 'WN'
            break;
        case 'M'
            break;
        case 'RM'
            break
        otherwise
            disp('非法，请重新输入')
    end
end
while 1                             %手动选择方法
    TFMethod = input('传递函数辨识方法(HK/DE)：', 's');     %Hankel矩阵法/Difference Equation差分方程法
    switch TFMethod
        case 'HK'
            break;
        case 'DE'
            break;
        otherwise
            disp('非法，请重新输入')
    end
end
disp('请输入起始拍，第0拍对应m=1，建议：m = 2')
while 1
    m = input('起始拍：m=');  %手动输入起始拍(第0拍，m=1)
    if m < 1 || mod(m,1) ~= 0
        disp('非法起始拍，请重新输入')
    else
        break;
    end
end
% k=5, T0=0.2效果好

%% 频率响应分析
num = [b0(k) b1(k)];
den = [a0(k) a1(k) a2(k) a3(k)];
sys = tf(num, den);  %sys为实际的闭环传递函数
sys_ol = minreal(sys / (1 - sys));  %sys为实际的开环传递函数
[~, ~, ~, Wpm] = margin(sys_ol);  %截止频率
Wl = Wpm * 0.03; %最低工作频率
figure('Name','开环传函频率分析'), margin(sys_ol)
disp(['截止频率：', num2str(Wpm), ' rad/s'])
disp(['最低工频：', num2str(Wl), ' rad/s'])

%% 阶跃响应分析
SRinfo = stepinfo(sys);
Ts = SRinfo.SettlingTime;   %过渡时间
figure('Name','闭环传函阶跃响应分析'), step(sys)
disp(['过渡时间：', num2str(Ts), ' s'])

%% 采样周期（时钟脉冲周期）、寄存器个数输入
while 1
    T0_req = (3 * Wpm)^(-1);
    disp(['请输入采样周期，建议：T0 <= ', num2str(T0_req)])
    while 1
        T0 = input('采样周期：T0=');      %手动输入采样周期
        if T0 <= 0 || T0 > 3
            disp('非法采样周期，请重新输入')
        else
            del = T0;  %时钟脉冲周期
            break;
        end
    end
    N_req = (del * Wl)^(-1);
    disp(['请输入生成PRBS的寄存器个数，建议：序列周期N = 2^n-1 >= ', num2str(N_req)])
    while 1
        n = input('寄存器个数：n=');      %手动输入寄存器个数
        N = 2^n - 1;  %M序列的周期
        disp(['序列周期为：N=', num2str(N)])
        if n < 3 || n > 10 || mod(n,1) ~= 0
            disp('非法寄存器个数，请重新输入')
        else
            break;
        end
    end
    if (N - 1) * del > Ts
        disp(['过渡时间要求满足：(N - 1) * del=', num2str((N - 1) * del), ' > Ts=', num2str(Ts)])
        input('按回车继续')
        break
    else
        disp(['(N - 1) * del <= ', num2str(Ts)])
    end
end

%% 生成系统响应y(t)
a = 2;  %二位式伪随机序列PRBS的电平 -a为1，+a为0
total = 4 * N;

[WNseq, Mseq, RMseq] = genIN(n, a, total);

tfin = total * del;
tim = 0:del:tfin-del;
TSim = tfin;

figure('Name','输入激励信号')
xlabel('Time（seconds)'), ylabel('Amplitude')
switch gMethod
    case 'WN'
        stairs(WNseq), grid on
        title('白噪声')
        axis tight
        disp(['白噪声', ' 均值:', num2str(mean(WNseq)), ' 方差:', num2str(var(WNseq))])
        IN = WNseq;
    case 'M'
        stairs(Mseq), grid on
        axis([0 total -2.5 2.5])
        title('M序列')
        disp(['M序列', ' 均值:', num2str(mean(Mseq)), ' 方差:', num2str(var(Mseq))])
        IN = Mseq;
    case 'RM'
        stairs(RMseq),grid on
        axis([0 total -2.5 2.5])
        title('逆M序列')
        disp(['逆M序列', ' 均值:', num2str(mean(RMseq)), ' 方差:', num2str(var(RMseq))])
        IN = RMseq;
end
y = lsim(sys, IN, tim);

% 绘制系统输入输出
figure('Name','系统激励输入和输出曲线')
stairs(tim, IN, 'r')
axis tight
hold on
plot(tim, y, 'b')
grid on
title('系统激励输入和输出')
xlabel('Time（seconds)'), ylabel('Amplitude')
legend({'in', 'out'}, 'Location', 'best')
hold off

% 计算互相关函数Rxy(ii*del)
sumUY = 0.0;
Rxy = [];
iDel_vec = [];
for ii = 1:N
    mu = ii - 1;
    iDel_vec = [iDel_vec; mu * del];
    if gMethod == 'RM'
        for k = N+1:N+2*N
            sumUY = sumUY + IN(k) * y(k+mu);
        end
        Rxy_mu = (1 / (2*N)) * sumUY;   
    else      
        for k = N+1:N+N
            sumUY = sumUY + IN(k) * y(k+mu);
        end
        Rxy_mu = (1 / N) * sumUY;   
    end
    sumUY = 0.0;
    Rxy = [Rxy; mu Rxy_mu];
end

% 计算g_hat和g
Rxy_iDel = Rxy(:, 2);
if gMethod == 'RM' 
    C = 0;
else
    C = sum(Rxy_iDel);
end
g_hat = N / (a^2 * del * (N + 1)) * (Rxy_iDel + C);
if gMethod == 'WN'
    g_hat = g_hat - g_hat(end);  %白噪声需要减去稳态值
end
[g, gt] = impulse(sys, iDel_vec(1):del:iDel_vec(end));
Result = [Rxy g_hat g];

% disp('-------------------------------------------')
% disp('i     Rxy(iDel)       g_hat            g')
% disp('-------------------------------------------')
% disp(num2str(Result))
% disp('-------------------------------------------')

figure('Name','系统脉冲响应曲线')
plot(iDel_vec, g, 'r-'), hold on
plot(iDel_vec, g_hat, 'b-.'), grid on
title('Impulse Response')
legend({'g', 'ghat'}, 'Location', 'best')
xlabel('Time（seconds)'), ylabel('Amplitude')

%% Hankel矩阵算法 
if TFMethod == 'HK'  
    sysd = c2d(sys, T0, 'zoh');            %传递函数离散化
    disp('实际传递函数为：')
    sys
%     disp('系统极点为：')
%     pole(sys)
    disp('离散化实际传函为：')
    sysd
    H = [g_hat(m) g_hat(m+1) g_hat(m+2)
        g_hat(m+1) g_hat(m+2) g_hat(m+3)
        g_hat(m+2) g_hat(m+3) g_hat(m+4)];
    if det(H) == 0
        error('Hankel矩阵奇异')
    else
        A = H^(-1) * [-g_hat(m+3); -g_hat(m+4); -g_hat(m+5)];
        B = [1 0 0; A(3) 1 0; A(2) A(3) 1] * [g_hat(m); g_hat(m+1); g_hat(m+2)];
        numd = B'*T0;   %乘以T0补偿求脉冲传递函数G(z^(-1))=c*(zI-A)^(-1)*b过程中，由采样时间引起的误差
        dend = [1 A(3) A(2) A(1)];
        sysd_identi = tf(numd, dend, T0);   %创建1个采样时间为T0的离散传递函数
    end
    sys_identi = d2c(sysd_identi, 'zoh');  %sys_identi为辨识出的传递函数
    disp('辨识所得脉冲传函为：')
    sysd_identi
    disp('辨识所得连续传函为：')
    sys_identi
%     figure('Name','单位脉冲响应曲线')
%     impulse(sysd, 'r-', TSim), hold on
%     impulse(sysd_identi, 'b-.', TSim)
%     legend({'实际脉冲传递函数', '辨识脉冲传递函数'}, 'Location', 'best')
    figure('Name','单位阶跃响应曲线')
    step(sys, 'r-', TSim), hold on
    step(sys_identi, 'b-.', TSim)
    legend({'实际传递函数', '辨识传递函数'}, 'Location', 'best')   
end

%% 差分方程算法
if TFMethod == 'DE'
    disp('实际传递函数为：')
    sys
%     disp('系统极点为：')
%     pole(sys)
%     [g_hat, gt] = impulse(sys, 0:T0:TSim);
    % m为起始拍，第0拍，m=1
    A = [g_hat(m+1) g_hat(m+2) g_hat(m+3); g_hat(m+2) g_hat(m+3) g_hat(m+4); g_hat(m+3) g_hat(m+4) g_hat(m+5)];
    B = [-g_hat(m); -g_hat(m+1); -g_hat(m+2)];
    a = A^(-1)*B;                   %待定系数a
    p = [a(3) a(2) a(1) 1];
    x = roots(p);                   %求解x为特征方程的单根
    s = log(x)/T0;                  %求解s为传递函数的极点
    c = ([(x.^(m-1)).'; (x.^m).'; (x.^(m+1)).'])^(-1)*[g_hat(m); g_hat(m+1); g_hat(m+2)]; %部分分式分子
    H = [];
    [num1, den1] = residue(c, s, H); %将部分分式展开式转换
    num2 = real(num1);              %返回复数阵列每个元素的实部
    den2 = real(den1);
    sys_identi = tf(num2, den2);    %sys_identi为辨识出的传递函数
    disp('辨识所得连续传函为：')
    sys_identi
%     figure('Name','单位脉冲响应曲线')
%     impulse(sys, 'r-', TSim), hold on     %脉冲响应
%     plot(gt, g_hat, 'bo')
%     legend({'单位脉冲响应', '单位脉冲响应点'}, 'Location', 'best')
    figure('Name','单位阶跃响应曲线')
    step(sys, 'r-', TSim), hold on
    step(sys_identi, 'b-.', TSim)
    legend({'实际传递函数', '辨识传递函数'}, 'Location', 'best')
end

%% 误差分析
[h, ht] = step(sys, 0:T0:TSim+10);
[h_identi, ht_identi] = step(sys_identi, 0:T0:TSim+10);
eh = (h_identi - h) ./ h * 100;
figure('Name','误差曲线')
plot(ht, eh), grid on
axis tight
xlabel('Time（seconds)'), ylabel('Error Percent (%)')
switch gMethod
    case 'M'
        title('M序列激励下，辨识传函单位阶跃响应误差百分比')
    case 'RM'
        title('逆M序列激励下，辨识传函单位阶跃响应误差百分比')
    case 'WN'
        title('白噪声激励下，辨识传函单位阶跃响应误差百分比')
end

%% 用到的函数
function [WNseq, Mseq, RMseq] = genIN(n, a, total)  %产生PRBS信号
%参数是寄存器个数n，M序列的电平峰值a，需要的M序列长度total

WNseq_mean = 0.1;  %正态随机数法，初始白噪声均值
%有限长度白噪声，希望生成的白噪声均值比M序列小
while WNseq_mean > 0.05
    WNseq = randn(1, total) * a;  %产生白噪声序列
    WNseq_mean = mean(WNseq);   %更新白噪声均值
end

% x0 = 1;  %乘同余法
% A = 5^13; M = 10^36;
% for k = 1:total
%     x2 = A * x0;
%     x1 = mod(x2, M);
%     v1 = x1/M;
%     v(:, k) = v1;
%     x0 = x1;
% end
% WNseq = (v - 0.5)*2*a;

% x0 = 1;  %混合同余法
% A = 2^n + 1; M = 2^35;
% c = (0.5+sqrt(3)/6) * 2^35;
% for k = 1:total
%     x2 = A * x0 + c;
%     x1 = mod(x2, M);
%     v1 = x1/M;
%     v(:, k) = v1;
%     x0 = x1;
% end
% WNseq = (v - 0.5)*2*a;

Mseq = [];  % 初始空M序列
RMseq = [];  %初始空逆M序列
% 初始化n个寄存器
for ii = 1:n
    R(ii) = 1;
end

S = 1;  %方波初值
for ii = 1:total
    temp = R(1);
    
    switch n  %不同级数寄存器生成M序列的反馈寄存器不同
        case {3, 4, 6, 7}
            R(1) = xor(R(n-1), R(n));  %模二和（异或）
        case 5
            R(1) = xor(R(n-2), R(n));  %模二和（异或）
        case 8
            tempxor = xor(R(n-2), R(n));  %模二和（异或）
            tempxor = xor(R(n-3), tempxor);
            R(1) = xor(R(n-4), tempxor);
        case 9
            R(1) = xor(R(n-4), R(n));  %模二和（异或）
        case 10
            R(1) = xor(R(n-3), R(n));  %模二和（异或）
        otherwise
            error('不支持的寄存器级数')
    end
    
    jj = 2;
    while jj <= n  %寄存器移位
        temp1 = R(jj);
        R(jj) = temp;
        jj = jj + 1;
        temp = temp1;
    end
    if R(n) == 1
        Mseq(ii) = -a;
    else
        Mseq(ii) = a;
    end
    Q(ii) = xor(R(n), S);
    if Q(ii) == 1
        RMseq(ii) = -a;
    else
        RMseq(ii) = a;
    end
    S = not(S);  %产生方波
end
end
