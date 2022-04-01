%%-------------------------------------------------------------------------
% 作者：       赵敏琨
% 日期：       2021年6月
% 说明：       相关分析法有关编程练习
% 版本：       MATLAB R2018a
% 注意：       分小节运行
%%-------------------------------------------------------------------------
%% 乘同余法产生白噪声
% ①公式：x_i = A * x_{i-1}(mod(M))
% 含义：下一个随机数是上一个随机数乘A对M取余
% ②公式：{xi_i}=x_i/M
% 含义：{xi_i}即为(0, 1)分布的随机序列
clc 
close all
clear

x0 = 1; 
% A = 7^10; M = 1^10;
A = 5^13; M = 10^36;
% A = 5^17; M = 2^42;

N = 100;
for k = 1:N
    x2 = A * x0;
    x1 = mod(x2, M);
    v1 = x1/M;
    v(:, k) = v1;
    x0 = x1;
end

xi = (v - 0.5)*2;
stairs(xi)
xlabel('k'), ylabel('噪声幅值'), title('白噪声序列')
disp(['白噪声', ' 均值：', num2str(mean(xi)), ' 方差：', num2str(var(xi))])

%% 混合同余法产生白噪声
% ①公式：x_i = (A * x_{i-1} + c)(mod(M))
% 含义：下一个随机数是上一个随机数乘A对M取余
% ②公式：{xi_i}=x_i/M
% 含义：{xi_i}即为(0, 1)分布的随机序列
clc 
close all
clear

x0 = 1; 
n = 7;
A = 2^n + 1; M = 2^35;
c = (0.5+sqrt(3)/6) * 2^35;

N = 100;
for k = 1:N
    x2 = A * x0 + c;
    x1 = mod(x2, M);
    v1 = x1/M;
    v(:, k) = v1;
    x0 = x1;
end

xi = (v - 0.5)*2;
stairs(xi)
xlabel('k'), ylabel('噪声幅值'), title('白噪声序列')
disp(['白噪声', ' 均值：', num2str(mean(xi)), ' 方差：', num2str(var(xi))])

%% 产生白噪声和有色噪声
clc 
close all
clear

L = 500;    %仿真长度
d = [1 -1.5 0.7 0.1];  c = [1 0.5 0.2];  % D、C多项式的系数（可用roots命令求其根）
nd = length(d) - 1; nc = length(c) - 1;  %nd、nc为d、c的阶次
xik = zeros(nc, 1);  %白噪声初值
ek = zeros(nd, 1);  %有色噪声初值
xi = randn(L, 1);  %randn产生均值为0，方差为1的高斯随机序列（白噪声序列）

for k = 1:L
    e(k) = -d(2:nd+1) * ek + c * [xi(k);xik];  %产生有色噪声
    
    %数据更新
    for i = nd:-1:2
        ek(i) = ek(i-1);
    end
    ek(1) = e(k);
    
    for i = nc:-1:2
        xik(i) = xik(i-1);
    end
    xik(1) = xi(k);
end
subplot(2,1,1)
stairs(xi)
xlabel('k'), ylabel('噪声幅值'), title('白噪声序列')
subplot(2,1,2)
stairs(e)
xlabel('k'), ylabel('噪声幅值'), title('有色噪声序列')

disp(['白噪声', ' 均值：', num2str(mean(xi)), ' 方差：', num2str(var(xi))])
disp(['有色噪声', ' 均值：', num2str(mean(e)), ' 方差：', num2str(var(e))])

%% 产生M序列和逆M序列
clc 
close all
clear

L = 60;  %M序列长度
x1 = 1; x2 = 1; x3 = 1; x4 = 0;  %移位寄存器初值
S = 1;  %方波初值

for k = 1:L
    IM = xor(S, x4);  %进行异或运算，产生逆M序列
    if IM == 0
        u(k) = -1;
    else
        u(k) = 1;
    end
    S = not(S);  %产生方波
    M(k) = xor(x3, x4);  %进行异或运算，产生M序列
    x4 = x3; x3 = x2; x2 = x1; x1 = M(k);  %寄存器移位
end
subplot(2,1,1)
stairs(M), grid on
axis([0 L/2 -0.5 1.5])
xlabel('k'), ylabel('M序列幅值'), title('M序列')
subplot(2,1,2)
stairs(u), grid on
axis([0 L -1.5 1.5])
xlabel('k'), ylabel('逆M序列幅值'), title('逆M序列')

disp(['M序列', ' 均值：', num2str(mean(M)), ' 方差：', num2str(var(M))])
disp(['逆M序列', ' 均值：', num2str(mean(u)), ' 方差：', num2str(var(u))])


%% 相关分析法综合练习
clc 
close all
clear

R =100e3;  %100千欧姆   
C = 1e-6;  %1微法
tc = R*C;  %时常数

% 产生M序列
n = 5;
a = 2;  %二位式伪随机序列PRBS的电平 -a为1，+a为0
del = 15e-3;  %时钟脉冲周期
N = 2^n-1;  %M序列的周期
total = 3 * N;

% 用函数genPRBS产生M序列
Out = genPRBS(n, a, del, total);
% 生成系统响应y(t)
s = tf('s');
G = 1/(tc*s+1)
tfin = total * del;
tim = 0:del:tfin-del;
y = lsim(G, Out, tim);

% 绘制系统输入输出
figure
stairs(tim, Out)
axis([0 1.0 -2.5 2.5])
hold on
plot(tim, y, 'r')
legend({'in', 'out'})
hold off

% 计算互相关函数Rxy(ii*del)
summ = 0.0;
Rxy = [];
iDel_vec = [];
for ii = 1:N
    tau = ii - 1;
    iDel_vec = [iDel_vec; tau * del];
    for jj = N+1:N+N
        summ = summ + Out(jj) * y(jj+tau);
    end
    Rxy_i = (1 / N) * summ;
    summ = 0.0;
    Rxy = [Rxy; tau Rxy_i];
end

% 计算g_hat和g
Rxy_iDel = Rxy(:, 2);
C = sum(Rxy_iDel);
% C = -Rxy(length(Rxy), 2);
S = (N+1) * a^2 * del / N;
g_hat = (Rxy_iDel + C) / S;
% g_hat(1) = 2 * g_hat(1);
g_hat = g_hat - g_hat(end);
g = 10 * exp(-10.*iDel_vec);
Result = [Rxy g_hat g];

disp('-------------------------------------------')
disp('i     Rxy(iDel)       ghat            g')
disp('-------------------------------------------')
disp(num2str(Result))
disp('-------------------------------------------')

figure
plot(iDel_vec, g), hold on
plot(iDel_vec, g_hat), grid on
legend({'g', 'ghat'})


%% 用到的函数
function Out = genPRBS(n, a, del, total)  %产生PRBS信号
%参数是寄存器个数n，M序列的电平峰值a，时钟脉冲周期del，需要的M序列长度total
Out = [];  % 初始空列表
% 初始化n个寄存器
for ii = 1:n
    R(ii) = 1;
end
if R(n) == 1
    Out(1) = -a;
elseif R(n) == 0
    Out(1) = a;
end

for ii = 2:total
    tempL = R(1);
    R(1) = xor(R(n-2), R(n));  %模二和（异或）
    for jj = 2:n  %寄存器移位
        tempR = R(jj);
        R(jj) = tempL;
        tempL = tempR;
    end
    if R(n) == 1
        Out(ii) = -a;
    elseif R(n) == 0
        Out(ii) = a;
    end
end
end
