%%-------------------------------------------------------------------------
% ���ߣ�       ������
% ���ڣ�       2021��6��
% ˵����       ��ط������йر����ϰ
% �汾��       MATLAB R2018a
% ע�⣺       ��С������
%%-------------------------------------------------------------------------
%% ��ͬ�෨����������
% �ٹ�ʽ��x_i = A * x_{i-1}(mod(M))
% ���壺��һ�����������һ���������A��Mȡ��
% �ڹ�ʽ��{xi_i}=x_i/M
% ���壺{xi_i}��Ϊ(0, 1)�ֲ����������
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
xlabel('k'), ylabel('������ֵ'), title('����������')
disp(['������', ' ��ֵ��', num2str(mean(xi)), ' ���', num2str(var(xi))])

%% ���ͬ�෨����������
% �ٹ�ʽ��x_i = (A * x_{i-1} + c)(mod(M))
% ���壺��һ�����������һ���������A��Mȡ��
% �ڹ�ʽ��{xi_i}=x_i/M
% ���壺{xi_i}��Ϊ(0, 1)�ֲ����������
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
xlabel('k'), ylabel('������ֵ'), title('����������')
disp(['������', ' ��ֵ��', num2str(mean(xi)), ' ���', num2str(var(xi))])

%% ��������������ɫ����
clc 
close all
clear

L = 500;    %���泤��
d = [1 -1.5 0.7 0.1];  c = [1 0.5 0.2];  % D��C����ʽ��ϵ��������roots�����������
nd = length(d) - 1; nc = length(c) - 1;  %nd��ncΪd��c�Ľ״�
xik = zeros(nc, 1);  %��������ֵ
ek = zeros(nd, 1);  %��ɫ������ֵ
xi = randn(L, 1);  %randn������ֵΪ0������Ϊ1�ĸ�˹������У����������У�

for k = 1:L
    e(k) = -d(2:nd+1) * ek + c * [xi(k);xik];  %������ɫ����
    
    %���ݸ���
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
xlabel('k'), ylabel('������ֵ'), title('����������')
subplot(2,1,2)
stairs(e)
xlabel('k'), ylabel('������ֵ'), title('��ɫ��������')

disp(['������', ' ��ֵ��', num2str(mean(xi)), ' ���', num2str(var(xi))])
disp(['��ɫ����', ' ��ֵ��', num2str(mean(e)), ' ���', num2str(var(e))])

%% ����M���к���M����
clc 
close all
clear

L = 60;  %M���г���
x1 = 1; x2 = 1; x3 = 1; x4 = 0;  %��λ�Ĵ�����ֵ
S = 1;  %������ֵ

for k = 1:L
    IM = xor(S, x4);  %����������㣬������M����
    if IM == 0
        u(k) = -1;
    else
        u(k) = 1;
    end
    S = not(S);  %��������
    M(k) = xor(x3, x4);  %����������㣬����M����
    x4 = x3; x3 = x2; x2 = x1; x1 = M(k);  %�Ĵ�����λ
end
subplot(2,1,1)
stairs(M), grid on
axis([0 L/2 -0.5 1.5])
xlabel('k'), ylabel('M���з�ֵ'), title('M����')
subplot(2,1,2)
stairs(u), grid on
axis([0 L -1.5 1.5])
xlabel('k'), ylabel('��M���з�ֵ'), title('��M����')

disp(['M����', ' ��ֵ��', num2str(mean(M)), ' ���', num2str(var(M))])
disp(['��M����', ' ��ֵ��', num2str(mean(u)), ' ���', num2str(var(u))])


%% ��ط������ۺ���ϰ
clc 
close all
clear

R =100e3;  %100ǧŷķ   
C = 1e-6;  %1΢��
tc = R*C;  %ʱ����

% ����M����
n = 5;
a = 2;  %��λʽα�������PRBS�ĵ�ƽ -aΪ1��+aΪ0
del = 15e-3;  %ʱ����������
N = 2^n-1;  %M���е�����
total = 3 * N;

% �ú���genPRBS����M����
Out = genPRBS(n, a, del, total);
% ����ϵͳ��Ӧy(t)
s = tf('s');
G = 1/(tc*s+1)
tfin = total * del;
tim = 0:del:tfin-del;
y = lsim(G, Out, tim);

% ����ϵͳ�������
figure
stairs(tim, Out)
axis([0 1.0 -2.5 2.5])
hold on
plot(tim, y, 'r')
legend({'in', 'out'})
hold off

% ���㻥��غ���Rxy(ii*del)
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

% ����g_hat��g
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


%% �õ��ĺ���
function Out = genPRBS(n, a, del, total)  %����PRBS�ź�
%�����ǼĴ�������n��M���еĵ�ƽ��ֵa��ʱ����������del����Ҫ��M���г���total
Out = [];  % ��ʼ���б�
% ��ʼ��n���Ĵ���
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
    R(1) = xor(R(n-2), R(n));  %ģ���ͣ����
    for jj = 2:n  %�Ĵ�����λ
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
