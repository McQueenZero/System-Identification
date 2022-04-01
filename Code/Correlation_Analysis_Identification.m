%%-------------------------------------------------------------------------
% ���ߣ�       ������
% ���ڣ�       2021��6��
% ˵����       ��ط�������ʶ������Ӧ��������ʶϵͳ���ݺ���
% �汾��       MATLAB R2018a
% ��ʶ������
% ������Ӧ��   1)������;  2)M����;  3)��M����  
% ���ݺ�����   1)����Hankel����;  2)����������Ӧ�Ĳ�ַ�.
% Ҫ��   ��������T0�ֱ�Ϊ0.2��0.5��0.8��
% ���ݺ�����ʽΪ��
%              b0 * s + b1
%   ---------------------------------
%   a0 * s^3 + a1 * s^2 + a2 * s + a3
% ���ݼ����£�
%     b0   b1   a0   a1   a2   a3
% 1    1   10    1    3    2   10
% 2    1   15    1    5    4   15
% 3    1   20    1    3    5   20
% 4    1   25    1    6    8   25
% 5    1   30    1    7    9   30
% 6    1   40    1    5    3   40
%%-------------------------------------------------------------------------
%% �������롢����ѡ�񡢷���ѡ��
clc 
close all
clear
% д�����ݼ�
DataSet = [ 1   10    1    3    2   10
            1   15    1    5    4   15
            1   20    1    3    5   20
            1   25    1    6    8   25
            1   30    1    7    9   30
            1   40    1    5    3   40  ];
b0 = DataSet(:, 1); b1 = DataSet(:, 2);
a0 = DataSet(:, 3); a1 = DataSet(:, 4); a2 = DataSet(:, 5); a3 = DataSet(:, 6);

while 1
    k = input('ѡ�����ݱ�ţ�k=');    %�ֶ�����ѡ�����ݱ��
    if k < 1 || k > 6
        disp('�Ƿ���ţ�����������')
    else
        break;
    end
end
while 1                             %�ֶ�ѡ�񷽷�
    gMethod = input('������Ӧ��ʶ����(WN/M/RM)��', 's');     %������/M����/��M����
    switch gMethod
        case 'WN'
            break;
        case 'M'
            break;
        case 'RM'
            break
        otherwise
            disp('�Ƿ�������������')
    end
end
while 1                             %�ֶ�ѡ�񷽷�
    TFMethod = input('���ݺ�����ʶ����(HK/DE)��', 's');     %Hankel����/Difference Equation��ַ��̷�
    switch TFMethod
        case 'HK'
            break;
        case 'DE'
            break;
        otherwise
            disp('�Ƿ�������������')
    end
end
disp('��������ʼ�ģ���0�Ķ�Ӧm=1�����飺m = 2')
while 1
    m = input('��ʼ�ģ�m=');  %�ֶ�������ʼ��(��0�ģ�m=1)
    if m < 1 || mod(m,1) ~= 0
        disp('�Ƿ���ʼ�ģ�����������')
    else
        break;
    end
end
% k=5, T0=0.2Ч����

%% Ƶ����Ӧ����
num = [b0(k) b1(k)];
den = [a0(k) a1(k) a2(k) a3(k)];
sys = tf(num, den);  %sysΪʵ�ʵıջ����ݺ���
sys_ol = minreal(sys / (1 - sys));  %sysΪʵ�ʵĿ������ݺ���
[~, ~, ~, Wpm] = margin(sys_ol);  %��ֹƵ��
Wl = Wpm * 0.03; %��͹���Ƶ��
figure('Name','��������Ƶ�ʷ���'), margin(sys_ol)
disp(['��ֹƵ�ʣ�', num2str(Wpm), ' rad/s'])
disp(['��͹�Ƶ��', num2str(Wl), ' rad/s'])

%% ��Ծ��Ӧ����
SRinfo = stepinfo(sys);
Ts = SRinfo.SettlingTime;   %����ʱ��
figure('Name','�ջ�������Ծ��Ӧ����'), step(sys)
disp(['����ʱ�䣺', num2str(Ts), ' s'])

%% �������ڣ�ʱ���������ڣ����Ĵ�����������
while 1
    T0_req = (3 * Wpm)^(-1);
    disp(['������������ڣ����飺T0 <= ', num2str(T0_req)])
    while 1
        T0 = input('�������ڣ�T0=');      %�ֶ������������
        if T0 <= 0 || T0 > 3
            disp('�Ƿ��������ڣ�����������')
        else
            del = T0;  %ʱ����������
            break;
        end
    end
    N_req = (del * Wl)^(-1);
    disp(['����������PRBS�ļĴ������������飺��������N = 2^n-1 >= ', num2str(N_req)])
    while 1
        n = input('�Ĵ���������n=');      %�ֶ�����Ĵ�������
        N = 2^n - 1;  %M���е�����
        disp(['��������Ϊ��N=', num2str(N)])
        if n < 3 || n > 10 || mod(n,1) ~= 0
            disp('�Ƿ��Ĵ�������������������')
        else
            break;
        end
    end
    if (N - 1) * del > Ts
        disp(['����ʱ��Ҫ�����㣺(N - 1) * del=', num2str((N - 1) * del), ' > Ts=', num2str(Ts)])
        input('���س�����')
        break
    else
        disp(['(N - 1) * del <= ', num2str(Ts)])
    end
end

%% ����ϵͳ��Ӧy(t)
a = 2;  %��λʽα�������PRBS�ĵ�ƽ -aΪ1��+aΪ0
total = 4 * N;

[WNseq, Mseq, RMseq] = genIN(n, a, total);

tfin = total * del;
tim = 0:del:tfin-del;
TSim = tfin;

figure('Name','���뼤���ź�')
xlabel('Time��seconds)'), ylabel('Amplitude')
switch gMethod
    case 'WN'
        stairs(WNseq), grid on
        title('������')
        axis tight
        disp(['������', ' ��ֵ:', num2str(mean(WNseq)), ' ����:', num2str(var(WNseq))])
        IN = WNseq;
    case 'M'
        stairs(Mseq), grid on
        axis([0 total -2.5 2.5])
        title('M����')
        disp(['M����', ' ��ֵ:', num2str(mean(Mseq)), ' ����:', num2str(var(Mseq))])
        IN = Mseq;
    case 'RM'
        stairs(RMseq),grid on
        axis([0 total -2.5 2.5])
        title('��M����')
        disp(['��M����', ' ��ֵ:', num2str(mean(RMseq)), ' ����:', num2str(var(RMseq))])
        IN = RMseq;
end
y = lsim(sys, IN, tim);

% ����ϵͳ�������
figure('Name','ϵͳ����������������')
stairs(tim, IN, 'r')
axis tight
hold on
plot(tim, y, 'b')
grid on
title('ϵͳ������������')
xlabel('Time��seconds)'), ylabel('Amplitude')
legend({'in', 'out'}, 'Location', 'best')
hold off

% ���㻥��غ���Rxy(ii*del)
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

% ����g_hat��g
Rxy_iDel = Rxy(:, 2);
if gMethod == 'RM' 
    C = 0;
else
    C = sum(Rxy_iDel);
end
g_hat = N / (a^2 * del * (N + 1)) * (Rxy_iDel + C);
if gMethod == 'WN'
    g_hat = g_hat - g_hat(end);  %��������Ҫ��ȥ��ֵ̬
end
[g, gt] = impulse(sys, iDel_vec(1):del:iDel_vec(end));
Result = [Rxy g_hat g];

% disp('-------------------------------------------')
% disp('i     Rxy(iDel)       g_hat            g')
% disp('-------------------------------------------')
% disp(num2str(Result))
% disp('-------------------------------------------')

figure('Name','ϵͳ������Ӧ����')
plot(iDel_vec, g, 'r-'), hold on
plot(iDel_vec, g_hat, 'b-.'), grid on
title('Impulse Response')
legend({'g', 'ghat'}, 'Location', 'best')
xlabel('Time��seconds)'), ylabel('Amplitude')

%% Hankel�����㷨 
if TFMethod == 'HK'  
    sysd = c2d(sys, T0, 'zoh');            %���ݺ�����ɢ��
    disp('ʵ�ʴ��ݺ���Ϊ��')
    sys
%     disp('ϵͳ����Ϊ��')
%     pole(sys)
    disp('��ɢ��ʵ�ʴ���Ϊ��')
    sysd
    H = [g_hat(m) g_hat(m+1) g_hat(m+2)
        g_hat(m+1) g_hat(m+2) g_hat(m+3)
        g_hat(m+2) g_hat(m+3) g_hat(m+4)];
    if det(H) == 0
        error('Hankel��������')
    else
        A = H^(-1) * [-g_hat(m+3); -g_hat(m+4); -g_hat(m+5)];
        B = [1 0 0; A(3) 1 0; A(2) A(3) 1] * [g_hat(m); g_hat(m+1); g_hat(m+2)];
        numd = B'*T0;   %����T0���������崫�ݺ���G(z^(-1))=c*(zI-A)^(-1)*b�����У��ɲ���ʱ����������
        dend = [1 A(3) A(2) A(1)];
        sysd_identi = tf(numd, dend, T0);   %����1������ʱ��ΪT0����ɢ���ݺ���
    end
    sys_identi = d2c(sysd_identi, 'zoh');  %sys_identiΪ��ʶ���Ĵ��ݺ���
    disp('��ʶ�������崫��Ϊ��')
    sysd_identi
    disp('��ʶ������������Ϊ��')
    sys_identi
%     figure('Name','��λ������Ӧ����')
%     impulse(sysd, 'r-', TSim), hold on
%     impulse(sysd_identi, 'b-.', TSim)
%     legend({'ʵ�����崫�ݺ���', '��ʶ���崫�ݺ���'}, 'Location', 'best')
    figure('Name','��λ��Ծ��Ӧ����')
    step(sys, 'r-', TSim), hold on
    step(sys_identi, 'b-.', TSim)
    legend({'ʵ�ʴ��ݺ���', '��ʶ���ݺ���'}, 'Location', 'best')   
end

%% ��ַ����㷨
if TFMethod == 'DE'
    disp('ʵ�ʴ��ݺ���Ϊ��')
    sys
%     disp('ϵͳ����Ϊ��')
%     pole(sys)
%     [g_hat, gt] = impulse(sys, 0:T0:TSim);
    % mΪ��ʼ�ģ���0�ģ�m=1
    A = [g_hat(m+1) g_hat(m+2) g_hat(m+3); g_hat(m+2) g_hat(m+3) g_hat(m+4); g_hat(m+3) g_hat(m+4) g_hat(m+5)];
    B = [-g_hat(m); -g_hat(m+1); -g_hat(m+2)];
    a = A^(-1)*B;                   %����ϵ��a
    p = [a(3) a(2) a(1) 1];
    x = roots(p);                   %���xΪ�������̵ĵ���
    s = log(x)/T0;                  %���sΪ���ݺ����ļ���
    c = ([(x.^(m-1)).'; (x.^m).'; (x.^(m+1)).'])^(-1)*[g_hat(m); g_hat(m+1); g_hat(m+2)]; %���ַ�ʽ����
    H = [];
    [num1, den1] = residue(c, s, H); %�����ַ�ʽչ��ʽת��
    num2 = real(num1);              %���ظ�������ÿ��Ԫ�ص�ʵ��
    den2 = real(den1);
    sys_identi = tf(num2, den2);    %sys_identiΪ��ʶ���Ĵ��ݺ���
    disp('��ʶ������������Ϊ��')
    sys_identi
%     figure('Name','��λ������Ӧ����')
%     impulse(sys, 'r-', TSim), hold on     %������Ӧ
%     plot(gt, g_hat, 'bo')
%     legend({'��λ������Ӧ', '��λ������Ӧ��'}, 'Location', 'best')
    figure('Name','��λ��Ծ��Ӧ����')
    step(sys, 'r-', TSim), hold on
    step(sys_identi, 'b-.', TSim)
    legend({'ʵ�ʴ��ݺ���', '��ʶ���ݺ���'}, 'Location', 'best')
end

%% ������
[h, ht] = step(sys, 0:T0:TSim+10);
[h_identi, ht_identi] = step(sys_identi, 0:T0:TSim+10);
eh = (h_identi - h) ./ h * 100;
figure('Name','�������')
plot(ht, eh), grid on
axis tight
xlabel('Time��seconds)'), ylabel('Error Percent (%)')
switch gMethod
    case 'M'
        title('M���м����£���ʶ������λ��Ծ��Ӧ���ٷֱ�')
    case 'RM'
        title('��M���м����£���ʶ������λ��Ծ��Ӧ���ٷֱ�')
    case 'WN'
        title('�����������£���ʶ������λ��Ծ��Ӧ���ٷֱ�')
end

%% �õ��ĺ���
function [WNseq, Mseq, RMseq] = genIN(n, a, total)  %����PRBS�ź�
%�����ǼĴ�������n��M���еĵ�ƽ��ֵa����Ҫ��M���г���total

WNseq_mean = 0.1;  %��̬�����������ʼ��������ֵ
%���޳��Ȱ�������ϣ�����ɵİ�������ֵ��M����С
while WNseq_mean > 0.05
    WNseq = randn(1, total) * a;  %��������������
    WNseq_mean = mean(WNseq);   %���°�������ֵ
end

% x0 = 1;  %��ͬ�෨
% A = 5^13; M = 10^36;
% for k = 1:total
%     x2 = A * x0;
%     x1 = mod(x2, M);
%     v1 = x1/M;
%     v(:, k) = v1;
%     x0 = x1;
% end
% WNseq = (v - 0.5)*2*a;

% x0 = 1;  %���ͬ�෨
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

Mseq = [];  % ��ʼ��M����
RMseq = [];  %��ʼ����M����
% ��ʼ��n���Ĵ���
for ii = 1:n
    R(ii) = 1;
end

S = 1;  %������ֵ
for ii = 1:total
    temp = R(1);
    
    switch n  %��ͬ�����Ĵ�������M���еķ����Ĵ�����ͬ
        case {3, 4, 6, 7}
            R(1) = xor(R(n-1), R(n));  %ģ���ͣ����
        case 5
            R(1) = xor(R(n-2), R(n));  %ģ���ͣ����
        case 8
            tempxor = xor(R(n-2), R(n));  %ģ���ͣ����
            tempxor = xor(R(n-3), tempxor);
            R(1) = xor(R(n-4), tempxor);
        case 9
            R(1) = xor(R(n-4), R(n));  %ģ���ͣ����
        case 10
            R(1) = xor(R(n-3), R(n));  %ģ���ͣ����
        otherwise
            error('��֧�ֵļĴ�������')
    end
    
    jj = 2;
    while jj <= n  %�Ĵ�����λ
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
    S = not(S);  %��������
end
end
