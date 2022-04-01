%%-------------------------------------------------------------------------
% ���ߣ�       ������
% ���ڣ�       2021��5��
% ˵����       ��ʶϵͳ���ݺ���
% �汾��       MATLAB R2018a
% ������   1)����Hankel����;  2)����������Ӧ�Ĳ�ַ�.
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
%% ����
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
while 1
    T0 = input('�������ڣ�T0=');      %�ֶ������������
    if T0 < 0 || T0 > 3
        disp('�Ƿ��������ڣ�����������')
    else
        break;
    end
end
while 1
    TSim = input('����ʱ�䣺TSim=');  %�ֶ��������ʱ��
    if TSim < 10 || TSim > 100
        disp('�Ƿ�����ʱ�䣬����������')
    else
        break;
    end
end
while 1                             %�ֶ�ѡ�񷽷�
    Method = input('������Ӧ��ʶ������HK/DE����', 's');     %Hankel����/Difference Equation��ַ��̷�
    switch Method
        case 'HK'
            break;
        case 'DE'
            while 1
                m = input('��ʼ��(��0��,m=1)��m=');  %�ֶ�������ʼ��(��0�ģ�m=1)
                if m < 1 || mod(m,1) ~= 0
                    disp('�Ƿ���ʼ�ģ�����������')
                else
                    break;
                end
            end
            break;
        otherwise
            disp('�Ƿ�������������')
    end
end
% k=5, T0=0.2Ч����

%% Hankel�����㷨 
if Method == 'HK'  
    num = [b0(k) b1(k)];
    den = [a0(k) a1(k) a2(k) a3(k)];
    sys = tf(num, den);             %sysΪʵ�ʵĴ��ݺ���
    sysd = c2d(sys, T0, 'zoh');            %���ݺ�����ɢ��
    disp('ʵ�ʴ��ݺ���Ϊ��')
    sys
%     disp('ϵͳ����Ϊ��')
%     pole(sys)
%     disp('��ɢ��ʵ�ʴ���Ϊ��')
%     sysd
    [g, gt] = impulse(sysd);
    H = [g(1+1) g(2+1) g(3+1)
        g(2+1) g(3+1) g(4+1)
        g(3+1) g(4+1) g(5+1)];
    if det(H) == 0
        error('Hankel��������')
    else
        A = H^(-1) * [-g(4+1); -g(5+1); -g(6+1)];
        B = [1 0 0; A(3) 1 0; A(2) A(3) 1] * [g(1+1); g(2+1); g(3+1)];
        numd = B'*T0;   %����T0���������崫�ݺ���G(z^(-1))=c*(zI-A)^(-1)*b�����У��ɲ���ʱ����������
        dend = [1 A(3) A(2) A(1)];
        sysd_identi = tf(numd, dend, T0);   %����1������ʱ��ΪT0����ɢ���ݺ���
    end
    sys_identi = d2c(sysd_identi, 'zoh');  %sys_identiΪ��ʶ���Ĵ��ݺ���
%     disp('��ʶ������ɢ����Ϊ��')
%     sysd_identi
    disp('��ʶ������������Ϊ��')
    sys_identi
    figure('Name','��λ������Ӧ����')
    impulse(sysd, 'r-', TSim), hold on
    impulse(sysd_identi, 'b-.', TSim)
    legend({'ʵ�����崫�ݺ���', '��ʶ���崫�ݺ���'}, 'Location', 'best')
    figure('Name','��λ��Ծ��Ӧ����')
    step(sys, 'r-', TSim), hold on
    step(sys_identi, 'b-.', TSim)
    legend({'ʵ�ʴ��ݺ���', '��ʶ���ݺ���'}, 'Location', 'best')
    
    % ׼���ж�������
    [h, ht] = step(sys, TSim+10);
    [h_identi, ht_identi] = step(sys_identi, TSim+10);
    figure('Name','�������')
    plot(ht, h_identi - h, 'g-'), grid on
    axis tight
    xlabel('Time��seconds)'), ylabel('Amplitude')
    eval(['title(''����ʱ�� T0=' num2str(T0) 's ʱ�ı�ʶ�������'')'])    
end

%% ��ַ����㷨
if Method == 'DE'
    num = [b0(k) b1(k)];
    den = [a0(k) a1(k) a2(k) a3(k)];
    sys = tf(num, den);             %sysΪʵ�ʵĴ��ݺ���
    disp('ʵ�ʴ��ݺ���Ϊ��')
    sys
%     disp('ϵͳ����Ϊ��')
%     pole(sys)
    [g, gt] = impulse(sys, 0:T0:TSim);
    % mΪ��ʼ�ģ���0�ģ�m=1
    A = [g(m+1) g(m+2) g(m+3); g(m+2) g(m+3) g(m+4); g(m+3) g(m+4) g(m+5)];
    B = [-g(m); -g(m+1); -g(m+2)];
    a = A^(-1)*B;                   %����ϵ��a
    p = [a(3) a(2) a(1) 1];
    x = roots(p);                   %���xΪ�������̵ĵ���
    s = log(x)/T0;                  %���sΪ���ݺ����ļ���
    c = ([(x.^(m-1)).'; (x.^m).'; (x.^(m+1)).'])^(-1)*[g(m); g(m+1); g(m+2)]; %���ַ�ʽ����
    H = [];
    [num1, den1] = residue(c, s, H); %�����ַ�ʽչ��ʽת��
    num2 = real(num1);              %���ظ�������ÿ��Ԫ�ص�ʵ��
    den2 = real(den1);
    sys_identi = tf(num2, den2);    %sys_identiΪ��ʶ���Ĵ��ݺ���
    disp('��ʶ������������Ϊ��')
    sys_identi
    figure('Name','��λ������Ӧ����')
    impulse(sys, 'r-', TSim), hold on     %������Ӧ
    plot(gt, g, 'bo')
    legend({'��λ������Ӧ', '��λ������Ӧ��'}, 'Location', 'best')
    figure('Name','��λ��Ծ��Ӧ����')
    step(sys, 'r-', TSim), hold on
    step(sys_identi, 'b-.', TSim)
    legend({'ʵ�ʴ��ݺ���', '��ʶ���ݺ���'}, 'Location', 'best')
     
    % ׼���ж�������
    [h, ht] = step(sys, TSim+10);
    [h_identi, ht_identi] = step(sys_identi, TSim+10);
    figure('Name','�������')
    plot(ht, h_identi - h, 'g-'), grid on
    axis tight
    xlabel('Time��seconds)'), ylabel('Amplitude')
    eval(['title(''����ʱ�� T0=' num2str(T0) 's ʱ�ı�ʶ�������'')'])    
end