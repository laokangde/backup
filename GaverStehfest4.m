clear all
clc
close all
% G-S变换的四个算例以及测试双边的laplace变换

%%
%%
% %G-S变换的四个算例
% N = 12;
% m = 1:N;
% K = [
%     -1.66666666666667e-02  1.60166666666667e+01 ...
%     -1.24700000000000e+03 2.75543333333333e+04 ...
%     -2.63280833333333e+05 1.32413870000000e+06 ...
%     -3.89170553333333e+06 7.05328633333333e+06 ...
%     -8.00533650000000e+06 5.55283050000000e+06 ...
%     -2.15550720000000e+06 3.59251200000000e+05
%     ];
% u = 4.0*pi*1e-7;
% sigma = 0.01;
% rho = 100;
% alpha = sqrt( u*sigma )*rho;
% t = 0.1:0.1:100;
% s = m' * (log(2)./t);
% 
% F1 = 1./s;
% f1 =K * F1 ./ t *log(2);
% f1_t = ones(1,length(t));
% err1=abs(f1-f1_t);
% fprintf('F1(s)\n')
% fprintf('最大相对误差：%5.4e \n',max(abs(err1)));
% fprintf('均方根误差：   %5.4e \n',sqrt(sum(abs(err1).^2)/length(t)));
% 
% F2 = 1./(s.^ 2);
% f2 = K * F2 ./ t *log(2);
% f2_t = t;
% err2=abs(f2-f2_t);
% fprintf('F2(s)\n')
% fprintf('最大相对误差：%5.4e \n',max(err2));
% fprintf('均方根误差：   %5.4e \n',sqrt(sum(err2.^2)/length(t)));
% 
% F3 = exp( -alpha .* sqrt( s ) ) ./ s;
% f3 = K * F3 ./ t *log(2);
% f3_t = erfc(0.5 * alpha ./ sqrt(t));
% err3=abs(f3-f3_t);
% fprintf('F3(s)\n')
% fprintf('最大相对误差：%5.4e \n',max(abs(err3)));
% fprintf('均方根误差：   %5.4e \n',sqrt(sum(abs(err3).^2)/length(t)));
% 
% F4 = exp(-alpha * sqrt(s)) ./ sqrt(s);
% f4 = K * F4 ./t * log(2);
% f4_t = exp(-alpha^2 ./ (4*t)) ./ sqrt(pi * t);
% err4=abs(f4-f4_t);
% fprintf('F4(s)\n')
% fprintf('最大相对误差：%5.4e \n',max(abs(err4)));
% fprintf('均方根误差：   %5.4e \n',sqrt(sum(abs(err4).^2)/length(t)));

%%
%%
%测试双边的laplace变换
N = 12;
m = 1:N;
K = [
    -1.66666666666667e-02  1.60166666666667e+01 ...
    -1.24700000000000e+03 2.75543333333333e+04 ...
    -2.63280833333333e+05 1.32413870000000e+06 ...
    -3.89170553333333e+06 7.05328633333333e+06 ...
    -8.00533650000000e+06 5.55283050000000e+06 ...
    -2.15550720000000e+06 3.59251200000000e+05
    ];

t = [-50:0.1:-0.1 0 0.1:0.1:50];
% t = -50:0.1:50;
s = m' * (log(2)./t);
t0 = 5;%脉冲持续的时间是2*t0

F5 = 1./s.*(exp(t0*s)-exp(-t0*s));
f5 =K * F5 ./ t *log(2);
f5_t = stepfun(t,-t0)-stepfun(t,t0);
% err5=abs(f5-f5_t);
% fprintf('F5(s)\n')
% fprintf('最大相对误差：%5.4e \n',max(abs(err5)));
% fprintf('均方根误差：   %5.4e \n',sqrt(sum(abs(err5).^2)/length(t)));
figure
plot(t,f5_t,'r',t,f5,'g')
% plot(t,f5,'g')
axis([-50 50 -0.1 1.1])

F6 = 1./s.*(1-exp(-t0*s));
f6 =K * F6 ./ t *log(2);
f6_t = stepfun(t,0)-stepfun(t,t0);
% err6=abs(f6-f6_t);
% fprintf('F6(s)\n')
% fprintf('最大相对误差：%5.4e \n',max(abs(err6)));
% fprintf('均方根误差：   %5.4e \n',sqrt(sum(abs(err6).^2)/length(t)));
figure
plot(t,f6_t,'r',t,f6,'g')
% figure
% plot(t,f6,'g')
axis([-50 50 -0.1 1.1])

F7 = 1./s.*(exp(t0*s)-1);
f7 =K * F7 ./ t *log(2);
f7_t = stepfun(t,-t0)-stepfun(t,0);
figure
plot(t,f7_t,'r',t,f7,'g')
% figure
% plot(t,f7,'g')
axis([-50 50 -0.1 1.1])