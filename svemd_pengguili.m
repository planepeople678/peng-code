%D=xlsread('channel22.xls');
D=[dataX,dataY,dataZ];
D_rank=rank(D);
[U,S,V]=svd(D);
function imf = emd(x)
x=transpose(x(:));
%x2=xlsread('channel22.xls');
imf = [];
while ~ismonotonic(x)
   x1 = x;
   sd = Inf;
   while (sd > 0.1) || ~isimf(x1)
      s1 = getspline(x1);
      s2 = -getspline(-x1);
      x2 = x1-(s1+s2)/2;
      sd = sum((x1-x2).^2)/sum(x1.^2);
      x1 = x2;
   end
   imf{end+1} = x1;
   x          = x-x1;
end
imf{end+1} = x;
function u = ismonotonic(x) 
u1 = length(findpeaks(x))*length(findpeaks(-x));
if u1 > 0, u = 0;
else      u = 1; 
end
function u = isimf(x)
N  = length(x);
u1 = sum(x(1:N-1).*x(2:N) < 0);
u2 = length(findpeaks(x))+length(findpeaks(-x));
if abs(u1-u2) > 1, u = 0;
else  u = 1; 
end
function s = getspline(x)
N = length(x);
p = findpeaks(x);
s = spline([0 p N+1],[0 x(p) 0],1:N);

z=ans{1,1};
z1=hilbert(z);
[A,fa,tt]=hhspectrum1(z);
[E,tt1]=toimage(A,fa,tt,length(tt));
figure()
disp_hhs(E,tt1)  
E=flipud(E);
for k=1:size(E,1)
bjp(k)=sum(E(k,:))*1/Fs; 
end
t=0:0.01:4-1/Fs;
N=2499;deta=t(2)-t(1);fs=1/deta;
f=(1:N-2)/N*(fs/2);
figure(5)
plot(f,bjp);
xlabel('feq/ Hz');
ylabel('single amplitude');
title('single marginal spectru ') 


o=xlsread('channel3-6000.xls');
x1=(o(:,2));
%y1=(o(:,3));
%z1=(o(:,4));
imf=emd(x1);
for i=1:11
    z=(imf(i,:));
    fn=size(z,2);
    for j=1:fn
        v=z(:,j);
        v1=hilbert(v);
        V=v1*v1;
        E(j)=sum(V);
    end
    frameTime=fram(fn,wlen,fs);
    figure(i);
    plot(frameTime,E,'b');
    title('imf',num2str(i),'energy');
    ylable('amplitude');xlable('time/s');
end

figure(1)
hold on
subplot 311;
plot(ans{1, 1});
xlabel('采样时间/ms');
ylabel('信号幅值');
title('imf1');
subplot 312;
plot(ans{1,2});
xlabel('采样时间/ms');
ylabel('信号幅值');
title('imf2');
subplot 313;
plot(ans{1, 3});
xlabel('采样时间/ms');
ylabel('信号幅值');
title('imf3');

figure(2)
hold on
subplot 311;
plot(ans{1, 4});
xlabel('采样时间/ms');
ylabel('信号幅值');
title('imf4');
subplot 312;
plot(ans{1,5});
xlabel('采样时间/ms');
ylabel('信号幅值');
title('imf5');
subplot 313;
plot(ans{1, 6});
xlabel('采样时间/ms');
ylabel('信号幅值');
title('imf6');

figure(3)
hold on
subplot 311;
plot(ans{1, 7});
xlabel('采样时间/ms');
ylabel('信号幅值');
title('imf7');
subplot 312;
plot(ans{1, 8});
xlabel('采样时间/ms');
ylabel('信号幅值');
title('imf8');
subplot 313;
plot(ans{1, 9});
xlabel('采样时间/ms');
ylabel('信号幅值');
title('res');


