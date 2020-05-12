%% LV=1, Xe-N only

successes=csvread('success.txt');


eta=successes(:,3);
a=successes(:,1);
b=successes(:,2);
Xe=successes(:,4);
N=successes(:,5);

figure(14)
data=(eta)*1E10;
pdAge = fitdist(data,'Kernel','BandWidth',.003);
x = 7.6:.001:8;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'k-','LineWidth',2)
title('Processing Rate KDE')
axis([6 8 0 50])
xlabel('\eta * 10^1^0')
ylabel('KDE')
hold off

figure(11)
data=b./1E9;
pdAge = fitdist(data,'Kernel','BandWidth',.15);
x = 0:.1:10;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'k-','LineWidth',2)
title('Inflection Point KDE')
xlabel('Time (Gyr)')
ylabel('KDE')
hold off

figure(12)
data=log(a);
pdAge = fitdist(data,'Kernel','BandWidth',.2);
x = -22:.1:-15;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'k-','LineWidth',2)
title('Growth Rate KDE')
xlabel('Log \alpha')
ylabel('KDE')
hold off


figure(13)
data1=log10(Xe);
data2=log10(N);
pdAge1 = fitdist(data1,'Kernel','BandWidth',.2);
pdAge2 = fitdist(data2,'Kernel','BandWidth',.2);
x = 5:.1:25;
ySix1 = pdf(pdAge1,x);
ySix2 = pdf(pdAge2,x);
hold on
plot(x,ySix1,'r-','LineWidth',2)
plot(x,ySix2,'g-','LineWidth',2)
title('Recycling KDE')
xlabel('log_1_0 Recycling')
ylabel('KDE')
legend('Xe','N')
hold off

%% LV=1

successes=csvread('IndXeNHsuccess.txt');


eta=successes(:,1);
Hd=successes(:,2);
Ha=successes(:,3);
Hb=successes(:,4);
Nd=successes(:,5);
Na=successes(:,6);
Nb=successes(:,7);
Xd=successes(:,8);
Xa=successes(:,9);
Xb=successes(:,10);

figure(4)
data=(eta)*1E10;
pdAge = fitdist(data,'Kernel','BandWidth',.003);
x = 7.6:.001:8;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'k-','LineWidth',2)
title('Processing Rate KDE')
axis([6 8 0 50])
xlabel('\eta * 10^1^0')
ylabel('KDE')
hold off

figure(1)
data=Xb./1E9;
pdAge = fitdist(data,'Kernel','BandWidth',.15);
x = 0:.1:10;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'r-','LineWidth',2)
title('Inflection Point KDE')
xlabel('Time (Gyr)')
ylabel('KDE')
hold off

figure(2)
data=log(Xa);
pdAge = fitdist(data,'Kernel','BandWidth',.2);
x = -22:.1:-15;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'r-','LineWidth',2)
title('Growth Rate KDE')
xlabel('Log \alpha')
ylabel('KDE')
hold off

figure(1)
data=Nb./1E9;
pdAge = fitdist(data,'Kernel','BandWidth',.15);
x = 0:.1:10;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'g-','LineWidth',2)
hold off

figure(2)
data=log(Na);
pdAge = fitdist(data,'Kernel','BandWidth',.2);
x = -25:.1:-15;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'g-','LineWidth',2)
hold off

figure(1)
data=Hb./1E9;
pdAge = fitdist(data,'Kernel','BandWidth',.15);
x = 0:.1:10;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'b-','LineWidth',2)
legend('Xe','N','H')
hold off

figure(2)
data=log(Ha);
pdAge = fitdist(data,'Kernel','BandWidth',.2);
x = -25:.1:-15;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'b-','LineWidth',2)
legend('Xe','N','H')
hold off


figure(3)
data1=log10(Xd);
data2=log10(Nd);
data3=log10(Hd);
pdAge1 = fitdist(data1,'Kernel','BandWidth',.2);
pdAge2 = fitdist(data2,'Kernel','BandWidth',.2);
pdAge3 = fitdist(data3,'Kernel','BandWidth',.2);
x = 5:.1:25;
ySix1 = pdf(pdAge1,x);
ySix2 = pdf(pdAge2,x);
ySix3 = pdf(pdAge3,x);
hold on
plot(x,ySix1,'r-','LineWidth',2)
plot(x,ySix2,'g-','LineWidth',2)
plot(x,ySix3,'b-','LineWidth',2)
title('Recycling KDE')
xlabel('log_1_0 Recycling')
ylabel('KDE')
legend('Xe','N','H')
hold off


%% LV=0.05

successes=csvread('LowLVIndXeNHsuccess.txt');

eta=successes(:,1);
Hd=successes(:,2);
Ha=successes(:,3);
Hb=successes(:,4);
Nd=successes(:,5);
Na=successes(:,6);
Nb=successes(:,7);
Xd=successes(:,8);
Xa=successes(:,9);
Xb=successes(:,10);

figure(4)
data=(eta)*1E10;
pdAge = fitdist(data,'Kernel','BandWidth',.003);
x = 6.25:.001:6.75;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'k--','LineWidth',2)
title('Processing Rate KDE')
axis([6 8 0 50])
xlabel('\eta * 10^1^0')
ylabel('KDE')
hold off

figure(1)
data=Xb./1E9;
pdAge = fitdist(data,'Kernel','BandWidth',.15);
x = 0:.1:10;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'r--','LineWidth',2)
title('Inflection Point KDE')
xlabel('Time (Gyr)')
ylabel('KDE')
hold off

figure(2)
data=log(Xa);
pdAge = fitdist(data,'Kernel','BandWidth',.2);
x = -22:.1:-15;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'r--','LineWidth',2)
title('Growth Rate KDE')
xlabel('Log \alpha')
ylabel('KDE')
hold off

figure(1)
data=Nb./1E9;
pdAge = fitdist(data,'Kernel','BandWidth',.15);
x = 0:.1:10;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'g--','LineWidth',2)
hold off

figure(2)
data=log(Na);
pdAge = fitdist(data,'Kernel','BandWidth',.2);
x = -25:.1:-15;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'g--','LineWidth',2)
hold off

figure(1)
data=Hb./1E9;
pdAge = fitdist(data,'Kernel','BandWidth',.15);
x = 0:.1:10;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'b--','LineWidth',2)
legend('Xe','N','H')
hold off

figure(2)
data=log(Ha);
pdAge = fitdist(data,'Kernel','BandWidth',.2);
x = -25:.1:-15;
ySix = pdf(pdAge,x);
hold on
plot(x,ySix,'b--','LineWidth',2)
legend('Xe','N','H')
hold off


figure(3)
data1=log10(Xd);
data2=log10(Nd);
data3=log10(Hd);
pdAge1 = fitdist(data1,'Kernel','BandWidth',.2);
pdAge2 = fitdist(data2,'Kernel','BandWidth',.2);
pdAge3 = fitdist(data3,'Kernel','BandWidth',.2);
x = 5:.1:25;
ySix1 = pdf(pdAge1,x);
ySix2 = pdf(pdAge2,x);
ySix3 = pdf(pdAge3,x);
hold on
plot(x,ySix1,'r--','LineWidth',2)
plot(x,ySix2,'g--','LineWidth',2)
plot(x,ySix3,'b--','LineWidth',2)
title('Recycling KDE')
xlabel('log_1_0 Recycling')
ylabel('KDE')
legend('Xe','N','H')
hold off
%% Histograms
% 
% figure(101)
% hold on
% histHa=histogram(successes(:,3),10);
% histNa=histogram(successes(:,6),10);
% histXa=histogram(successes(:,9),10);
% histHa.Normalization = 'probability';
% histNa.Normalization = 'probability';
% histXa.Normalization = 'probability';
% histHa.BinWidth=2.1E-9;
% histNa.BinWidth=2.1E-9;
% histXa.BinWidth=2.1E-9;
% histHa.FaceColor='b';
% histNa.FaceColor='g';
% histXa.FaceColor='r';
% title('Alpha values')
% legend('H','N','Xe')
% hold off
% 
% 
% Hbmean=mean(successes(:,4));
% Hstd=std(successes(:,4));
% Nbmean=mean(successes(:,7));
% Nstd=std(successes(:,7));
% Xbmean=mean(successes(:,10));
% Xstd=std(successes(:,10));
% table([Hbmean;Nbmean;Xbmean],[Hstd;Nstd;Xstd],'VariableNames',{'Mean','STD'},'RowNames',{'H','N','Xe'})
% 
% 
% 
% figure(102)
% 
% subplot(3,1,1)
% hold on
% histHb=histogram(successes(:,4),10);
% histHb.Normalization = 'probability';
% histHb.BinWidth=120000000;
% histHb.FaceColor='b';
% y = 2E9:0.1E9:6E9;
% mu = Hbmean;
% sigma = Hstd;
% f = 1E8*exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',1.5)
% axis([2E9 6E9 0 0.3])
% ylabel('H')
% title('Beta values')
% hold off
% 
% subplot(3,1,2)
% hold on
% histNb=histogram(successes(:,7),10);
% histNb.Normalization = 'probability';
% histNb.BinWidth=120000000;
% histNb.FaceColor='g';
% axis([2E9 6E9 0 0.3])
% ylabel('N')
% y = 2E9:0.1E9:6E9;
% mu = Nbmean;
% sigma = Nstd;
% f = 1E8*exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',1.5)
% hold off
% 
% subplot(3,1,3)
% hold on
% histXb=histogram(successes(:,10),10);
% histXb.BinWidth=120000000;
% histXb.Normalization = 'probability';
% histXb.FaceColor='r';
% axis([2E9 6E9 0 0.3])
% xlabel('Time')
% ylabel('Xe')
% y = 2E9:0.1E9:6E9;
% mu = Xbmean;
% sigma = Xstd;
% f = 1E8*exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',1.5)
% hold off