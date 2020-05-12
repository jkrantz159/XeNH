function [XeSucc]=ParallelXeModel(XeCap,alpha,beta,eta,Resfrac,LVfrac,PLOTCHECK,atm,t,T,deltaXe,count)
%% XeModel
% After Parai and Mukhopadhyay 2018

%% Downwelling

% t=[0:0.1E6:200E6,201E6:1E6:3300E6,3305E6:5E6:4.565E9,4.568E9];
% T=4.568E9; % Age of Earth in yrs
% Above passed as globals from XenonModel
% alpha: growth rate (10E-10 - 10E-8 /Gyr)
% beta: sigmoid inflection point (0.08 - 10 Gyr)
% tau: time constant (10E-11 - 5E-8 /Gyr)
% XeCap: Carrying capacity/present day downwelling (0-5E8 atoms/gram)
Xed=zeros(1,length(t));
Xe130M=zeros(1,length(t));
Xe128M=zeros(1,length(t));
Mres=3.6E27*Resfrac;
ME=5.972E27;

% Initial Mantle
% AVCC at 1% of Mearth
% Concentration of 130Xe in CC = 5.38E-14 mol/g (Marty 2012)
% 128/130 AVCC = .5073 (Pepin 2000: isotopic composition of primordial...)


mol130M=(5.38E-14*ME)*(LVfrac/100); % mol 130Xe in mantle
Initial130=mol130M*6.022E23/Mres; % atomic concentration of 130Xe in mantle reservoir
Initial128=Initial130*.5073;

for j = 1:length(t)

    % Sigmoidal
    Xed(j) = XeCap/(1+exp(-alpha*(t(j)-beta)));
end
    
for i = 2:length(t)
    %% Degassing

    % Qp: Present day mantle processing rate (based on He flux at ridges) (6.1E17 g/yr)
    % Eta: 1.6E-10 - 9.9E-10


    % Mantle processing rate
    % Q(t) = Qp*exp(eta*(T-t));

    % Mass of mantle processed
    Qp = 6.1E17; % g/yr, Present day processing rate
    dM = (Qp/eta)*(exp(eta*(T-t(i-1)))-exp(eta*(T-t(i))));

    % Normalized dM
    % dM/Mres (Mres=mass of reservoir)
    % assume 90% of mantle is convecting mantle -> Mres = 3.6E27 grams
    % i.e. Resfrac == 0.9
    



    if i==2
        Xe130Mlast = Initial130;
        Xe128Mlast = Initial128;
        atmlast = atm(1);
        Xedlast=Xed(i);
    else
        Xe130Mlast = Xe130M(i-1);
        Xe128Mlast = Xe128M(i-1);
        atmlast = atm(i-1);
        Xedlast=Xed(i-1);
    end

    Xe130M (i) = Xe130Mlast + (dM/Mres)*(Xedlast - Xe130Mlast);
    Xe128M (i) = Xe128Mlast + (dM/Mres)*(Xedlast*(atmlast)-Xe128Mlast);

end


%% Criterion for success

if (Xe130M(end) >= 4.3E5) && (Xe130M(end) <= 9.2E5) % atomic concentration of 130Xe in mantle
    Succ1=1;
else
    Succ1=0;
end

if (Xe128M(end)/Xe130M(end) >= 0.475) && (Xe128M(end)/Xe130M(end) <= 0.478)
    Succ2=1;
else
    Succ2=0;
end

if Succ1==1 && Succ2==1
    XeSucc=1;

else
    XeSucc=0;
end


if PLOTCHECK==1
    figure(3)
    hold on
    plot(t,deltaXe(Xe128M./Xe130M),'-k')
    xlabel('Myr')
    ylabel('delta Xe 128/130 in mantle')
    hold off
end
    



    


%     figure(11)
%     hold on
%     plot(t,Xe128M,'-r')
%     plot(t,Xe130M,'-k')
%     xlabel('Myr')
%     ylabel('Xe in Mantle')
%     legend('128','130')
%     hold off
%     
%     figure(12)
%     hold on
%     plot(t,Xed,'k')
%     plot(t,Xed.*atm,'r')
%     xlabel('Myr')
%     ylabel('Downwelling Xenon')
%     legend('130','128')
%     hold off
%     
%     figure(13)
%     hold on
%     plot(t,log10(Xedlast./Xe130M),'k')
%     plot(t,log10(Xedlast.*(atmlast)./Xe128M),'r')
%     xlabel('Myr')
%     ylabel('Downwelling Xe Relative to Mantle Xe')
%     legend('130','128')
%     hold off

end
