function [NSucc]=ParallelNewNModel(NCap,alpha,beta,eta,Resfrac,LVfrac,PLOTCHECK,t,T,deltaN,count)
%% New N Model
% After Parai and Mukhopadhyay 2018 and Barry and Hilton 2016

%% Downwelling

% t=[0:0.1E6:200E6,201E6:1E6:3300E6,3305E6:5E6:4.565E9,4.568E9];
% T=4.568E9; % Age of Earth in yrs
% Above passed as globals from XenonModel
% alpha: growth rate (10E-10 - 10E-8 /Gyr)
% beta: sigmoid inflection point (0.08 - 10 Gyr)
% tau: time constant (10E-11 - 5E-8 /Gyr)
% XeCap: Carrying capacity/present day downwelling (0-5E8 atoms/gram)

Nd=zeros(1,length(t));
N14M=zeros(1,length(t));
N15M=zeros(1,length(t));
Mres=3.6E27*Resfrac;
ME=5.972E27;

% Initial Mantle
% AVCC at 1% of Mearth
% 1% of Earth mass as AVCC (Murchison and Orgueil)
% Concentration of N in CC = 0.1519 wt% (Sephton et al. 2003)


g14M=(.001519*ME*LVfrac/100); % g of 14N in mantle
Initial14=(g14M/14*.99636*6.022E23)/(Mres); % atomic concentration of 14N in mantle reservoir
Initial15=Initial14*2.3E-3;           %15/14 ratio=2.3E-3 from Owen et al. 2001

for j = 1:length(t)

    % Sigmoidal
    Nd(j) = NCap/(1+exp(-alpha*(t(j)-beta)));
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
        N14Mlast = Initial14;
        N15Mlast = Initial15;
        Ndlast=Nd(i);
    else
        N14Mlast = N14M(i-1);
        N15Mlast = N15M(i-1);
        Ndlast=Nd(i-1);
    end

    N14M (i) = N14Mlast + (dM/Mres)*(Ndlast - N14Mlast);
    N15M (i) = N15Mlast + (dM/Mres)*(Ndlast*(0.003671565)-N15Mlast); % Sed 15/14 from d15=5

end

%% Criterion for success

if (N14M(end) >= 7.06E19*.99636*6.022E23/Mres) && (N14M(end) <= 9.78E21*.99636*6.022E23/Mres) % atomic concentration of 14N in mantle reservoir
                    %1.3074E16                                         %1.8111E18
    Succ1=1;
else
    Succ1=0;
end
    

if (N15M(end)/N14M(end) >= 0.0036275) && (N15M(end)/N14M(end) <= 0.0036425)
    Succ2=1;
else
    Succ2=0;
end

if Succ1==1 && Succ2==1
    NSucc=1;

else
    NSucc=0;
end
    
if PLOTCHECK==1
    figure(1)
    hold on
    plot(t,deltaN(N15M./N14M),'-k')
    xlabel('Myr')
    ylabel('Delta N')
    hold off
end

    




%     figure(2)
%     hold on
%     plot(t,N15M,'-r')
%     xlabel('Myr')
%     ylabel('N15 in Mantle')
%     hold off
%     
%     figure(3)
%     hold on
%     plot(t,N14M,'-k')
%     xlabel('Myr')
%     ylabel('N14 in Mantle')
%     hold off
%     
%     figure(4)
%     hold on
%     plot(t,Nd)
%     xlabel('Myr')
%     ylabel('Downwelling N14')
%     hold off
%     
%     figure(5)
%     hold on
%     plot(t,log10(Nd./N14M))
%     xlabel('Myr')
%     ylabel('Downwelling N14 Relative to Mantle N14')
%     hold off
%     
%     figure(6)
%     hold on
%     plot(t,Nd*(0.003671565))
%     xlabel('Myr')
%     ylabel('Downwelling N15')
%     hold off
%     
%     figure(7)
%     hold on
%     plot(t,log10(Nd*(0.003671565)./N15M))
%     xlabel('Myr')
%     ylabel('Downwelling N15 Relative to Mantle N15')
%     hold off
end
