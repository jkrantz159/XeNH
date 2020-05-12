function [NeSucc]=ParallelNeModel(NeCap,alpha,beta,eta,Resfrac,LVfrac,PLOTCHECK,t,T,count)
%% Ne Model
% After Parai and Mukhopadhyay 2018

%% Downwelling

% t=[0:0.1E6:200E6,201E6:1E6:3300E6,3305E6:5E6:4.565E9,4.568E9];
% T=4.568E9; % Age of Earth in yrs
% Above passed as globals from XenonModel
% alpha: growth rate (10E-10 - 10E-8 /Gyr)
% beta: sigmoid inflection point (0.08 - 10 Gyr)
% tau: time constant (10E-11 - 5E-8 /Gyr)
% XeCap: Carrying capacity/present day downwelling (0-5E8 atoms/gram)

Ned=zeros(1,length(t));
Ne20M=zeros(1,length(t));
Ne22M=zeros(1,length(t));
Mres=3.6E27*Resfrac;
ME=5.972E27;

% Initial Mantle
% AVCC at 1% of Mearth
% Concentration of 22Ne in CC = 1.62E-12 mol/g (Marty 2012)
% 20/22 AVCC(PSN) = 13.36 (Williams and Mukhopadhyay 2018)


mol22M=(1.62E-12*ME)*(LVfrac/100); % mol 22Ne in mantle
Initial22=mol22M*6.022E23/Mres; % atomic concentration of 22Ne in mantle reservoir
Initial20=Initial22.*13.36;

for j = 1:length(t)

    % Sigmoidal
    Ned(j) = NeCap/(1+exp(-alpha*(t(j)-beta)));
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
        Ne20Mlast = Initial20;
        Ne22Mlast = Initial22;
        Nedlast=Ned(i);
    else
        Ne20Mlast = Ne20M(i-1);
        Ne22Mlast = Ne22M(i-1);
        Nedlast=Ned(i-1);
    end

    Ne20M (i) = Ne20Mlast + (dM/Mres)*(Nedlast - Ne20Mlast);
    Ne22M (i) = Ne22Mlast + (dM/Mres)*(Nedlast./(9.80)-Ne22Mlast);

end

%% Criterion for success

if (Ne22M(end) >= (5.8E-15-3.2E-15)*6.022E23/Resfrac) && (Ne22M(end) <= (5.8E-15+3.2E-15)*6.022E23/Resfrac) % atomic concentration of 22Ne in mantle
    Succ1=1;
else
    Succ1=0;
end
    

if (Ne20M(end)/Ne22M(end) >= 12-.5) && (Ne20M(end)/Ne22M(end) <= 12+.5)
    Succ2=1;
else
    Succ2=0;
end

if Succ1==1 && Succ2==1
    NeSucc=1;
else
    NeSucc=0;
end
    
if PLOTCHECK==1
    fprintf('Success Found (Ne) ')
    figure(2)
    hold on
    plot(t,Ne20M./Ne22M,'-k')
    xlabel('Myr')
    ylabel('^2^0Ne/^2^2Ne')
    hold off
    saveas(gcf,[pwd '\NeFigs\temp' num2str(count) '.jpg']);
elseif PLOTCHECK==3
    fprintf('Success Found (Ne) ')
    figure(2)
    hold on
    plot(t,Ne20M./Ne22M,'-g')
    xlabel('Myr')
    ylabel('^2^0Ne/^2^2Ne')
    hold off
elseif PLOTCHECK==4
    fprintf('Success Found (Ne) ')
    figure(2)
    hold on
    plot(t,Ne20M./Ne22M,'-b')
    xlabel('Myr')
    ylabel('^2^0Ne/^2^2Ne')
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
