%% Parallel Combined Random Model

%% Initial Conditions
fprintf('Start time: %s.\n',datestr(now))
% Late Veneer of carbonaceous chondrites (1%, consider also 0.5%, 0.1%)
% 1% of Earth mass as AVCC (Murchison and Orgueil)

% Initial d15 as PSN = -40
% Subduction of sediments at d15 = +5

% Atmosphere starting at 39 per mille frac, linear to modern at 2 Ga

t=[0:0.1E6:200E6,201E6:1E6:3300E6,3305E6:5E6:4.565E9,4.568E9];
T=4.568E9; % Age of Earth in yrs
atm=zeros(1,length(t));

% Earth Atmosphere (Porcelli, Ballentine, and Wieler, 2002)
% 128/132 = 0.0714
% 130/132 = 0.1514
% 128/130 today = 0.0714/0.1514 -= 0.4716

% -39 per mille fractionation at t=0
% 128/130 t0 = 0.49
deltaXe= @(x) (((x/(0.0714/0.1514))-1)*1000);       % Converts 128/130 to dXe
deltaN= @(x) (((x/((1-.99636)/.99636))-1)*1000);

for i=1:length(t)
    time=t(i);
        if time < 2.568E9
            a=0.5178-(time/2.568E9)*(0.5178-0.4716);
            atm(i)=a;            
        else
            atm(i)=0.4716;                          % Atomic ratio of 128/130
        end
end


run=0; Resfrac=0.9; LVfrac=1; FullSuccess=0; DoubleCount=0;


runs=1E8;

parfor count=1:runs
        
%     i = alpha                                    
%     j = beta                                     
%     m = eta                                      
%     n = Xed aka Xe carrying capacity             
%     p = Nd aka N carrying capacity

    i=10.^(-10 + (-7+10).*rand(1));                                         % 1E-7 TO 1E-10
    j=(10.*rand(1))*1E9;                                                    % 0 to 12
    m=7.5E-10+(rand(1).*.5E-10);                                            % 7E-10 to 8E-10
    p=4*10.^(17.*rand(1));                                                  % 4 to 4E17
    n=5*10.^(8.*rand(1));                                                   % 5 to 5E8
    
    [NSucc]=ParallelNewNModel(p,i,j,m,Resfrac,LVfrac,0,t,T,deltaN,count);
    [XeSucc]=ParallelXeModel(n,i,j,m,Resfrac,LVfrac,0,atm,t,T,deltaXe,count);


    if NSucc==1 && XeSucc==1
        ParallelXeModel(n,i,j,m,Resfrac,LVfrac,1,atm,t,T,deltaXe,count);
        ParallelNewNModel(p,i,j,m,Resfrac,LVfrac,1,t,T,deltaN,count);
        fprintf('Success Found (N) ')
        fprintf('Success Found (Xe) ')
        fprintf('\n')
        parsave(i,j,m,n,p,Resfrac,LVfrac);
    end
    
    
    

    
end

fprintf('End time: %s.\n',datestr(now))