%% Mass Calculation tomfoolery

successes=csvread('LowLVIndXeNHsuccess.txt');
eta=successes(:,1);
eta_min=min(eta);
eta_max=max(eta);
Resfrac_min=min(successes(:,11));
Resfrac_max=max(successes(:,11));

t=[0:0.1E6:200E6,201E6:1E6:3300E6,3305E6:5E6:4.565E9,4.568E9];
T=4.568E9; % Age of Earth in yrs
TotalMassProcessed_Min = 0;
TotalMassProcessed_Max = 0;
Mres_min=3.6E27*Resfrac_min;
Mres_max=3.6E27*Resfrac_max;
ME=5.972E27;


for i = 2:length(t)
    Qp = 6.1E17; % g/yr, Present day processing rate
    dM = (Qp/eta_min)*(exp(eta_min*(T-t(i-1)))-exp(eta_min*(T-t(i))));
    
    TotalMassProcessed_Min = TotalMassProcessed_Min + dM;
end

for i = 2:length(t)
    Qp = 6.1E17; % g/yr, Present day processing rate
    dM = (Qp/eta_max)*(exp(eta_max*(T-t(i-1)))-exp(eta_max*(T-t(i))));
    
    TotalMassProcessed_Max = TotalMassProcessed_Max + dM;
end


TotalMassProcessed_Min./Mres_min
TotalMassProcessed_Max./Mres_max