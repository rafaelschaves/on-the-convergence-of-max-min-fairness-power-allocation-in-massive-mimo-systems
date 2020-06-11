addpath('./functions/');

% Cheking deirectory

dir_save    = './results/power_allocation/downlink/';

if ~exist(dir_save,'dir')
    mkdir(dir_save);
end

root_save_1 = [dir_save 'results_error_prob_'];                            % Results necessary to plot Fig.2
root_save_2 = [dir_save 'results_lower_bound_prob_'];                      % Results necessary to plot Fig.3
root_save_3 = [dir_save 'results_iterations_'];                            % Results necessary to plot Fig.4

% Checking variables

if ~exist('MC','var')
    MC = 1000;                                                             % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
end

if ~exist('M','var')
    M = 64;                                                                % Number of antennas at the base station
end

if ~exist('K','var')
    K = 8;                                                                 % Number of users at the cell
end

if ~exist('bs_power','var')
    bs_power = 1;                                                          % BS power
end

if ~exist('radius','var')
    radius = 1000;                                                         % Cell's radius
end

if ~exist('channel_type','var')
    channel_type = 'ur-los';
end

N_ALG = 3;

commcell.nAntennas       = M;                                              % Number of Antennas
commcell.nUsers          = K;                                              % Number of Users
commcell.radius          = radius;                                         % Cell's raidus (circumradius) in meters
commcell.bsHeight        = 32;                                             % Height of base station in meters
commcell.userHeight      = [1 2];                                          % Height of user terminals in meters ([min max])
commcell.frequency       = 1.9e9;                                          % Carrier frequency in Hz
commcell.meanShadowFad   = 0;                                              % Shadow fading mean in dB
commcell.stdDevShadowFad = 8;                                              % Shadow fading standard deviation in dB
commcell.city            = 'large';                                        % Type of city

linkprop.bsPower         = bs_power;                                       % in Watts
linkprop.userPower       = 0.2;                                            % in Watts
linkprop.AntennaGainBS   = 0;                                              % in dBi
linkprop.AntennaGainUser = 0;                                              % in dBi
linkprop.noiseFigureBS   = 9;                                              % in dB
linkprop.noiseFigureUser = 9 ;                                             % in dB
linkprop.bandwidth       = 20e6;                                           % in Hz

[~,snr_db] = linkBudgetCalculation(linkprop);                              % SNR in dB
snr        = 10.^(snr_db/10);

gamma      = zeros(MC,N_ALG);
n_it       = zeros(MC,N_ALG);
lambda_max = zeros(MC,1);

for mc = 1:MC
    mc
    
    [G,beta]   = massiveMIMOChannel(commcell,channel_type);
    
    [gamma(mc,1),~,~,n_it(mc,1),R] = maxMinFairness(G,beta,snr,'algorithm 1');
    [gamma(mc,2),~,~,n_it(mc,2)]   = maxMinFairness(G,beta,snr,'algorithm 2');
    [gamma(mc,3),~,~,n_it(mc,3)]   = maxMinFairness(G,beta,snr,'algorithm 3');
    
    lambda_max(mc) = max(sum(R,2));
end

error_prob  = sum(abs(gamma(:,1) - gamma(:,3)) > 1e-6)/MC;                 % Compute error probability
n_low_bound = sum(gamma > 1./lambda_max);                                  % Compute the number of times that the right endpoint would work as a left endpoint

% Results necessary to plot Fig.2

save([root_save_1 strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' ...
      num2str(K) '_cell_radius_' num2str(radius) '_m_BS_power_' ...
      num2str(bs_power) '_W_MC_' num2str(MC) '.mat'],'error_prob');
  
% Results necessary to plot Fig.3

save([root_save_2 strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' ...
      num2str(K) '_cell_radius_' num2str(radius) '_m_BS_power_' ...
      num2str(bs_power) '_W_MC_' num2str(MC) '.mat'],'n_low_bound');

% Results necessary to plot Fig.4

save([root_save_3 strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' ...
      num2str(K) '_cell_radius_' num2str(radius) '_m_BS_power_' ...
      num2str(bs_power) '_W_MC_' num2str(MC) '.mat'],'n_it');
