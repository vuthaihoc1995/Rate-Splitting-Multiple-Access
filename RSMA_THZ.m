clc, clear; close all;
% Network
K = 2;                % Number of users
N = 16;              % Number of antenna setting at BS
BW = 100*1e6;       % Bandwidth in MHz (100 MHz )
NFdB = 6.5;         % From 280 GHz - 330 GHz
PTdB = -5:2.5:20;
PTdB1 = min(PTdB):2:max(PTdB);
% Noise power in dBm
sig2dBm = -174 + 10*log10(BW) + NFdB;
% Carrier frequency in GHz (300 GHz)
fc = 300*1e9;
% Absorbtion Loss Coefficient Measured at 300 GHz
kabs = 0.0033;
lambda = physconst('LightSpeed')/fc;
% Antenna gain of Tx
G_Rx= db2pow(3);
% power allocation for common + private stream
alpha_c = 0.4;
alpha_k = (1-alpha_c)/K;
% rate splitting for common + private stream
zeta_c = 0.6;
zeta_k = (1-zeta_c)/K;
% Target rate for common + private stream
r_c = 0.5;
r_k = 0.25;
% impefect SIC coefficient
psik = 0.01;
% Threshold for decoding common + private stream
threshold_c = 2^(r_c/zeta_c) - 1 ;
threshold_k = 2^(r_k/zeta_k)  - 1 ;
% User-k location
R_cell = 10;
d_BS_Uk       = (1:1:K)*R_cell/K;
% LSF * MAP part of channel gk
tau_k = exp(-d_BS_Uk*kabs)./d_BS_Uk.^2*(lambda/(4*pi))^2*G_Rx;
%fading coefficient
m_k = K+ 1 - (1:1:K); Omega_k = 1;

%% Simulation for RSMA-THz only
% channel realization
trial = 1e4;
for iter = 1:trial
    
    % Channel matrix Nakagami-m: G (N*K);
    G = zeros(N,K);
    for k=1:K
        angle_k = -pi + pi*rand(N,1);
        G(:,k) = sqrt(gamrnd(m_k(k),Omega_k/m_k(k),[N,1])).*exp(1i*angle_k);
    end

    % Generate precoding weight p_k and p_c based on hybrid ZF and MRT:
    P_k = G*inv(G'*G)*diag(vecnorm(G));
    p_c = P_k*ones(K,1); % or sum(P_k,2);
    % Channel gain = Channel and weight multiplication
    for k=1:K
        gkh_pc(k,iter)         = abs(p_c'*G(:,k)).^2;
        gkh_pk(k,iter)         = abs(P_k(:,k)'*G(:,k)).^2;
        G_j = G; G_j(:,k) = [];
        for j=1:K-1
            gkh_pj(k,iter)          = abs(P_k(:,k)'*G_j(:,j)).^2;
        end
    end
end
% compute outage for RSMA-THz only
CDF_Uk_sim = zeros(K,length(PTdB));
Rate_Uk_sim = zeros(K,length(PTdB));
parfor idx1 =1:length(PTdB)
    snr = db2pow(PTdB(idx1) - sig2dBm);

    for k=1:K
        % SINR at user k;
        gamma_kc     = alpha_c*snr*tau_k(k)*gkh_pc(k,:)./(alpha_k*snr*tau_k(k)*gkh_pk(k,:)  + alpha_k*snr*tau_k(k)*gkh_pj(k,:) + 1);
        gamma_kp     = alpha_k*snr*tau_k(k)*gkh_pk(k,:)./(psik*alpha_c*snr*tau_k(k)*gkh_pc(k,:) + alpha_k*snr*tau_k(k)*gkh_pj(k,:) + 1);
        
        count_nonop = 0; Rate_Uk_common = 0; Rate_Uk_private = 0
        for iter = 1:trial
            % Count Non-Outage Event
            if gamma_kc(iter) > threshold_c && gamma_kp(iter) > threshold_k
                count_nonop = count_nonop + 1 ;
            end
            % Sum Rate for each trial 
            Rate_Uk_common = Rate_Uk_common + zeta_c*log2(1+ gamma_kc(iter));
            Rate_Uk_private = Rate_Uk_private + zeta_k*log2(1+ gamma_kc(iter));
        end
        % Count Outage Event
        
        CDF_Uk_sim(k,idx1) = 1- count_nonop/trial;
        
        % Count Achievable Rate 
        
        Rate_Uk_sim_common(k,idx1) = Rate_Uk_common/trial;
        Rate_Uk_sim_private(k,idx1) = Rate_Uk_private/trial
    end
end

%% Analytical for Outage Probability only
CDF_Uk_ana = zeros(K,length(PTdB1));
for idx = 1:length(PTdB1)
    snr = db2pow(PTdB1(idx) - sig2dBm);
    for k=1:K
        if alpha_c - threshold_c*alpha_k > 0 && alpha_k - threshold_k*alpha_c*psik > 0
            temp =   max(threshold_c/(alpha_c - threshold_c*alpha_k),threshold_k/(alpha_k - threshold_k*alpha_c*psik));

            CDF_Uk_ana(k,idx) =  1 - 1/gamma(N*m_k(k))*igamma(N*m_k(k),m_k(k)/Omega_k/snr/tau_k(k)*temp);
        else
            CDF_Uk_ana(k,idx) = 1;
        end
    end

end

%% % Define color
cyan        = [0.00,0.45,0.74];
red         = [1, 0, 0];
green       = [0, 1, 0];
black       = [0 0 0];
blue        = [0 0 1];
yellow      = [0.93,0.69,0.13];
orange      = [0.85,0.33,0.10];
purple      = [0.49,0.18,0.56];
magenta = [1 0 1];
%
figure(1);
semilogy(PTdB,CDF_Uk_sim(1,:),'o','Color',red,'MarkerSize',8,'LineWidth',1.1); hold on;
semilogy(PTdB,CDF_Uk_sim(2,:),'s','Color',blue,'MarkerSize',8,'LineWidth',1.1); hold on;
semilogy(PTdB1,CDF_Uk_ana,'-','Color',black,'MarkerSize',8,'LineWidth',1.1); hold on;
legend('$U_1$','$U_2$',...
    'interpreter','latex','FontSize',13,'FontName','Times New Roman')
ylim([1e-6 1]); hold on;
xlim([min(PTdB) max(PTdB)]); hold on;
xlabel('Transmit power [dBm]','FontSize',14,'FontName','Times New Roman');
ylabel('Outage Probability','FontSize',14,'FontName','Times New Roman');

%
figure(2);
plot(PTdB,Rate_Uk_sim_common(1,:),'o-','Color',red,'MarkerSize',8,'LineWidth',1.1); hold on;
plot(PTdB,Rate_Uk_sim_common(2,:),'+-','Color',blue,'MarkerSize',8,'LineWidth',1.1); hold on;

plot(PTdB,Rate_Uk_sim_private(1,:),'d-','Color',green,'MarkerSize',8,'LineWidth',1.1); hold on;
plot(PTdB,Rate_Uk_sim_private(2,:),'s-','Color',magenta,'MarkerSize',8,'LineWidth',1.1); hold on;

legend('$U_1$ - common','$U_2$  - common','$U_1$  - private','$U_2$  - private',...
    'interpreter','latex','FontSize',13,'FontName','Times New Roman')
ylim([0 1]); hold on;
xlim([min(PTdB) max(PTdB)]); hold on;
xlabel('Transmit power [dBm]','FontSize',14,'FontName','Times New Roman');
ylabel('Capacity [bps/Hz]','FontSize',14,'FontName','Times New Roman');
