clear all;
clc;

Nt = 4;
Nr = 2;
fc = 28*10^9;
d = 50;
P_tx_dB = [0:0.5:10]; %transmit power(dB)
P_tx = 10.^(P_tx_dB/10); %transmit power(W)
P_noise = 1; %noise power
N_beam_t = 513; %transmitter에서 고를 수 있는 beam 개수
N_beam_r = 513; %receiver에서 고를 수 있는 beam 개수
    
snr_table = [];
snr_table_raw = [];
mean_count = [];
for signal_power_ind = 1:numel(P_tx)
    s = P_tx(signal_power_ind); %pilot signal = 1*transmit power
   
    snr = 0;
    snr_raw = 0;
    count_table = [];
    threshold = 8*s;
    
    for mc = 1:1000 %이만큼 반복에서 평균을 낼 예정이다
        
        count = 0;
        rho = randn + 1i*randn;
        tau = rand * 20*1e-9;
        theta_t = -pi + rand*2*pi;
        theta_r = -pi + rand*2*pi;
        a_t = [];
        a_r = [];

        for k = 1:Nt
            a_t = [a_t; exp(1i*pi*(k-1)*sin(theta_t))]; %aod
        end
        for k= 1:Nr
            a_r = [a_r; exp(1i*pi*(k-1)*sin(theta_r))]; %aoa
        end
 
        H = rho*exp(-1i*2*pi*fc*tau)*a_r*a_t'; %channel matrix
        
        
        phi = [];
        for k = 1:N_beam_t
            phi = [phi, -pi/2+(k-1)*pi/(N_beam_t-1)];
        end

        F = [];
        for k = 1:Nt
            F = [F; 1/sqrt(Nt)*exp(1i*pi*(k-1)*sin(phi))];
        end

        r = [];
        r_norm = [];
        r_max = 0;
        
        center = (N_beam_t - 1) / 2;
        range = (N_beam_t - 1) / 2;
        for k = 1:log2(N_beam_t - 1)
            count = count + 1;
            noise = sqrt(P_noise) * (randn(Nr,1) + 1i*randn(Nr,1));
            left = ceil(center - range/2);
            right = ceil(center + range/2);
            range = ceil(range/2);
            left_norm = norm(H*F(:,left)*s + noise);
            right_norm = norm(H*F(:,right)*s + noise);
            
            if left_norm >= threshold
                r_max= left_norm;
                localmax_t = left;
                break;
            elseif right_norm >= threshold
                r_max = right_norm;
                localmax_t = right;
                break;
            elseif left_norm >= 0.65*threshold
                range = ceil(range/8);
            elseif right_norm >= 0.65*threshold
                range = ceil(range/8);
            end
            
            if left_norm > r_max
                r_max = left_norm;
                localmax_t = left;
                center = left;
            end
            if right_norm > r_max
                r_max = right_norm;
                localmax_t = right;
                center = right;
            end
            if range <= 1
                break;
            end
        end
% Exhaustive search   
%         for k = 1:N_beam_t
%             noise = sqrt(P_noise) * (randn(Nr,1) + 1i*randn(Nr,1));
%             r = [r, H*F(:,k)*s + noise];
%             rk_norm = norm(r(:,k)) ;
%             r_norm = [r_norm, rk_norm];      
%             if rk_norm > r_max
%                 r_max = rk_norm;
%                 localmax_t = k;
%             end
%         end

        psi = [];
        for k = 1:N_beam_r
            psi = [psi, -pi/2+(k-1)*pi/(N_beam_r-1)];
        end

        W = [];
        for k = 1:Nr
            W = [W; 1/sqrt(Nr)*exp(1i*pi*(k-1)*sin(psi))];
        end


        y = [];
        y_norm = [];
        y_max = 0;
        
        center = (N_beam_r - 1) / 2;
        range = (N_beam_r - 1) / 2;
        for k = 1:log2(N_beam_r - 1)
            count = count + 1;
            noise = sqrt(P_noise) * (randn(Nr,1) + 1i*randn(Nr,1));
            left = ceil(center - range/2);
            right = ceil(center + range/2);
            range = ceil(range/2);
            left_norm = norm(W(:,left)'*H*F(:,localmax_t)*s + W(:,left)'*noise);
            right_norm = norm(W(:,left)'*H*F(:,localmax_t)*s + W(:,left)'*noise);
            
            if left_norm >= threshold
                y_max= left_norm;
                localmax_r = left;
                break;
            elseif right_norm >= threshold
                y_max = right_norm;
                localmax_r = right;
                break;
            elseif left_norm >= 0.65*threshold
                range = ceil(range/8);
            elseif right_norm >= 0.65*threshold
                range = ceil(range/8);
            end
            
            if left_norm > y_max
                y_max = left_norm;
                localmax_r = left;
                center = left;
            end
            if right_norm > y_max
                y_max = right_norm;
                localmax_r = right;
                center = right;
            end
            if range <= 1
                break;
            end
        end
        
% Exhaustive search
%         for k = 1:N_beam_r
%             noise = sqrt(P_noise) * (randn(Nr,1) + 1i*randn(Nr,1));
%             y = [y, W(:,k)'*H*F(:,localmax_t)*s + W(:,k)'*noise];
%             yk_norm = norm(y(:,k));
%             y_norm = [y_norm, yk_norm];
%             if yk_norm > y_max
%                 y_max = yk_norm;
%                 localmax_r = k;
%             end
%         end     

        snr = snr + P_tx(signal_power_ind) * norm(W(:,localmax_r)'*H*F(:,localmax_t))^2 / P_noise^2;
        snr_raw = snr_raw + P_tx(signal_power_ind) * norm(W(:,1)'*H*F(:,1))^2 / P_noise^2;
        count_table = [count_table, count];
    end  %MC end 
    mean_count = [mean_count, mean(count_table)];
    snr = snr / mc;
    snr_raw = snr_raw / mc;
    snr_table = [snr_table, log2(1+(snr))];
    snr_table_raw = [snr_table_raw, log2(1+(snr_raw))];
end
mean = mean(mean_count)
plot(P_tx_dB, snr_table);
hold on
plot(P_tx_dB, snr_table_raw);
xlabel('transmit power(dB)');
ylabel('SNR');
legend('SNR, trained', 'SNR, untrained');
hold on

% snr = [];
% snr_raw = [];
% for ind = 1:numel(P_tx)
%   
%     rho = randn + 1i*randn;
%     tau = rand * 20*1e-9;
%     theta_t = -pi + rand*2*pi;
%     theta_r = -pi + rand*2*pi;
%     a_t = [];
%     a_r = [];
% 
%     for k = 1:Nt
%         a_t = [a_t; exp(1i*pi*(k-1)*sin(theta_t))];
%     end
%     for k= 1:Nr
%         a_r = [a_r; exp(1i*pi*(k-1)*sin(theta_r))];
%     end
% 
%     H = rho*exp(-1i*2*pi*fc*tau)*a_r*a_t';
%     
%     F = [];
%     phi = -pi/2+(beam_t(ind)-1)*pi/(N_beam_t-1);
%     for k = 1:Nt
%         F = [F; 1/sqrt(Nt)*exp(1i*pi*(k-1)*sin(phi))];
%     end
%     
%     W = [];
%     psi = -pi/2+(beam_r(ind)-1)*pi/(N_beam_r-1);
%     for k = 1:Nr
%         W = [W; 1/sqrt(Nr)*exp(1i*pi*(k-1)*sin(psi))];
%     end
%     
%     F_raw = 1/sqrt(Nt)*ones(Nt, 1);
%     W_raw = 1/sqrt(Nr)*ones(Nr, 1);
%     
%     snr = [snr, P_tx(ind)*norm(W'*H*F)^2/(P_noise^2)];
%     snr_raw = [snr_raw, P_tx(ind)*norm(W_raw'*H*F_raw)^2/(P_noise^2)];
% end
% plot(P_tx_dB, snr)
% plot(P_tx_dB, log2(snr));
%hold on
%legend('snr','snr raw')

% number = 1000;
% data = sign(-0.5+rand(1, number));
% signal = P_tx * data;
% correct = 0;
% correct_BT = 0;
% for k = 1:number
%     noise = sqrt(P_noise) * (randn(Nr,1) + 1i*randn(Nr,1));
%     rcv = H*repmat(signal(k),Nt,1) + noise;
%     rcv_BT = W(:,argmax_r)'*H*F(:,argmax_t)*signal(k) + W(:,argmax_r)'*noise;
%     if sign(real(mean(rcv))) == data(k)
%         correct = correct + 1;
%     end
%     if sign(real(mean(rcv_BT))) == data(k)
%         correct_BT = correct_BT + 1;
%     end
% end
% correct_rate = correct / number;
% correct_rate_BT = correct_BT / number;
