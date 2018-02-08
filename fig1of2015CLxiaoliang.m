% simulation for cl2015
%% 'anti-jamming transmission stackelberg game with observation errors'
% liangxiao2015
%fig.1(a) (b) (c)

clc
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%parameter setting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cs = [0.4 0.6 0.9];%C_s:cost parameter of SU
M=size(Cs,2);
Cj= 0.2;%C_j:cost parameter of JA

sigma = 0.1; % noise in paper (\sigma)
P_max =5;
J_max =5;
hs = 0.5;% h_s
hj = 0.5;% h_j
epsilon=0:0.05:0.9;% observation error of the jammer
N=size(epsilon,2);

P_powerMat= zeros(M,N);
P_powerMatErr=zeros(M,N);
J_powerMat =  zeros(M,N);
SU_utility=zeros(M,N);
SU_SINR=zeros(M,N);
JA_utility=zeros(M,N);
%% primary user is active

for ii=1:3 % three Cs iteration
    M=0;
    for   jj=1:N % 19 epsilon iteration
        M=M+1;
        %% leader game: secondary user
        %condition \prod_1 in(8)
        if Cs(ii)> hs/sigma;
            P_power=0;
            %condition \prod_2 in (8)
        elseif (  hs/(J_max * hj+sigma) <= Cs(ii) && Cs(ii) < hs/(2*sigma)...
                && Cj < 4 * hj * Cs(ii)^2 * P_max  )||...
                (  Cs(ii) < min(hs/sigma/2,hs/(J_max*hj+sigma))...
                && Cj < 4 * hj * Cs(ii)^2 * P_max/hs &&...
                Cj > 4 * P_max *hj * Cs(ii) / hs * ( hs/( P_max * hj + sigma ) -Cs(ii) )  );
            
            
            P_power=hs * Cj  / (4 * hj * Cs(ii)^2);
            
            %this step seems wrong in original paper, maybe it is a good choice to delete it 
%         elseif Cs(ii) < min(hs/sigma/2,hs/(J_max*hj+sigma))...
%                 && Cj < 4 * hj * Cs(ii)^2 * P_max/hs 
            
%             P_power=hs * Cj  / (4 * hj * Cs(ii)^2);
            
            %condition \prod_3 in (8)
        elseif Cs(ii)>=max(hs/2/sigma,hs/(J_max * hj + sigma)) &&...
                Cs(ii)<=hs/sigma && Cj<P_max * hs * hj/sigma/sigma;
            
            P_power=sigma^2 * Cj/hs/hj;
            
        elseif hs/2/sigma<Cs(ii) && Cs(ii) < hs/ (J_jammer *hj + sigma) &&...
                Cj >= P_max * hs * hj / (hs * sigma - sigma^2 * Cs(ii))...
                *(hs/(J_jammer *hj + sigma)-Cs(ii))&& Cj <P_max * hs * hj/sigma/sigma;
            
            P_power=sigma^2 * Cj/hs/hj;
            
        else
            P_power=P_max;
        end
        
        % this step is missed at first, which lead to utility of SU is
        % affected by this matrix because this is the power matrix of SU.
        % without this step, SU have to utilize the estimated
        % P_power, i.e., P_powerMatErr to compute the utility  by equation
        % (4) in the CL2015 paper.
        P_powerMat(ii,M) = P_power;
        %% follower game: jammer

%         P_powerMatErr(ii,M) = P_power+epsilon*P_power;  %  observed transmit power
        P_powerMatErr(ii,M) = P_power-epsilon(jj)*P_power;% another observed transmit power
        if P_powerMatErr(ii,M) <= sigma^2 * Cj /hs/hj;
            J_power=0;
        elseif P_powerMatErr(ii,M)>=Cj * (J_max * hj +sigma)^2/hs/hj;
            J_power=J_max;
        else
            J_power=( (hs*hj*P_powerMatErr(ii,M)/Cj)^0.5 - sigma )/hj;
        end
        J_powerMat(ii,M)=J_power;
        SU_SINR(ii,M)=hs * P_powerMat(ii,M) / (sigma+hj * J_power);
        SU_utility(ii,M)=hs * P_powerMat(ii,M) /...
            (sigma+hj * J_power)-Cs(ii) * P_powerMat(ii,M);
        JA_utility(ii,M)=-hs * P_powerMat(ii,M) /...
            (sigma+hj * J_power)+Cs(ii) * P_powerMat(ii,M)-Cj * J_powerMat(ii,M);
        
    end
end
figure(2)
plot(epsilon,SU_SINR(1,:),'k-d',epsilon,SU_SINR(2,:),'g-o',epsilon,SU_SINR(3,:),'b-*')
xlabel('\epsilon');
ylabel('SINR of the SU');
legend('C_s=0.4','C_s=0.6','C_s=0.9');
grid on;
figure(3)
plot(epsilon,SU_utility(1,:),'k-d',epsilon,SU_utility(2,:),'g-o',epsilon,SU_utility(3,:),'b-*')
xlabel('\epsilon');
ylabel('Utility of the SU');
legend('C_s=0.4','C_s=0.6','C_s=0.9');
grid on;
figure(4)
plot(epsilon,JA_utility(1,:),'k-d',epsilon,JA_utility(2,:),'g-o',epsilon,JA_utility(3,:),'b-*')
xlabel('\epsilon');
ylabel('Utility of the jammer');
legend('C_s=0.4','C_s=0.6','C_s=0.9');
grid on;



