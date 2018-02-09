% appear in IEEE ICCE2018 
% simulation for the proposed CICG  compared with SICG of
%% 'anti-jamming transmission stackelberg game with observation errors' 2015xiaoliang IEEE communication letters

%just for different Cs, we compare the utility of SU, under different PI£¨PI=2P_max£©
%need to Note: we especially set P_max=J_max=PI/2, which is reasonable
%without interrupting PU's transmission from an optimization viewpoint, readers can refer to
%figure 1 of our paper IEEE ICCE2018 "Anti-Jamming CRN Transmission Stackelberg Game under Cumulative
%Interference Constraint"

clc
clear all;
close all;

%% the SICG in liang xiao CL2015: anti-jamming transmission stackelberg game with observation errors'
%%%%%%%%%%%%%%%%%%%%%%parameter setting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PI = [4 6 8];% maximum interfering power gain for PU
M=size(PI,2);
Cj= 0.1;%C_j:cost parameter of JA
epsilon=0;% without loss of generality, we set observation error =0
sigma = 0.1; % noise in paper (\sigma)

hs = 0.5;% h_s
hj = 0.5;% h_j
Cs=0:0.1:6;% trasmmistion power cost of SU
N=size(Cs,2);

P_powerMat= zeros(M,N);
P_powerMatErr=zeros(M,N);
J_powerMat =  zeros(M,N);
SU_utility=zeros(M,N);
SU_SINR=zeros(M,N);
JA_utility=zeros(M,N);
%% primary user is active

for ii=1:3 % three PI iteration
    P_max =PI(ii);
    J_max =PI(ii);
    M=0;
    for   jj=1:N % 19 Cs iteration
        M=M+1;
        %% leader game: secondary user
        %condition \prod_1 in(8)
        if Cs(jj)> hs/sigma;
            P_power=0;
            %condition \prod_2 in (8)
        elseif (  hs/(J_max * hj+sigma) <= Cs(jj) && Cs(jj) < hs/(2*sigma)...
                && Cj < 4 * hj * Cs(jj)^2 * P_max  )||...
                (  Cs(jj) < min(hs/sigma/2,hs/(J_max*hj+sigma))...
                && Cj < 4 * hj * Cs(jj)^2 * P_max/hs &&...
                Cj > 4 * P_max *hj * Cs(jj) / hs * ( hs/( P_max * hj + sigma ) -Cs(jj) )  );
            
            
            P_power=hs * Cj  / (4 * hj * Cs(jj)^2);
            
            %this step seems wrong in original paper, maybe it is a good choice to delete it 
%         elseif Cs(jj) < min(hs/sigma/2,hs/(J_max*hj+sigma))...
%                 && Cj < 4 * hj * Cs(jj)^2 * P_max/hs 
            
%             P_power=hs * Cj  / (4 * hj * Cs(jj)^2);
            
            %condition \prod_3 in (8)
        elseif Cs(jj)>=max(hs/2/sigma,hs/(J_max * hj + sigma)) &&...
                Cs(jj)<=hs/sigma && Cj<P_max * hs * hj/sigma/sigma;
            
            P_power=sigma^2 * Cj/hs/hj;
            
        elseif hs/2/sigma<Cs(jj) && Cs(jj) < hs/ (J_jammer *hj + sigma) &&...
                Cj >= P_max * hs * hj / (hs * sigma - sigma^2 * Cs(jj))...
                *(hs/(J_jammer *hj + sigma)-Cs(jj))&& Cj <P_max * hs * hj/sigma/sigma;
            
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
        P_powerMatErr(ii,M) = P_power-epsilon*P_power;% another observed transmit power
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
            (sigma+hj * J_power)-Cs(jj) * P_powerMat(ii,M);
        JA_utility(ii,M)=-hs * P_powerMat(ii,M) /...
            (sigma+hj * J_power)+Cs(jj) * P_powerMat(ii,M)-Cj * J_powerMat(ii,M);
        
    end
end
% figure(2)
% plot(Cs,SU_SINR(1,:),'k-d',Cs,SU_SINR(2,:),'g-o',Cs,SU_SINR(3,:),'b-*')
% xlabel('C_s');
% ylabel('SINR of the SU');
% legend('P_I=4','P_I=6','P_I=8');
% grid on;
% figure(3)
% plot(Cs,SU_utility(1,:),'k-d',Cs,SU_utility(2,:),'g-o',Cs,SU_utility(3,:),'b-*')
% xlabel('C_s');
% ylabel('Utility of the SU');
% legend('P_I=4','P_I=6','P_I=8');
% grid on;
% figure(4)
% plot(Cs,JA_utility(1,:),'k-d',Cs,JA_utility(2,:),'g-o',Cs,JA_utility(3,:),'b-*')
% xlabel('C_s');
% ylabel('Utility of the jammer');
% legend('P_I=4','P_I=6','P_I=8');
% grid on;

%% the proposed CICG
%%%%%%%%%%%%%%%%%%%%%%parameter setting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon=0;% we assume the observation error is zero
Cj= 0.1;%C_j:cost parameter of JA

sigma = 0.1; % noise in paper (\sigma)
% P_max =5;
% J_max =5;
PI=[4 6 8];% for clear comparison we denote PI(ii)=P_max+J_max, where P_max and J_max are P_s^M and P_j^M in CL2015xiaoliang, respectively
M=size(PI,2);
hs = 0.5;% h_s
hj = 0.5;% h_j
gs=1;
gj=1;
opt.gs=gs;opt.gj=gj;opt.hj=hj;opt.hs=hs;opt.sigma=sigma;opt.Cj=Cj;
% Cs=0:0.1:6;% cost of SU
% N=size(Cs,2);

P_powerMat= zeros(M,N);
P_powerMatErr=zeros(M,N);
J_powerMat =  zeros(M,N);
SU_utility2=zeros(M,N);
SU_SINR2=zeros(M,N);
JA_utility2=zeros(M,N);
P_o=zeros(1,M);
P_tilde=zeros(1,M);
%% primary user is active



for ii=1:3 % three PI(ii) iteration
    M=0;
    
    for   jj=1:N % 19 Cs iteration
        P_diamond=(-gj * (hs*hj/Cj)^.5 + (gj^2*hs*hj/Cj+4*gs*hj* ...
            (sigma*gj+hj*PI(ii)) )^.5)^2 / ( 4*gs^2*hj^2 );
        P_o(ii)= ( PI(ii) + sigma * gj - ( (hj*hs*gj*PI(ii)+hs*gj*gj*sigma)/Cs(jj)  )^.5 ) /gs/gj;
        P_tilde(ii)=hs*Cj/(4* hj * Cs(jj)^2);
        M=M+1;
        %% leader game: secondary user
        %condition \prod_0 in(7)
        if ( Cs(jj)> hs/sigma && P_diamond>=PI(ii)/gs)
            P_power=0;
        elseif       (Cs(jj)> hs && P_diamond<PI(ii)/gs && P_o(ii)>=PI(ii)/gs)
            P_power=0;
        elseif   (Cs(jj)> hs && P_diamond<PI(ii)/gs && P_o(ii)<PI(ii)/gs && hfunc(PI(ii)/gs,Cs(jj),opt,PI(ii))<=0);
            P_power=0;
            %condition \prod_1 in (7)
        elseif ( hs/sigma >=Cs(jj)>=hs/2/sigma && P_diamond>=PI(ii)/gs &&  Cj<hj*hs*PI(ii)/sigma/sigma/gs  )
            P_power=sigma^2*Cj/hs/hj;
        elseif ( hs/sigma >=Cs(jj)>=hs/2/sigma &&  P_diamond<PI(ii)/gs && P_o(ii)>=PI(ii)/gs  )
            P_power=sigma^2*Cj/hs/hj;
        elseif   ( hs/sigma >=Cs(jj)>=hs/2/sigma &&  P_diamond<PI(ii)/gs && P_o(ii)<PI(ii)/gs &&...
                lfunc(sigma^2*Cj/hs/hj,Cs(jj),opt) >=hfunc(PI(ii)/gs,Cs(jj),opt,PI(ii)) );
            P_power=sigma^2*Cj/hs/hj;
            
            
            %condition \prod_2 in (7)
        elseif P_diamond>=PI(ii)/gs && Cs(jj)< hs/2/sigma && cj<4*hj*Cs(jj)^2 * PI(ii)/hs/gs;
            P_power=P_tilde(ii);
        elseif P_diamond<PI(ii)/gs && P_o(ii)>=PI(ii)/gs && P_tilde(ii) <= P_diamond && Cs(jj)< hs/2/sigma;
            P_power=P_tilde(ii);
        elseif      (P_diamond<PI(ii)/gs && P_o(ii)<PI(ii)/gs &&...
                P_tilde(ii) <= P_diamond && Cs(jj)< hs/2/sigma &&...
                lfunc(P_tilde(ii),Cs(jj),opt)>=hfunc(PI(ii)/gs,Cs(jj),opt,PI(ii)) );
            P_power=P_tilde(ii);
            
            %condition \prod_3 in (7)
        elseif (Cs(jj)< hs/2/sigma && P_diamond <PI(ii)/gs &&...
                P_tilde(ii)>P_diamond && P_o(ii)>=PI(ii)/gs)
            P_power=P_diamond;
        elseif   (Cs(jj)< hs/2/sigma && P_diamond <PI(ii)/gs &&...
                P_tilde(ii)>P_diamond &&P_o(ii)<PI(ii)/gs &&...
                lfunc(P_diamond,Cs(jj),opt)>=hfunc(PI(ii)/gs,Cs(jj),opt,PI(ii)) );
            P_power=P_diamond;
            
        else
            P_power=PI(ii)/gs;
        end
        
        % this step is missed at first, which lead to utility of SU is
        % affected by this matrix because this is the power matrix of SU.
        % without this step, SU have to utilize the estimated
        % P_power, i.e., P_powerMatErr to compute the utility  by equation
        % (4) in the CL2015 paper.
        P_powerMat(ii,M) = P_power;
        %% follower game: jammer
        
        %         P_powerMatErr(ii,M) = P_power+epsilon*P_power;  %  observed transmit power
        P_powerMatErr(ii,M) = P_power-epsilon*P_power;% another observed transmit power
        if P_powerMatErr(ii,M) <= sigma^2 * Cj /hs/hj;
            J_power=0;
        elseif P_diamond>=P_powerMatErr(ii,M)>=sigma^2* Cj /hs/hj;
            J_power=( (hs*hj*P_powerMatErr(ii,M)/Cj)^0.5 - sigma )/hj;
        else
            J_power=( PI(ii)-gs*P_powerMatErr(ii,M) )/gj;
        end
        J_powerMat(ii,M)=J_power;
        SU_SINR2(ii,M)=hs * P_powerMat(ii,M) / (sigma+hj * J_power);
        SU_utility2(ii,M)=hs * P_powerMat(ii,M) /...
            (sigma+hj * J_power)-Cs(jj) * P_powerMat(ii,M);
        JA_utility2(ii,M)=-hs * P_powerMat(ii,M) /...
            (sigma+hj * J_power)+Cs(jj) * P_powerMat(ii,M)-Cj * J_powerMat(ii,M);
        
    end
end

figure(3)
plot(Cs,SU_utility(1,:),'k-.d',Cs,SU_utility(3,:),'b-.*');
hold on;
plot(Cs,SU_utility2(1,:),'k-o',Cs,SU_utility2(3,:),'b-p')
xlabel('C_s');
ylabel('Utility of the SU');
legend('SICG with P_I=4','SICG with P_I=8',...
    'CICG with P_I=4','CICG with P_I=8');
grid on;





