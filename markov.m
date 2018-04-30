function []=markov(run_sim,state_num,run)

tic
close all

if run_sim==1
    
    % setting up parameters of system to be simulated
    N=100;
    x=linspace(-4*pi,4*pi,N);
    y1=-2*sin((x-pi)/2);
    y2=x/(x(end)-x(1));
    y=y1+y2; % consists of symmetric y1 plus an assymetric y2 term
    y=y-min(y); % shift the potential so the minimum value is zero
    A=10;
    KbT=0.596;
    % construct rate matrix using Arhennius rates
    for i=1:N-1
        K(i,i+1)=A*exp((y(i+1)-y(i))/2/KbT);
        K(i+1,i)=A*exp((y(i)-y(i+1))/2/KbT);
    end
    for i=1:N
        K(i,i)=0;
        K(i,i)=-sum(K(:,i));
    end
    [~,d]=eig(K);
    [k,~]=sort(diag(d),'descend');
    rel_exact=-1/k(2); % exact relaxation time

    %% Simulation of Markov chain system (using Gillespie gives almost identical results)
    steps=40000;
    M=expm(K'*0.25); %timestep of 0.25
    s_traj=zeros(steps,N);
    for j=1:N
        s_traj(1,j)=j; %initialise a trajectory from each state
        
        for i=1:steps
            rand_num=rand(1);
            for k=1:N
                if rand_num<sum(squeeze(M(s_traj(i,j),1:k)))
                    s_traj(i+1,j)=k;
                    break                
                end
            end
        end
    end
    save(['trajectories' num2str(state_num) '_' num2str(run) '_shorter_res_asym.mat'])
    

end
if run_sim==0
    state_num=3;
    load(['trajectories' num2str(state_num) '_' num2str(run) '_shorter_res_asym.mat'])
end

%% 2 state coarse grain of trajectory
if state_num==2
    for k=1:N
        for i=1:steps % we have NN transitions and we have divided the timescale in 100*NN steps
            if s_traj(i,k)<(0.51*N)
                state(i,k)=1;
            else
                state(i,k)=2;
            end
        end
    end
    state=state';
    dlmwrite(['coor_bin_analytic_long_' num2str(state_num) '_' num2str(run)  '.txt'],state,'delimiter','\t','precision',6)
    
end
%% 3 state coarse grain of trajectory
if state_num==3
    for k=1:N
        for i=1:steps % we have NN transitions and we have divided the timescale in 100*NN steps
            if s_traj(i,k)<0.48*N
                state(i,k)=1;
            elseif s_traj(i,k)>=0.48*N && s_traj(i,k)<0.54*N
                state(i,k)=2;
            else
                state(i,k)=3;
            end
        end
    end
    state=state';
    dlmwrite(['coor_bin_analytic_long_' num2str(state_num) '_' num2str(run)  'asym_downhill.txt'],state,'delimiter','\t','precision',6)
     
end

toc

end