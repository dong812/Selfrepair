

% function [] =  Frontiers_Partial_Fault()
clc
clear
close all
tic
addpath ('Functions')
addpath ('Setup_Files')
% for loop = 1:1
%Set the mode of Ca2+ modulation (AM = 1, FM = 2, AM-FM = 3)
CaMode = 1;
% I = 5250;
I = 300;
SIC = 50;

ATP = 0;

switch CaMode
    case 1
        CaModeSTR = 'AM_Mode';
    case 2
        CaModeSTR = 'FM_Mode';
    case 3
        CaModeSTR = 'AM-FM_Mode';
end
deltaCA2Total = 0;  % Change in total Calcium
deltaAGtot = 0;     % Change in total AG
deltaIP3 = 0;       % Change in IP3
deltaDse = 0;       % Change in DSE

[Par, AstPar,AstVar, SynPar, SynVar, NeuPar, NeuVar, Flag, Silence, Input]= Frontiers_Partial_Fault_Setup(CaMode);%Sets up all variables, parameters etc.
[Input_Spikes] = Self_Repair_Input_Generation(Par,Input, Silence, Flag);%Generates input to the network.

% if a fault is to be simulated then set the time at which it is to occur
if Flag.Fault == 1
    Fault_Time = Par.Fault_Time;
end


for j = 2:Par.T
    % if a fault at current time step then set the probability of the
    % faulty syanpses and also change the initial synaptic probabilities as
    % this is used in the calculation of of synaptic activity
    if Flag.Fault && j == Fault_Time
        ATP = 10;
        for faultsyns = Par.FaultSyns(1):Par.FaultSyns(2)
        SynVar.PrStore(faultsyns,j) = SynVar.Fault(faultsyns);
        SynVar.PrInit(faultsyns) = SynVar.Fault(faultsyns);
        AstVar.IP3ATPSpikes(faultsyns,1) = 1;
        end
    end
    %1 : Number of Astrocytes
    for k = 1:Par.Num_Astro_OP
        deltaCA2Total = 0; %Temporary storage for total change in Calcium;
        for l = 1:size(Par.Astro_map,2)
            % gets each neuron number that the astrocyte is connected to
            n = Par.Astro_map(k,l);
            %Calculates the range of input synapses that correspond to that
            %neuron
            beg_Syn = (n * Par.Num_Inputs_OP) - (Par.Num_Inputs_OP - 1);
            end_Syn = (n * Par.Num_Inputs_OP);
            %% Self Repair Module
%             if Flag.SelfRepair
                deltaGlu = Par.Ts.*((-AstVar.Glu(k)/AstPar.tau_Glu) + (AstVar.r_Glu.*AstVar.CA2_Spike(k,j-1)));
                AstVar.Glu(k) = AstVar.Glu(k) + deltaGlu;
                AstVar.GluStore(k,j) = AstVar.Glu; %store new glutamate value
                
                                             
               
                for h = beg_Syn : end_Syn
                    AG_Release = NeuVar.Spikes(n,j-1)>0; %Will be 1 if there is a spike
                    %Calculate change in 2-AG dependant on ouput spikes from
                    %the post synaptic neuron.(equation 1, Wade et al 2012)
                    deltaAG = Par.Ts .*((-SynVar.AG(h)/SynPar.tau_AG) + ((SynPar.r_AG.*AG_Release)));
                    SynVar.AG(h) = SynVar.AG(h) + deltaAG;
                    SynVar.AGStore(h,j) = SynVar.AG(h); % Stores AG levels
                    
                    deltalocalGlu = Par.Ts.*(-AstVar.localGluStore(h,j-1)/AstPar.tau_Glu + (AstVar.r_Glu.*AstVar.localCA2_Spike(h,j-1)));
                    AstVar.localGluStore(h,j) = AstVar.localGluStore(h,j-1) + deltalocalGlu;
                    if Flag.Esp
                 %Change in Glutamate (equation 14, Wade et al 2012)
                       

                      %Calculate Esp which results from Glu release of astrocyte
                      %as a result of Astrocytic Ca2+ transients (equation 15, Wade et al 2012)
                      deltaEsp = Par.Ts .* ((-AstVar.Esp(h)/AstPar.tau_Esp) + ((AstVar.Esp_weight.*(AstVar.Glu(k))/AstPar.tau_Esp)));
                      AstVar.Esp(h) = AstVar.Esp(h) + deltaEsp;
                      AstVar.EspStore(h,j) = AstVar.Esp(h); % Store new Esp value
                    end

                    
                    if Flag.Dse
                        % Change in DSE (equation 13, Wade et al 2012)
                        deltaDse = (deltaAG .* SynPar.Kag) * -1;
                        SynVar.Dse(h) = SynVar.Dse(h) + deltaDse;
                        SynVar.DseStore(h,j) = SynVar.Dse(h);
                    end
                    % Calculate synaptic probability increase due to Esp
                    PrIncrease = (SynVar.PrInit(h)/100).*AstVar.EspStore(h,j-1);
                    % Calculate synaptic probability decrease due to DSE
                    PrDecrease = (SynVar.PrInit(h)/100).*SynVar.DseStore(h,j-1);
                    % Change in synaptic probability
                    deltaPr = PrIncrease + PrDecrease;
                    % Store synaptic probability
                    SynVar.PrStore(h,j) = SynVar.PrInit(h) + deltaPr;
                end
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Individual synapse calculations
            for h = beg_Syn : end_Syn
                Spike = Input_Spikes(h,j-1)>0; % Will be 1 if there is a spike
                %SynVar.Itot(n,j) = SynVar.Itot(n,j) + I.*Spike.*SynVar.PrStore(h,j-1);
            % Probabilistic Synapse
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if Spike
                    RandNum = rand; % Generate random number between 0 and 1
                    % if > synaptic probability current is injected into
                    % neuron
                    if RandNum <= SynVar.PrStore(h,j-1)
                        SynVar.SynapseOut(h,j) = 1; %Stores the activity of the synapse for that time step
                        SynVar.Itot(n,j) = SynVar.Itot(n,j) + I ; % adds current to the total neuronal current
                    end
                end
                SynVar.Itot(n,j) = SynVar.Itot(n,j) +SIC.*AstVar.localCA2_Spike(h,j-1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Li-Rinzel Calcium Dynamics
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (Flag.Astro ==1)

                    deltaATPIP3 = Par.Ts.*(-AstVar.IP3ATPStore(h,j-1)/AstPar.tau_IP3ATP+AstPar.r_IP3ATP.*ATP.*AstVar.IP3ATPSpikes(h,1));
                    AstVar.IP3ATPStore(h,j) = AstVar.IP3ATPStore(h,j-1) + deltaATPIP3;
                    IP3 = AstVar.IP3Store(h,j-1) + AstVar.IP3ATPStore(h,j-1) ;%IP3 at previous Ts
                    %Calculate m_inf for the current time step (equation 9, Wade et al 2012)
                    m_inf = IP3 / (IP3 + AstPar.d1);
                    %Calculate n_inf, (equation 10, Wade et al 2012)
                    n_inf = AstVar.CA2(h) / (AstVar.CA2(h) + AstPar.d5);
                    % Calculate Jpump (equation 12, Wade et al 2012)
                    Jpump = (AstPar.vER .* AstVar.CA2(h)^2)/(AstPar.KER^2 + AstVar.CA2(h)^2);
                    %Calculate Jleak (equation 11, Wade et al 2012)
                    Jleak = AstPar.v2.*(AstPar.c0 - ((1+AstPar.c1).*AstVar.CA2(h)));
                    %Calculate JChan (equation 8, Wade et al 2012)
                    Jchan = AstPar.v1 .* m_inf^3 .* n_inf^3 .* AstVar.h(h)^3 .* (AstPar.c0-(1+AstPar.c1).*AstVar.CA2(h));
                    %Calculate the change in CA2 (equation 3, Wade et al 2012)
                    deltaCA2 = Par.Ts .* (Jchan + Jleak - Jpump );% Calculates the cange in Ca2+;
                    deltaCA2Total =  deltaCA2Total + deltaCA2; %total change in calcium for this astrocyte
                    %Calculate Q2 (equation 7, Wade et al 2012)
                    Q2 = AstPar.d2 .* ((IP3 + AstPar.d1)/(IP3 + AstPar.d3));
                    %calculate  h_inf (equation 5, Wade et al 2012)
                    h_inf = Q2/(Q2+AstVar.CA2(h));
                    %(equation 6, Wade et al 2012)
                    tauh = 1/(AstPar.a2.*(Q2+AstVar.CA2(h)));
                    %Calcul;ate change in h (equation 4, Wade et al 2012)
                    deltah = Par.Ts .* ((h_inf - AstVar.h(h))/ tauh);
                    %Calculate change in IP3 (equation 2, Wade et al 2012)
                    deltaIP3 = Par.Ts .*(((AstPar.IP_star_3 - AstVar.IP3(h))/AstPar.tau_IP3) + (AstPar.r_IP3 .* SynVar.AG(h)));
                    %Upadates CA2+, h, IP3, IP3Store, AG and AGStore for this time step
                    AstVar.CA2(h) = AstVar.CA2(h) + deltaCA2;
                    AstVar.CA2Store(h,j) = AstVar.CA2(h);

                    if (AstVar.CA2(h) > AstPar.localCA2Thresh) && (AstVar.localCA2_Spike_Timer(h)  == AstVar.CA2_Spike_Timer_Reset)
                         AstVar.localCA2_Spike(h,j) = 1;
                         AstVar.localCA2_Spike_Timer(h) =  AstVar.localCA2_Spike_Timer(h) -1;
                    elseif (AstVar.CA2(h) > AstPar.localCA2Thresh) && (AstVar.localCA2_Spike_Timer(h) < AstVar.CA2_Spike_Timer_Reset)
                           AstVar.localCA2_Spike_Timer(h) = AstVar.localCA2_Spike_Timer(h) - 1;
                    end
                    if (AstVar.CA2(h) < AstPar.localCA2Thresh) || (AstVar.localCA2_Spike_Timer(h) == 0)
                         AstVar.localCA2_Spike_Timer(h) = AstVar.CA2_Spike_Timer_Reset;
                    end



                    
                    AstVar.h(h) = AstVar.h(h) + deltah;
                    AstVar.IP3(h) = AstVar.IP3(h) + deltaIP3;
                    AstVar.IP3Store(h,j) = AstVar.IP3(h);
                    AstVar.CA2TotalStore(k,j) = AstVar.CA2TotalStore(k,j-1) + deltaCA2Total;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end %if (Flag.Astro ==1)

                %Updates Variables for next time step and storage
                % Variable storage for plotting, only active if flag is set
%                 AstVar.CA2TotalStore(k,j) = AstVar.CA2TotalStore(k,j-1) + deltaCA2Total;
            end %end for h = beg_syn : end_syn

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% LIF Neuron Model
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %check if neuron fired in the previous time step
            if NeuVar.v(n) > NeuPar.Vth;
                NeuVar.Spikes(n,j) = 1; %records spike
                NeuVar.Refract(n) = NeuPar.RefPeriod; % puts neuron into a refractory period
                NeuVar.v(n) = NeuPar.Reset; % resets neuronal voltage
                %if not then checks if neuron is in refractory period
            elseif NeuVar.Refract(n) > 0
                NeuVar.Refract(n) = NeuVar.Refract(n) - Par.Ts;%subtracts 1 time step from Refractory period
                NeuVar.v(n) = NeuPar.Reset; %Sets mem voltage for this time step to 0
                %if neither of these then it calculates the new membrane voltage
            else
                % equation 16, Wade et al. 2012
                deltav = Par.Ts .*((-1 .* (NeuVar.v(n)/NeuPar.tau_mem)) + ((NeuPar.R_mem .* SynVar.Itot(n,j-1))/NeuPar.tau_mem));
                NeuVar.v(n) = NeuVar.v(n) + deltav;
            end

            % records membrane voltage if needed for plot
            if Flag.vPlot == 1
                %Stores Membrane voltage of neuron n for this time step
                NeuVar.vStore(n,j) = NeuVar.v(n);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end %end for l = 1:size(Par.Astro_map,2)

        if Flag.Astro ==1
            % if total CA > Ca threshold and has not occured in the time
            % alotted between calcium spikes then record the time of the CA
            % spike and decrement the timee
             if (AstVar.CA2TotalStore(1,j) > AstPar.CA2Thresh) && (AstVar.CA2_Spike_Timer == AstVar.CA2_Spike_Timer_Reset)
                AstVar.CA2_Spike(k,j) = 1;
                AstVar.CA2_Spike_Timer =  AstVar.CA2_Spike_Timer -1;
             % else if if total CA > Ca threshold and the timer is < reset
             % then decrement the timer.
             elseif (AstVar.CA2TotalStore(1,j) > AstPar.CA2Thresh) && (AstVar.CA2_Spike_Timer < AstVar.CA2_Spike_Timer_Reset)
                 AstVar.CA2_Spike_Timer = AstVar.CA2_Spike_Timer - 1;
             end
             
             % if total CA > Ca threshold or the timer = 0 set the timer to
             % the reset value.
             if (AstVar.CA2TotalStore(1,j) < AstPar.CA2Thresh) || (AstVar.CA2_Spike_Timer == 0)
                 AstVar.CA2_Spike_Timer = AstVar.CA2_Spike_Timer_Reset;
             end
        end %end if Flag.Astro ==1
    end %end for k = 1:Par.NumAstro_IP
end %end for j = 2:Par.T

%%

%calculate the frequency of each neuron at each time step
%The frequency is calculated by counting the number of spikes fired by the
%neuron over a 1 second window.
freq = zeros(Par.Num_Neurons_OP,Par.T);
for i = 1:Par.Num_Neurons_OP
    for x = 1000 : Par.T
        freq(i,x) =  sum(NeuVar.Spikes(i,x-999:x));
    end
end

figure
plot(Par.t,freq(2,:))
xlabel('Time ( s )')
ylabel('Freq (Hz)')


figure
plot(Par.t,SIC.*AstVar.localCA2_Spike(11,:))
xlabel('Time ( s )')
ylabel('SIC(pA)')



%Plots and saves the results
Plot_Results(Flag, Par, AstVar, SynVar, freq)
%Timer(toc);
%disp('Simulation Complete')