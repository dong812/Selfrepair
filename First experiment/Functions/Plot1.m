function[]=Plot1(Flag, Par, AstVar, SynVar, freq)
% plots ESP, DSE,  and PR of the first synapse of each neuron

figure
subplot(5,1,1)
plot(Par.t,AstVar.EspStore(11,:))
ylabel('eSP (%)')
subplot(5,1,2)
plot(Par.t, SynVar.DseStore(11,:))
ylabel('DSE (%)')
subplot(5,1,3)
plot(Par.t, SynVar.PrStore(11,:))
ylabel('PR')
subplot(5,1,4)
plot(Par.t, SynVar.PrStore(20,:))
ylabel('PR')
subplot(5,1,5)
plot(Par.t,freq(2,:))
xlabel('Time ( s )')
ylabel('Freq (Hz)')
    

    
    

    
if Flag.Astro
    %plots the total calcium of the Astrocyte
    figure
    subplot(2,1,1)
    plot(Par.t,AstVar.CA2Store(11,:),Par.t,0.2)
    ylabel('Ca^2^+(\muM)')
    xlabel('Time (s)')
    subplot(2,1,2)
    plot(Par.t,AstVar.CA2Store(20,:),Par.t,0.2)
    ylabel('Ca^2^+(\muM)')
    xlabel('Time (s)')
    % sname = strcat(STR1,CaModeSTR,STR3);
    % saveas(gcf,sname)
   
end

%Plots the probability of all OP neuron synapses
z = size(SynVar.PrStore,1); %Works out the total number of synapses in the OP layer
counter = 0; %counter used to monitor the number of synapses plotted
while counter < z
    figure
    for i = 1:5
        subplot(5,1,i)
        plot(Par.t, SynVar.PrStore(i+counter,:))
        if i == 1
            strng1 = int2str(counter+1);
            strng2 = int2str(counter+5);
            name = strcat('Synapses',strng1,'-',strng2);
            title(name)
        end
        ylabel('PR')
        if i == 5
            xlabel('Time (ms)');
        end
    end
    strng1 = int2str(counter+1);
    strng2 = int2str(counter+5);
    sname = strcat('Results\',Par.ExpName,'\PR_Synapses',strng1,'-',strng2);
    saveas(gcf,sname)
    counter = counter+5;
end

% if the No Fault experiment is run then this produces Fig. 3 of Wade et
% al. 2012.

