function[]=Plot_Results(Flag, Par, AstVar, SynVar, freq)
% plots ESP, DSE,  and PR of the first synapse of each neuron
for i = 1:Par.Num_Neurons_OP
    strng = int2str(i);
    strng = strcat('Neuron',strng);
    figure
    subplot(4,1,1)
    plot(Par.t,AstVar.EspStore(1,:))
    ylabel('Esp (%)')
    title(strng)
    subplot(4,1,2)
    plot(Par.t, SynVar.DseStore((10*i)-9,:))
    ylabel('Dse (%)')
    subplot(4,1,3)
    plot(Par.t, SynVar.PrStore((10*i)-9,:))
    ylabel('Pr')
    subplot(4,1,4)
    plot(Par.t,freq(i,:))
    xlabel('Time ( s )')
    ylabel('Freq (Hz)')
    sname = strcat('Results\',Par.ExpName, '\',strng);
    saveas(gcf,sname)
end
if Flag.Astro
    %plots the total calcium of the Astrocyte
    figure
    plot(Par.t,AstVar.CA2TotalStore)
    ylabel('Ca^2^+, IP_3 (\muM)')
    xlabel('Time (s)')
    % sname = strcat(STR1,CaModeSTR,STR3);
    % saveas(gcf,sname)
    sname = strcat('Results\',Par.ExpName,'\TotalCalcium');
    saveas(gcf,sname)
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
    sname = strcat('Results\',Par.ExpName,'\Pr_Synapses',strng1,'-',strng2);
    saveas(gcf,sname)
    counter = counter+5;
end

% if the No Fault experiment is run then this produces Fig. 3 of Wade et
% al. 2012.
if strcmp('No_Fault', Par.ExpName);
    PlotColour = [0 0 1; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 0 1 0; 1 1 1 ]; % sets up a coloor array so that each plot of the subplot is a different colour
    figure
    subplot(4,1,1) %Plot ESP
    for x = 1:Par.Num_Neurons_OP
        plot(Par.t,AstVar.EspStore(1,:),'color',PlotColour(x,:))
        ylabel('Esp (%)')
        hold on
    end
    subplot(4,1,2) % plot DSE of the first synapse of each OP neuron
    for x = 1:Par.Num_Neurons_OP
        plot(Par.t, SynVar.DseStore((10*x)-9,:),'color',PlotColour(x,:))
        ylabel('Dse (%)')
        hold on
    end
    subplot(4,1,3)
    for x = 1:Par.Num_Neurons_OP % plot the probability of the first synapse of each OP neuron
        plot(Par.t, SynVar.PrStore((10*x)-9,:),'color',PlotColour(x,:))
        ylabel('Pr')
        hold on
    end
    subplot(4,1,4)
    for x = 1:Par.Num_Neurons_OP % plot the average frequency of each OP neuron
        plot(Par.t,freq(x,:),'color',PlotColour(x,:))
        xlabel('Time ( s )')
        ylabel('Freq (Hz)')
        hold on
    end
    sname = strcat('Results\',Par.ExpName, '\No_Fault');
    saveas(gcf,sname)
end
