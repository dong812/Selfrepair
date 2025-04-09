function [Input_Spikes, cellSpikes] = Self_Repair_Input_Generation(Par,Input, Silence, Flag)
% Generates the input spike trains based on the Input.Frequency array and
% the Silence flag (Flag.Silence). The type of input is chosen by entering
% an integer value between 0 and 3:
% 0 = Linear without any silence,
% 1 = Linear with silent periods,
% 2 = Poisson without any silence, and
% 3 = Poisson with Silent Periods
%
%See Linear_Input_Train and Poisson_Input_Train for further help on how the
%trains are generated

Input_Spikes = (zeros(Par.Num_Neurons_OP.*Par.Num_Inputs_OP, Par.T));% 1x1x100 matrix
Sec = 1/Par.Ts; %This is the number of time steps in a second
if Flag.InputType == 1
    switch Flag.ManIn
        case 1
            disp('0 = Linear Complete');
            disp('1 = Linear Silence');
            disp('2 = Poisson Complete');
            disp('3 = Poisson Silence')
            Flag.InType = input('Please Enter Input Type: [0]', 's');
            if isempty(Flag.InType)
                Flag.InType = '0';
            end

            while ((Flag.InType ~= '0') && (Flag.InType ~= '1')&& (Flag.InType ~= '2')&& (Flag.InType ~= '3'))
                disp('0 = Linear Complete');
                disp('1 = Linear Silence');
                disp('2 = Poisson Complete');
                disp('3 = Poisson Silence')
                Flag.InType = input('Please enter a valid Input Type: [0]', 's');
                if isempty(Flag.InType)
                    Flag.InType = 0;
                end
            end
    end
% else
%     Flag.InType = int2str(InputMode);
end
switch Flag.InType
    case  '0'
        Flag.Silence = 0;
        Linear(Par,Input, Silence, Flag);
    case  '1'
        Flag.Silence = 1;
        Linear(Par,Input, Silence, Flag);
    case  '2'
        Flag.Silence = 0;
        Poisson(Par,Input, Silence, Flag);
    case  '3'
        Flag.Silence = 1;
        Poisson(Par,Input, Silence, Flag);
    otherwise
        disp('0 = Linear Complete');
        disp('1 = Linear Silence');
        disp('2 = Poisson Complete');
        disp('3 = Poisson Silence')
        error('Input_Type is invalid. Please enter a value between 0-3');
end

        %%
    function Linear(Par,Input, Silence, Flag)

        %function [Input_Spikes] = Linear_Input_Train(Num_Inputs,frequency, T, Silence, Silence_Flag, Train_Plot)
        %Generates a set of linear input spike trains based on the frequencies
        %passed to it. Where Num_Inputs is the number of inputs to generate,
        %Frequency is a matrix of train frequencies, T is length of the
        %sample window, Silence is the portions of silence to be added, Silence Flag
        %determines weither to add the silence or not, and Train_Plot determines
        %weither to plot the Input trains or not.

        if Flag.Silence == 1
            disp('Generating Linear trains with silence')
        else
            disp('Generating Linear trains')
        end
        %% generates linear spike trains
        % Input_Spikes = zeros(Par.Num_Inputs_OP, Par.T+1);% 1x1x100 matrix
        lambda = zeros(1, Par.Num_Inputs_OP);
        counter = 1;
        for z = 1:Par.Num_Neurons_OP
            for s = 1:Par.Num_Inputs_OP
                
                lambda(1,s) = Sec/Input.frequency(z,s);        %ISI time in ms
                t = 0;
                while t < Par.T
                    x = lambda(1,s);
                    t = t + round(x); %rounds x to the nearest integer
                    if t > Par.T
                        break;
                    end;
                    Input_Spikes(counter,t) = 1;
                    %             spike(1,:)= linearSpikeTrain(1,1,:);
                end
                counter = counter+1;
            end
        end

        % Insert Silence periods as specified by 'Silence'
        if Flag.Silence == 1
                for j = 1:size(Silence)
                    index = 1;
                    for i = 1:size(Silence,2)/2
                        a = Silence(j,index);
                        b = Silence(j,index+1);
                        if a == 0
                           a = 1;
                        end
                        Input_Spikes(j,a:b) = 0;
                        index = index+2;
                    end
                end
        end
    end

        %%
    function Poisson(Par,Input, Silence, Flag)

        %function [Input_Spikes] = Poisson_Input_Train(Num_Inputs,frequency, T, Silence, Silence_Flag, Train_Plot)
        %Generates a set of Poisson input spike trains based on the frequencies
        %passed to it. Where Num_Inputs is the number of inputs to generate,
        %Frequency is a matrix of train frequencies, T is length of the
        %sample window, Silence is the portions of silence to be added, Silence Flag
        %determines weither to add the silence or not, and Train_Plot determines
        %weither to plot the Input trains or not.
        
        if Flag.Silence == 1
            disp('Generating Poisson trains with silence')
        else
            disp('Generating Poisson trains')
        end
        %% generates Poisson spike trains
        counter = 1;
%         Input_Spikes = zeros(Par.Num_Inputs_OP, Par.T);% 
        for z = 1:Par.Num_Neurons_OP
            for s = 1:Par.Num_Inputs_OP
                lamda = Sec/Input.frequency(z,s);
                num_spikes = round ((Par.T+1)/Sec)*Input.frequency(z,s);
                spikes = poissrnd(lamda,1,num_spikes);%ISI time in ms
                total = 0;
                for i = 1 : size(spikes, 2)
                    total = total + spikes(:,i);
                    if total >=Par.T+1;
                        break
                    else
                        Input_Spikes(counter,total) = 1;
                    end
                end
                counter = counter + 1;
            end
        end
        
        % Insert Silence periods as specified by 'Silence'
        if Flag.Silence == 1
            
            for j = 1:size(Silence)
                index = 1;
                for i = 1:size(Silence,2)/2
                    a = Silence(j,index);
                    b = Silence(j,index+1);
                    if a == 0
                        a = 1;
                    end
                    Input_Spikes(j,a:b) = 0;
                    index = index+2;
                end
            end
            
        end

       
    end
if Flag.InSpikePlot == 1
    Num_Inputs = Par.Num_Inputs_OP .* Par.Num_Neurons_OP;
    [cellSpikes] = layerPlot(Num_Inputs,Par.T,Input_Spikes,Flag.InSpikePlot);
    saveas(gcf,'C:\Matlab\SelfRepair\TrainingResults\Input_Spikes.fig')
%     close
else
    cellSpikes = [];
end
disp('Input Spike Trains Generated')
clear Num_Inputs
end
