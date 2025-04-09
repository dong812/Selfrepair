function [] = Timer(TrainTime)
%% Calculate Training and Testing Times
TrainTimeh = floor(TrainTime/3600);
Trainhours = int2str(TrainTimeh);
TrainTimem = floor((TrainTime - TrainTimeh*3600)/60);
Trainmins = int2str(TrainTimem);
TrainTimes = (TrainTime - (TrainTimeh*3600) - (TrainTimem*60));
Trainsecs = int2str(TrainTimes);
TrainingTimeCaption = ['Training Time = ',Trainhours,':',Trainmins,'.',Trainsecs,' (hh:mm:ss)'];
caption = TrainingTimeCaption;
h = msgbox(caption, 'Finished', 'none', 'modal');