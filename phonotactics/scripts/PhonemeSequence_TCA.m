global TASK_DIR

global DUKEDIR
global BOX_DIR

RECONDIR=[BOX_DIR '\ECoG_Recon'];

fDown = 100; %Downsampled Sampling Frequency
timeExtract = [-2 2];

Task=[];

Task.Name='Phoneme_Sequencing';
Subject = popTaskSubjectData(Task);

Task.Base.Name='Start';
Task.Base.Epoch='Start';
Task.Base.Time=[-0.5 0];
TASK_DIR=([BOX_DIR '/CoganLab/D_Data/' Task.Name]);
DUKEDIR=TASK_DIR;

Task.Conds(1).Name='All';
Task.Conds(1).Field(1).Name='Response';
Task.Conds(1).Field(1).Epoch='ResponseStart';
Task.Conds(1).Field(1).Time=[-1 1];
