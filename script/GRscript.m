


%% Step 1
% Read Mvnx
clc 
clear
% load data
tree = load_mvnx_v2('gaitreporters1.mvnx');
% read some basic data from the file
% mvnxVersion = tree.metaData.mvnx_version;
% if (isfield(tree.metaData, 'comment'))
%     fileComments = tree.metaData.comment;
% end
% %read some basic properties of the subject;
% frameRate = tree.metaData.subject_frameRate;
% suitLabel = tree.metaData.subject_label;
% originalFilename = tree.metaData.subject_originalFilename;
% recDate = tree.metaData.subject_recDate;
% segmentCount = tree.metaData.subject_segmentCount;
% %retrieve sensor labels
% %creates a struct with sensor data
% if isfield(tree,'sensorData') && isstruct(tree.sensorData)
%     sensorData = tree.sensorData;
% end
% %retrieve segment labels
% %creates a struct with segment definitions
% if isfield(tree,'segmentData') && isstruct(tree.segmentData)
%     segmentData = tree.segmentData;
% end
% 
% %read the data from the structure e.g. segment 1
% if isfield(tree.segmentData,'position')
%     %Plot
%     figure('name','Position of first segment')
%     plot(tree.segmentData(1).position)
%     figure('name','Position of first segment in 3D')
%     plot3(tree.segmentData(1).position(:,1),tree.segmentData(1).position(:,2),tree.segmentData(1).position(:,3));
% end



%% Step 2
% load excel file (not necessary)& find the locs of rhs 
[subjectzxy,FTXT,FRAW] = xlsread('gaitreporters1.xlsx','Joint Angles ZXY', 'A:BO' );
a = (diff(tree.footContact(3).footContacts));
rhs= find(a==1)+1;

count_rhs= length(rhs);
number_of_cycles = count_rhs-1;



%% Step 3
%create cycles


for i= 1: number_of_cycles
    
    temp_var = strcat( 'gait_',num2str(i) );
    eval(sprintf('%s = subjectzxy(rhs(i):rhs(i+1),:)',temp_var));
    
      
end



%% Step 4
%Plot cycles for Knee flexion/extension angle

% 
% for j= 1:number_of_cycles
%        figure (1);hold on
%    plot(eval(['gait_', num2str(j),'(:,49)' ]));
%    legendInfo{j} = ['gait' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)

 
 
%% Step 5
%remove unmwanted cycles for Knee flexion/extension angle

for j= 1:number_of_cycles
    [peaks, peaks_ind] = findpeaks(eval(['gait_', num2str(j),'(:,49)' ]));
    if max(eval(['gait_', num2str(j),'(:,49)' ])) < 58
         clearvars(['gait_', num2str(j) ] );
    elseif eval(['gait_', num2str(j),'(1,49)' ])>18
        clearvars(['gait_', num2str(j) ] );
    
    elseif peaks(1)< 20
        clearvars(['gait_', num2str(j) ] );
    end
end



%% Step 6 
%plot new cycles 

% 
% hold on;
% for j= 4:10 %number_of_cycles
%        figure (1);
%    plot(eval(['gait_', num2str(j),'(:,49)' ]));
%    legendInfo{j} = ['gait' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)



%% Step 7
%Normalize Cycles



for i = 4 : 10 % number_of_cycles
    
    y=(0:100)';
    x=(linspace (0,100,length(eval(['gait_', num2str(i) ]))))';
    normalized_cycle_var = strcat( ['normalized_gait_', num2str(i) ]);
    a = interp1(x,eval(['gait_', num2str(i) ]), y,'spline');
    
    eval(sprintf('%s = a',normalized_cycle_var));
    
   
   
end




%% Step 8
%Plot normalized cycles for Knee flexion/extension angle

% 
% for j= 4:10 %number_of_cycles
%        figure (1);
%    plot(eval(['normalized_gait_', num2str(j),'(:,49)' ]));hold on; 
% end
%    
% hold off;

%% Step 9   knee flex/ext
%combine

gait_knee_fe = [];


for j= 4:10 %number_of_cycles
      
   B =  eval(['normalized_gait_', num2str(j),'(:,49)'] ) ;
   
    gait_knee_fe = [gait_knee_fe B];
    
    
end

       

gait_knee_fe_mean = mean(gait_knee_fe,2);
gait_knee_fe_std = std(gait_knee_fe')';
gait_cycle=0:1:100;

%% HIP DEGREE

% Step 1
% Read Mvnx
% clc 
% clear
% load data
tree = load_mvnx_v2('gaitreporters1.mvnx');
%% Step 2
% load excel file (not necessary)& find the locs of rhs 
[subjecthipzxy,hipTXT,hipRAW] = xlsread('gaitreporters1.xlsx','Joint Angles ZXY', 'A:BO' );
a = (diff(tree.footContact(3).footContacts));
rhs= find(a==1)+1;

count_rhs= length(rhs);
number_of_cycles = count_rhs-1;



%% Step 3
%create cycles


for i= 1: number_of_cycles
    
    temp_var = strcat( 'gait_hip_',num2str(i) );
    eval(sprintf('%s = subjecthipzxy(rhs(i):rhs(i+1),:)',temp_var));
    
      
end



%% Step 4
%Plot cycles for hip angle
% 
% hold on;
% for j= 1:number_of_cycles
%        figure (1);
%    plot(eval(['gait_hip_', num2str(j),'(:,46)' ]));
%    legendInfo{j} = ['gait_hip' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)

 
 
%% Step 5
%remove unmwanted cycles for hip angle

for j= 1:number_of_cycles
    [peaks, peaks_ind] = findpeaks(eval(['gait_hip_', num2str(j),'(:,46)' ]));
    %if max(eval(['gait_hip', num2str(j),'(:,46)' ])) < 58
     %    clearvars(['gait_hip', num2str(j) ] );
    if eval(['gait_hip_', num2str(j),'(1,46)' ])> 30
        clearvars(['gait_hip_', num2str(j) ] );
     elseif eval(['gait_hip_', num2str(j),'(1,46)' ])< 10
         clearvars(['gait_hip_', num2str(j) ] );
%     elseif eval(['gait_hip_', num2str(j),'(36,46)' ])< 0.4
%         clearvars(['gait_hip_', num2str(j) ] );
%     
    %elseif peaks(1)< 20
    %    clearvars(['gait_', num2str(j) ] );
    end
end



%% Step 6 
%plot new cycles 

% 
% hold on;
% for j= 4:11  %number_of_cycles
%        figure (1);
%    plot(eval(['gait_hip_', num2str(j),'(:,46)' ]));
%    legendInfo{j} = ['gait_hip' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)



%% Step 7
%Normalize Cycles



for i = 4:11 % number_of_cycles
    
    y=(0:100)';
    x=(linspace (0,100,length(eval(['gait_hip_', num2str(i) ]))))';
    normalized_cycle_var = strcat( ['normalized_gait_hip_', num2str(i) ]);
    a = interp1(x,eval(['gait_hip_', num2str(i) ]), y,'spline');
    
    eval(sprintf('%s = a',normalized_cycle_var));
    
   
   
end




%% Step 8
%Plot normalized cycles for hip angle

% 
% for j= 4:11 %number_of_cycles
%        figure (1);
%    plot(eval(['normalized_gait_hip_', num2str(j),'(:,46)' ]));hold on; 
% end
%    
% hold off;

%% Step 9         hip
%combine

gait_hip_fe = [];


for j= 4:11 %number_of_cycles
      
   B =  eval(['normalized_gait_hip_', num2str(j),'(:,46)'] ) ;
   
    gait_hip_fe = [gait_hip_fe B];
    
    
end

       

gait_hip_fe_mean = mean(gait_hip_fe,2);
gait_hip_fe_std = std(gait_hip_fe')';
gait_cycle=0:1:100;


%% ANKLE DEGREE

% Step 1
% Read Mvnx
% clc 
% clear
% load data
tree = load_mvnx_v2('gaitreporters1.mvnx');
%% Step 2
% load excel file (not necessary)& find the locs of rhs 
[subjectankzxy,ankTXT,ankRAW] = xlsread('gaitreporters1.xlsx','Joint Angles ZXY', 'A:BO' );
a = (diff(tree.footContact(3).footContacts));
rhs= find(a==1)+1;

count_rhs= length(rhs);
number_of_cycles = count_rhs-1;



%% Step 3
%create cycles


for i= 1: number_of_cycles
    
    temp_var = strcat( 'gait_ank_',num2str(i) );
    eval(sprintf('%s = subjectankzxy(rhs(i):rhs(i+1),:)',temp_var));
    
      
end



%% Step 4
%Plot cycles for ankle angle

% hold on;
% for j= 1:number_of_cycles
%        figure (1);
%    plot(eval(['gait_ank_', num2str(j),'(:,52)' ]));
%    legendInfo{j} = ['gait_ank' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)

 
 
%% Step 5
%remove unmwanted cycles for ankle

for j= 1:number_of_cycles
    [peaks, peaks_ind] = findpeaks(eval(['gait_ank_', num2str(j),'(:,52)' ]));
    %if max(eval(['gait_hip', num2str(j),'(:,46)' ])) < 58
     %    clearvars(['gait_hip', num2str(j) ] );
    if eval(['gait_ank_', num2str(j),'(1,52)' ])> 10
        clearvars(['gait_ank_', num2str(j) ] );
    elseif eval(['gait_ank_', num2str(j),'(68,52)' ])< -3
        clearvars(['gait_ank_', num2str(j) ] );
    elseif eval(['gait_ank_', num2str(j),'(39,52)' ])> 17
        clearvars(['gait_ank_', num2str(j) ] );
    elseif eval(['gait_ank_', num2str(j),'(47,52)' ])< -30
        clearvars(['gait_ank_', num2str(j) ] );
%    
    %elseif peaks(1)< 20
    %    clearvars(['gait_', num2str(j) ] );
    end
end



%% Step 6 
%plot new cycles 

% 
% hold on;
% for j= 4:9  %number_of_cycles
%        figure (1);
%    plot(eval(['gait_ank_', num2str(j),'(:,52)' ]));
%    legendInfo{j} = ['gait_ank' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)



%% Step 7
%Normalize Cycles



for i = 4:9 % number_of_cycles
    
    y=(0:100)';
    x=(linspace (0,100,length(eval(['gait_ank_', num2str(i) ]))))';
    normalized_cycle_var = strcat( ['normalized_gait_ank_', num2str(i) ]);
    a = interp1(x,eval(['gait_ank_', num2str(i) ]), y,'spline');
    
    eval(sprintf('%s = a',normalized_cycle_var));
    
   
   
end




%% Step 8
%Plot normalized cycles for ankle angle

% 
% for j= 4:9 %number_of_cycles
%        %figure (1);
%    plot(eval(['normalized_gait_ank_', num2str(j),'(:,52)' ]));hold on; 
% end
%    
% hold off;

%% Step 9
%combine

gait_ank_fe = [];


for j= 4:9 %number_of_cycles
      
   B =  eval(['normalized_gait_ank_', num2str(j),'(:,52)'] ) ;
   
    gait_ank_fe = [gait_ank_fe B];
    
    
end

       

gait_ank_fe_mean = mean(gait_ank_fe,2);
gait_ank_fe_std = std(gait_ank_fe')';
gait_cycle=0:1:100;
%% Plot all
figure(2)
subplot(3,2,4)
errorbar(gait_cycle, gait_ank_fe_mean,gait_ank_fe_std );
title('Right ankle FLEX/EXT plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)') 
subplot(3,2,2)
errorbar(gait_cycle, gait_knee_fe_mean,gait_knee_fe_std );
title('Right knee FLEX/EXT plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)') 
subplot(3,2,6)
errorbar(gait_cycle, gait_hip_fe_mean,gait_hip_fe_std );
title('Right hip FLEX/EXT plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)') 

%% LEFT SIDE OF THE BODY




%% Step 1
% Read Mvnx
%clc 
%clear
% load data
tree = load_mvnx_v2('gaitreporters1.mvnx');
% read some basic data from the file
% mvnxVersion = tree.metaData.mvnx_version;
% if (isfield(tree.metaData, 'comment'))
%     fileComments = tree.metaData.comment;
% end
% %read some basic properties of the subject;
% frameRate = tree.metaData.subject_frameRate;
% suitLabel = tree.metaData.subject_label;
% originalFilename = tree.metaData.subject_originalFilename;
% recDate = tree.metaData.subject_recDate;
% segmentCount = tree.metaData.subject_segmentCount;
% %retrieve sensor labels
% %creates a struct with sensor data
% if isfield(tree,'sensorData') && isstruct(tree.sensorData)
%     sensorData = tree.sensorData;
% end
% %retrieve segment labels
% %creates a struct with segment definitions
% if isfield(tree,'segmentData') && isstruct(tree.segmentData)
%     segmentData = tree.segmentData;
% end
% 
% %read the data from the structure e.g. segment 1
% if isfield(tree.segmentData,'position')
%     %Plot
%     figure('name','Position of first segment')
%     plot(tree.segmentData(1).position)
%     figure('name','Position of first segment in 3D')
%     plot3(tree.segmentData(1).position(:,1),tree.segmentData(1).position(:,2),tree.segmentData(1).position(:,3));
% end



%% Step 2
% load excel file (not necessary)& find the locs of lhs 
[lsubjectzxy,FTXT,FRAW] = xlsread('gaitreporters1.xlsx','Joint Angles ZXY', 'A:BO' );
b = (diff(tree.footContact(1).footContacts));
lhs= find(b==1)+1;

count_lhs= length(lhs);
number_of_cycles = count_lhs-1;



%% Step 3
%create cycles


for i= 1: number_of_cycles
    
    temp_var = strcat( 'lgait_',num2str(i) );
    eval(sprintf('%s = lsubjectzxy(lhs(i):lhs(i+1),:)',temp_var));
    
      
end

%% Step 4
%Plot cycles for Knee flexion/extension angle

% 
% for j= 1:number_of_cycles
%        figure (1);hold on;
%    plot(eval(['lgait_', num2str(j),'(:,61)' ]));
%    legendInfo{j} = ['lgait' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)

 
 
%% Step 5
%remove unmwanted cycles for Knee flexion/extension angle

for j= 1:number_of_cycles
%     [peaks, peaks_ind] = findpeaks(eval(['lgait_', num2str(j),'(:,61)' ]));
%     if max(eval(['lgait_', num2str(j),'(:,61)' ])) < 58
%          clearvars(['lgait_', num2str(j) ] );
%     if eval(['lgait_', num2str(j),'(35,61)' ]) > 22
%         clearvars(['lgait_', num2str(j) ] );
     if eval(['lgait_', num2str(j),'(1,61)' ]) < 6
         clearvars(['lgait_', num2str(j) ] );
     elseif eval(['lgait_', num2str(j),'(1,61)' ]) > 20
        clearvars(['lgait_', num2str(j) ] );
%     elseif peaks(1)< 20
%         clearvars(['lgait_', num2str(j) ] );
     end
end



%% Step 6 
%plot new cycles 



% for j= 3:10 %number_of_cycles
%        figure (1);hold on;
%    plot(eval(['lgait_', num2str(j),'(:,61)' ]));
%    legendInfo{j} = ['lgait' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)
% 


%% Step 7
%Normalize Cycles



for i = 3 : 10 % number_of_cycles
    
    y=(0:100)';
    x=(linspace (0,100,length(eval(['lgait_', num2str(i) ]))))';
    normalized_cycle_var = strcat( ['normalized_lgait_', num2str(i) ]);
    b = interp1(x,eval(['lgait_', num2str(i) ]), y,'spline');
    
    eval(sprintf('%s = b',normalized_cycle_var));
    
   
   
end




%% Step 8
%Plot normalized cycles for Knee flexion/extension angle


% for j= 4:10 %number_of_cycles
%        figure (1);
%    plot(eval(['normalized_lgait_', num2str(j),'(:,61)' ]));hold on; 
% end
%    
% hold off;

%% Step 9
%combine

lgait_knee_fe = [];


for j= 4:10 %number_of_cycles
      
   B =  eval(['normalized_lgait_', num2str(j),'(:,61)'] ) ;
   
    lgait_knee_fe = [lgait_knee_fe B];
    
    
end

       

lgait_knee_fe_mean = mean(lgait_knee_fe,2);
lgait_knee_fe_std = std(lgait_knee_fe')';
lgait_cycle=0:1:100;

%% HIP DEGREE

%% Step 1
% Read Mvnx
% clc 
% clear
% load data
tree = load_mvnx_v2('gaitreporters1.mvnx');
%% Step 2
% load excel file (not necessary)& find the locs of lhs 
[subjectlhipzxy,hipTXT,hipRAW] = xlsread('gaitreporters1.xlsx','Joint Angles ZXY', 'A:BO' );
b = (diff(tree.footContact(1).footContacts));
lhs= find(b==1)+1;

count_lhs= length(lhs);
number_of_cycles = count_lhs-1;



%% Step 3
%create cycles


for i= 1: number_of_cycles
    
    temp_var = strcat( 'lgait_hip_',num2str(i) );
    eval(sprintf('%s = subjectlhipzxy(lhs(i):lhs(i+1),:)',temp_var));
    
      
end



%% Step 4
%Plot cycles for Left Hip angle
% 
% for j= 1:number_of_cycles
%        figure (1);hold on;
%    plot(eval(['lgait_hip_', num2str(j),'(:,58)' ]));
%    legendInfo{j} = ['lgait_hip' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)

 
 
%% Step 5
%remove unmwanted cycles for Left Hip angle

for j= 1:number_of_cycles
   % [peaks, peaks_ind] = findpeaks(eval(['lgait_hip_', num2str(j),'(:,58)' ]));
    %if max(eval(['lgait_hip', num2str(j),'(:,58)' ])) < 58
     %    clearvars(['lgait_hip', num2str(j) ] );
    if eval(['lgait_hip_', num2str(j),'(1,58)' ]) < 15
        clearvars(['lgait_hip_', num2str(j) ] );
     elseif eval(['lgait_hip_', num2str(j),'(30,58)' ])> 5
         clearvars(['lgait_hip_', num2str(j) ] );
%     elseif eval(['lgait_hip_', num2str(j),'(36,58)' ])< 0.4
%         clearvars(['lgait_hip_', num2str(j) ] );
%     
    %elseif peaks(1)< 20
    %    clearvars(['lgait_', num2str(j) ] );
    end
end



%% Step 6 
%plot new cycles 


% 
% for j= 3:10  %number_of_cycles
%        figure (1);hold on;
%    plot(eval(['lgait_hip_', num2str(j),'(:,58)' ]));
%    legendInfo{j} = ['lgait_hip' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)



%% Step 7
%Normalize Cycles



for i = 3:10 % number_of_cycles
    
    y=(0:100)';
    x=(linspace (0,100,length(eval(['lgait_hip_', num2str(i) ]))))';
    normalized_cycle_var = strcat( ['normalized_lgait_hip_', num2str(i) ]);
    b = interp1(x,eval(['lgait_hip_', num2str(i) ]), y,'spline');
    
    eval(sprintf('%s = b',normalized_cycle_var));
    
   
   
end




%% Step 8
%Plot normalized cycles for Left Hip angle


% for j= 3:10 %number_of_cycles
%        figure (1);
%    plot(eval(['normalized_lgait_hip_', num2str(j),'(:,58)' ]));hold on; 
% end
%    
% hold off;

%% Step 9
%combine

lgait_hip_fe = [];


for j= 3:10 %number_of_cycles
      
   B =  eval(['normalized_lgait_hip_', num2str(j),'(:,58)'] ) ;
   
    lgait_hip_fe = [lgait_hip_fe B];
    
    
end

       

lgait_hip_fe_mean = mean(lgait_hip_fe,2);
lgait_hip_fe_std = std(lgait_hip_fe')';
lgait_cycle=0:1:100;


%% ANKLE DEGREE

%% Step 1
% Read Mvnx
% clc 
% clear
% load data
tree = load_mvnx_v2('gaitreporters1.mvnx');
%% Step 2
% load excel file (not necessary)& find the locs of lhs 
[subjectlankzxy,ankTXT,ankRAW] = xlsread('gaitreporters1.xlsx','Joint Angles ZXY', 'A:BO' );
b = (diff(tree.footContact(1).footContacts));
lhs= find(b==1)+1;

count_lhs= length(lhs);
number_of_cycles = count_lhs-1;



%% Step 3
%create cycles


for i= 1: number_of_cycles
    
    temp_var = strcat( 'lgait_ank_',num2str(i) );
    eval(sprintf('%s = subjectlankzxy(lhs(i):lhs(i+1),:)',temp_var));
    
      
end



%% Step 4
%Plot cycles for ANKLE angle
% 
% for j= 1:number_of_cycles
%        figure (1);hold on;
%    plot(eval(['lgait_ank_', num2str(j),'(:,64)' ]));
%    legendInfo{j} = ['lgait_ank' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)
% 
%  
 
%% Step 5
%remove unmwanted cycles for ANKLE angle

for j= 1:number_of_cycles
   % [peaks, peaks_ind] = findpeaks(eval(['lgait_ank_', num2str(j),'(:,64)' ]));
    %if max(eval(['lgait_hip', num2str(j),'(:,58)' ])) < 58
     %    clearvars(['lgait_hip', num2str(j) ] );
    if eval(['lgait_ank_', num2str(j),'(1,64)' ]) > 6
        clearvars(['lgait_ank_', num2str(j) ] );
     elseif eval(['lgait_ank_', num2str(j),'(48,64)' ])> -7
         clearvars(['lgait_ank_', num2str(j) ] );
%     elseif eval(['lgait_ank_', num2str(j),'(39,64)' ])> 17
%         clearvars(['lgait_ank_', num2str(j) ] );
%     elseif eval(['lgait_ank_', num2str(j),'(47,64)' ])< -30
%         clearvars(['lgait_ank_', num2str(j) ] );
%    
    %elseif peaks(1)< 20
    %    clearvars(['lgait_', num2str(j) ] );
    end
end



%% Step 6 
%plot new cycles 


% 
% for j= 3:10  %number_of_cycles
%        figure (1);hold on;
%    plot(eval(['lgait_ank_', num2str(j),'(:,64)' ]));
%    legendInfo{j} = ['lgait_ank' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)



%% Step 7
%Normalize Cycles



for i = 3:10 % number_of_cycles
    
    y=(0:100)';
    x=(linspace (0,100,length(eval(['lgait_ank_', num2str(i) ]))))';
    normalized_cycle_var = strcat( ['normalized_lgait_ank_', num2str(i) ]);
    b = interp1(x,eval(['lgait_ank_', num2str(i) ]), y,'spline');
    
    eval(sprintf('%s = b',normalized_cycle_var));
    
   
   
end




%% Step 8
%Plot normalized cycles for ANKLE angle

% 
% for j= 3:10 %number_of_cycles
%        figure (1);
%    plot(eval(['normalized_lgait_ank_', num2str(j),'(:,64)' ]));hold on; 
% end
%    
% hold off;

%% Step 9
%combine

lgait_ank_fe = [];

for j= 3:10 %number_of_cycles
      
   B =  eval(['normalized_lgait_ank_', num2str(j),'(:,64)'] ) ;
   
    lgait_ank_fe = [lgait_ank_fe B];
    
    
end

lgait_ank_fe_mean = mean(lgait_ank_fe,2);
lgait_ank_fe_std = std(lgait_ank_fe')';
lgait_cycle=0:1:100;
%% Plot all
figure(2)
subplot(3,2,3)
errorbar(lgait_cycle, lgait_ank_fe_mean,lgait_ank_fe_std );
title('Left ankle FLEX/EXT plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)') 
subplot(3,2,1)
errorbar(lgait_cycle, lgait_knee_fe_mean,lgait_knee_fe_std );
title('Left knee FLEX/EXT plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)') 
subplot(3,2,5)
errorbar(lgait_cycle, lgait_hip_fe_mean,lgait_hip_fe_std );
title('Left hip FLEX/EXT plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)') 

%% PLOT MORE
%% INT/EXT
%% Right KNEE INT/EXT
clear legendInfo
for j= 4:10 %number_of_cycles
       figure (3);
       subplot(3,2,2)
       hold on;
   plot(eval(['normalized_gait_', num2str(j),'(:,48)' ]));hold on; 
   legendInfo{j} = ['Right knee_' (num2str(j))];
end
title('Right knee INT/EXT plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')

hold off;
legend(legendInfo(4:10))
%% LEFT KNEE INT/EXT
clear legendInfo
for j= 3:10 %number_of_cycles
       figure (3);
       subplot(3,2,1)
       hold on;
   plot(eval(['normalized_lgait_', num2str(j),'(:,60)' ]));hold on; 
   legendInfo{j} = ['Left knee_' (num2str(j))];
end

title('Left knee INT/EXT plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo(3:10))


 
%% RIGHT ANKLE INT/EXT
clear legendInfo
for j= 4:9 %number_of_cycles
       figure (3);
       subplot(3,2,4)
       hold on;
   plot(eval(['normalized_gait_ank_', num2str(j),'(:,51)' ]));hold on; 
   legendInfo{j} = ['Right ankle_' (num2str(j))];
end
title(' Right ankle INT/EXT plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo(4:9))
 %% LEFT ANKLE INT/EXT
clear legendInfo
for j= 3:10 %number_of_cycles
       figure (3);
       subplot(3,2,3)
       hold on;
   plot(eval(['normalized_lgait_ank_', num2str(j),'(:,63)' ]));hold on; 
   legendInfo{j} = ['Left ankle_' (num2str(j))];
end
title('Left ankle INT/EXT plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo(3:10))
 

%% RIGHT HIP INT/EXT
clear legendInfo
for j= 4:11 %number_of_cycles
       figure (3);
       subplot(3,2,6)
       hold on;
   plot(eval(['normalized_gait_hip_', num2str(j),'(:,45)' ]));hold on; 
   legendInfo{j} = ['Right hip_' (num2str(j))];
end
title('Right hip INT/EXT plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo(4:11))
 %% LEFT HIP INT/EXT
clear legendInfo
for j= 3:10 %number_of_cycles
       figure (3);
       subplot(3,2,5)
       hold on;
   plot(eval(['normalized_lgait_hip_', num2str(j),'(:,57)' ]));hold on; 
   legendInfo{j} = ['Left hip_' (num2str(j))];
end
title('Left hip INT/EXT plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo(3:10))

%% left ballfoot flexion and extension
clear legendInfo
for j= [4,5,7:10] %number_of_cycles
       figure (4);
       subplot(3,2,1)
       hold on;
   plot(eval(['normalized_lgait_ank_', num2str(j),'(:,67)' ]),'LineWidth', 1 ); 
   legendInfo{j} = ['Left Ballfoot_' (num2str(j))];
end
title('Left Ballfoot Flex/Exten plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([4,5,7:10]))
%% Left Ballfoot Internal and External Rotation
clear legendInfo
for j= [4,5,8,10] %number_of_cycles
       figure (4);
       subplot(3,2,3)
       hold on;
   plot(eval(['normalized_lgait_ank_', num2str(j),'(:,66)' ]),'LineWidth', 1 ); 
   legendInfo{j} = ['Left Ballfoot_' (num2str(j))];
end
title('Left Ballfoot Int/Ext plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([4,5,8,10]))
%% Left Ballfoot Abduction and Adduction
clear legendInfo
for j= [5,6,8,10] %number_of_cycles
       figure (4);
       subplot(3,2,5)
       hold on;
   plot(eval(['normalized_lgait_ank_', num2str(j),'(:,65)' ]),'LineWidth', 1 ); 
   legendInfo{j} = ['Left Ballfoot_' (num2str(j))];
end
title('Left Ballfoot Abd/Add plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([5,6,8,10]))
%% Right Ballfoot Flexion and Extension
clear legendInfo
for j= [6:9] %number_of_cycles
       figure (4);
       subplot(3,2,2)
       hold on;
   plot(eval(['normalized_gait_ank_', num2str(j),'(:,55)' ]),'LineWidth', 1 ); 
   legendInfo{j} = ['Right Ballfoot_' (num2str(j))];
end
title('Right Ballfoot Flex/Exten plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([6:9]))
%% Right Ballfoot Internal and External Rotation
clear legendInfo
for j= [4,6:8] %number_of_cycles
       figure (4);
       subplot(3,2,4)
       hold on;
   plot(eval(['normalized_gait_ank_', num2str(j),'(:,54)' ]),'LineWidth', 1 ); 
   legendInfo{j} = ['Right Ballfoot_' (num2str(j))];
end
title('Right Ballfoot Int/Ext plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([4,6:8]))
%% Right Ballfoot Abduction and Adduction
clear legendInfo
for j= [4,6:9] %number_of_cycles
       figure (4);
       subplot(3,2,6)
       hold on;
   plot(eval(['normalized_gait_ank_', num2str(j),'(:,53)' ]),'LineWidth', 1 ); 
   legendInfo{j} = ['Right Ballfoot_' (num2str(j))];
end
title('Right Ballfoot Abd/Add plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([4,6:9]))
%% Right Hip Abduction and Adduction
clear legendInfo
for j= [4:10] %number_of_cycles
       figure (5);
       subplot(3,2,6)
       hold on;
   plot(eval(['normalized_gait_hip_', num2str(j),'(:,44)' ]),'LineWidth', 2 ); 
   legendInfo{j} = ['Right hip_' (num2str(j))];
end
title('Right hip ABD/ADD plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([4:10]))
%% Right Knee Abduction and Adduction
clear legendInfo
for j= [4,5,7,9,10] %number_of_cycles
       figure (5);
       subplot(3,2,2)
       hold on;
   plot(eval(['normalized_gait_', num2str(j),'(:,47)' ]),'LineWidth', 2 ); 
   legendInfo{j} = ['Right Knee_' (num2str(j))];
end
title('Right Knee ABD/ADD plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([4,5,7,9,10]))

%% Right Ankle Abduction and Adduction 
clear legendInfo
for j= [4,6,7,9] %number_of_cycles
       figure (5);
       subplot(3,2,4)
       hold on;
   plot(eval(['normalized_gait_ank_', num2str(j),'(:,50)' ]),'LineWidth', 1 ); 
   legendInfo{j} = ['Right Ankle_' (num2str(j))];
end
title('Right Ankle ABD/ADD plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([4,6,7,9]))
%% Left Hip Abduction and Adduction
clear legendInfo
for j= [3:8] %number_of_cycles
       figure (5);
       subplot(3,2,5)
       hold on;
   plot(eval(['normalized_lgait_hip_', num2str(j),'(:,56)' ]),'LineWidth', 1 ); 
   legendInfo{j} = ['Left Hip_' (num2str(j))];
end
title('Left Hip Abd/Add plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([3:8]))
%% Left Knee Abduction and Adduction
clear legendInfo
for j= [4,5,7,8,10] %number_of_cycles
       figure (5);
       subplot(3,2,1)
       hold on;
   plot(eval(['normalized_lgait_', num2str(j),'(:,59)' ]),'LineWidth', 1 ); 
   legendInfo{j} = ['Left Knee_' (num2str(j))];
end
title('Left Knee Abd/Add plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([4,5,7,8,10]))
%% Left Ankle Abduction and Adduction
clear legendInfo
for j= [4,5,7:10] %number_of_cycles
       figure (5);
       subplot(3,2,3)
       hold on;
   plot(eval(['normalized_lgait_ank_', num2str(j),'(:,62)' ]),'LineWidth', 1 ); 
   legendInfo{j} = ['Left Ankle_' (num2str(j))];
end
title('Left Ankle Abd/Add plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([4,5,7:10]))
%%%% Step 2
% load excel file (not necessary)& find the locs of lhs
[subjectlankzxy,ankTXT,ankRAW] = xlsread('gaitreporters1.xlsx','Ergonomic Joint Angles ZXY', 'A:BO' );
b = (diff(tree.footContact(1).footContacts));
lhs= find(b==1)+1;

count_lhs= length(lhs);
number_of_cycles = count_lhs-1;

%% Step 3
%create cycles
for i= 1: number_of_cycles
    
    temp_var = strcat( 'gait_pelvis_',num2str(i) );
    eval(sprintf('%s = subjectlankzxy(lhs(i):lhs(i+1),:)',temp_var));
    
      
end

%% Step 4
%Plot cycles for Vertical Axial Pelvis Bending Plot
for j= 1:number_of_cycles
       figure (6);hold on;
   plot(eval(['gait_pelvis_', num2str(j),'(:,15)' ]));
   legendInfo{j} = ['gait_pelvis_' num2str(j)];
  
end
 hold off;  
 legend(legendInfo)
%% Step 5 plotting new cycles for Vertical Axial Pelvis Bending Plot 
clear legendInfo
for j= [3:5,9,10] %number_of_cycles
       figure (6);
       %subplot(3,2,5)
       hold on;
   plot(eval(['gait_pelvis_', num2str(j),'(:,15)' ]),'LineWidth', 0.5 ); 
   legendInfo{j} = ['gait_pelvis_' (num2str(j))];
end
title('Vertical Axial Pelvis Bending Plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([3:5,9,10]))
%% Step 6
%Normalize Cycles

for i = [3:5,9,10] % number_of_cycles
    
    y=(0:100)';
    x=(linspace (0,100,length(eval(['gait_pelvis_', num2str(i) ]))))';
    normalized_cycle_var = strcat( ['normalized_gait_pelvis_', num2str(i) ]);
    b = interp1(x,eval(['gait_pelvis_', num2str(i) ]), y,'spline');
    
    eval(sprintf('%s = b',normalized_cycle_var));
    
end
%% Step 7
%Plot normalized cycles for Vertical Axial Pelvis Bending Plot

for j= [3:5,9,10] %number_of_cycles
       figure (6);
   plot(eval(['normalized_gait_pelvis_', num2str(j),'(:,15)' ]));hold on; 
end
   
hold off;
title('Vertical Axial Pelvis Bending Plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([3:5,9,10]))
%% Vertical Pelvis Flexion/Extension



% Step 1
% Read Mvnx
% clc 
% clear
% load data
tree = load_mvnx_v2('gaitreporters1.mvnx');
%% Step 2
% load excel file (not necessary)& find the locs of lhs
[subjectlankzxy,ankTXT,ankRAW] = xlsread('gaitreporters1.xlsx','Ergonomic Joint Angles ZXY', 'A:BO' );
b = (diff(tree.footContact(1).footContacts));
lhs= find(b==1)+1;

count_lhs= length(lhs);
number_of_cycles = count_lhs-1;

%% Step 3
%create cycles
for i= 1: number_of_cycles
    
    temp_var = strcat( 'gait_pelvis_',num2str(i) );
    eval(sprintf('%s = subjectlankzxy(lhs(i):lhs(i+1),:)',temp_var));   
end

%% Step 4
%Plot cycles for Vertical Pelvis Flexion/Extension Plot
for j= 1:number_of_cycles
       figure (7);hold on;
   plot(eval(['gait_pelvis_', num2str(j),'(:,16)' ]));
   legendInfo{j} = ['gait_pelvis_' num2str(j)];
  
end
 hold off;  
 legend(legendInfo)
%% Step 5 plotting new cycles for Vertical Pelvis Flexion/Extension Plot 
clear legendInfo
for j= [4:10] %number_of_cycles
       figure (7);
       %subplot(3,2,5)
       hold on;
   plot(eval(['gait_pelvis_', num2str(j),'(:,16)' ]),'LineWidth', 0.5 ); 
   legendInfo{j} = ['gait_pelvis_' (num2str(j))];
end
title('Vertical Pelvis Flexion/Extension Plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([4:10]))
%% Step 6
%Normalize Cycles

for i = [4:10] % number_of_cycles
    
    y=(0:100)';
    x=(linspace (0,100,length(eval(['gait_pelvis_', num2str(i) ]))))';
    normalized_cycle_var = strcat( ['normalized_gait_pelvis_', num2str(i) ]);
    b = interp1(x,eval(['gait_pelvis_', num2str(i) ]), y,'spline');
    
    eval(sprintf('%s = b',normalized_cycle_var));
    
end
%% Step 7
%Plot normalized cycles for Vertical Pelvis Flexion/Extension Plot

for j= [4:10] %number_of_cycles
       figure (7);
   plot(eval(['normalized_gait_pelvis_', num2str(j),'(:,16)' ]));hold on; 
end
   
hold off;
title('Vertical Pelvis Flexion/Extension Plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([4:10]))

%% Vertical Pelvis Lateral Bending 
% Step 1
tree = load_mvnx_v2('gaitreporters1.mvnx');

%% Step 2 for Vertical Pelvis Lateral Bending
% load excel file (not necessary)& find the locs of lhs
[subjectlankzxy,ankTXT,ankRAW] = xlsread('gaitreporters1.xlsx','Ergonomic Joint Angles ZXY', 'A:BO' );
b = (diff(tree.footContact(1).footContacts));
lhs= find(b==1)+1;

count_lhs= length(lhs);
number_of_cycles = count_lhs-1;

%% Step 3 for Vertical Pelvis Lateral Bending
%create cycles
for i= 1: number_of_cycles
    
    temp_var = strcat( 'gait_pelvis_',num2str(i) );
    eval(sprintf('%s = subjectlankzxy(lhs(i):lhs(i+1),:)',temp_var));   
end

%% Step 4
%Plot cycles for Vertical Pelvis Lateral Bending 
for j= 1:number_of_cycles
       figure (8);hold on;
   plot(eval(['gait_pelvis_', num2str(j),'(:,14)' ]));
   legendInfo{j} = ['gait_pelvis_' num2str(j)];
  
end
 hold off;  
 legend(legendInfo)
%% Step 5 plotting new cycles for Vertical Pelvis Lateral Bending
clear legendInfo
for j= [3:10] %number_of_cycles
       figure (8);
       %subplot(3,2,5)
       hold on;
   plot(eval(['gait_pelvis_', num2str(j),'(:,14)' ]),'LineWidth', 0.5 ); 
   legendInfo{j} = ['gait_pelvis_' (num2str(j))];
end
title('Vertical Pelvis Lateral Bending Plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([3:10]))
%% Step 6
%Normalize Cycles for Vertical Pelvis Lateral Bending Plot

for i = [3:10] % number_of_cycles
    
    y=(0:100)';
    x=(linspace (0,100,length(eval(['gait_pelvis_', num2str(i) ]))))';
    normalized_cycle_var = strcat( ['normalized_gait_pelvis_', num2str(i) ]);
    b = interp1(x,eval(['gait_pelvis_', num2str(i) ]), y,'spline');
    
    eval(sprintf('%s = b',normalized_cycle_var));
    
end
%% Step 7
%Plot normalized cycles for Vertical Pelvis Lateral Bending Plot

for j= [3:10] %number_of_cycles
       figure (8);
   plot(eval(['normalized_gait_pelvis_', num2str(j),'(:,14)' ]));hold on; 
end
   
hold off;
title('Vertical Pelvis Lateral Bending Plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
hold off;
legend(legendInfo([3:10]))

%% Step 9   RIGHT KNEE FLEX/EXT (Not Normalized Data)

gait_knee_fe = [];


for j= 4:10 %number_of_cycles
      
   B =  eval(['gait_', num2str(j),'(:,49)'] ) ;
   
    xgait_knee_fe = [gait_knee_fe B]; % added x before names to change
    
end
xgait_knee_fe_mean = mean(xgait_knee_fe,2);
xgait_knee_fe_std = std(xgait_knee_fe')';
%gait_cycle=0:1:100;

%% Step 9       RIGHT HIP          (Not Normalized Data)

gait_hip_fe = [];
for j= 4:11 %number_of_cycles
      
   B =  eval(['gait_hip_', num2str(j),'(:,46)'] ) ;
   
    xgait_hip_fe = [gait_hip_fe B];
    
    
end
xgait_hip_fe_mean = mean(xgait_hip_fe,2);
xgait_hip_fe_std = std(xgait_hip_fe')';
%gait_cycle=0:1:100;

%% Step 9   RIGHT ANKLE        (Not Normalized Data)
%combine

gait_ank_fe = [];
for j= 4:9 %number_of_cycles
      
   B =  eval(['gait_ank_', num2str(j),'(:,52)'] ) ;
   
    xgait_ank_fe = [gait_ank_fe B];
    
    
end
xgait_ank_fe_mean = mean(xgait_ank_fe,2);
xgait_ank_fe_std = std(xgait_ank_fe')';
%gait_cycle=0:1:100;

%% Step 9  LEFT KNEE FLEX/EXT         (Not Normalized Data)
%combine

lgait_knee_fe = [];
for j= 4:10 %number_of_cycles
      
   B =  eval(['lgait_', num2str(j),'(:,61)'] ) ;
   
    xlgait_knee_fe = [lgait_knee_fe B];
    
    
end
xlgait_knee_fe_mean = mean(xlgait_knee_fe,2);
xlgait_knee_fe_std = std(xlgait_knee_fe')';
%gait_cycle=0:1:100;

%% Step 9    LEFT HIP           (Not Normalized Data)
%combine

lgait_hip_fe = [];
for j= 3:10 %number_of_cycles
      
   B =  eval(['lgait_hip_', num2str(j),'(:,58)'] ) ;
   
    xlgait_hip_fe = [lgait_hip_fe B];
    
    
end
xlgait_hip_fe_mean = mean(xlgait_hip_fe,2);
xlgait_hip_fe_std = std(xlgait_hip_fe')';
%gait_cycle=0:1:100;

%% Step 9     LEFT ANKLE      (Not Normalized Data)
%combine

lgait_ank_fe = [];
for j= 3:10 %number_of_cycles
      
   B =  eval(['lgait_ank_', num2str(j),'(:,64)'] ) ;
   
    xlgait_ank_fe = [lgait_ank_fe B];
    
    
end
xlgait_ank_fe_mean = mean(xlgait_ank_fe,2);
xlgait_ank_fe_std = std(xlgait_ank_fe')';
%xlgait_cycle=0:1:100;

%% Step 9 for Vertical Pelvis Flexion/Extension  (Not Normalized Data)
gait_pelvis_fe = [];
for j= [4:10] %number_of_cycles
      
   B =  eval(['gait_pelvis_', num2str(j),'(:,16)'] ) ;
   
    xgait_pelvis_fe = [gait_pelvis_fe B];

end
xflex_pelv_fe_mean = mean(xgait_pelvis_fe,2);
xflex_pelv_fe_std = std(xgait_pelvis_fe')';
%xlgait_cycle=0:1:100;

%% Step 9 for Vertical Axial Pelvis Bending  (Not Normalized Data)
gait2_pelvis_fe = [];
for j = [3:5,9,10] %number_of_cycles
      
   B =  eval(['gait_pelvis_', num2str(j),'(:,15)'] ) ;
   
    xgait_axpelvis_fe = [gait2_pelvis_fe B];

end
xaxial_pelv_fe_mean = mean(xgait_axpelvis_fe,2);
xaxial_apelv_fe_std = std(xgait_axpelvis_fe')';
%xlgait_cycle=0:1:100;

%% Step 9 for Vertical Pelvis Lateral Bending Plot  (Not Normalized Data)
gait3_pelvis_fe = [];
for j= [3:10] %number_of_cycles
      
   B =  eval(['gait_pelvis_', num2str(j),'(:,14)'] ) ;
   
    xgait_lapelvis_fe = [gait3_pelvis_fe B];

end
xla_pelv_fe_mean = mean(xgait_lapelvis_fe,2);
xla_apelv_fe_std = std(xgait_lapelvis_fe')';
%xlgait_cycle=0:1:100;

%% RIGHT KNEE INT/EXT ROTATION MEAN   (normalized)
gait_knee_fe = [];

for j= 4:10 %number_of_cycles
      
   B =  eval(['normalized_gait_', num2str(j),'(:,48)'] ) ;
   
    gait_knee_fe = [gait_knee_fe B]; 
    
end
int_knee_mean = mean(gait_knee_fe,2);
int_knee_std = std(gait_knee_fe')';
gait_cycle=0:1:100;

figure(9)
 subplot(2,2,4)
errorbar(gait_cycle, int_knee_mean,int_knee_std );
title('Right Knee Int/Ext Mean plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')

%% RIGHT KNEE ABD/ADD MEAN    (normalized)
gait_knee_fe = [];

for j= [4,5,7,9,10] %number_of_cycles
      
   B =  eval(['normalized_gait_', num2str(j),'(:,47)'] ) ;
   
    gait_knee_fe = [gait_knee_fe B]; 
    
end
abd_knee_mean = mean(gait_knee_fe,2);
abd_knee_std = std(gait_knee_fe')';
gait_cycle=0:1:100;

figure(9)
  subplot(2,2,2)

errorbar(gait_cycle, abd_knee_mean,abd_knee_std );
title('Right Knee Abd/Add Mean plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
%% LEFT KNEE ABD/ADD MEAN   (normalized)
lgait_knee_fe = [];

for j= [4,5,7,8,10] %number_of_cycles
      
   B =  eval(['normalized_lgait_', num2str(j),'(:,59)'] ) ;
   
    lgait_knee_fe = [lgait_knee_fe B]; 
    
end
abd_lknee_mean = mean(lgait_knee_fe,2);
abd_lknee_std = std(lgait_knee_fe')';
gait_cycle=0:1:100;

figure(9)
  subplot(2,2,1)
 
errorbar(gait_cycle, abd_lknee_mean,abd_lknee_std );
title('Left Knee Abd/Add Mean plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
%% LEFT KNEE INT/EXT MEAN   (normalized)
lgait_knee_fe = [];

for j= 3:10 %number_of_cycles
      
   B =  eval(['normalized_lgait_', num2str(j),'(:,60)'] ) ;
   
    lgait_knee_fe = [lgait_knee_fe B]; 
    
end
int_lknee_mean = mean(lgait_knee_fe,2);
int_lknee_std = std(lgait_knee_fe')';
gait_cycle=0:1:100;

figure(9)
  subplot(2,2,3) 
errorbar(gait_cycle, int_lknee_mean, int_lknee_std );
title('Left Knee Int/Ext Mean plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
%% RIGHT HIP INT/EXT ROTATION MEAN   (normalized)
gait_hip_fe = [];

for j= 4:11 %number_of_cycles
      
   B =  eval(['normalized_gait_hip_', num2str(j),'(:,45)'] ) ;
   
    gait_hip_fe = [gait_hip_fe B];
    
    
end

int_hip_fe_mean = mean(gait_hip_fe,2);
int_hip_fe_std = std(gait_hip_fe')';
gait_cycle=0:1:100;

figure(10)
  subplot(2,2,4)

errorbar(gait_cycle, int_hip_fe_mean, int_hip_fe_std );
title('Right Hip Int/Ext Mean plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
%% RIGHT HIP ABD/ADD MEAN   (normalized)
gait_hip_fe = [];

for j= [4:10] %number_of_cycles
      
   B =  eval(['normalized_gait_hip_', num2str(j),'(:,44)'] ) ;
   
    gait_hip_fe = [gait_hip_fe B];
    
    
end

abd_hip_fe_mean = mean(gait_hip_fe,2);
abd_hip_fe_std = std(gait_hip_fe')';
gait_cycle=0:1:100;

figure(10)
  subplot(2,2,2)

errorbar(gait_cycle, abd_hip_fe_mean, abd_hip_fe_std );
title('Right Hip Abd/Add Mean plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
%% LEFT HIP INT/EXT ROTATION MEAN    (normalized)

lgait_hip_fe = [];

for j= 3:10 %number_of_cycles
      
   B =  eval(['normalized_lgait_hip_', num2str(j),'(:,57)'] ) ;
   
    lgait_hip_fe = [lgait_hip_fe B];
    
    
end
int_lhip_fe_mean = mean(lgait_hip_fe,2);
int_lhip_fe_std = std(lgait_hip_fe')';
gait_cycle=0:1:100;

figure(10)
  subplot(2,2,3)

errorbar(gait_cycle, int_lhip_fe_mean, int_lhip_fe_std );
title('Left Hip Int/Ext Mean plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
%% LEFT HIP ABD/ADD MEAN      (normalized)

lgait_hip_fe = [];

for j= 3:8 %number_of_cycles
      
   B =  eval(['normalized_lgait_hip_', num2str(j),'(:,56)'] ) ;
   
    lgait_hip_fe = [lgait_hip_fe B];
    
    
end
abd_lhip_fe_mean = mean(lgait_hip_fe,2);
abd_lhip_fe_std = std(lgait_hip_fe')';
gait_cycle=0:1:100;

figure(10)
  subplot(2,2,1)

errorbar(gait_cycle, abd_lhip_fe_mean, abd_lhip_fe_std );
title('Left Hip Abd/Add Mean plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
%% Vertical Pelvis Flexion/Extension Mean Plot (normalized)
gait_pelvis_fe = [];

for j= [4:10] %number_of_cycles
      
   B =  eval(['normalized_gait_pelvis_', num2str(j),'(:,16)'] ) ;
   
    gait_pelvis_fe = [gait_pelvis_fe B];
    
    
end
flex_pelv_fe_mean = mean(gait_pelvis_fe,2);
flex_pelv_fe_std = std(gait_pelvis_fe')';
gait_cycle=0:1:100;

figure(17)
 
errorbar(gait_cycle, flex_pelv_fe_mean, flex_pelv_fe_std );
title('Pelvis Flex/Exten Mean plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
%% Vertical Axial Pelvis Bending Mean Plot       (normalized)
gait_axpelvis_fe = [];

for j= [3:5,9,10] %number_of_cycles
      
   B =  eval(['normalized_gait_pelvis_', num2str(j),'(:,15)'] ) ;
   
    gait_axpelvis_fe = [gait_axpelvis_fe B];
    
    
end
axial_pelv_fe_mean = mean(gait_axpelvis_fe,2);
axial_pelv_fe_std = std(gait_axpelvis_fe')';
gait_cycle=0:1:100;

figure(18)
 
errorbar(gait_cycle, axial_pelv_fe_mean, axial_pelv_fe_std );
title('Vertical Axial Pelvis Bending Mean Plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')
%% Vertical Pelvis Lateral Bending Mean Plot  (normalized)
gait_lapelvis_fe = [];

for j= [3:10] %number_of_cycles
      
   B =  eval(['normalized_gait_pelvis_', num2str(j),'(:,14)'] ) ;
   
    gait_lapelvis_fe = [gait_lapelvis_fe B];
    
    
end
lat_pelv_fe_mean = mean(gait_lapelvis_fe,2);
lat_pelv_fe_std = std(gait_lapelvis_fe')';
gait_cycle=0:1:100;

figure(19)
 
errorbar(gait_cycle, lat_pelv_fe_mean, lat_pelv_fe_std );
title('Vertical Pelvis Lateral Bending Mean Plot')
xlabel('GC (%)') 
ylabel('joint angle (deg)')

%% ROM Right Ankle
ma_ank = max(xgait_ank_fe_mean);
mi_ank = min(xgait_ank_fe_mean);
rom_ankle = ma_ank - mi_ank
%% ROM Right Hip
ma_hip = max(xgait_hip_fe_mean);
mi_hip = min(xgait_hip_fe_mean);
rom_hip = ma_hip - mi_hip
%% ROM Right Knee
ma_knee = max(xgait_knee_fe_mean);
mi_knee = min(xgait_knee_fe_mean);
rom_knee = ma_knee - mi_knee
%% ROM Left Ankle
ma_lank = max(xlgait_ank_fe_mean);
mi_lank = min(xlgait_ank_fe_mean);
rom_lankle = ma_lank - mi_lank
%% ROM Left Hip
ma_lhip = max(xlgait_hip_fe_mean);
mi_lhip = min(xlgait_hip_fe_mean);
rom_lhip = ma_lhip - mi_lhip
%% ROM Left Knee
ma_lknee = max(xlgait_knee_fe_mean);
mi_lknee = min(xlgait_knee_fe_mean);
rom_lknee = ma_lknee - mi_lknee
%% ROM Vertical Pelvis Flexion/Extension
ma_flex_pelv = max(xflex_pelv_fe_mean);
mi_flex_pelv = min(xflex_pelv_fe_mean);
rom_flex_pelv = ma_flex_pelv - mi_flex_pelv
%% ROM Vertical Axial Pelvis Bending
ma_ax_pelv =  max(xaxial_pelv_fe_mean);
mi_ax_pelv = min(xaxial_pelv_fe_mean);
rom_ax_pelv = ma_ax_pelv - mi_ax_pelv
%% Vertical Pelvis Lateral Bending
ma_la_pelv = max(xla_pelv_fe_mean);
mi_la_pelv = min(xla_pelv_fe_mean);
rom_la_pelv = ma_la_pelv - mi_la_pelv
%% ROM Figure Code
figure;
T1 = uicontrol('Style', 'edit', 'Units', 'normalized', ...
               'Position', [0, 0, 1, 1], ...
               'Min', 0, 'Max', 2);
set(T1, 'String', {'ROM Right Ankle',rom_ankle,'', 'ROM Right Hip' , rom_hip,'','ROM Right Knee' , rom_knee,'', 'ROM Left Ankle' , rom_lankle,'','ROM Left Hip' , rom_lhip,'','ROM Left Knee', rom_lknee,'','ROM Vertical Pelvis Flexion/Extension', rom_flex_pelv,'','ROM Vertical Axial Pelvis Bending', rom_ax_pelv,'','Vertical Pelvis Lateral Bending', rom_la_pelv});


%% right initial contact and peaks

for j= 4:10 %number_of_cycles
       
   right_ic_kne(j - 3,1) = eval(['gait_', num2str(j),'(1,49)' ]);
   
   left_kne_max(j - 3,1) = max(eval(['gait_', num2str(j),'(:,49)' ]))
   left_kne_min(j - 3,1) = min(eval(['gait_', num2str(j),'(:,49)' ]))
  
end


for j= 4:9 %number_of_cycles
       
   right_ic_ank(j - 3,1) = eval(['gait_ank_', num2str(j),'(1,52)' ]);
   right_ank_max(j - 3,1) = max(eval(['gait_ank_', num2str(j),'(:,52)' ]))
   right_ank_min(j - 3,1) = min(eval(['gait_ank_', num2str(j),'(:,52)' ]))

  
end


for j= 4:11 %number_of_cycles
       
   right_ic_hip(j - 3,1) = eval(['gait_hip_', num2str(j),'(1,46)' ]);
   right_hip_max(j - 3,1) = max(eval(['gait_hip_', num2str(j),'(:,46)' ]))
   right_hip_min(j - 3,1) = min(eval(['gait_hip_', num2str(j),'(:,46)' ]))

  
end

%% left initial contact and peaks

for j= 3:10 %number_of_cycles
       
   left_ic_kne(j - 2,1) = eval(['lgait_', num2str(j),'(1,61)' ]);
   left_kne_max(j - 2,1) = max(eval(['lgait_', num2str(j),'(:,61)' ]))
   left_kne_min(j - 2,1) = min(eval(['lgait_', num2str(j),'(:,61)' ]))

  
end


for j= 3:10 %number_of_cycles
       
   left_ic_ank(j - 2,1) = eval(['lgait_ank_', num2str(j),'(1,64)' ]);
   left_ank_max(j - 2,1) = max(eval(['lgait_ank_', num2str(j),'(:,64)' ]))
   left_ank_min(j - 2,1) = min(eval(['lgait_ank_', num2str(j),'(:,64)' ]))

  
end


for j= 3:10 %number_of_cycles
       
   left_ic_hip(j - 2,1) = eval(['lgait_hip_', num2str(j),'(1,58)' ]);
   left_hip_max(j - 2,1) = max(eval(['lgait_hip_', num2str(j),'(:,58)' ]))
   left_hip_min(j - 2,1) = min(eval(['lgait_hip_', num2str(j),'(:,58)' ]))

  
end

%% temporospatial parameters

%% normal gait cycles

a = (diff(tree.footContact(3).footContacts));
rhs= find(a==1)+1;
time_right_gait = rhs * 1/60;
c = (diff(tree.footContact(1).footContacts));
lhs= find(c==1)+1;
time_left_gait = lhs * 1/60;

GCR = diff(time_right_gait);
GCL = diff(time_left_gait);
%% Left toe off
a = (diff(tree.footContact(2).footContacts));
ltoff= find(a==-1)+1;

count_ltoff= length(ltoff);
number_of_cycles_ltoff = count_ltoff-1;

%% Right toe off
b = (diff(tree.footContact(4).footContacts));
rtoff= find(b==-1)+1;

count_rtoff= length(rtoff);
number_of_cycles_rtoff = count_rtoff-1;

%% Left heel strike
c = (diff(tree.footContact(1).footContacts));
lhs= find(c==1)+1;

count_lhs= length(lhs);
number_of_cycles_lhs = count_lhs-1;

%% Right heel strike
d = (diff(tree.footContact(3).footContacts));
rhs= find(d==1)+1;

count_rhs= length(rhs);
number_of_cycles_rhs = count_rhs-1;

%% timings

time_rhs = rhs * 1/60;
time_lhs = lhs * 1/60;
time_ltoff = ltoff * 1/60;
time_rtoff = rtoff * 1/60;

%% calculating temporospatial parameters
%% Stride TIME
rstrdur = diff(time_rhs);
lstrdur = diff(time_lhs);

%% duration of step

rduration = [];
lduration = time_lhs - time_rhs;
i=1:11;
j=1:11;
rduration = time_rhs(i+1,:) - time_lhs(i,:);

%% Stance phase time and percent

lstance_phase_time = (time_ltoff(i+2,:) - time_lhs(i,:));

lstance_phase_perc = (lstance_phase_time / GCL(i,:))*100;


rstance_phase_time = (time_rtoff(j,:) - time_rhs(j,:));

rstance_phase_perc = (rstance_phase_time / GCR(j,:))*100;

%% Swing phase time and percent


lswing_phase_time = (time_lhs(i,:) - time_ltoff(i+1,:));

lswing_phase_perc = (lswing_phase_time / GCL(i,:))*100;


rswing_phase_time = (time_rhs(j+1,:) - time_rtoff(j,:));

rswing_phase_perc = (rswing_phase_time / GCR(j,:))*100;

%% Double support time and percent


indoublesu = (time_rtoff + 0.01666) - time_lhs(1:11,:);

fidoublesu = (time_ltoff(2:13,:) + 0.01666) - time_rhs;

totdoublesu = indoublesu + fidoublesu(2:12);

totdoublesu_perc = (totdoublesu / GCR(j,:))*100;


%% Cadence



s = max(time_ltoff) < max(time_rtoff);
if s == 1
    big = max(time_rtoff);
else 
    big = max(time_ltoff);
end
n = max(time_lhs) < max(time_rhs);
if n == 1 
    alsobig = max(time_rhs);
else
    alsobig = max(time_lhs);
end
w = big < alsobig;
if w == 1
    biggest = alsobig;
else
biggest = big;
end

cadence = ((count_lhs + count_rhs)/biggest)*60;

%% right stride length 
[subjectzxy,FTXT,FRAW] = xlsread('gaitreporters1.xlsx','Segment Velocity', 'A:BR' );



for i= 1: number_of_cycles_rhs;
    
    temp_var = strcat( 'gait_vel_',num2str(i) );
    eval(sprintf('%s = subjectzxy(rhs(i):rhs(i+1),:)',temp_var));
    
      
end

% 
% for j= 1:number_of_cycles_rhs
%        figure (1);hold on
%    plot(eval(['gait_vel_', num2str(j),'(:,20)' ]));
%    legendInfo{j} = ['gait' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)

clear j
for j= 1:number_of_cycles_rhs
%     [peaks, peaks_ind] = findpeaks(eval(['gait_vel_', num2str(j),'(:,49)' ]));
%     if max(eval(['gait_vel_', num2str(j),'(:,49)' ])) < 58
%          clearvars(['gait_vel_', num2str(j) ] );
    if eval(['gait_vel_', num2str(j),'(1,20)' ]) < 1
        clearvars(['gait_vel_', num2str(j) ] );
    elseif eval(['gait_vel_', num2str(j),'(70,20)' ]) < 1
        clearvars(['gait_vel_', num2str(j) ] );
%     elseif peaks(1)< 20
%         clearvars(['gait_', num2str(j) ] );
    end
end

% 
% hold on;
% for j= 4:10 %number_of_cycles
%        figure (1);
%    plot(eval(['gait_vel_', num2str(j),'(:,20)' ]));
%    legendInfo{j} = ['gait' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)

 
 gait_vel = [];


for j= 4:10 %number_of_cycles
      
   B =  eval(['gait_vel_', num2str(j),'(1:70,20)'] ) ;
   
    gait_vel = [gait_vel B];
    
    
end

       

gait_vel_mean = mean(gait_vel,2);
% 
% figure(20)
% plot(gait_vel_mean);
% title('Velocity mean')
% xlabel('GC (%)') 
% ylabel('joint angle (deg)')
vel_mean = mean(gait_vel_mean)

%rstridelength = (gait_vel_mean*120)/cadence


rstridelength = vel_mean * GCR

%% Left stride length
%lstridelength = (gait_vel_mean*120)/cadence

lstridelength = vel_mean * GCL


%% right swing velocity


for i= 1: number_of_cycles_rhs;
    
    temp_var = strcat( 'gait_vel_swi_',num2str(i) );
    eval(sprintf('%s = subjectzxy(rhs(i):rhs(i+1),:)',temp_var));
    
      
end

% 
% for j= 1:number_of_cycles_rhs
%        figure (1);hold on
%    plot(eval(['gait_vel_swi_', num2str(j),'(:,56)' ]));
%    legendInfo{j} = ['gait' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)

 
clear j
for j= 1:number_of_cycles_rhs
%     [peaks, peaks_ind] = findpeaks(eval(['gait_vel_', num2str(j),'(:,49)' ]));
%     if max(eval(['gait_vel_', num2str(j),'(:,49)' ])) < 58
%          clearvars(['gait_vel_', num2str(j) ] );
    if eval(['gait_vel_swi_', num2str(j),'(60,56)' ]) < 3
        clearvars(['gait_vel_swi_', num2str(j) ] );
%     elseif eval(['gait_vel_', num2str(j),'(70,20)' ]) < 1
%         clearvars(['gait_vel_', num2str(j) ] );
%     elseif peaks(1)< 20
%         clearvars(['gait_', num2str(j) ] );
    end
end


hold on;
for j= 3:10 %number_of_cycles
       figure (21);
       
   plot(eval(['gait_vel_swi_', num2str(j),'(30:70,56)' ]));
   legendInfo{j} = ['gait' num2str(j)];
  
end
 hold off;  
 title('Right velocity swing Plot')
xlabel('Right Gait Cycle') 
ylabel('Velocity')
 legend(legendInfo)
 %% left swing velocity 
 
 
for i= 1: number_of_cycles_lhs;
    
    temp_var = strcat( 'lgait_vel_swi_',num2str(i) );
    eval(sprintf('%s = subjectzxy(lhs(i):lhs(i+1),:)',temp_var));
    
      
end

% 
% for j= 1:number_of_cycles_lhs
%        figure (1);hold on
%    plot(eval(['lgait_vel_swi_', num2str(j),'(:,68)' ]));
%    legendInfo{j} = ['gait' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)

 
clear j
for j= 1:number_of_cycles_lhs
%     [peaks, peaks_ind] = findpeaks(eval(['gait_vel_', num2str(j),'(:,49)' ]));
%     if max(eval(['gait_vel_', num2str(j),'(:,49)' ])) < 58
%          clearvars(['gait_vel_', num2str(j) ] );
    if eval(['lgait_vel_swi_', num2str(j),'(60,68)' ]) < 2.5
        clearvars(['lgait_vel_swi_', num2str(j) ] );
%     elseif eval(['gait_vel_', num2str(j),'(70,20)' ]) < 1
%         clearvars(['gait_vel_', num2str(j) ] );
%     elseif peaks(1)< 20
%         clearvars(['gait_', num2str(j) ] );
    end
end


hold on;
for j= 2:10 %number_of_cycles
       figure (1);
   plot(eval(['lgait_vel_swi_', num2str(j),'(30:71,68)' ]));
   legendInfo{j} = ['gait' num2str(j)];
  
end
title('Left velocity swing Plot')
xlabel('Left Gait Cycle') 
ylabel('Velocity')
 hold off;  
 legend(legendInfo)
 
 %% Left step width
 
[subjectzxy,FTXT,FRAW] = xlsread('gaitreporters1.xlsx','Segment Velocity', 'A:BR' );

for i= 1: number_of_cycles_lhs;
    
    temp_var = strcat( 'lwgait_vel_',num2str(i) );
    eval(sprintf('%s = subjectzxy(lhs(i):lhs(i+1),:)',temp_var));
    
      
end

% 
% for j= 1:number_of_cycles_lhs
%        figure (1);hold on
%    plot(eval(['lwgait_vel_', num2str(j),'(:,66)' ]));
%    legendInfo{j} = ['gait' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)
  
clear j
for j= 1:number_of_cycles_lhs
%     [peaks, peaks_ind] = findpeaks(eval(['gait_vel_', num2str(j),'(:,49)' ]));
%     if max(eval(['gait_vel_', num2str(j),'(:,49)' ])) < 58
%          clearvars(['gait_vel_', num2str(j) ] );
    if eval(['lwgait_vel_', num2str(j),'(60,66)' ]) > 0.4
        clearvars(['lwgait_vel_', num2str(j) ] );
    elseif eval(['lwgait_vel_', num2str(j),'(43,66)' ]) > 0.5
        clearvars(['lwgait_vel_', num2str(j) ] );
    elseif eval(['lwgait_vel_', num2str(j),'(63,66)' ]) < 0.1
        clearvars(['lwgait_vel_', num2str(j) ] );
%     elseif eval(['gait_vel_', num2str(j),'(70,20)' ]) < 1
%         clearvars(['gait_vel_', num2str(j) ] );
%     elseif peaks(1)< 20
%         clearvars(['gait_', num2str(j) ] );
    end
end
% 
% hold on;
% for j= [2,3,5,7:10] %number_of_cycles
%        figure (1);
%    plot(eval(['lwgait_vel_', num2str(j),'(:,66)' ]));
%    legendInfo{j} = ['gait' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)

 
 lwgait_vel = [];


for j= [2,3,5,7:10] %number_of_cycles
      
   B =  eval(['lwgait_vel_', num2str(j),'(1:70,66)'] ) ;
   
    lwgait_vel = [lwgait_vel B];
    
    
end

       

lwgait_vel_mean = mean(lwgait_vel,2);
% 
% figure(2)
% plot(lwgait_vel_mean);
lw_vel_mean = mean(lwgait_vel_mean)

lstepwidth = lw_vel_mean * GCL



%% Right step width




[subjectzxy,FTXT,FRAW] = xlsread('gaitreporters1.xlsx','Segment Velocity', 'A:BR' );

for i= 1: number_of_cycles_rhs;
    
    temp_var = strcat( 'rwgait_vel_',num2str(i) );
    eval(sprintf('%s = subjectzxy(rhs(i):rhs(i+1),:)',temp_var));
    
      
end

% 
% for j= 1:number_of_cycles_rhs
%        figure (1);hold on
%    plot(eval(['rwgait_vel_', num2str(j),'(:,54)' ]));
%    legendInfo{j} = ['gait' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)
  
clear j
for j= 1:number_of_cycles_rhs
%     [peaks, peaks_ind] = findpeaks(eval(['gait_vel_', num2str(j),'(:,49)' ]));
%     if max(eval(['gait_vel_', num2str(j),'(:,49)' ])) < 58
%          clearvars(['gait_vel_', num2str(j) ] );
    if eval(['rwgait_vel_', num2str(j),'(1,54)' ]) > 0.3
        clearvars(['rwgait_vel_', num2str(j) ] );
    elseif eval(['rwgait_vel_', num2str(j),'(59,54)' ]) > 0.3
        clearvars(['rwgait_vel_', num2str(j) ] );
    elseif eval(['rwgait_vel_', num2str(j),'(59,54)' ]) < 0
        clearvars(['rwgait_vel_', num2str(j) ] );
%     elseif eval(['gait_vel_', num2str(j),'(70,20)' ]) < 1
%         clearvars(['gait_vel_', num2str(j) ] );
%     elseif peaks(1)< 20
%         clearvars(['gait_', num2str(j) ] );
    end
end
% 
% hold on;
% for j= [2:7,10,11] %number_of_cycles
%        figure (1);
%    plot(eval(['rwgait_vel_', num2str(j),'(:,54)' ]));
%    legendInfo{j} = ['gait' num2str(j)];
%   
% end
%  hold off;  
%  legend(legendInfo)

 
 rwgait_vel = [];


for j= [2:7,10,11] %number_of_cycles
      
   B =  eval(['rwgait_vel_', num2str(j),'(1:70,54)'] ) ;
   
    rwgait_vel = [rwgait_vel B];
    
    
end

       

rwgait_vel_mean = mean(rwgait_vel,2);
% 
% figure(2)
% plot(rwgait_vel_mean);
rw_vel_mean = mean(rwgait_vel_mean)

rstepwidth = rw_vel_mean * GCR
