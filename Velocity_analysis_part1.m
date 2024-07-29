% Velocity calculation code - Part 1 - Optogenetic control of kinesins -1, -2, -3 and dynein reveals their specific roles in vesicular transport

%% NOTE: Save the "velocity_1", "Slid_time" manually from the Workspace! After using this code, proceed to "Velo_processing_intact.m" code to plot the velocity time series plot.

clear; close all; clc;
set(0,'DefaultFigureWindowStyle','docked');

%% A) Sourcing the code file path:
% set(0,'DefaultFigureWindowStyle','docked');
% %addpath('C:\Users\HendricksLab\Desktop\GitHub (Karthik)\Karthik_Codes\Final codes\Imp_Off-on-off - With Multiple Trajectory analysis\') % Directory containing codes
% addpath('C:\Users\karth\Dropbox\My PC (Karthikeyan-PC)\Desktop\Code and Sure internship project\Karthik_Codes\Karthik_Code')
% addpath('C:\Users\karth\Dropbox\My PC (Karthikeyan-PC)\Desktop\Code and Sure internship project\Karthik_Codes\Karthik_Code\STAC\'); % Single trajectory
% addpath('C:\Users\karth\Dropbox\My PC (Karthikeyan-PC)\Desktop\Code and Sure internship project\Karthik_Codes\Karthik_Code\MTAC\'); % Single cell - multitrajectory
% addpath('C:\Users\karth\Dropbox\My PC (Karthikeyan-PC)\Desktop\Code and Sure internship project\Karthik_Codes\Karthik_Code\MCMTAC\'); % Multi cell - multitrajectory


%% B) Variables (Change as per your conditions):
dT = 0.302; % Exposure time
WindowSize = 15; % Should be adjusted

%% C) Uploading data and initializing the variables and determining the number of trajectories
cd('/Users/sahilnagpal/Desktop/Test/') % Location of the data
File = 'Spots in tracks statistics_full cell' % Name of the .csv file generated from Trackmate (old version) - with x-coordinate of cell center in V1, and y-coordinate of cell center in W1
State = readmatrix(File);
LastValueTrackID = State(end,3)+1 % No. of trajectories of the cell in the chosen file
Delay= 1:ceil(WindowSize/2); % For MSD, Delay time: Gives the delay time from 1 to half of window size

% Initializing the variables:
X =[];
Y =[];
POSITION=[];
TIME=[];
velocity =[];

tic

%% Enter the timeseries values:
Quest = input('Do you want to manually time series for the conditions? 1 for yes, 2 for no : '); % Enter 1 to add the time series
% Manually entering the time series if you haven't mentioned the time
% series values in the xlsx file (for example, you are interested in
% analysing the trajectory from 0 to 120 seconds for condition 1/dark-state,
% and then from 140 to 260 seconds for condition 2/lit-state).

if Quest == 1 % For this program, put the lower time series as 0
    disp('How many transistions do you want to consider?');
    Transition_no = input('Enter the number of conditions (either 1, 2 or 3): ');
    for tn = 1:Transition_no
        Time_S_Trans1 = input('Enter lower time series: ');
        Time_S_Trans2 = input('Enter upper time series: ');
        Time_Pul_Trans{tn} = [Time_S_Trans1; Time_S_Trans2];
        times_pulse{tn} = Time_Pul_Trans{tn};
    end
else
    Time_Pul_Trans = 0; % If the user already has mentioned the transition series values in the xlsx file, then they can mention 0 in the command prompt.
end

Origin_X = State(1,22); % Cell Origin in X Coordinate
Origin_Y = State(1,23); % Cell Origin in Y Coordinate

New_OriX = 0; % Initializing the variable for the trajectory's X Origin
New_OriY = 0; % Initializing the variable for the trajectory's Y Origin


%% Calculating the X, Y, Position and time:
% figure; hold on;
for ks = 1:LastValueTrackID %6800:6900 %1:LastValueTrackID %LastValueTrackID

    %% Determining X, Y and Time from the trackmate data:

    jkd=find(State(:,3)==(ks-1)); % TrackID
    X = State(jkd,5); % X Coordinates of that particular TrackID
    Y = State(jkd,6); % Y Coordinates of that particular TrackID
    Time = State(jkd,8); % Time Coordinates of that particular TrackID

    %% Ignoring the trajectories whose data points are less than 5 values:
    if numel(X)<=5 && X(1)~=X(2) && Y(1)~=Y(2) % NOTE: Considering only those trajectories having datapoints greater than 5.
        Velocity = [];
        Position_1= {};
        time_1= {};
        X_1={};
        Y_1={};
        X_n={};
        Y_n={};
    end


    %% Calculating the position (Magnitude - Vector length):
    dx_origin=(X-Origin_X);
    dy_origin=(Y-Origin_Y);
    mag_vect=sqrt(dx_origin.^2 + dy_origin.^2); % Finding the magnitude (from origin to trajectory first)

    indx = min(find(mag_vect == min(mag_vect))); % Finding the index which shows less magnitude indicating that it is closer to the origin

    %Sometimes there is a trajectory with same magnitude for different indices. So used min function to choose the first point closest to cell centre

    New_OriX = X(indx); % Using the index, determine the trajectory X Coordinate which is the closest to the Origin X Coordinate
    New_OriY = Y(indx); % Same with Y coordinates as well.
    g = numel(X);
    b = numel(New_OriX);


    % Finding the new distance from the trajectory origin point to the other points in trajectory.
    New_dx=(X - New_OriX); 
    New_dy=(Y - New_OriY);

    Position = sqrt(New_dx.^2 + New_dy.^2); % Vector length

    %% Identifying X, Y, Position and Time based on the time series the user inputted:
    kpr = 0;
    for tn = 1:numel(times_pulse) % Checking number of columns
        kpr=kpr+1;
        tp1=times_pulse{tn}; % i has 2 row values (as it is in cell format) --> 1st row: Time in which the condition started (eg: 0, 135, 270) and 2nd Row: Time in which the condition ended (eg: 120, 255, 390)
        Position_1{kpr}=Position(Time>=tp1(1) & Time<tp1(2)); % Finding positions that are between the condition starting state and condition ending state
        time_1{kpr}=Time(Time>=tp1(1) & Time<tp1(2)); % Finding time frames that are between the condition starting state and condition ending state
        X_1{kpr}=X(Time>=tp1(1) & Time<tp1(2)); % Finding positions that are between the condition starting state and condition ending state
        Y_1{kpr}=Y(Time>=tp1(1) & Time<tp1(2));
    end

%% If you want to divide the trajectory into 3 conditions, follow this:
    Position_2{ks} = Position_1;
    X_2{ks} = X_1;
    Y_2{ks} = Y_1;
    Time_2{ks} = time_1;

%% If you want to have the trajectory intact without dividing, then follow this:
    Position_3{ks} = Position;
    Time_3{ks} = Time;
    
end


%% Calculating velocity (intact trajectories):
for mu = 1:LastValueTrackID % 1 to total # of trajectories
    mu
    if isempty(Position_3{mu}) % if the code finds any position matrix as empty, then the velocity for that trajectory is also empty
        velocity_1{mu} =  [];
        velocity_trans{mu} = [];
        Slid_time{mu} = [];

    else
        % Each trajectory has 3 conditions: Off1, On and Off2.
            Posit = Position_3{mu};
            T = Time_3{mu};
            NofWindows = 0;
            jk = 0;
            NofWindows = numel(Posit) - WindowSize;
            N_Total = NofWindows;
            if isempty(Posit) 
                isempty(Posit);
                Velocity{mu} ={};
                Sliding_Time{mu} ={};
                newPosit= {};
                newTime = {};
            elseif NofWindows <= 1
                Velocity{mu} ={};
                Sliding_Time{mu} ={};
                newPosit= {};
                newTime = {};
            else
                NofWindows;
                for WinSz = 1:NofWindows
                    WinSz;
                    jk = WinSz:(WinSz+WindowSize);
                    newPosit = Posit(jk);
                    newTime = T(jk);
                    meanTime(WinSz) = mean(newTime);
                    V = polyfit(newTime, newPosit, 1);
                    Velocity_2(WinSz) = V(1); % Identifying velocity using the slope
                end
                velocity_1{mu}=Velocity_2'; % storing separate velocities for the 3 conditions 
                Slid_time{mu}=meanTime'; % storing separate times for the 3 conditions
                clear Velocity_2 meanTime;
        end
    end


end


toc