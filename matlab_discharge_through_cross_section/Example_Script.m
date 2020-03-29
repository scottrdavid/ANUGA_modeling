%% Example Script showsing how to call functions to calculate discharge 
%%% through a cross-section for ANUGA model outputs
clear all;close all; clc
%% Input data from ANUGA--this simulation Q=100m^3/s
swwfile=[pwd,'\topo.sww'];

%% Example 1--Make the cross sections using the cross section generator tool
%%% This will bring up a figure that you can click out cross-sections on 
%%% the code will clip the cross-section to the domain extent and densify
%%% the line. 

%%% Cross-section options
maxspacing=5; %% max spacing of points 
num_xs=2; %% number of cross-sections

%%% Plotting Options
%line thickness mesh when plotting. Can set to 0 to not show the mesh
mesh_size=0.5;

%%% Generate cross-section: left click to pick points; right click to
%%% choose last point on a cross-section
XS=cross_section_generator(swwfile,maxspacing,num_xs,mesh_size);

%%% Compute the time varying discharge for each cross-section
doplotxs=true;%% show remapping of input cross-section to the grid
T_Q_ex1=anuga_xsection_discharge(swwfile,XS,doplotxs);

%%% Plot up the results
field_id=fieldnames(T_Q_ex1); %% get names of x-sections
figure; hold on

for i=1:numel(field_id) %% loop through cross-sections
    time_discharge=T_Q_ex1.(field_id{i}); %% get time and discharge
    plot(time_discharge(:,1),time_discharge(:,2)) %% plot time vs discharge
end
xlabel('Time (s)')
ylabel('Discharge (m^3/s)' )
legend(field_id) %% add legend with cross section ids

%% Example 2-- Load in cross-sections from shapefile and densify the lines
%%% The shp2xs will also clip the line to the domain extent

maxspacing=10; %%% point density
%%% note do NOT include the extension on the file name for the shapefile
shp=[pwd,'\GIS_data\xs']; %% shapefile that hasnt had a point density set
shp_xs=shp2xs(swwfile,shp,maxspacing);

%%% Compute the time varying discharge for each cross-section
doplotxs=false;%% show remapping of input cross-section to the grid
T_Q_ex2=anuga_xsection_discharge(swwfile,shp_xs,doplotxs);

%% Example 3-- Similar to example 2 but now set maxspacing to 0 as to not
%%% densify the line and use the density I set in ArcMap

maxspacing=0; %%% point density
shp=[pwd,'\GIS_data\xsdensified']; %% shapefile where point density was set to 10m in Arcmap
shp_xs=shp2xs(swwfile,shp,maxspacing);

%%% Compute the time varying discharge for each cross-section
doplotxs=false;%% show remapping of input cross-section to the grid
T_Q_ex3=anuga_xsection_discharge(swwfile,shp_xs,doplotxs);


%% Example 4 Inputing your own cooridates for cross-sections
%%% write some cross-sections to a structure
xs.xs1=[ones(1,200)*800;1:200]';
xs.xs2=[ones(1,200)*400;1:200]';
xs.xs3=[ones(1,200)*450;1:200]';
doplotxs=false;
Q_t_ex4=anuga_xsection_discharge(swwfile,xs,doplotxs);
