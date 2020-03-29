function output=cross_section_generator(swwfile,maxspacing,num_xs,mesh_size)
%%% ------------------------Descriptions and notes-------------------------
%%% Inputs
%%% swwfile=file name of anuga output
%%% maxspacing=maximum distance between points 
%%% num_xs=number of cross-sections you would like to make
%%% mesh_size=thickness of mesh lines; 0 will shut the mesh off

%%% Outputs
%%% output=sturcture array containing the xy coorindates of the densified
%%% line clipped to the extent of the model domain

%%% Description
%%% This code will bring up a map view of elevation in the model domain.
%%% You can then click out cross-sections of interest. Left click adds
%%% points and right click indicates the last point in the cross-section.
%%% The code will then clip the data to the extent of the domain. If the
%%% code is clipping to much of the data try
%%%------------------------------------------------------------------------

%% Set up names for output stuctures
xs_name=cell(1,num_xs);
for j=1:num_xs
    xs_name{j}=['xs',num2str(j)];
end
clear j

%% Begin Analysis
%%%look into nc file and see what variables there are to extract
swwinfo=ncinfo(swwfile); 

%%% extract non temporal dependent protions of sww file
x=double(ncread(swwfile,'x'));%% x coordinates
y=double(ncread(swwfile,'y'));%% y  coordinates
elev=ncread(swwfile,'elevation'); %% time

xref=swwinfo.Attributes(8).Value; %% get x corner reference point
yref=swwinfo.Attributes(9).Value; %% get y corner reference point
x=x+xref; %% adjust x values to reference point
y=y+yref; %% adjust y values to reference point

%% Find the domain extent
bid=boundary(x,y);%% get index of outter most points
%%% Build polygon of the domain
bx=x(bid);%% outer x extent of domain
by=y(bid);%% outer y extent of domain

%% Plot up model domain and grid
if mesh_size==0 %% if line thickness is 0 set to none
    mesh_size=1;
    meshstyle='none';
else
    meshstyle='-';
end

figure
tri=delaunay([x y]);
trisurf(tri,x,y,elev,'linestyle',meshstyle,'linewidth',mesh_size);
view(2)
axis equal
hold on
%% Build the Crossections
for j=1:num_xs
    %%% Click out Cross-section
    title(['Draw cross-section ',num2str(j)])
    xy=ginput2(elev);  %%% see below function
    
    %%% Densify the cross-sections
    dX=diff(xy(:,1)); %% get x directed distance
    dY=diff(xy(:,2)); %% get y directed distance 
    d=sqrt(dX.^2+dY.^2); %% get distance between points
    ns=ceil(d./maxspacing)+1; %% find number of segments need

    xyout=[]; %% create a empty matrix
    for i=1:length(dX)
        xt=linspace(xy(i,1),xy(i+1,1),ns(i));%%linearly densify x direction
        yt=linspace(xy(i,2),xy(i+1,2),ns(i));%%linearly densify y direction
        xyout=[xyout;[xt' yt']]; %% build temp output
    end

%% Clip the data to the extent of the domain
    %%% find the points that fall in or edge of the domain
    [in,on]=inpolygon(xyout(:,1),xyout(:,2),bx,by);
    in=logical(in+on);%% comine inside and edge of domain
    xyout=xyout(in,:); %%keep only points in domain
    %%% plot up the new line
    plot3(xyout(:,1),xyout(:,2),ones(size(xyout(:,1))).*max(elev),'k','linewidth',2)

%% Write data to the ouput variable
    output.(xs_name{j})=xyout;
end



%% Additoanl functions 
function [xy]=ginput2(elev) %% function to click out cross sections
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')
but = 1; %% button command for ginput
ind=1; %% counter for index
while but == 1
    [xi,yi,but] = ginput(1); %%% records xy for each click left click but=1 right click but=0
    plot3(xi,yi,max(elev),'k.','markersize',20) %%% plots the location of your click
    hold on
    xy(ind,1) = xi; %% records x data 
    xy(ind,2) = yi; %% records y data
    ind=ind+1; %% counter for indexing

end
end

end