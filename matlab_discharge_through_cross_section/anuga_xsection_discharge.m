function output=anuga_xsection_discharge(swwfile,xs,doplotxs)

%%% ------------------------Descriptions and notes-------------------------
%%% Inputs
%%% swwfile=file name of anuga output
%%% xs=number of cross-sections you would like to make
%%% doplotxs=True of False: Plots up input xs-section and remapped xsection

%%% Outputs
%%% output=sturcture array containing the Discharge and time for each cross
%%% section column 1 is time and column 2 is discharge 

%%% Description
%%% Code computes the discharge through cross-sections for anuga ouput
%%% files. It can handle multiple cross sections and outputs a structure
%%% containing the time and discharge for each cross-section. The code
%%% remaps the cross-sections to the model grid and computes the flux
%%% perpendicular to each cell face. 
%%%------------------------------------------------------------------------

field_id=fieldnames(xs); %% get names of x-sections
for i=1:numel(field_id)
    %% Extract the data for each x-section
    xspoints=xs.(field_id{i});%% get xy points for given cross-section
   
    %% Load in coordiantes for .sww file
    %%% extract non temporal dependent protions of sww file
    swwinfo=ncinfo(swwfile); %% look into nc file and see what variables there are to extract
    x=double(ncread(swwfile,'x'));%% x coordinates
    y=double(ncread(swwfile,'y'));%% y  coordinates
    t=ncread(swwfile,'time'); %% time
    xref=swwinfo.Attributes(8).Value; %% get x corner reference point
    yref=swwinfo.Attributes(9).Value; %% get y corner reference point
    x=x+xref; %% adjust x values to reference point
    y=y+yref; %% adjust y values to reference point


    %% find corresponding mesh locations to nearest xsection
    kn=dsearchn([x,y],xspoints);%% find points closest in xsection
    %Optional plotting
    if doplotxs == true
        warning off
        figure
        tri=delaunay([x y]);
        triplot(tri,x,y);
        axis equal
        hold on
        plot(xspoints(:,1),xspoints(:,2),'k','linewidth',2)
        plot(x(kn),y(kn),'-or','linewidth',2)
        title('Remapped x-section to cell faces')
        legend('mesh','xsection','remapped xsection')
    end

    x=x(kn); %%keep only x coordinates of interest
    y=y(kn); %%keep only y coordinates of interst

    %% Extract x and y directed data for points of interest
    %%% pre-allocate mememory for x&y directed specific discharge
    qxt=ones(length(kn),length(t)).*NaN;
    qyt=ones(length(kn),length(t)).*NaN;

    %%% Extract data from .sww file
    for j=1:length(t)%% stepping through each timestep to help with huge .sww files
        %%get x momentum data-->x directed specific discharge or h*u    
        xmomt=double(ncread(swwfile,'xmomentum',[1,j], [inf,1]));
        %%get y momentum data-->y directed specific discharge h*v
        ymomt=double(ncread(swwfile,'ymomentum',[1,j], [inf,1])); 
        %%% write specific discharge to a temp variable will be edited in the
        %%% next section
        qxt(:,j)=xmomt(kn,:);%% combine x dir specific Q of interest
        qyt(:,j)=ymomt(kn,:); %% combine y dir speific Q  of interest
        clear xmomt ymomt    
    end

    %%% Possibly explore how to use indexing to extract the data instead of loop 
    %%% This will increase spead and reduce the amount of ram needed 

    %% Finding unique data points
    %%% Setting store vertices uniquely causes duplicate points here I am
    %%% finding those locations taking the median qx and qy vaues

    xyt=[x,y];%% temp variable to house data will be empty after next loop
    xy=ones(size(xyt)).*NaN; %%% variable to store all the 
    qx=ones(size(qxt)).*NaN; %%% variable to store all the 
    qy=ones(size(qyt)).*NaN; %%% variable to store all the 

    c=1;%% counter

    while ~isempty(xyt)
        z=xyt(1,1:2)-xyt(:,1:2); %%% reference everything to a point of interst
        idxy=find(abs(z(:,1))+abs(z(:,2))==0); %%%% 0 indcates same point

        %%% Get the unique points and median specific discharge at each
        xy(c,:)=xyt(idxy(1),:);
        qx(c,:)=median(qxt(idxy,:));
        qy(c,:)=median(qyt(idxy,:));

        %%% Remove processed data form temp variables
        xyt(idxy,:)=[];    
        qxt(idxy,:)=[];
        qyt(idxy,:)=[];
        c=c+1;%% update counter
    end
    xy(any(isnan(xy),2),:)=[]; %%% delete leftover nans from memory allocation
    qx(any(isnan(qx),2),:)=[]; %%% delete leftover nans from memory allocation
    qy(any(isnan(qy),2),:)=[]; %%% delete leftover nans from memory allocation

    clear xyt qxt qyt %%% clear out the now empty values
    %% Vector decomposition to get x and y directed discharge accross each face

    dxy=diff(xy); %% get difference between x&y point locations
    d=sqrt(dxy(:,1).^2+dxy(:,2).^2); %% get length of each face


    %%%Set average xy specific discharge for each line segment
    qxavg=(qx(1:end-1,:)+qx(2:end,:))./2;
    qyavg=(qy(1:end-1,:)+qy(2:end,:))./2;

    %%%Decomposition of unit vector along parallel to cross section (transverse to flow)         
    parx=dxy(:,1)./d;
    pary=dxy(:,2)./d;

    %%%Decomposition of unit vector in streamwise direction or perpendicular of cross section
    perx=pary;
    pery=-parx;

    %%%Solve for fluxes perpendicular to the cross-section        
    % solve for qx  perpendicular to line segment
    perqx=((qxavg.*perx)+(qyavg.*pery))./(perx.^2+pery.^2).*perx;
    % solve for qy perpendiculat to line segment
    perqy=((qxavg.*perx)+(qyavg.*pery))./(perx.^2+pery.^2).*pery;
    % solve for magnitude
    qmag=sqrt(perqx.^2+perqy.^2);
    Q_seg=qmag.*d; %% convert specific discharge to discharge for each seg.
    Q=sum(Q_seg); %% sum up segments to get a totale discharge

%% Write data to output file
    output.(field_id{i})=[t,Q'];
end

