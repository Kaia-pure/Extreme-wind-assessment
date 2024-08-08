%% Read data and extract the monthly maximum
clc;clear;
n=0;
year=2023;
% for year=1940:1958
    for month=1:12
        n=n+1;
        file1{n,:}=char(['G:\ERA5\v_wind100m\',num2str(year),num2str(month,'%02d'),'v_wind.nc']);
        file2{n,:}=char(['G:\ERA5\u_wind100m\',num2str(year),num2str(month,'%02d'),'u_wind.nc']);
    end
% end
A = nan(1440,681,12);
for y1=21:701
    y1
    %read the full WS array to do Weibull fitting.
    %read also maximum values to compute Annual maxima
    parfor n=1:12
        A(:,y1,n)=max(sqrt(double(ncread(char(file1{n,:}),'v100',[1 y1 1],[1440 1 Inf])).^2+ ....
            +double(ncread(char(file2{n,:}),'u100',[1 y1 1],[1440 1 Inf])).^2),[],3);
    end
end  
    save('E:\ExtremeWind\Max_2023.mat','A')

%% Evaluate U50 and the significance
clc;clear;
A = load('E:\ExtremeWind\Ex19402022\Max_19402023.mat'); A = A.A;
latSize=681;
NorAnnmax=NaN(83,1440,latSize);

for year=1:84
    NorAnnmax(year,:,:)=max(A(:,:,(((year-1)*12+1):year*12)),[],3);
end
R50=NaN(1440,latSize,45,3);
U50=NaN(1440,latSize,3);
k=1;
for y1=1:latSize
    for x1=1:1440
%         %GEV-Gumbel paramaters from maximum likelihood fitting:
%         pd=fitdist(NorAnnmax(1:84,x1,y1),'ExtremeValue');
%         U50(x1,y1,1)=pd.ParameterValues(2)*log(-1*log(0.02))+pd.ParameterValues(1);
%         pci1=paramci(pd,'Parameter','mu');
%         pci2=paramci(pd,'Parameter','sigma');
%         U50(x1,y1,2)=pci2(1)*log(-1*log(0.02))+pci1(1);
%         U50(x1,y1,3)=pci2(2)*log(-1*log(0.02))+pci1(2);
        %estimates for 20-year blocks;
%         for k=1:34
%             pd=fitdist(NorAnnmax((1+(k-1):30+(k-1)),x1,y1),'ExtremeVa2lue');
            pd=fitdist(NorAnnmax(65:84,x1,y1),'ExtremeValue'); %65:84
            R50(x1,y1,k,1)=pd.ParameterValues(2)*log(-1*log(0.02))+pd.ParameterValues(1);
            pci3=paramci(pd,'Parameter','mu');
            pci4=paramci(pd,'Parameter','sigma');
            R50(x1,y1,k,2)=pci4(1)*log(-1*log(0.02))+pci3(1);
            R50(x1,y1,k,3)=pci4(2)*log(-1*log(0.02))+pci3(2);
%             zalpha=1.96;
%             sigma=(pd.ParameterValues(1)-pci3(1))/sqrt(20)*zalpha;
%         end
    end
    disp(y1)
end
% save('E:\ExtremeWind\U50estimates19402023.mat','U50');
save('E:\ExtremeWind\Ex19402022\R50estimates20yr200402023.mat','R50'); %20042023

%% The trend of U50 and the significance
clc;clear;
latSize=681;
R1=load('E:\ExtremeWind\Ex19402022\R50estimates30yr19402023.mat');R1=R1.R50;
R1=R1(:,:,:,1);


parameter = zeros(1440,latSize,4); 
yr = (1940:1:1994)';

for i = 1:1440
    for j= 1:latSize
        stats = regstats(squeeze(R1(i,j,:,1)),yr);
        parameter(i,j,1) = stats.beta(2);
        parameter(i,j,2) = stats.beta(1);
        parameter(i,j,3) = stats.rsquare;
        parameter(i,j,4) = stats.tstat.pval(2);
    end
end
save('E:\ExtremeWind\Ex19402022\Parameter30yr19402023.mat','parameter');

%% Data visualization of U50 
clc;clear;
U50=load('E:\ExtremeWind\Ex19402022\U50estimates19402022.mat');U50=U50.U50;
LAT=load('E:\ExtremeWind\Ex19402022\LAT.mat');LAT=LAT.LAT;
LON=load('E:\ExtremeWind\Ex19402022\LON.mat');LON=LON.LON;
land=readgeotable('E:\ExtremeWind\land\land_merge.shp');
farm=readgeotable('E:\ExtremeWind\China\China.shp');
LON1=nan(1440,681);
LAT1=nan(1440,681);
LON2=nan(1440,681);
LAT2=nan(1440,681);
LON3=nan(1440,681);
LAT3=nan(1440,681);
for i=1:1440
    for j=1:681
        if U50(i,j,1)>37.5
            LON1(i,j)=LON(i,j);
            LAT1(i,j)=LAT(i,j);
        end
        if U50(i,j,1)>42.5
            LON2(i,j)=LON(i,j);
            LAT2(i,j)=LAT(i,j);
        end
        if U50(i,j,1)>50
            LON3(i,j)=LON(i,j);
            LAT3(i,j)=LAT(i,j);
        end
    end
end

sz=15; 
U50=U50(:,:,1);
fig = figure('Units','centimeters','Position',[5,5,20,12.5]);
% axesm('MapProjection','robinson','MapLatLimit',[-85 85],'MapLonLimit',[-180 180],'Frame','off');
axesm('MapProjection','miller','MapLatLimit',[3 45],'MapLonLimit',[105 150],'Frame','off');

h=framem; 
set(h,'LineWidth',2);
pcolorm(LON,LAT,U50(:,:));
% set(j,'alphadata',~isnan(P20yr(:,:)))
load coast; 
plotm(lat,long,'-','Color','k','LineWidth',1);
cb=colorbar();
set(get(cb,'title'),'string','m s^{-1}','FontSize',14);
clim([0 45])
plotm(lat,long,'-','Color',[0.4 0.4 0.4],'LineWidth',1)

h1=scatterm(LON1,LAT1,7,[145, 0, 72]/255,'filled');
h2=scatterm(LON2,LAT2,7,[151, 130, 75]/255,'filled');
h3=scatterm(LON3,LAT3,7,[102, 102, 102]/255,'filled');
% legend([h1,h2,h3],'U_{50} > 37.5','U_{50} > 42.5','U_{50} > 50')

%set(gca,'visible','off');
% title('Trend of U50 (every 20 years)',FontSize=sz)
axis off 
geoshow(land,'FaceColor',[0, 53, 102]/255)
 geoshow(farm,'DisplayType','point')

%  print(fig,'E:\ExtremeWind\U50 621.png','-dpng','-r600')

%% Data visualization of U50 trend and p value
clc;clear;
P20yr=load('E:\ExtremeWind\Ex19402022\Parameter30yr.mat'); P20yr=P20yr.parameter;
LAT=load('E:\ExtremeWind\Ex19402022\LAT.mat');LAT=LAT.LAT;
LON=load('E:\ExtremeWind\Ex19402022\LON.mat');LON=LON.LON;
index=(P20yr(:,:,4)>0.05);
land=readgeotable('E:\ExtremeWind\land\land_merge.shp');
farm=readgeotable('E:\ExtremeWind\China\China.shp');


% map = shaperead('E:\ExtremeWind\land\land_merge.shp','UseGeoCoords',true);
% [in, on] = inpolygon(lon(i,1),lat(i,1),[map(j).Lon],[map(j).Lat]); %判断网格点是否在多边形内

sz=15;
LON1=nan(1440,681,3);
LAT1=nan(1440,681,3);
for i=1:1440
    for j=1:681
        if index(i,j)==1
            LON1(i,j)=LON(i,j);
            LAT1(i,j)=LAT(i,j);
        end
    end
end
P20yr=P20yr(:,:,1)*10;
a=sum(sum(P20yr>0));
fig = figure('Units','centimeters','Position',[5,5,35,20]);
kk=1; ax(kk) = axes; 
axesm('MapProjection','robinson','MapLatLimit',[-85 85],'MapLonLimit',[0 360],'Frame','off');
% axesm('MapProjection','miller','MapLatLimit',[3 45],'MapLonLimit',[105 150],'Frame','off');

h=framem; 
set(h,'LineWidth',2);
pcolorm(LON,LAT,P20yr(:,:));
% set(j,'alphadata',~isnan(P20yr(:,:)))
load coast; 
plotm(lat,long,'-','Color','k','LineWidth',1);
scatterm(LON1(:,:,1),LAT1(:,:,1),0.1,[0.5 0.5 0.5])
%set(gca,'visible','off');
% title('Trend of U50 (every 20 years)',FontSize=sz)
cb=colorbarpwn(ax(1),-3.5,3.5,'colorP',[218, 85, 82]/255,'colorN',[0, 106, 163]/255);
set(get(cb,'title'),'string','m s^{-1} decade^{-1}','FontSize',14);
axis off 
plotm(lat,long,'-','Color',[0.4 0.4 0.4],'LineWidth',1);
geoshow(land,'FaceColor',[0, 53, 102]/255,'EdgeColor',[0.5 0.5 0.5])
loc_here = [0 90 100 180 ; 0 90 180 -85; 0 90 -90 30; 0 -90 -70 10]; 
loc_here = [0 90 100 180 ; 0 90 180 -85; 0 90 -90 30; 0 -90 -70 10]; 
for k1 = 1:2
    lat_na = [linspace(loc_here(k1,1),loc_here(k1,2),50), loc_here(k1,2),linspace(loc_here(k1,2),loc_here(k1,1),50), loc_here(k1,1)];
    lon_na = [ones(1,50)*loc_here(k1,3), loc_here(k1,4), ones(1,50)*loc_here(k1,4), loc_here(k1,3)];
    plotm(lat_na,lon_na,'k-','linewidth',1.5,'color',[118 90 124]/255)
end
%  geoshow(farm,'DisplayType','point')

%  print(fig,'E:\ExtremeWind\Figure\U50 Trend every 30.png','-dpng','-r600')
%% Extract the location of wind farms and Visualization
clc;clear;
windfarm=readgeoraster('E:\ExtremeWind\offshore__windfarm\assemble_OWF_proj.tif');
limit=[67.1375038413,-130.148471662,162.817181591,-43.0078918004];
psz=0.039072506;

latindex_farm=zeros(size(windfarm,1),1);
lonindex_farm=zeros(size(windfarm,2),1);
for i=1:size(windfarm,1)
    latindex_farm(i,1)=-(limit(1)-i*psz)/0.25+341;
end
for j=1:size(windfarm,2)
    if j*psz+limit(2)<0
            lonindex_farm(j,1)=(j*psz+limit(2)+360)/0.25;
    end
    if j*psz+limit(2)>0
    lonindex_farm(j,1)=(j*psz+limit(2))/0.25;
    end
end
latindex_farm=ceil(latindex_farm)-1;
lonindex_farm=ceil(lonindex_farm);
[farmlon,farmlat]=meshgrid(latindex_farm,lonindex_farm);%farmlon是纬度信息,farmlat是经度信息

windfarm=windfarm';
farmlat(windfarm==0)=nan;
farmlon(windfarm==0)=nan;
farmlat=reshape(farmlat,[size(farmlat,1)*size(farmlat,2),1]);farmlat(isnan(farmlat))=[];
farmlon=reshape(farmlon,[size(farmlon,1)*size(farmlon,2),1]);farmlon(isnan(farmlon))=[];
loc(:,1)=farmlat;
loc(:,2)=farmlon;

P20yr=load('E:\ExtremeWind\Ex19402022\Parameter20yr.mat'); P20yr=P20yr.parameter;
LAT=load('E:\ExtremeWind\Ex19402022\LAT.mat');LAT=LAT.LAT;
LON=load('E:\ExtremeWind\Ex19402022\LON.mat');LON=LON.LON;
farmtrend=zeros(1440,681);
for j=1:size(loc,1)
    farmtrend(loc(j,1),loc(j,2))=P20yr(loc(j,1),loc(j,2),1);
end

% Europe
fig = figure('Units','centimeters','Position',[5,5,20,15]);
kk=1; ax(kk) = axes; 
axesm('MapProjection','robinson','MapLatLimit',[49 70],'MapLonLimit',[-17 25],'Frame','off');
h=framem; 
set(h,'LineWidth',2);
pcolorm(LON,LAT,farmtrend(:,:));
load coast; 
plotm(lat,long,'-','Color','k','LineWidth',1);
colorbarpwn(ax(1),-0.1,0.1,'colorP',[218, 85, 82]/255,'colorN',[0, 106, 163]/255);
axis off %关闭外部坐标轴,外部坐标轴不同于map axes
plotm(lat,long,'-','Color',[0.4 0.4 0.4],'LineWidth',1);

%Asia
fig = figure('Units','centimeters','Position',[5,5,20,15]);
kk=1; ax(kk) = axes; 
axesm('MapProjection','miller','MapLatLimit',[10 45],'MapLonLimit',[105 150],'Frame','off');h=framem;
set(h,'LineWidth',2);
pcolorm(LON,LAT,farmtrend(:,:));
load coast; 
plotm(lat,long,'-','Color','k','LineWidth',1);
colorbarpwn(ax(1),-0.1,0.11,'colorP',[218, 85, 82]/255,'colorN',[0, 106, 163]/255);
axis off 
plotm(lat,long,'-','Color',[0.4 0.4 0.4],'LineWidth',1);

% America
fig = figure('Units','centimeters','Position',[5,5,20,15]);
kk=1; ax(kk) = axes; 
axesm('MapProjection','robinson','MapLatLimit',[32 45],'MapLonLimit',[-83 -65],'Frame','off');
h=framem; 
set(h,'LineWidth',2);
pcolorm(LON,LAT,farmtrend(:,:));
load coast; 
plotm(lat,long,'-','Color','k','LineWidth',1);
colorbarpwn(ax(1),-0.1,0.1,'colorP',[218, 85, 82]/255,'colorN',[0, 106, 163]/255);
axis off 
plotm(lat,long,'-','Color',[0.4 0.4 0.4],'LineWidth',1);

% Australia
fig = figure('Units','centimeters','Position',[5,5,20,15]);
kk=1; ax(kk) = axes; 
axesm('MapProjection','robinson','MapLatLimit',[-46 -31],'MapLonLimit',[140 154],'Frame','off');
h=framem; 
set(h,'LineWidth',2);
pcolorm(LON,LAT,farmtrend(:,:));
load coast; 
plotm(lat,long,'-','Color','k','LineWidth',1);
colorbarpwn(ax(1),-0.1,0.1,'colorP',[218, 85, 82]/255,'colorN',[0, 106, 163]/255);
axis off 
plotm(lat,long,'-','Color',[0.4 0.4 0.4],'LineWidth',1);

% histogram for all wind farms
farmtrend=nan(1440,681);
for j=1:size(loc,1)
    farmtrend(loc(j,1),loc(j,2))=P20yr(loc(j,1),loc(j,2),1);
end
farmtrend=reshape(farmtrend,[1440*681,1]);
farmtrend(isnan(farmtrend))=[];
histogram(farmtrend,'Normalization','probability')
xlabel('Trend of U50 (m/s)')
title('Trend of U50 in wind farms')
