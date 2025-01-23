%
% Take the C/C++ simulation output files to generate kymograph of traction
% stress and doublet length over time
%
%
%  Created on 01 Jan 2025
%          by Ying Zhang (ying1.zhang@northeastern.edu)
%
%

clf; 
close all;

% Read in simulation data
FileDir = '/YOUR/DATA/SAVING/DIRECTORY/';

filename1 = [FileDir,'morse_cell1_all_1_1_1_0.16_0.008_10.5_1.04_6_20_0.75_4.txt'];
filename2 = [FileDir,'morse_cell1_surface_1_1_1_0.16_0.008_10.5_1.04_6_20_0.75_4.txt'];
filename3 = [FileDir,'morse_cell2_all_1_1_1_0.16_0.008_10.5_1.04_6_20_0.75_4.txt'];
filename4 = [FileDir,'morse_cell2_surface_1_1_1_0.16_0.008_10.5_1.04_6_20_0.75_4.txt'];

c1_all = dlmread(filename1,'',6);
c1_wall = dlmread(filename2,'',6);
c2_all = dlmread(filename3,'',6);
c2_wall = dlmread(filename4,'',6);

% Number of discrete points per cell boundary and total number of frames
M = 162;
Nframes = length(c1_all)/M;

% Start generating the figure
figure(1)
subplot(3,2,[1 2 3 4])
hold all
for( i = 1:Nframes )

    startWPos = (i-1)*M + 1;
    endWPos = i*M;

    c1_x = c1_all(startWPos:endWPos,2);
    c1_vx = c1_all(startWPos:endWPos,24);
    c1_length(i) = max(c1_x)-min(c1_x);
    c1_xfront(i) = max(c1_x);
    c1_centroid(i) = min(c1_x)+0.5*(max(c1_x)-min(c1_x));

    c2_x = c2_all(startWPos:endWPos,2);
    c2_vx = c2_all(startWPos:endWPos,24);
    c2_length(i) = max(c2_x)-min(c2_x);
    c2_xback(i) = min(c2_x);
    c2_centroid(i) = min(c2_x)+0.5*(max(c2_x)-min(c2_x));

    ctotal_length(i) = max(c1_x)-min(c2_x);
    pair_centroid(i) = min(c2_x)+0.5*(max(c1_x)-min(c2_x));

    c1_xw = c1_wall(startWPos:endWPos,2);
    c1_fxw = c1_wall(startWPos:endWPos,4);
    c2_xw = c2_wall(startWPos:endWPos,2);
    c2_fxw = c2_wall(startWPos:endWPos,4);

    idxw = find(c2_fxw~=0);
    xw_min = min(c2_xw(idxw));
    if(isempty(xw_min))
        dx = 0;
    else
        dx = abs(xw_min-c2_xback(i));
    end
    c1_xw = c1_xw+1;
    c2_xw = c2_xw+1;

    c1_a = [c1_xw c1_fxw];
    c1_a = sortrows(c1_a);
    c1_F = repmat(c1_a(:,2),[1,2]);
    c1_XW = c1_a(:,1);
    c1_kkk = find(c1_F(:)==0);
    c1_FF = c1_F;
    c1_FF(c1_kkk) = NaN;

    c2_a = [c2_xw c2_fxw];
    c2_a = sortrows(c2_a);
    c2_F = repmat(c2_a(:,2),[1,2]);
    c2_XW = c2_a(:,1);
    c2_kkk = find(c2_F(:)==0);
    c2_FF = c2_F;
    c2_FF(c2_kkk) = NaN;
    
    % Making kymographs
    pcolor(c2_XW,(60/33.55705)/60*[i-1 i],c2_FF');
    hold on
    pcolor(c1_XW,(60/33.55705)/60*[i-1 i],c1_FF'); 
    hold on   
    
end

% Making kymographs
plot(c1_xfront,((1:1:Nframes-1)/33.55705)*60/60, '-w','linewidth',1);
hold on
plot(c2_xback,((1:1:Nframes-1)/33.55705)*60/60,'-w','linewidth',1);
hold on
plot(pair_centroid,((1:1:Nframes-1)/33.55705)*60/60,'-w','linewidth',1);
xlim([30 inf]);
ylim([200/60 (Nframes/33.55705)*60/60])

bipolarMap = generate_bipolar_colormap(256, 0, 'cubic');
colormap(bipolarMap)
clim([-0.2 0.2])
set(gca, 'color', [10 10 10]/255)
cb = colorbar;
cb.Position = [0.819727891156463,0.430437418522362,0.056698988052932,0.089463078990076];
shading interp;
xlabel('Distance ($\mu$m)','FontSize',20, 'Interpreter', 'latex');
view(-90,90); 
axis ij;
box on
set(gca,'FontSize',20);
set(gca,'FontName','Times New Roman')
set

% Plotting doublet length over time
subplot(3,2,[5 6])
plot(((1:1:Nframes)/33.55705)*60/60,ctotal_length(1:1:Nframes),'-k','linewidth',2);
ylim([35 75])
Tf = (Nframes/33.55705)*60/60;
xlim([200/60 (Nframes/33.55705)*60/60])
grid on
ylabel('Length ($\mu$m)', 'Interpreter', 'latex', 'FontSize', 20);
xlabel('Time (min)','Interpreter','latex','FontSize', 20);
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',20);