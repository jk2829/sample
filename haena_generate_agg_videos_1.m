% Indexing Many Circles Test Video 
% For ML analysis
% 11/18/2020

%%%%%%%%% Plot circles
% Output should be a .avi file that contains circles moving over a preset
% amount of frames

%%% overlapping
% a = -10; b = 10; %interval (a,b)
% r = a + (b-a).*rand(100,1); %radii
% 
% c = -100; d = 100;
% cx = c + (d-c).*rand(100,1); %center(s)
% cy = c + (d-c).*rand(100,1); %center(s)
% c = [cx cy];
% plot(cx,cy,'*')
% figure
%%%

clear all
close all
% blahblahbalh
% non-overlapping
numCirclesMax = 25;
[centers,radii] = tonsNonOverlappingCirlces(numCirclesMax+1);
% viscircles(centers, radii,'Color','k');
% axis equal

figure
for i = 1:length(radii)
    set(gca,'Color','k')
    plot_circles(radii(i),centers(i,:),'w')
    hold on
    axis ([0 70 0 20])
end

%%%%%%%%%%

% Move circles
frames = 20; 

% Write video
video = VideoWriter('barely_move_50_try2.avi');
video.FrameRate = 1;
open(video)

aggreg_remaining = numCirclesMax;
figure
for j=1:frames
    clf
    for i = 1:length(radii)
        set(gca,'Color','k')
        plot_circles(radii(i),centers(i,:),'w')
        hold on
        axis ([0 70 0 20]) %video frame size
    end
    hold on
    %plot(centers(50,1), centers(50,2),'r*')
    
    %%%%% RULE HERE:
    %   horizontal movement -- change centers(:,1)
    %   vertical movement -- change centers(:,2)
    %   random movement -- add rand(1) to both centers  
    centers(:,1) = centers(:,1) + 1; 
    
    %%%%%
    remove_row = randi(aggreg_remaining);
    centers(remove_row,:) = [];
    radii(remove_row) = [];
    aggreg_remaining=aggreg_remaining-1;
    
    % Capture current figure and convert to intensity
    %set(gca,'XColor', 'none','YColor','none')
    H = getframe(gca);
    %imwrite(H.cdata, 'myfilename.png')
    writeVideo(video,H);

    f = H.cdata;
    g=f(:,:,2); % 0 is black, 255 is white
    %remove first and last column
    g(:,1) = []; g(:,end) = [];
    g1=double(g); 
    h1=mat2gray(g1); %convert matrix to intensity image
    h2 = imbinarize(h1);
    [L1,N]=bwlabel(h1);
    disp(j)

end
close(video)

%% Extra
% only necessary if you want circles to blink, converge or diverge
blink_flag = 0;
conv_flag = 0;
div_flag = 1;


% fh = figure('Menu','none','ToolBar','none'); 
% ah = axes('Units','Normalize','Position',[0 0 1 1]);
for i=1:frames
    clf
    % Plot circles 
    %%%%% Blink circles
    if blink_flag == 1
        if i ~= frame_change
            set(gca,'Color','k')
            plot_circles(r,c(1,:),'w')
            plot_circles(r,c(2,:),'w')
            plot_circles(r,c(3,:),'w')
            %viscircles([0 0], 5, 'Color','g')
            hold on
        else %only plot 1 circle
            hold on
            set(gca,'Color','k')
            plot_circles(r,c(1,:),'w')
            plot_circles(r,c(3,:),'w')
        end
        axis ([-100 100 -100 100]) 
        axis square
    end
    
    %%%%% Converge circles
    if conv_flag == 1
        if i < frame_change
            set(gca,'Color','k')
            plot_circles(r,c(1,:),'w')
            plot_circles(r,c(2,:),'w')
            plot_circles(r,c(3,:),'w')
        elseif i >= frame_change
            set(gca,'Color','k')
            plot_circles(r+10,c(2,:),'w')
        end
        axis ([-100 100 -100 100])
        axis square       
    end
    
    %%%%% Diverge into multiple circles
    if div_flag == 1
        if i < frame_change
            set(gca,'Color','k')
            plot_circles(r,c(1,:),'w')
            plot_circles(r,c(2,:),'w')
            plot_circles(r,c(3,:),'w')
        elseif i >= frame_change %plot big circle
            set(gca,'Color','k')
            plot_circles(r,[c(1,1) c(1,2)-20],'w')
            plot_circles(r,c(2,:),'w')
            plot_circles(r,[c(3,1) c(3,2)+20],'w')
        end
        axis ([-100 100 -100 100])
        axis square       
    end
    
    
    
    %%%%%%%%
    % Capture current figure and convert to intensity
    %set(gca,'XColor', 'none','YColor','none')
    H = getframe(gca);
    %imwrite(H.cdata, 'myfilename.png')
    writeVideo(video,H);

    f = H.cdata;
    g=f(:,:,2); % 0 is black, 255 is white
    %remove first and last column
    g(:,1) = []; g(:,end) = [];
    g1=double(g); 
    h1=mat2gray(g1); %convert matrix to intensity image
    
    [L1,N]=bwlabel(h1);
    
%     % Visualize the sort order of the circles
%     %Found on https://blogs.mathworks.com/steve/2008/03/25/bwlabel-search-order/
%     s = regionprops(L1,'BoundingBox', 'Extrema', 'Centroid');
%     
%     %Store the x- and y-coordinates of the centroids into a two-column matrix.
%     centroids = cat(1,s.Centroid);
%     l=length(centroids(:,1));
%     centx(1:l,i)=centroids(:,1);
%     centy(1:l,i)=centroids(:,2);
%     
%     boxes = cat(1, s.BoundingBox);
%     left_edge = boxes(:,1);
%     [sorted, sort_order] = sort(left_edge);
%     s2 = s(sort_order);
%     imshow(f, 'InitialMag', 'fit')
%     hold on
%     for k = 1:numel(s2)
%        centroid = s2(k).Centroid;
%        text(centroid(1), centroid(2), sprintf('%d', k));
%     end
%     hold off
    
    %Color circles
    %rgb = label2rgb(L1);
    %imshow(rgb)
    %find location of circles
    %[r, c] = find(L==n);
    
    disp(i)
    
    % Repeat movement
    pause(.5)
    c(:,1) = c(:,1)+20; % move over along x axis
end
close(video)


%%%%%%% Helper functions
function plot_circles(r,c,colors)
pos = [c-r 2*r 2*r];
rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', colors, 'Edgecolor','none')
end

function [centers,radii] = tonsNonOverlappingCirlces(numCirclesMax)
numCircles = 1;
radius = 0.7;

maxIterations = 100 * numCirclesMax; % Fail Safe.
iteration = 1;
while  numCircles < numCirclesMax && iteration < maxIterations
  xTrial = 40*rand;
  yTrial = 20*rand;
  iteration = iteration + 1; % Fail safe.
  % Only need to check for overlap for second and later circles.
  if numCircles > 1
    % Find distance from other, prior circles.
    distances = sqrt((xTrial - x) .^ 2 + (yTrial - y) .^ 2);
    if min(distances) < 2 * radius + 2
      % It's overlapping at least one of the prior ones
      continue; % Skip to end of loop and continue with loop.
    end
  end
  x(numCircles) = xTrial;
  y(numCircles) = yTrial;
  numCircles = numCircles + 1;
end

radii = radius * ones(1, length(x));
centers = [x', y'];
end
