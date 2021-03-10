% Indexing_Kirby_Suggestions
% Most recent update: 12/11/20
% Andrew Edit 2/25/21

% GOAL: Track the location of aggregates 
%   INPUT: .avi video file with moving aggregates
%   OUTPUT: indexed centroid and area matrix 
%           (r,c) = (aggregate #, frame #)

% Andrew EDIT 2/25/21: addressing removed aggregates per frame: 
%%%%%%%%%

% PART 1:
%   -Input video file with aggregates
%   -Extract centroids and areas into centx, centy, area

clear all
clf
shg
% Import video 
cd 'C:\Users\Andrew''s Laptop\Documents\KirbyLab\MATLAB\Indexing Videos\'
%cd 'C:\Users\Andrew Kang\Documents\KirbyLab\MATLAB\Indexing\'
v = VideoReader('barely_move_50_try2.avi');
frames = v.NumFrame;
%frames = 10;
centx = []; centy = []; area = [];


video = VideoWriter('Remove_per_frame_color.avi');
video.FrameRate = 2;

%imtool
%imtool(new_L1) % Use this tool to zoom in on ellipses and identify aggregate number
% index 39 and 42 are so close that they become joint. 
%%
open(video)
num = [];
L1matrix = [];
for i = 1:10
    f=read(v,i); %read data in the i frame
    g=sum(double(f),3); %pull out all rows and columns in 2nd dimension
    %remove first and last column
    g(:,1) = []; g(:,end) = [];
    g(1,:) = []; g(end,:) = [];
    
%     rowsadded = zeros(20, size(g,2)); 
%     g = cat(1, rowsadded, g); g = cat(1 ,g, rowsadded);
%     colsadded = zeros(size(g,1), 20);
%     g = cat(2, colsadded, g); g = cat(2, g, colsadded);
    %imshow(g)    
    %axis on
    g1=double(g);
    h1=mat2gray(g1); %convert matrix to intensity image
    h1 = imbinarize(h1, .2); %binarize, remove tickmarks in image
    
    [L1,N] =bwlabel(h1); %labels each aggregate with a unique number in the matrix
                        %N is the total number of objects
    num = [num; N]
    L1matrix(:,:,i) = L1;
    mycolormap = jet(100);
%     for p = 1:size(L1,1) 
%         for q = 1:size(L1,2)
%             if mod(L1(p,q), 10) == 1
%                 L1(p,q) = 30;
%             elseif mod(L1(p,q), 10) == 2
%                 L1(p,q) = 37;
%             elseif mod(L1(p,q), 10) == 3
%                 L1(p,q) = 52;
%             elseif mod(L1(p,q), 10) == 4
%                 L1(p,q) = 94;
%             elseif mod(L1(p,q), 10) == 5
%                 L1(p,q) = 63;
%         	elseif mod(L1(p,q), 10) == 6
%                 L1(p,q) = 26;
%             elseif mod(L1(p,q), 10) == 7
%                 L1(p,q) = 81;
%             elseif mod(L1(p,q), 10) == 8
%                 L1(p,q) = 76;
%             elseif mod(L1(p,q), 10) == 9
%                 L1(p,q) = 44;
%             elseif mod(L1(p,q), 10) == 0 && L1(p,q) ~= 0
%                 L1(p,q) = 88;
%             end
%             
%         end
%     end
    %L1 = mat2gray(L1);
    imshow(L1)
    colormap(jet(100))
    title(sprintf('Frame = %0.0f%', i), 'FontSize', 16);
    axis on
   
    [new_L1,new_N]=bwlabel(h1); % 

    s1 = regionprops(new_L1,'centroid');
    a1 = regionprops(new_L1,'Area');
    bb1 = regionprops(new_L1, 'BoundingBox');
    centroids = cat(1,s1.Centroid);

   
%     if isempty(centroids) %account for blinking aggs
%         centx(1,i)=0;
%         centy(1,i)=0;
%         area(1,i)=0;
%     else
        l=length(centroids(:,1));
        centx(1:l,i)=centroids(:,1);
        centy(1:l,i)=centroids(:,2);
        area(1:l,i)=cat(1,a1.Area);
        
%     end                                            
    
%     hold on
%     plot(centx(1,i), centy(1,i),'r*')
%     axis on
%     set(gca,'YDir','normal')
%     hold off
    
    pause(.3)
    disp(i)

      H = getframe(gcf);
      writeVideo(video,H);

end
close(video)

writematrix(centx,'Centx.csv')
writematrix(centy,'Centy.csv')
writematrix(area,'Area.csv')
%%
%%%%%%%%% Andrew Code
% Index circles 
[m,n]=size(centx); %size of centroid matrix

%new matrices with reassigned indices
centx_new=zeros(m,n); 
centy_new=zeros(m,n); 
area_new=zeros(m,n);  

%the coordinates from the first frame remain the same
centx_new(:,1)=centx(:,1); 
centy_new(:,1)=centy(:,1); 
area_new(:,1)=area(:,1); 
totalagg = m;

d = zeros(m, m, 1);
overlapmatrix = [];
%%
L11 = L1matrix(:,:,5);
L12 = L1matrix(:,:,6);
imtool(L11)
imtool(L12)
%%
idxlost = []; 
greatercount = 1;
%
for i = 2:10
    % Code for discerning aggregate to remove
    % 
    for agg = 1:num(i-1);
        [r0, c0] = find(L1matrix(:,:,i-1) == agg) % image at frame i - 1
        [r1, c1] = find(L1matrix(:,:,i) == agg) % image at frame i 
        overlappix = 0;
        for rowidx = 1:length(r0);
            for next_coord = 1:length(r1);
                if r0(rowidx) == r1(next_coord) && c0(rowidx) == c1(next_coord)
                    overlappix = overlappix + 1; 
                end
            end
        end
        overlapmatrix(agg, i-1) = overlappix
    end
    
    %
    idx = find(overlapmatrix(:, i-1) == 0); 
    idx = idx(1);
    if i-1 >= 2;
        if idx >= idxlost(end);
            idx = idx + greatercount;
            greatercount = greatercount + 1;
            idxlost = [idxlost;idx];
        else;
            %idx = idx + greatercount -1;
            idxlost = [idxlost;idx];
            
        end
    else
        idxlost = [idxlost;idx]
    end
    
    %
    for agg = 1:num;
    xprev = centx_new(agg,(i-1)); %centx_new(agg, (i-1)); %
    yprev = centy_new(agg,(i-1)); %centy_new(agg, (i-1)); % 
        for j = 1:size(centx_new, 1)
            xnext = centx(j, i); % 
            ynext = centy(j, i); % 
            d(agg, j, i) = sqrt((xprev-xnext).^(2)+(yprev-ynext).^(2));
        end
        d(idx, :, i) = 1000;
    % 
        
%         for e = 1:size(d,1)
%             for f = 1:size(d,2)
%                 if d(e,f,i) == 0
%                     d(e,f,i) = 100;
%                 end
%             end
%         end

        [val,ind] = min(d(agg,:,i));
        centx_new(agg,i) = centx(ind,i); 
        centy_new(agg,i) = centy(ind,i);
        for z = 1:size(idxlost,1);
            remove = idxlost(z);
            centx_new(remove,i) = -1;
            centy_new(remove,i) = -1;
        end
        area_new(agg,i)=area(ind,i);
        %centx(agg,i) = centx(ind,i); % update previous frame
        %centy(agg,i) = centy(ind,i); % update previous frame
    end

    disp(i)
% Check proper indexing
%xrearr = centx_new(31,:) >= centx_new(30,:)
%yrearr = centy_new(30,:) >= centy_new(31,:)
end
% Tracking indexing: 
%
idxmatrix = zeros(m,n);
idxmatrix(:,1) = linspace(1, m, m);

for frame = 2:size(d,3);
    for row = 1:size(d,1);
        [val,ind] = min(d(row,:,frame));
        idxmatrix(row,frame) = ind;
        idxmatrix(idxlost(frame-1), frame) = -1;
        if idxmatrix(row,frame-1) == -1
            idxmatrix(row,frame) = -1;
        end
        
    end
end
%writematrix(idxmatrix,'IndexTrack.csv')

%%

%centy
%centy_new

%writematrix(centx_new,'Centx_new.csv')
%writematrix(centy_new,'Centy_new.csv')
%writematrix(area_new,'Area.csv')

%% Recreation of new video with consistent color coding of aggregates 
% L1 is the old image with the matrices 
% We need to take this image, switch indices based on index matrix and then
% recreate image; call this new color coded image L2. 
% Finally, take L2 and replace entities as 0 whenever L2 is a positive
% integer. 

L2matrix = []; 

for i = 1:10 %frames
    numaggs = max(idxmatrix(:,i));
    beforeimage = L1matrix(:,:,i);
    [beforeimage, nbefore] = bwlabel(beforeimage);
    frameone_idx = 1; 
    L2matrix(:,:,i) = zeros(size(beforeimage,1), size(beforeimage,2));
    for j = 1:num
        framei_idx =idxmatrix(j,i);
        for row = 1: size(beforeimage,1)
            for col = 1:size(beforeimage,2)
                if beforeimage(row,col) == framei_idx
                    L2matrix(row, col, i) = frameone_idx;
                end
            end
        end
        frameone_idx = frameone_idx + 1;
    end
end
%% Corrected Indexing Video Generation
v2 = VideoWriter('remove_1_corrected.avi');
v2.FrameRate = 2;
open(v2)

for i = 1:10
    L2 = L2matrix(:,:,i);
    mycolormap = jet(100);
    for p = 1:size(L2,1) 
        for q = 1:size(L2,2)
            if mod(L2(p,q), 10) == 1
                L2(p,q) = 30;
            elseif mod(L2(p,q), 10) == 2
                L2(p,q) = 37;
            elseif mod(L2(p,q), 10) == 3
                L2(p,q) = 52;
            elseif mod(L2(p,q), 10) == 4
                L2(p,q) = 94;
            elseif mod(L2(p,q), 10) == 5
                L2(p,q) = 63;
        	elseif mod(L2(p,q), 10) == 6
                L2(p,q) = 26;
            elseif mod(L2(p,q), 10) == 7
                L2(p,q) = 81;
            elseif mod(L2(p,q), 10) == 8
                L2(p,q) = 76;
            elseif mod(L2(p,q), 10) == 9
                L2(p,q) = 44;
            elseif mod(L2(p,q), 10) == 0 && L2(p,q) ~= 0
                L2(p,q) = 88;
            end
            
        end
    end
    L2 = mat2gray(L2);
    imshow(L2)
    colormap(mycolormap)
    title(sprintf('Corrected Indexing, Frame = %0.0f%', i), 'FontSize', 16);
    axis on
    frame = getframe(gcf)
    writeVideo(v2,frame);
    pause(0.5)
end
close(v2);
%%
% PART 2
%   -Generate empty matrices to be filled based on aggregate index
%   -Apply closest distance algorithm to sort aggregates

% Index circles 
[m,n]=size(centx); %size of centroid matrix

%new matrices with reassigned indices
centx_new=zeros(m,n); 
centy_new=zeros(m,n); 
area_new=zeros(m,n);  

%the coordinates from the first frame remain the same
centx_new(:,1)=centx(:,1); 
centy_new(:,1)=centy(:,1); 
area_new(:,1)=area(:,1); 

% Nearest distance algorithm
row_counter = m;

%
d = zeros(m, m, 1);
cell_size = 40;
dist_away = 40*sqrt(2);

%%
for i=2:50%frames %n-1 %frame number
    for z = 1:m % 1:m %each aggregate in i frame
        x_current = centx(z,i);
        y_current = centy(z,i);
       for j=1:size(centy_new,1) % i-1 frame  
           %%% issue here with no aggregates in frame
           % Frame i-1 aggregate centroid
           x_prev = centx(j,i-1);
           y_prev = centy(j,i-1);
           if x_current <= (x_prev + cell_size) && x_current >= (x_prev - cell_size) && ... 
              y_current <= (y_prev + cell_size) && y_current >= (y_prev - cell_size) 
              %x_new_pos(z,1) ~= 0 && y_new_pos(z,1) ~= 0
              %calculate the distance between aggs
              d(z, j, i) = sqrt((x_prev-x_current).^(2)+(y_prev-y_current).^(2));
              sqrt((x_prev-x_current).^(2)+(y_prev-y_current).^(2))
           else
              d(z, j, i) = 1000; %might slow down code
           end
           [x_current, x_prev; y_current y_prev;z j; ...
               sqrt((x_prev-x_current).^(2)+(y_prev-y_current).^(2)) sqrt((x_prev-x_current).^(2)+(y_prev-y_current).^(2))]
       end
       

       % Choose smallest d and give it those agg properties
       % if distance is greater than some value, throw that distance away
       
       % if a distance is the smallest below a threshold d, say that
       % aggregate is the same as the previous aggregate
       
       %%%%%%%%%%%% HERE
       if min(d(z,:,i)) <= dist_away %&& centy_new(z,i) == 0
           %%% enters here bc min(d) = 96
           [val,ind] = min(d(z,:,i));
                      
           centx_new(ind,i)=centx(z,i); 
           centy_new(ind,i)=centy(z,i);
           area_new(ind,i)=area(z,i);
           
       elseif min(d(1,:)) == 1000 && centx(z,i) ~= 0 
           % assume that a new aggregate has appeared, give it a new row in
           % cent_new matrices
           row_counter = row_counter+1; % add a new row, check this works
           centx_new(row_counter, i) = centx(z,i); 
           centy_new(row_counter, i) = centy(z,i); 
           area_new(row_counter, i) = area(z,i); 
           disp('new row')
        end
    end 
    disp(i)
end
%%
% combine forward and backward looking algo
% send summary of problem

% For checking purposes:
centy
centy_new

writematrix(centx_new,'Centx_new.csv')
writematrix(centy_new,'Centy_new.csv')
writematrix(area_new,'Area.csv')

% plot(centx_new(1,:), centy_new(1,:),'r*','LineWidth',8)
% hold on
% plot(centx_new(2,:), centy_new(2,:),'g*','LineWidth',8)
% hold on
% plot(centx_new(3,:), centy_new(3,:),'b*','LineWidth',8)
% hold on
% plot(centx_new(4,:), centy_new(4,:),'m*','LineWidth',8)
% 
% %axis on
% axis([0 700 0 700])
% legend('Agg 1','Agg 2', 'Agg 3','Agg 4')

