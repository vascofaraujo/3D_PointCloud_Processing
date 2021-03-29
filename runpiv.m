%% -- CRIAR LISTA IMAGENS
%DATA DIRECTORY - PATH (with / at the end)
base_data_dir=[pwd '/newpiv/'];
d1=dir([base_data_dir 'depth*']);
r1=dir([base_data_dir 'rgb*']);
if exist('im1'),
    clear im1;
end
for i=1:length(d1),
    im1(i).rgb=[base_data_dir r1(i).name];
    im1(i).depth=[base_data_dir d1(i).name];
end
%load calibration data
load cameraparametersAsus;
maxnpts=500000;

%%
imgseq = im1
wr = 1
cal = cam_params

%
nr_imgs = length(imgseq);
merge = 0;
for i=1:nr_imgs
    files(i).rgbs = imread(imgseq(i).rgb);
    files(i).depths = imread(imgseq(i).depth);
end

ransac_its = round(1.5 * log10(1-0.999)/log10(1-0.5^4));

FromCam2W(1).R = eye(3,3);
FromCam2W(1).T = zeros(3,1);

check_vec = zeros(1,nr_imgs);
check_vec(wr) = 1;

i = 0;
tok = 0;
while(~all(check_vec))
    a = (wr+1);
    b = 1;
    c = nr_imgs; 
    if i == nr_imgs
        a = (wr-1);
        b = -1;
        c = 1;
        tok = 1;
    end
    h = waitbar(0,'Please wait...');
    for i = a:b:c
        if tok
            dep1 = files(i+1).depths;
            dep2 = files(i).depths;
            rgb1 = files(i+1).rgbs;
            rgb2 = files(i).rgbs;
        else
            dep1 = files(i-1).depths;
            dep2 = files(i).depths;
            rgb1 = files(i-1).rgbs;
            rgb2 = files(i).rgbs;
        end

        %Load Colors
        color1 = reshape(rgb1, [307200, 3]);
        color1 = double(color1)./double(255);
        color2 = reshape(rgb2, [307200, 3]);
        color2 = double(color2)./double(255);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        k = 160;

        % Finds RT and Inliers

       cam_R = cal.R;
cam_T = cal.T;
K_rgb = cal.Krgb;
K_dep = cal.Kdepth;

xyz_array1das = get_xyzasus(dep1(:), [size(rgb1,1) size(rgb1,2)], 1:size(rgb1,1)*size(rgb1,2), K_dep, 1, 0);
rgbd1 = get_rgbd(xyz_array1das, rgb1, cam_R, cam_T, K_rgb);

xyz_array2das = get_xyzasus(dep2(:), [size(rgb2,1) size(rgb2,2)], 1:size(rgb2,1)*size(rgb2,2), K_dep, 1, 0);
rgbd2 = get_rgbd(xyz_array2das, rgb2, cam_R, cam_T, K_rgb);

% % Removes 0 and depth > 4m depths and converts to XYZ
% idx = [];
% for n = 1:size(xyz_array1das,1)
%     A = xyz_array1das(n,:);
%     B = xyz_array2das(n,:);
%     if ~(A(3) && B(3) && A(3) < 4 && A(3) > -4 && B(3) < 4 && B(3) > -4)% is not 0 and below 4 meters
%         idx = [idx n];
%     end
% end
% xyz_array1das(idx,:) = [];
% xyz_array2das(idx,:) = [];

%Reshape matrices
xyz_array1 = reshape(xyz_array1das,size(rgb1));
xyz_array2 = reshape(xyz_array2das,size(rgb2));

%Find Matches between figures
I1 = single(rgb2gray(rgbd1));
I2 = single(rgb2gray(rgbd2));
[f1,d1] = vl_sift(I1);
[f2,d2] = vl_sift(I2);
[matches, ~] = vl_ubcmatch(d1, d2) ;

%Round coordinates
match1_cord = round(f1(1:2, matches(1,:)));
match2_cord = round(f2(1:2, matches(2,:)));

%Remove black points
idx = [];
for n = 1:size(match1_cord,2)
    A = ~any(rgbd1(match1_cord(2,n),match1_cord(1,n),:));
    B = ~any(rgbd2(match2_cord(2,n),match2_cord(1,n),:));
   if (A || B)
      idx = [idx n];
   end
end
match1_cord(:, idx) = [];
match2_cord(:, idx) = [];

% Removes 0 and depth > 4m depths and converts to XYZ
for n = 1:size(match1_cord,2)
    A = reshape((xyz_array1(match1_cord(2,n),match1_cord(1,n),:)), [1 3]);
    B = reshape((xyz_array2(match2_cord(2,n),match2_cord(1,n),:)), [1 3]);
    if A(3) && B(3) && A(3) < 4 && A(3) > -4 && B(3) < 4 && B(3) > -4% is not 0 and below 4 meters
        match1_cord3(n,:) = A;
        match2_cord3(n,:) = B;
    end
end

% for n = 1:size(match1_cord,2)
%     match1_cord3(n,:) = xyz_array1(match1_cord(2,n),match1_cord(1,n),:);
%     match2_cord3(n,:) = xyz_array2(match2_cord(2,n),match2_cord(1,n),:);
% end

% figure(60);
% showMatchedFeatures(I1,I2,round(match1_cord'),round(match2_cord'),'montage')

best = 0;
inliers = [0 0 0 0 0 0];
inliers_px = [0 0 0 0];
sumofins = 0;

% Ransac
for n = 1:k
    % Select 4 random pairs
    for i = 1:4
        x = randi(size(match1_cord,2));
        r_pairs_a(i,:) = [match1_cord(1,x), match1_cord(2,x)];
        r_pairs_b(i,:) = [match2_cord(1,x), match2_cord(2,x)];
    end
    
    %showMatchedFeatures(I1,I2,r_pairs_a(:,1:2),r_pairs_b(:,1:2),'montage')
    
    % Convert 4 pairs to XYZ 
    for n = 1:4
        r_pa3(n,:) = xyz_array1(r_pairs_a(n,2),r_pairs_a(n,1),:);
        r_pb3(n,:) = xyz_array2(r_pairs_b(n,2),r_pairs_b(n,1),:);
    end
    
    %Find R and T
    [~, ~, transform] = procrustes(double(r_pa3), double(r_pb3));
    R = transform.T;
    T = transform.c;

    % Apply R and T to all matches from 2nd image
    change_matches_2 = (double(match2_cord3)*R);
    for i = 1:size(change_matches_2,1)
       change_matches_2(i,:) = change_matches_2(i,:) + T(1,:);
    end
%     change_matches_2(:,1) = change_matches_2(:,1) + T(1,1);
%     change_matches_2(:,2) = change_matches_2(:,2) + T(1,2);
%     change_matches_2(:,3) = change_matches_2(:,3) + T(1,3);

    Pi = change_matches_2;
    Pi_p = match1_cord3;
    
    euclidean = vecnorm(Pi_p - Pi,2,2);
    checks = euclidean < 0.02;
    num_right = sum(checks);
    sumofins = sumofins + num_right;
    if num_right > best
        best = num_right;
        chosen_R = R;
        chosen_T = T(1,:);
        %save inliers
        for m = 1:length(checks)
           if checks(m) == 1 %Is inlier
               inliers = [inliers; match1_cord3(m,:), match2_cord3(m,:)];
               inliers_px = [inliers_px; match1_cord(:,m)', match2_cord(:,m)'];
           end
        end
    end
    
end

inliers_px = inliers_px(2:size(inliers_px,1),:);
inliers = inliers(2:size(inliers,1),:);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~, ~, transform] = procrustes(inliers(:,1:3), inliers(:,4:6), 'reflection', false);
        R = transform.T;
        T = transform.c(1,:); 

        RTs.R = R;
        RTs.T = T;
        if tok
            RTgraph(i, i+1) = RTs;    
        else
            RTgraph(i, i-1) = RTs;
        end

        %Transform 2nd Image
        d = i;
        e = -1;
        f = wr+1;
        if tok
           d = i;
           e = 1;
           f = wr-1;
        end 
        t_xyz2 = xyz2;
        for m = d:e:f
            row = m;
            col = m-1;
            if tok
                row = m;
                col = m+1;
            end
            t_xyz2 = t_xyz2*RTgraph(row,col).R;
            for ct = 1:size(t_xyz2,1)
                t_xyz2(ct,:) = t_xyz2(ct,:) + RTgraph(row,col).T;
            end
        end
        
        check_vec(i) = 1;
        disp("-----------------------------------")
        disp(check_vec)

        % Find R_iw and T_iw
        [~, ~, FC2W] = procrustes(t_xyz2, xyz2, 'reflection', false);
        FromCam2W(i).R = FC2W.T;
        FromCam2W(i).T = FC2W.c(1,:)';
        
        %Save Reference
        if i-1 == wr
            refpc = pointCloud(xyz1);
            refpc.Color = uint8(color1*255);
        end

        % Save Point Cloud
        box(i).pcs = pointCloud(t_xyz2);
        box(i).pcs.Color = uint8(color2*255);

        waitbar(i/nr_imgs, h, sprintf('%d/%d',i, nr_imgs))
    end
end
close(h)

% Merge Point Clouds
merge = refpc; % World Frame
for m = (wr+1):size(box, 2) % Merge up wf
    merge = pcmerge(merge, box(m).pcs, 0.005);
end
for m = (wr-1):-1:1 % Merge down wf
    merge = pcmerge(merge, box(m).pcs, 0.005);
end

%Downsample
percentage = maxnpts/length(merge.Location);
merge = pcdownsample(merge, 'random', percentage);

%OUTPUTS
t2w = FromCam2W;
xyz = merge.Location;
rgb = merge.Color;




%%
if numel(xyz)>3*maxnpts || numel(rgb)>3*maxnpts,
    error(' too many points returned');
else,
    pc=pointCloud(xyz, 'Color',rgb);
end
figure(1);
pcshow(pc);
%figure(2);
%replace pointsbig for the point cloud of each image and check that the
%transformations are ok
% load auxfile;
% for i=2:length(t2w),
%     xyzc=pointsbig{i,1}.Location;
%     xyzw=[t2w(i).R t2w(i).T]*[xyzc';ones(1,length(xyzc))];
%     pcshow(pointCloud(xyzw','Color',pointsbig{i,2}));
%     pause;
% end


