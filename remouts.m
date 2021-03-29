function [xyz_array1das, xyz_array2das, inliers, inliers_px] = remouts(dep1,dep2,rgb1,rgb2,cal,k)
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

% figure
% showMatchedFeatures(I1,I2,inliers_px(:,1:2),inliers_px(:,3:4),'montage')

end

