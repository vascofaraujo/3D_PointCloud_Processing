function [t2w, xyz, rgb] = rigid_transforms(imgseq,wr,cal,maxnpts)
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

        [xyz1, xyz2, inliers, ~] = remouts(dep1, dep2, rgb1, rgb2, cal, 160);

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
end

