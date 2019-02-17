function [all_x,all_y] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots)
% given the video image "video", find the xy coordinates in each frame
% where the bucket is located (i.e. tracks the bucket over time)

% var_scale: controls the variance filter; looks for points with high
% variance and have variance > var_scale * mean_variance

% max_pixel_val: filters the pixels by color, looking for pixels with
%                grayscale color above this value.

% xrange & yrange: describe where to look for the bucket it in pixelspace
%                  chosen manually after viewing unfiltered points

% plots: logical vector for which describing which figures to show


numFrames = size(video, 4);
all_x = [];
all_y = [];

for k = 1 : numFrames
    mov(k).cdata = video(:,:,:,k); 
    mov(k).colormap = [];
end

all_Xg = zeros(480,640, numFrames);

% convert video to grayscale
for j=1:numFrames 
    X=frame2im(mov(j));
    Xg = rgb2gray(X);
    all_Xg(:,:, j) = Xg;
end

% calculate variances of each pixel over time and the mean variance
% among all pixels
all_variances = var(all_Xg, [], 3);
mean_var = mean(all_variances, 'all');


for j=1:numFrames
     
    D = uint8(all_Xg(:,:,j));
    
    % filter pixels by color (in this frame) and by their variance (over
    % course of the entire film)
    % points that meet the given criteria
    
    points = logical((all_variances > var_scale*mean_var) .* (D >= max_pixel_val));
    D(points) = 0;
    D(~points) = 255;
    
    [y_vals, x_vals] = find(points == 1);
    
    % indices of points within yrange
    filtered_y = (y_vals >= yrange(1)) .* (y_vals <= yrange(2));
    
    % indices of points within xrange
    filtered_x = (x_vals >= xrange(1)) .* (x_vals <= xrange(2));
    
    % only take points that are both within yrange AND xrange
    filtered_points = [x_vals(logical(filtered_y.*filtered_x)),  y_vals(logical(filtered_x.*filtered_y))];
    ave_x = mean(filtered_points(:, 1));
    ave_y = mean(filtered_points(:, 2));
    
    all_x = [all_x, ave_x];
    all_y = [all_y, ave_y];
    

    % X points of interest vs time
    
    if plots(1) == 1
        figure(1)
        subplot(211)
        title("X versus time")
        plot(j*ones(1, length(x_vals)), x_vals ,'r.');
        xlim([0, numFrames])
        ylim([0, size(points, 2)])
        hold on

        % Y points of interest vs time
        subplot(212)
        title("Y versus time")
        plot(j*ones(1, length(y_vals)), y_vals ,'r.');
        xlim([0, numFrames])
        ylim([0, size(points, 1)])
        hold on
    end
    
    if plots(2) == 1
        figure(2)
        imshow(D); 
    end
    
    if plots(3) == 1
        % All points that meet target conditions as a function of time
        figure(3)
        plot3(x_vals, y_vals, j*ones(1,length(y_vals)), 'r.'), grid on;
        xlabel("X")
        ylabel("Y")
        title("Unfiltered points")
        zlabel("Time (frame number)")
        xlim([0, 640]);
        ylim([0, 480]);
        zlim([0, numFrames]);
        hold on;
    end
    
    if plots(4) == 1
        % filtered paint can trajectory in 3D
        figure(4)
        plot3(filtered_points(:,1), filtered_points(:,2), j*ones(1,length(filtered_points)), 'r.'), grid on;
        hold on
        % also plot "averaged" point
        plot3(ave_x, ave_y, j, 'ko')
        xlabel("X")
        ylabel("Y")
        title("filtered points")
        zlabel("Time (frame number)")
        xlim([0, 640]);
        ylim([0, 480]);
        zlim([0, numFrames]);
        hold on;
    end
    
    if plots(5) == 1
        %averaged trajectory in 3D
        figure(5)
        plot3(ave_x, ave_y, j, 'ko'), grid on;
        xlabel("X")
        ylabel("Y")
        title("averaged trajectory")
        zlabel("Time (frame number)")
        xlim([0, 640]);
        ylim([0, 480]);
        zlim([0, numFrames]);
        hold on;
    end
    
    if plots(6) == 1
        figure(6)
        all_Xg(round(ave_y)-5:round(ave_y)+5, round(ave_x)-5:round(ave_x)+5, j) = 0;
        imshow(uint8(all_Xg(:,:, j)))
    end
    
    drawnow
end


end

