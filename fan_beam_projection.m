% ***This program is written by Kun Wang, whose student ID is P19104048***
% This program is to simulate the sinogram for shepp-logan using fan-beam
% CT. Fanbeam function in matlab, siddon algorithms, ... are used to
% implement forward projection of the phantom.

clear
clc
load Shepp_logan.mat
figure, imshow(Shepp_logan, [0.98 1.05])

%%parameters for reconstruction
% dso : the distance from the x-ray source to the object in mm
% dsd : the distance from the x-ray source to the detector in mm
% angle : the number of scanning angle
% pixel : the detector pixel size in mm
% detector_num : the number of detector bin
% voxel : the voxel size in mm
% N : image resolution in height
% Nx: the number of planes along x direction
% Ny: the number of planes along y direction
% fov : the diameter of the field of view

dso = 500;
dsd = 800;
theta = 0:1:359;
pixel = 1.25;
detector_num = 320;
voxel = pixel * dso / dsd;
N = 256;
Nx = 257;
Ny = 257;
fov = 200;

tstart = tic;

% tmp = zeros(N);
% for ii=1:N
%     for jj=1:N
%         tmp(ii,jj) = Shepp_logan(N+1-ii,jj);
%     end
% end

%%Forward projection using fanbeam function
% [F,Fpos,Fangles] = fanbeam(tmp', 500/voxel, 'FanRotationIncrement', 1, 'FanSensorGeometry', 'line');
% tellap = toc(tstart);
% disp(['time is ', num2str(tellap)]);
% 
% figure, imshow(F',[])
% xlabel('Detector Position')
% ylabel('Rotation Angle (degree)')

%%Forward projection using siddon algorithm
% p1_source : the original position of the ray source
p1_source = [-dso 0];
det_pos = (-(detector_num-1)/2:(detector_num-1)/2) * pixel;
p2_det = zeros(detector_num, 2);
for ii=1:length(det_pos)
    p2_det(ii,1) = dsd - dso;
    p2_det(ii,2) = det_pos(ii);
end
Xplanes = (-N/2:1:N/2) * voxel;
Yplanes = (-N/2:1:N/2) * voxel;

% figure;
% hold on;
% X0 = -Nx/2*voxel;
% Y0 = -Ny/2*voxel;
% 
% X_start = X0;
% X_end = X0 + N * voxel;
% Y_start = Y0;
% Y_end = Y0 + N * voxel;
% 
% Xb = X0 + N*voxel;
% Yb = Y0 + N*voxel;
% line([X0 X0], [Y_start Y_end]);
% line([Xb Xb], [Y_start Y_end]);
% line([X_start X_end], [Yb Yb]);
% line([X_start X_end], [Y0 Y0]);

% plot([p1_x p2_x], [p1_y p2_y], 'r');
% plot([p1_x p2_det(320, 1)], [p1_y p2_det(320, 2)], 'r');

sino = zeros(length(theta), detector_num);
beta = theta * pi / 180;

for ii=1:length(beta)
    
    p1_x = p1_source(1) * cos(beta(ii)) - p1_source(2) * sin(beta(ii));
    p1_y = p1_source(1) * sin(beta(ii)) + p1_source(2) * cos(beta(ii));
    
    for jj=1:detector_num
        % the two points rotate
        p2_x = p2_det(jj, 1) * cos(beta(ii)) - p2_det(jj, 2) * sin(beta(ii));
        p2_y = p2_det(jj, 1) * sin(beta(ii)) + p2_det(jj, 2) * cos(beta(ii));
        
%         if ii<=10
%             if jj == 1 || jj == detector_num
%                 plot([p1_x p2_x], [p1_y p2_y], 'r-');
%             end
%         end
        
        if p1_x~=p2_x
            alpha_x1 = (Xplanes(1) - p1_x) / (p2_x - p1_x);
            alpha_Nx = (Xplanes(Nx) - p1_x) / (p2_x - p1_x);
        end
        % p1_x == p2_x
        if p1_y~=p2_y
            alpha_y1 = (Yplanes(1) - p1_y) / (p2_y - p1_y);
            alpha_Ny = (Yplanes(Ny) - p1_y) / (p2_y - p1_y);
        end
        %p1_y == p2_y
        
        alpha_min = max([0, min(alpha_x1, alpha_Nx), min(alpha_y1, alpha_Ny)]);
        alpha_max = min([1, max(alpha_x1, alpha_Nx), max(alpha_y1, alpha_Ny)]);
        
        if alpha_min < alpha_max
            if p2_x >= p1_x
                i_min = Nx - floor((Xplanes(Nx) - alpha_min*(p2_x - p1_x) - p1_x) / voxel);
                i_max = 1 + floor((p1_x + alpha_max*(p2_x - p1_x) - Xplanes(1)) / voxel);
            else
                i_min = Nx - floor((Xplanes(Nx) - alpha_max*(p2_x - p1_x) - p1_x) / voxel);
                i_max = 1 + floor((p1_x + alpha_min*(p2_x - p1_x) - Xplanes(1)) / voxel);
            end

            if p2_y >= p1_y
                j_min = Ny - floor((Yplanes(Ny) - alpha_min*(p2_y - p1_y) - p1_y) / voxel);
                j_max = 1 + floor((p1_y + alpha_max*(p2_y - p1_y) - Yplanes(1)) / voxel);
            else
                j_min = Ny - floor((Yplanes(Ny) - alpha_max*(p2_y - p1_y) - p1_y) / voxel);
                j_max = 1 + floor((p1_y + alpha_min*(p2_y - p1_y) - Yplanes(1)) / voxel);
            end
        
            if p2_x > p1_x
                alpha_x = (Xplanes(i_min:1:i_max) - p1_x) / (p2_x - p1_x);
            elseif p2_x < p1_x
                alpha_x = (Xplanes(i_max:-1:i_min) - p1_x) / (p2_x - p1_x);
            end
        
            if p2_y > p1_y
                alpha_y = (Yplanes(j_min:1:j_max) - p1_y) / (p2_y - p1_y);
            elseif p2_y < p1_y
                alpha_y = (Yplanes(j_max:-1:j_min) - p1_y) / (p2_y - p1_y);
            end
        
            alpha = [alpha_min, alpha_x, alpha_y, alpha_max];
            alpha = unique(alpha);
            alpha = sort(alpha);
        
            d12 = sqrt((p1_x - p2_x)^2 + (p1_y - p2_y)^2);
            len = d12 * (diff(alpha));
        
            i_index = zeros(length(len), 1);
            j_index = zeros(length(len), 1);
            for kk=2:length(alpha)
                alpha_mid = (alpha(kk) + alpha(kk-1)) / 2;
                
                i_index(kk-1) = 1 + floor((p1_x + alpha_mid * (p2_x - p1_x) - Xplanes(1)) / voxel);
                i_index(kk-1) = min( i_index(kk-1), 256 );
                i_index(kk-1) = max( i_index(kk-1), 1 );
                    
                j_index(kk-1) = 1 + floor((p1_y + alpha_mid * (p2_y - p1_y) - Yplanes(1)) / voxel);
                j_index(kk-1) = min( j_index(kk-1), 256 );
                j_index(kk-1) = max( j_index(kk-1), 1 ); 
            end

            for mm=1:(length(len))
                sino(ii, jj) = sino(ii, jj) + len(mm) * Shepp_logan(N+1-j_index(mm), i_index(mm));
            end        
        end
    end
    disp(['processing...', num2str(ii)]);
end

figure, imshow(sino, [])
xlabel('Detector Position')
ylabel('Rotation Angle (degree)')
save('mSiddon.mat', 'sino')