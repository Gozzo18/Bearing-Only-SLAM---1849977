close all
clear
clc
warning('off', 'all');

#load dependencies
addpath "../tools/g2o_wrapper"
addpath "../tools/visualization"
addpath "./Scripts"
source "./Scripts/multi_ICP_3d.m"
source "./Scripts/triangulation.m"

global pose_dim=3;
global landmark_dim=2;

[_, robot_positions, robot_odometries, observations_guess] = loadG2o("../datasets/02-BearingOnlySLAM/slam2D_bearing_only_initial_guess.g2o");
retrieved_landmark_positions = triangulate(robot_positions, observations_guess);
retrieved_landmark_positions = filterPositions(retrieved_landmark_positions);
[Xr_guess, Xl_guess, Z_guess, associations_guess, num_poses, num_landmarks, landmark_ids_guess] = extractData(retrieved_landmark_positions, robot_positions, observations_guess);

#figure;
#plot(Xr_guess(1,3,:), Xr_guess(2,3,:), 'r');
#hold on;
#scatter(Xl_guess(1,:), Xl_guess(2,:), 'b');

num_iterations = 20;
damping = 0.01;
kernel_threshold = 1.0;
[XR, XL, chi_stats, _] = doMultiICP(Xr_guess, Xl_guess, landmark_ids_guess, Z_guess, 
							associations_guess,
              num_poses, 
							num_landmarks, 
						  num_iterations,
              damping,
              kernel_threshold);              
disp(chi_stats);

[landmark_positions_true, robot_positions_true, _, observations_true] = loadG2o("../datasets/02-BearingOnlySLAM/slam2D_bearing_only_ground_truth.g2o");
[Xr_true, Xl_true, _, _, _, _, _] = extractData(landmark_positions_true, robot_positions_true, observations_true);


figure;
plot(XR(1,3,:), XR(2,3,:), 'g');
hold on;
plot(Xr_true(1,3,:), Xr_true(2,3,:), 'k');
hold on;
scatter(XL(1,:), XL(2,:), 'b');
hold on;
scatter(Xl_true(1,:), Xl_true(2,:), 'r');
hold on;
lgd = legend("estimated robot map", "true robot map", "estimated landmarks","true landmarks", "location", "southeast");