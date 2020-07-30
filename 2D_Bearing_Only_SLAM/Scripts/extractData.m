function [X_r, X_l, Z, associations, num_poses, num_landmarks, landmark_ids] = extractData(landmark_positions, robot_positions, observations);
  
  num_poses = size(robot_positions,2)-1;
  num_landmarks = size(landmark_positions,2);
  
  X_r = zeros(3, 3, num_poses);
  X_l = zeros(2, num_landmarks);
  
  for(i=2:size(robot_positions,2))
    s = sin(robot_positions(i).theta);
    c = cos(robot_positions(i).theta);
    Rp = [c,-s;
          s,c];
    X_r(1:2,1:2,i-1) = Rp;
    X_r(1:2,3,i-1) = [robot_positions(i).x; robot_positions(i).y];
    X_r(3,3, i-1) = 1;
  end

  landmark_ids = zeros(1, num_landmarks);
  for(i=1:num_landmarks)
    landmark_ids(1,i) = landmark_positions(i).id;
    X_l(:,i) = [landmark_positions(i).x_pose; landmark_positions(i).y_pose];
  end

  measurement_num = 1;
  associations = [];
  Z = [];
  for(i=1:size(observations,2))
    observed_landmarks = observations(i).observation;
    for(j=1:size(observed_landmarks,2))
      obs = observed_landmarks(j);
      for (k=1:size(landmark_ids,2))
        if obs.id == landmark_ids(1,k)
          associations(:, measurement_num) = [i; obs.id];
          Z(1,measurement_num) = normalizeAngle(obs.bearing);
          measurement_num++;
        end
      end
    end
  end
  

  
end