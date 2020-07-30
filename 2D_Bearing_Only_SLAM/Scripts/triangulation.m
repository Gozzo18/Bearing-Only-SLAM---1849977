source "../tools/utilities/geometry_helpers_2d.m"

function filtered_positions = filterPositions(landmark_positions)
  
  num_of_landmark = 1;
  for(i=1:size(landmark_positions,2))
    if (landmark_positions(1,i) ~= 0 && landmark_positions(2,i) ~= 0 && landmark_positions(3,i) ~= 0)
      filtered_positions(num_of_landmark).id = landmark_positions(1,i)-1;
      filtered_positions(num_of_landmark).x_pose = landmark_positions(2,i);
      filtered_positions(num_of_landmark).y_pose = landmark_positions(3,i);
      num_of_landmark += 1;
    end
  end
end


function landmark_positions = triangulate(robot_positions, observations)
  
  landmark_positions = zeros(3,205);
  num_landmarks = 0;
  
  #Compute the number and retrieve the ids of the landmarks in the environment
  landmark_ids = zeros(1, num_landmarks);
  for (i=1:size(observations,2))
    robot_observation = observations(i).observation;
    for (j=1:size(robot_observation,2))
      observed_landmark = robot_observation(j).id + 1;
      if (~ismember(observed_landmark, landmark_ids))
        num_landmarks += 1;
        landmark_ids(1,num_landmarks) = observed_landmark;
      end
    end
  end

  #Remove the unused positions
  for (i=1:num_landmarks)
    if (landmark_ids(1,i) ~= 0)
      filtered_positions(1,i) = landmark_ids(1,i);
    end   
  end
  
  #For each landmark retrieve all the poses that observe it
  for (i=1:num_landmarks)
    landmark_id = landmark_ids(1,i);
    #Robot poses that observe the landmark
    robot_poses_id = zeros(1,15);
    #Robot bearing measurement
    global_bearings = zeros(1,15);
    #Index 
    counter_poses = 0;
    for (j=1:size(observations,2))
      robot_observation = observations(j).observation;
      for (k=1:size(robot_observation,2))
        observed_landmark = robot_observation(k).id + 1;
        if (landmark_id == observed_landmark)
          counter_poses += 1;
          robot_poses_id(1,counter_poses) = observations(j).pose_id;
          global_bearings(1, counter_poses) = robot_observation(k).bearing + robot_positions(j+1).theta ;
        end
      end
    end
    #With enough poses, try to initialize the landmark
    if (landmark_positions(1,landmark_id) == 0)
      landmark_positions(:,landmark_id) = triangulateLandmarks(landmark_positions,landmark_id, robot_poses_id, robot_positions, global_bearings);
    end
  end
end

function [updated_landmark_position] = triangulateLandmarks(landmark_positions, landmark_id, robot_poses_id, robot_positions, global_bearings)
  
  updated_landmark_position = zeros(3,1);
  notUpdated = 1;
  
  index1 = 1;
  for(index1=1:size(robot_poses_id,2))
    for(index2=index1+1:size(robot_poses_id,2))
      if (robot_poses_id(index2) ~= 0 && notUpdated)
        current_robot_id = robot_poses_id(index1) - 1199;
        current_robot_position = robot_positions(current_robot_id);
        next_robot_id = robot_poses_id(index2) - 1199;
        next_robot_position = robot_positions(next_robot_id);
        #Check if the two poses are close in time
        if ( (next_robot_id - current_robot_id) <= 10)
          bearing1 = normalizeAngle(global_bearings(index1));
          bearing2 = normalizeAngle(global_bearings(index2));
          #Check if parallax is big enough
          if ( rad2deg(abs(bearing2 - bearing1)) > 5 ) 
            m1 = tan(bearing1);
            b1 = -m1*current_robot_position.x + current_robot_position.y;
            #Retrieve the component of the line for the second robot pose
            m2 = tan(bearing2);
            b2 = -m2*next_robot_position.x + next_robot_position.y;
            #Compute x-y positions
            x = (b2-b1)/(m1 - m2);
            y = m1 * x + b1;
            updated_landmark_position = [landmark_id; x; y];
            notUpdated = 0;
          end
        end
      end
    end
  end 
end