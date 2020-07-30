source "../tools/utilities/geometry_helpers_2d.m"

%(minimal) size of pose and landmarks
global pose_dim=3;
global landmark_dim=2;


# retrieves the index in the perturbation vector, that corresponds to
# a certain pose
# input:
#   pose_index:     the index of the pose for which we want to compute the
#                   index
#   num_poses:      number of pose variables in the state
#   num_landmarks:  number of pose variables in the state
# output:
#   v_idx: the index of the sub-vector corrsponding to 
#          pose_index, in the array of perturbations  (-1 if error)
function v_idx=poseMatrixIndex(pose_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;

  if (pose_index>num_poses)
    v_idx=-1;
    return;
  endif;
  v_idx=1+(pose_index-1)*pose_dim;
endfunction;


# retrieves the index in the perturbation vector, that corresponds to
# a certain landmark
# input:
#   landmark_index:     the index of the landmark for which we want to compute the
#                   index
#   num_poses:      number of pose variables in the state
#   num_landmarks:  number of pose variables in the state
# output:
#   v_idx: the index of the perturnation corrsponding to the
#           landmark_index, in the array of perturbations
function v_idx=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;
  if (landmark_index>num_landmarks)
    v_idx=-1;
    return;
  endif;
  v_idx=1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;
endfunction;

# error and jacobian of a measured landmark
# input:
#   Xr: the robot pose
#   Xl: the landmark pose
#   z:  measured bearing angle of landmark w.r.t. robot
# output:
#   e: 1x1 is the difference between prediction and measurement
#   Jr: 1x3 derivative w.r.t a the error and a perturbation on the
#       pose
#   Jl: 1x2 derivative w.r.t a the error and a perturbation on the
#       landmark
function [e,Jr,Jl]=errorAndJacobian(Xr,Xl,z)
  
  R = Xr(1:2, 1:2);
  t = Xr(1:2,3);
  theta = atan2(R(2,1),R(1,1));

  p_hat = R'*(Xl - t);
  z_hat = atan2(p_hat(2),p_hat(1));
  e = normalizeAngle(z_hat - z);

  Jr=zeros(1,3);
  Jl=zeros(1,2);
  
  atan_der = J_atan2(p_hat);
  #atan_der = [-p_hat(2)/(pow2(p_hat(1))+pow2(p_hat(2))), p_hat(1)/(pow2(p_hat(1))+pow2(p_hat(2)))] 
  #inner_der = R'*[-eye(2), p_hat];
  Jr(1,1:3) = atan_der * R' * [-eye(2) [Xl(2);-Xl(1)]];
  Jl(1,1:2) = atan_der * R';

endfunction;

# implementation of the boxplus
# applies a perturbation to a set of landmarks and robot poses
# input:
#   XR: the robot poses 
#   XL: the landmark pose 
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   dx: the perturbation vector of appropriate dimensions
#       the poses come first, then the landmarks
# output:
#   XR: the robot poses obtained by applying the perturbation
#   XL: the landmarks obtained by applying the perturbation
function [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx)
  global pose_dim;
  global landmark_dim;
  for(pose_index=1:num_poses)
    pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
    dxr=dx(pose_matrix_index:pose_matrix_index+pose_dim-1);
    XR(:,:,pose_index)=v2t(dxr)*XR(:,:,pose_index);
  endfor;
  for(landmark_index=1:num_landmarks)
    landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);
    dxl=dx(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,:);
    XL(:,landmark_index)+=dxl;
  endfor;
endfunction;

# implementation of the optimization loop with robust kernel
# applies a perturbation to a set of landmarks and robot poses
# input:
#   XR: the initial robot poses 
#   XL: the initial landmark estimates
#   Z:  the measurements (1xnum_measurements)
#   associations: 2xnum_measurements. 
#                 associations(:,k)=[p_idx,l_idx]' means the kth measurement
#                 refers to an observation made from pose p_idx, that
#                 observed landmark l_idx
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   num_iterations: the number of iterations of least squares
#   damping:      damping factor (in case system not spd)
#   kernel_threshod: robust kernel threshold

# output:
#   XR: the robot poses after optimization
#   XL: the landmarks after optimization
#   chi_stats: array 1:num_iterations, containing evolution of chi2
#   num_inliers: array 1:num_iterations, containing evolution of inliers
function [XR, XL, chi_stats, num_inliers]=doMultiICP(XR, XL, landmark_ids, Z,
              associations, 
              num_poses, 
							num_landmarks, 
							num_iterations,
              damping, 
							kernel_threshold)
              
  global pose_dim;
  global landmark_dim;
  chi_stats=zeros(1,num_iterations);
  num_inliers=zeros(1,num_iterations);
  system_size=pose_dim*num_poses+landmark_dim*num_landmarks; 
  
  for (iteration=1:num_iterations)
    H=zeros(system_size, system_size);
    b=zeros(system_size,1);
    chi_stats(iteration)=0;
    for (measurement_num=1:size(Z,2))
      pose_index=associations(1,measurement_num);
      landmark_index=associations(2,measurement_num);
      z=Z(1,measurement_num);
      Xr=XR(:,:,pose_index);
      #Search for the landmark id with index in [1, num_landmarks] associated to landmark_index
      for (i=1:num_landmarks)
         if landmark_index == landmark_ids(i)
           landmark_index = i;
         endif
      endfor
      Xl=XL(:,landmark_index);
      [e,Jr,Jl] = errorAndJacobian(Xr, Xl, z);
      chi=e'*e;
      if (chi>kernel_threshold)
      	e*=sqrt(kernel_threshold/chi);
     	  chi=kernel_threshold;
      else
      	num_inliers(iteration)++;
      endif;
      chi_stats(iteration)+=chi;

      Hrr=Jr'*Jr;
      Hrl=Jr'*Jl;
      Hll=Jl'*Jl;
      br=Jr'*e;
      bl=Jl'*e;

      pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
      landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);
      H(pose_matrix_index:pose_matrix_index+pose_dim-1, pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hrr;

      H(pose_matrix_index:pose_matrix_index+pose_dim-1, landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hrl;

      H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1, landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hll;

      H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1, pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hrl';

      b(pose_matrix_index:pose_matrix_index+pose_dim-1)+=br;
      b(landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=bl;
    endfor
    H+=eye(system_size)*damping;
    dx=zeros(system_size,1);
    % we solve the linear system, blocking the first pose
    % this corresponds to "remove" from H and b the locks
    % of the 1st pose, while solving the system
    #dx(pose_dim+1:end)= -(H(pose_dim+1:end,pose_dim+1:end)\b(pose_dim+1:end,1));
    dx = -H\b;
    
    [XR, XL] = boxPlus(XR, XL, num_poses, num_landmarks, dx);
    #figure;
    #plot(XR(1,3,:), XR(2,3,:), 'r');
    #hold on;
    #scatter(XL(1,:), XL(2,:), 'b');
  endfor
endfunction