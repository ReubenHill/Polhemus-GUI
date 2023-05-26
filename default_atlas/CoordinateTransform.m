%This function outputs a vector to transform inion to origin
%also outputs rotation transformation matrix to allign the head model to
%intuitive coordinates.
%To apply to row vector, the vector should be multiplied by the transpose
%with the vector on the left
function [Matrix,vector] = CoordinateTransform(landmarks)
  %calculate lengths between vectors

  %these are the untransformed reference points, defined here in case I want
  %to use them
  Nasion = landmarks(1,:);
  Inion = landmarks(2,:);
  Ar = landmarks(3,:);
  Al = landmarks(4,:);
  Cz = landmarks(5,:);

  %------TRANSLATION-------

  %translate inion to origion
  vector = -Inion;

  %translate Al
  Al = Al + vector;

  %------ROTATE AL TO Y AXIS-------

  %calculate rotation to y axis
  AlToYAxisRot = vrrotvec(Al,[0,1,0]);

  %convert to rotation matrix
  AlToYAxisMatrix = vrrotvec2mat(AlToYAxisRot);

  %repmat
  % Apply translation and rotation to points
  for k = 1:5
      landmarks(k,:) = landmarks(k,:)+vector;
      landmarks(k,:) = landmarks(k,:)*AlToYAxisMatrix';
  end

  %------ROTATE AR TO XY PLANE ABOUT INION-AL AXIS-------

  % find angle to rotate nasion into XY plane about the new y axis
  [ArToXYRotAngle,~] = cart2pol(landmarks(3,1),landmarks(3,3));

  %Find second rotation matrix
  ArToXYMatrix = vrrotvec2mat([0,1,0,ArToXYRotAngle]);

  % Apply second rotation to points
  for k = 1:5
      landmarks(k,:) = landmarks(k,:)*ArToXYMatrix';
  end

  %------FINAL ROTATION ABOUT AL-AR AXIS TO ALLIGN INION AND NASION-------

  %find angle of nasion to xy plane
  [~,NasionToXYRotAngle,~] = cart2sph(landmarks(1,1),landmarks(1,2),landmarks(1,3));
  %define vector to rotate around (the line joining AL and AR)
  NasionRotVector = landmarks(4,:) - landmarks(3,:);
  %find rotation matrix
  NasionToXYRotMatrix = vrrotvec2mat([NasionRotVector, NasionToXYRotAngle]);

  %------OUTPUT FINAL MATRIX-------
  Matrix = NasionToXYRotMatrix*ArToXYMatrix*AlToYAxisMatrix;

end