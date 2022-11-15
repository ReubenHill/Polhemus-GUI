function [output] = affine_trans_RJC(input,A,B)

%Performs 3d transformation of nx3 'input' using affine transformation matrix
%formed from A (3X3) and B (3x1)

size_in = size(input);

if size_in(2) ~= 3;
    display('ERROR: Input must be of dimensions nx3...');
    return
end

output = ((A*input') + repmat(B,1,size_in(1)))';

