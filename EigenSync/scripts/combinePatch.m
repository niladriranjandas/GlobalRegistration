function [c] = combinePatch( Patch )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


c=Patch{1};

for i=2:length(Patch)

    c=vertcat(c,Patch{i});
end

end

