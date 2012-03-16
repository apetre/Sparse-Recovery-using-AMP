function [mse] = mse(x,y)
% mse calculates the Mean Squared Error between two vectors
mse = mean((x-y).^2);
end

