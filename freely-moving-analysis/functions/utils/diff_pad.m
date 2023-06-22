function x_prime = diff_pad(x, dim)
% Same as diff but pads with zero so that x and x_prime have the same
% dimensions

if dim == 2
   x_prime = [zeros(size(x,1), 1), diff(x, [], 2)]; 
elseif dim == 1
   x_prime = [zeros(1, size(x,2)); diff(x, [], 1)]; 
end

end

