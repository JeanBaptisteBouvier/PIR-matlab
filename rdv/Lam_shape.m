function Lam_shaped = Lam_shape( Lam )
% Lam_shaped change the order of the components of the STM saved in every
% line of Lam, from position 13 to 48.
% In fact ode78_... saved it columns by columns but the .m algorithm use it
% as if it was saved line by line
% 
% 2017
% JBB

[lines, ~] = size(Lam);
Lam_shaped = zeros(size(Lam));
for i=1:lines
    M = reshape(Lam(i,13:end), [6,6]);
    Lam_shaped(i, 1:12) = Lam(i, 1:12); 
    Lam_shaped(i, 13:end) = reshape(M', [1,36]);
end



end

