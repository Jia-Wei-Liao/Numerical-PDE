function F = FivePointF(F, U, dx, dy)

F = (circshift(F,-1) + circshift(F,1) + ...
     circshift(F,[0 -1]) + circshift(F,[0 1]) + 8*F)/12;

F = F - ( 4*( circshift(U,1)     + circshift(U,-1)     + circshift(U, [0 1])  + circshift(U, [0 -1])  ) ...
        + 1*( circshift(U, [1 1]) + circshift(U, [1 -1]) + circshift(U, [-1 1]) + circshift(U, [-1 -1]) ))/(6*dx*dy);

F = F(2:end-1, 2:end-1);
end