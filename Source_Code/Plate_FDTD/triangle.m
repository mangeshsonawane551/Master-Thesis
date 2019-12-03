function ANS = triangle(l,m,operator,bctype)
  I_y = speye(l); % for y axis
  I_x = speye(m); % for x axis
  ss = l*m;
  XX = spdiags([ones(m,1),-2*ones(m,1),ones(m,1)],-1:1,I_x);
  XX(1,2) = (bctype-1)*2; XX(m,m-1) = (bctype-1)*2;
  YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,I_y);
  YY(1,2) = (bctype-1)*2; YY(l,l-1) = (bctype-1)*2;
  % laplacian
  if operator== 1
  ANS = kron(XX,I_y) + kron(I_x,YY);
  %biharmonic
  elseif operator==2
          Lap = kron(XX,I_y) + kron(I_x,YY);
          ANS = Lap*Lap;
  end
  %update points
ANS(:,[1:l+1, ss-l:ss]) = 0;
ANS([1:l+1, ss-l:ss],:) = 0;
ANS(:,[2*l:l:ss, (2*l)+1:l:ss]) = 0; 
ANS([2*l:l:ss, (2*l)+1:l:ss],:) = 0;
end