function outM = Delta(arg1,arg2,arg3)
  % Arguement check
  if nargin<2
    error('Not enough input arguements')

  elseif nargin==2

    if ischar(arg2)
      ord = arg2;
      p = arg1;
      q = 1;
      bctype = 1;
    else
      bctype = arg2;
      ord = 'xx';
      p = arg1;
      q = 1;
    end

  elseif nargin==3

    if ischar(arg2)
      ord = arg2;
      p = arg1;
      q = 1;
      bctype = arg3;
   
   
    end

  elseif nargin>1
    error('Too many input arguements')

  end

  y = speye(p); % identity matrix for y axis
  ss = p*q;
 
%1D
  if q == 1

    switch ord
    case 'xx'

      XX = spdiags([ones(p,1),-2*ones(p,1),ones(p,1)],-1:1,speye(p));
      XX(1,2) = (bctype-1)*2; XX(p,p-1) = (bctype-1)*2;
      outM = XX;

    case 'xxxx'
      XX = spdiags([ones(p,1),-2*ones(p,1),ones(p,1)],-1:1,speye(p));
      XX(1,2) = (bctype-1)*2; XX(p,p-1) = (bctype-1)*2;
      outM = XX;
      outM = XX^2;

    case 'I'
      I = speye(ss);
      I([1 end],:) = 0;
      outM = I;

    otherwise
      error('check your arguements');

    end
%2D
  else
   switch ord
    case 'xx'
      XX = spdiags([ones(q,1),-2*ones(q,1),ones(q,1)],-1:1,speye(q));
      XX = spdiags([ones(q,1),-2*ones(q,1),ones(q,1)],-1:1,speye(q));
      XX(1,2) = (bctype-1)*2; XX(q,q-1) = (bctype-1)*2;
      outM = kron(XX,y);

    case 'xxxx'
      XX = spdiags([ones(q,1),-2*ones(q,1),ones(q,1)],-1:1,speye(q));
      XX = spdiags([ones(q,1),-2*ones(q,1),ones(q,1)],-1:1,speye(q));
      XX(1,2) = (bctype-1)*2; XX(q,q-1) = (bctype-1)*2;
      outM = kron(XX,y)^2;

    otherwise
      error('check your arguements');

    end

    % alter points to 0
 outM(:,[1:p+1, ss-p:ss]) = 0;outM([1:p+1, ss-p:ss],:) = 0;
 outM(:,[2*p:p:ss, (2*p)+1:p:ss]) = 0; outM([2*p:p:ss, (2*p)+1:p:ss],:)= 0;

  end
