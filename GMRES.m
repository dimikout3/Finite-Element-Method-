function x=GMRES(A,b,tol,N,restart,M)
    x=zeros(size(b));
    e1=[1;zeros(restart,1)];
    r=M*b;
    for i=1:N
    [V,H]=arnoldi(r,A,restart,M);
    [Q,R]=qr(H,0);
    y=R\Q'*norm(r)*e1;
    x=x+V*y;
    r=M*(b-A*x);
    if norm(r)<tol*norm(b)
        break;
    end
    end
    if i==N
    disp('GMRES: Max Iterations Reached');
    else
    disp(['GMRES: Iterations: ' num2str(i) ' Residual: ' num2str(norm(r))]);
    end