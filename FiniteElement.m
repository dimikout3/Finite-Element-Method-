clear;
clc;
%c=1; % Source coefficient
m=30; % Number of Squares per dim
h=2/m; % Mesh size
ne=m^2; % Number of elements
non=(m+1)^2; % Number of nodes
coo=zeros(non,2); % Initializing the matrix of coordinates
nn=zeros(ne,4); % Initializing the matrix of elements
t1=(-1:h:1)';
t2=(0:h:1)';
ax=1;
ay=1;
c=1;

bounds=[];
bounds=[bounds 1:m+1  (m+2):(m+1):((m+1)^2) ];

for i=1:m+1
    coo(((i-1)*(m+1)+1):i*(m+1),1)=t1;
    coo(((i-1)*(m+1)+1):i*(m+1),2)=(i-1)*h-1;
end

k=1;
for i=1:m
    for j=1:m
        nn(k,:)=[j j+1 j+m+2 j+m+1]+(i-1)*(m+1);
        k=k+1;
    end
end
K=spalloc(non,non,6*non);
b=zeros(non,1);


    J=(h^2)/4;
    M(1,1)=ax/3+ay/3;
    M(1,2)=-ax/3+ay/6;
    M(1,3)=-ax/6-ay/6;
    M(1,4)=ax/6-ay/3;
    M(2,1)=M(1,2);
    M(2,2)=ax/3+ay/3;
    M(2,3)=ax/6-ay/3;
    M(2,4)=-ax/6-ay/6;
    M(3,1)=M(1,3);
    M(3,2)=M(2,3);
    M(3,3)=ax/3+ay/3;
    M(3,4)=-ax/3+ay/6;
    M(4,1)=M(1,4);
    M(4,2)=M(2,4);
    M(4,3)=M(3,4);
    M(4,4)=ax/3+ay/3;
    
    T(1,1)=64*c*J/144;
    T(1,2)=32*c*J/144;
    T(1,3)=16*c*J/144;
    T(1,4)=32*c*J/144;
    T(2,1)=T(1,2);
    T(2,2)=64*c*J/144;
    T(2,3)=32*c*J/144;
    T(2,4)=16*c*J/144;
    T(3,1)=T(1,3);
    T(3,2)=T(2,3);
    T(3,3)=64*c*J/144;
    T(3,4)=32*c*J/144;
    T(4,1)=T(1,4);
    T(4,2)=T(2,4);
    T(4,3)=T(3,4);
    T(4,4)=64*c*J/144;
    
    F(1,1)=J;
    F(2,1)=J;
    F(3,1)=J;
    F(4,1)=J;

for i=1:ne
    x1=coo(nn(i,1),1);
    x2=coo(nn(i,2),1);
    x3=coo(nn(i,3),1);
    x4=coo(nn(i,4),1);
    y1=coo(nn(i,1),2);
    y2=coo(nn(i,2),2);
    y3=coo(nn(i,3),2);
    y4=coo(nn(i,4),2);
    

    for j=1:4
        for k=1:4
            K(nn(i,j),nn(i,k))=K(nn(i,j),nn(i,k))+M(j,k)+T(j,k);
        end
        b(nn(i,j))=b(nn(i,j))+F(j);
    end
end
K(bounds,:)=[];
K(:,bounds)=[];
b(bounds)=[];

%   M1=diag(diag(K)) % Jacobi precondition
%   M2=tril(K); % Gauss precondition

%    tic
%      x1=pcg(K,b,1e-10,7000,M2);
%    toc

% tic
% x1=GMRES(K,b,1e-10,1000,10,M1);
% toc
% tic
% x2=GMRES(K,b,1e-10,1000,20,M1);
% toc
% tic
% x3=GMRES(K,b,1e-10,1000,40,M1);
% toc
% tic
% x4=GMRES(K,b,1e-10,1000,80,M1);
% toc

% tic
% x5=GMRES(K,b,1e-10,1000,10,M2);
% toc
% tic
% x6=GMRES(K,b,1e-10,1000,20,M2);
% toc
% tic
% x7=GMRES(K,b,1e-10,1000,40,M2);
% toc
% tic
% x8=GMRES(K,b,1e-10,1000,80,M2);
% toc

%  tic
   x1=BicgstabUnpreco( K,b,1000); % no precondition
%  toc
%  tic
%   x2=bicgstabPreco( K,b,500,M1); % jacobi precondition
%  toc
%  tic
%  x3=BicgstabPreco( K,b,1000,M2); % Gauss precondition 
%  toc

 mesh(reshape(x1,m,m));

