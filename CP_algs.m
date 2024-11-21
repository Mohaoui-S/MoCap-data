function [Fseq, out] = CP_algs(mcsequence,Data,opts)

%% CP decomposition algorithms for MOCAP compeltion problem %%
  
warning off; 

%% Parameters  
if isfield(opts, 'tol')   tol = opts.tol;           else  tol = 1e-4;      end
if isfield(opts, 'eps')   eps = opts.eps;           else  eps =1e-3;       end
if isfield(opts, 'K')     K = opts.K;               else  K = 10;          end
if isfield(opts, 'Rmax')  maxR = opts.Rmax;         else  maxR = Inf;      end
if isfield(opts, 'maxiter') maxiter = opts.maxiter; else  maxiter = 1000;  end
if isfield(opts, 'maxiter') maxK = opts.maxiter;    else  maxK = 1000;     end
if isfield(opts, 'tol_cnv') tol_cnv = opts.tol_cnv; else  tol_cnv = 1e-6;  end
gama=opts.gama; rho=opts.rho;
CP = opts.model;
%% C3D to Tensor format
n_f = mcsequence.nFrames;
n_m = mcsequence.nMarkers;

incomplete = reshape(mcsequence.data, n_f, n_m,3);
Data=reshape(Data, n_f, n_m,3);

Omeg = ~isnan(incomplete);
F = zeros(size(incomplete));
F(Omeg) = incomplete(Omeg);


%% Main
R  = 1;  Len = length(Data(:));
N  = ndims(F);  Nt = size(F);  NN = sum_all(Omeg);
Er  = zeros(Nt);  X  = Er;  Y  = X;  Ar= Er;
X(Omeg) = F(Omeg);  X(~Omeg) = sum_all(X)/NN;  

for n = 1:N
  A{n} = randn(Nt(n),R);
  A{n} = A{n}/max(max(A{n}));
end
 
for n = 1:N
  dd  = eye(Nt(n)); dd(1,:) = [];
  D{n} = eye(Nt(n)-1,Nt(n)) - dd;
  DtD{n} = rho(n)*(D{n}'*D{n}) + eye(Nt(n));
  A{n} = (inv(rho(n)*DtD{n} + eye(Nt(n))))*randn(Nt(n),R);
  A{n} = A{n}/norm(A{n});
end
 
Lam = tr_allpr(X,A,1); 
[~, ind] = sort(abs(Lam));
for r = 1:R
 Y = Y + Lam(r)*outerprod(A,r);
end
  
obj = sum_all((F(Omeg) - Y(Omeg)).^2);
X(~Omeg) = Y(~Omeg); 
Er = X - Y; 
Xold=X;
Cherr=zeros(1,maxiter);
RMSE=zeros(1,maxiter);
ReE=zeros(1,maxiter);

for iter = 1:maxiter
 for r = ind(1:min(end,K))
     
  Ar = Lam(r)*outerprod(A,r);
  Y = Y - Ar; 
  Er = Er + Ar;
  
   for n = 1:N
      
    if strcmp(CP,'Cp')   
           a = innerprod_Ar_exc(Er,A,r,n);
           a = a(:); a=a/Lam(r);a = a/norm(a);
           A{n}(:,r) =a; b{n} = a;
          
    elseif strcmp(CP,'sparse')
           a = innerprod_Ar_exc(Er,A,r,n);
           a = a(:); a0 = a/norm(a); 
          for nn = 1:maxK
            DH = (rho(n)*sign(a0) - a + Lam(r)*a0)/Lam(r);
            a0 = a0 - gama*DH; a0 = a0/norm(a0);
          end   
          A{n}(:,r) = a0; b{n} = a0;
  
    elseif strcmp(CP,'smooth') 
           a = innerprod_Ar_exc(Er,A,r,n);
           a = a(:); a0 = a/norm(a);
           for nn = 1:maxK
             DH = (DtD{n}*a0 - a + Lam(r)*a0 )/Lam(r);
             a0 = (a0 - gama*DH)/norm(a0); a = a0/norm(a0);            
           end
           A{n}(:,r) = a0; b{n} = a0;
    end
    
   end
  
 % udpating lambda_r  
  Lam(r) = tr_allpr(Er,b,1);
  if Lam(r) < 0
    Lam(r) = -Lam(r);
    A{1}(:,r) = - A{1}(:,r);
  end
   
 % the residues
   Ar = Lam(r)*outerprod(A,r);    
   Er = Er - Ar; Er(~Omeg) = 0; Y = Y + Ar;
    
 % updating X
   X(~Omeg) = Y(~Omeg); X(Omeg)= F(Omeg);     
 end
 
  Rech= (norm(X(:) - Xold(:)))/norm(Xold(:));  
  Rerr= (norm(X(:) - Data(:)))/norm(Data(:));
  Xold=X;  
  % rank update   
  objec = sum_all(Er.^2);
  relative = abs(objec - obj)/abs(objec);
  % if relative < eps
  if relative < eps || abs(objec-obj)/NN < tol    
    if R ~= maxR
      R = R + 1;
      for n = 1:N
       a = randn(Nt(n),1);a = a/norm(a); 
       A{n}(:,R) = a; b{n} = a;
      end
       Lam(R) = tr_allpr(Er,b,1); 
       Er = Er - Lam(R)*outerprod(A,R); Er(~Omeg) = 0;      
    else
      break;
    end
  end 
 rmse = sqrt(sum((X(:) - Data(:)).^2)/ Len); obj = objec;
  if mod(iter,5)==0
    fprintf('%d:  %f   \n',iter,rmse);
    %fprintf('%d: \n',iter);
  end
  [~, ind] = sort(abs(Lam));   
    
   % Stoping criteria       
    if Rech < tol_cnv
       break;
    end
  
 RMSE(iter)=rmse;
 Cherr(iter)=Rech;
 ReE(iter) =Rerr;
end
out=[];out.X=X;
out.Cherr=Cherr;
out.RMSE=RMSE;
out.Rer=ReE;


%% Tensor to C3D format 
F_matrix = reshape(X, [size(X, 1) 3*size(X, 2)]);
mcsequence.data =double(F_matrix);
Fseq=mcsequence;  
end



