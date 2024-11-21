function u = innerprod_Ar_exc(E,U,k,n)

  N = length(U);
  for m = 1:N
    V{m} = U{m}(:,k);
  end
  u = tr_allpr_ex(E,V,1,n);
