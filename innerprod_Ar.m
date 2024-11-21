function u = innerprod_Ar(E,U,k)

  N = length(U);
  for m = 1:N
    V{m} = U{m}(:,k);
  end
  u = tr_allpr(E,V,1);
