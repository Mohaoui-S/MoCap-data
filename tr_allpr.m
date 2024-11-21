function Z = tr_allpr(L,A,tr)

  N = length(A);  Z = L;
  for n = 1:N
    if ~isempty(A{n})
      if tr == 0
        Z = tmult(Z,A{n},n);
      else
        Z = tmult(Z,A{n}',n);
      end
    end
  end
end
