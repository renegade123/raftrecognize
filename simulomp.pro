FUNCTION SimulOMP,S, Phi, sigma, T, normType
  ;  % Simultaneous OMP, specify the type OF norm used FOR selecting the atoms
  ;  % based on paper (normType = 1 in Tropp's algorithm)
  ;  % "Algorithms for Simultaneous Sparse Approximation Part I: Greedy  Pursuit"
  ;  % J. Tropp, A. Gilbert, and M. Strauss
  ;  % MIN ||Coeff||_{row,0}   subject to: S = Phi * Coeff
  ;  % Input:       S -- signal to be approximated, d x K matrix
  ;  %            Phi -- DICTIONARY
  ;  %          sigma -- error tolerance
  ;  %              T -- number OF iterations
  ;  %       normType -- type OF norm, 1 FOR L1, 2 FOR L2, 3 FOR infinity norm
  ;  % Ouptut:  Coeff -- coefficients, a sparse matrix with only T nonzero rows
  ;  %         indSet -- index set FOR selected atoms (common support)
  ;  %         Approx -- approximation OF S
  ;  %            Res -- residuals
  ;
  ;  % normalizing columns OF the DICTIONARY
  ;  %Phi0 = Phi;
  ;  %Phi = zeros(SIZE(Phi0));
  phiSize = SIZE(Phi)
  sSize = SIZE(S)
  N = phiSize[1]; % the number of atoms
  d = sSize[1];
  K = sSize[2];
  Coeff = MAKE_ARRAY(N,d,VALUE=0,/DOUBLE)
  Res = S;
  indSet = MAKE_ARRAY(T,1,VALUE=0,/DOUBLE)

  iter = 0;
  norm_res = MAKE_ARRAY(T,1,VALUE=0,/DOUBLE);
  WHILE ((norm(Res[*]) GT sigma) && (iter LT T)) DO BEGIN
    ;% compute the projection
    IF normType EQ 1 THEN BEGIN
      temp = SORT(TOTAL(ABS(TRANSPOSE(Phi)##Res),1))
      dim = reverse(temp)
    ENDIF ELSE BEGIN
      IF normType EQ 2 THEN temp = reverse(SORT(TOTAL(ABS(TRANSPOSE(Phi)##Res)^2,2))) ELSE temp = reverse(SORT(MAX(ABS(TRANSPOSE(Phi)##Res), DIMENSION=2)))
    ENDELSE
    ;% update the index set
    val = temp[0]
    idx = temp
    indSet[iter] = idx[0];
    Coeff[indSet[0:iter],*] =  PSEUDO_INVERSE(Phi[indSet[0:iter],*])##S;
    Approx = Phi##Coeff;
    Res = S - Approx;
    norm_res(iter) = norm(Res[*]);
    iter = iter + 1;
  END
  indSet = indSet[1:iter-1];
  RETURN,Coeff
END