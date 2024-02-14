% A DFO Method Using a New Model
% Codes for the paper entitled
% "A Derivative-free Method Using a New Under-determined Quadratic Interpolation Model"
% Copyright: Pengcheng Xie & Ya-xiang Yuan

% Connect: xpc@lsec.cc.ac.cn
% A DFO Method Using a New Model
% ----------------------------------------------------------
% License Information

% ----------------------------------------------------------
% This code is distributed under the MIT License.
% You should have received a copy of the MIT License along
% with this program. If not, see <https://opensource.org/licenses/MIT>.

% For further information or questions, contact the authors
% via the provided email address.
% ----------------------------------------------------------
% Code Version Information

% ----------------------------------------------------------
% Version: 1.0
% Changes: Initial release.
% ----------------------------------------------------------

% ----------------------------------------------------------
% References
% ----------------------------------------------------------
% For more information, refer to the paper:

% "A Derivative-free Method Using a New Under-determined Quadratic Interpolation Model"
% by Pengcheng Xie & Ya-xiang Yuan.
%
% If you use this code in your research, please cite the above paper.

% ----------------------------------------------------------
% ----------------------------------------------------------
% Contributors
% ----------------------------------------------------------

% This code was written by Pengcheng Xie & Ya-xiang Yuan.
% ----------------------------------------------------------
function [W, VLAG, BETA, S] = ...
    bigden (N, NPT, XOPT, XPT, BMAT, ZMAT, IDZ, NDIM, KOPT, ...
    KNEW, D, W, VLAG, BETA, S)
  % IMPLICIT REAL(8) (A-H,O-Z)
  % DIMENSION XOPT(*),XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),D(*), &
  %    W(*),VLAG(*),S(*),WVEC(NDIM,*),PROD(NDIM,*)
  % DIMENSION DEN(9),DENEX(9),PAR(9)
  [WVEC, PROD] = deal(0);
  DEN = zeros(1, 9);
  DENEX = zeros(1, 9);
  PAR = zeros(1, 9);
  %     N is the number of variables.
  %     NPT is the number of interpolation equations.
  %     XOPT is the best interpolation point so far.
  %     XPT contains the coordinates of the current interpolation points.
  %     BMAT provides the last N columns of H.
  %     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
  %     NDIM is the first dimension of BMAT and has the value NPT+N.
  %     KOPT is the index of the optimal interpolation point.
  %     KNEW is the index of the interpolation point that is going to be moved.
  %     D will be set to the step from XOPT to the new point, and on entry it
  %       should be the D that was calculated by the last call of biglag. The
  %       length of the initial D provides a trust region bound on the final D.
  %     W will be set to Wcheck for the final choice of D.
  %     VLAG will be set to Theta*Wcheck+e_b for the final choice of D.
  %     BETA will be set to the value that will occur in the updating formula
  %       when the KNEW-th interpolation point is moved to its new position.
  %     S, WVEC, PROD and the private arrays DEN, DENEX and PAR will be used
  %       for working space.

  %     D is calculated in a way that should provide a denominator with a large
  %     modulus in the updating formula when the KNEW-th interpolation point is
  %     shifted to the new position XOPT+D.

  %     Set some constants.

  HALF = 0.5e0;
  ONE = 1.0e0;
  QUART = 0.25e0;
  TWO = 2.0e0;
  ZERO = 0.0e0;
  TWOPI = 8.0e0 * atan(ONE);
  NPTM = NPT - N - 1;

  %     Store the first NPT elements of the KNEW-th column of H in W(N+1)
  %     to W(N+NPT).

  for K = 1:NPT
    W(N + K) = ZERO;
  end

  for J = 1:NPTM
    TEMP = ZMAT(KNEW, J);
    if (J < IDZ)
      TEMP = -TEMP;
    end
    for K = 1:NPT
      W(N + K) = W(N + K) + TEMP * ZMAT(K, J);
    end
  end
  ALPHA = W(N + KNEW);
  %     The initial search direction D is taken from the last call of biglag,
  %     and the initial S is set below, usually to the direction from X_OPT
  %     to X_KNEW, but a different direction to an interpolation point may
  %     be chosen, in order to prevent S from being nearly parallel to D.
  DD = ZERO;
  DS = ZERO;
  SS = ZERO;
  XOPTSQ = ZERO;
  for I = 1:N
    DD = DD + D(I)^2;
    S(I) = XPT(KNEW, I) - XOPT(I);
    DS = DS + D(I) * S(I);
    SS = SS + S(I)^2;
    XOPTSQ = XOPTSQ + XOPT(I)^2;
  end

  if (DS * DS > 0.99e0 * DD * SS)
    KSAV = KNEW;
    DTEST = DS * DS / SS;
    for K = 1:NPT
      if (K ~= KOPT)
        DSTEMP = ZERO;
        SSTEMP = ZERO;
        for I = 1:N
          DIFF = XPT(K, I) - XOPT(I);
          DSTEMP = DSTEMP + D(I) * DIFF;
          SSTEMP = SSTEMP + DIFF * DIFF;
        end
        if (DSTEMP * DSTEMP / SSTEMP < DTEST)
          KSAV = K;
          DTEST = DSTEMP * DSTEMP / SSTEMP;
          DS = DSTEMP;
          SS = SSTEMP;
        end
      end
    end

    for I = 1:N
      S(I) = XPT(KSAV, I) - XOPT(I);
    end
  end

  SSDEN = DD * SS - DS * DS;
  ITERC = 0;
  DENSAV = ZERO;
  %     Begin the iteration by overwriting S with a vector that has the
  %     required length and direction.
  while (1)
    ITERC = ITERC + 1;
    TEMP = ONE / sqrt(SSDEN);
    XOPTD = ZERO;
    XOPTS = ZERO;
    for I = 1:N
      S(I) = TEMP * (DD * S(I) - DS * D(I));
      XOPTD = XOPTD + XOPT(I) * D(I);
      XOPTS = XOPTS + XOPT(I) * S(I);
    end
    %     Set the coefficients of the first two terms of BETA.
    TEMPA = HALF * XOPTD * XOPTD;
    TEMPB = HALF * XOPTS * XOPTS;
    DEN(1) = DD * (XOPTSQ + HALF * DD) + TEMPA + TEMPB;
    DEN(2) = TWO * XOPTD * DD;
    DEN(3) = TWO * XOPTS * DD;
    DEN(4) = TEMPA - TEMPB;
    DEN(5) = XOPTD * XOPTS;
    for I = 6:9
      DEN(I) = ZERO;
    end
    %     Put the coefficients of Wcheck in WVEC.
    for K = 1:NPT
      TEMPA = ZERO;
      TEMPB = ZERO;
      TEMPC = ZERO;
      for I = 1:N
        TEMPA = TEMPA + XPT(K, I) * D(I);
        TEMPB = TEMPB + XPT(K, I) * S(I);
        TEMPC = TEMPC + XPT(K, I) * XOPT(I);
      end
      WVEC(K, 1) = QUART * (TEMPA * TEMPA + TEMPB * TEMPB);
      WVEC(K, 2) = TEMPA * TEMPC;
      WVEC(K, 3) = TEMPB * TEMPC;
      WVEC(K, 4) = QUART * (TEMPA * TEMPA - TEMPB * TEMPB);
      WVEC(K, 5) = HALF * TEMPA * TEMPB;
    end

    for I = 1:N
      IP = I + NPT;
      WVEC(IP, 1) = ZERO;
      WVEC(IP, 2) = D(I);
      WVEC(IP, 3) = S(I);
      WVEC(IP, 4) = ZERO;
      WVEC(IP, 5) = ZERO;
    end

    %     Put the coefficents of THETA*Wcheck in PROD.
    for JC = 1:5
      NW = NPT;
      if (JC == 2 || JC == 3)
        NW = NDIM;
      end
      for K = 1:NPT
        PROD(K, JC) = ZERO;
      end
      for J = 1:NPTM
        SUM = ZERO;
        for K = 1:NPT
          SUM = SUM + ZMAT(K, J) * WVEC(K, JC);
        end
        if (J < IDZ)
          SUM = -SUM;
        end
        for K = 1:NPT
          PROD(K, JC) = PROD(K, JC) + SUM * ZMAT(K, J);
        end
      end

      if (NW == NDIM)
        for K = 1:NPT
          SUM = ZERO;
          for J = 1:N
            SUM = SUM + BMAT(K, J) * WVEC(NPT + J, JC);
          end
          PROD(K, JC) = PROD(K, JC) + SUM;
        end
      end

      for J = 1:N
        SUM = ZERO;
        for I = 1:NW
          SUM = SUM + BMAT(I, J) * WVEC(I, JC);
        end
        PROD(NPT + J, JC) = SUM;
      end
    end
    %     Include in DEN the part of BETA that depends on THETA.
    for K = 1:NDIM
      SUM = ZERO;
      for I = 1:5
        PAR(I) = HALF * PROD(K, I) * WVEC(K, I);
        SUM = SUM + PAR(I);
      end
      DEN(1) = DEN(1) - PAR(1) - SUM;
      TEMPA = PROD(K, 1) * WVEC(K, 2) + PROD(K, 2) * WVEC(K, 1);
      TEMPB = PROD(K, 2) * WVEC(K, 4) + PROD(K, 4) * WVEC(K, 2);
      TEMPC = PROD(K, 3) * WVEC(K, 5) + PROD(K, 5) * WVEC(K, 3);
      DEN(2) = DEN(2) - TEMPA - HALF * (TEMPB + TEMPC);
      DEN(6) = DEN(6) - HALF * (TEMPB - TEMPC);
      TEMPA = PROD(K, 1) * WVEC(K, 3) + PROD(K, 3) * WVEC(K, 1);
      TEMPB = PROD(K, 2) * WVEC(K, 5) + PROD(K, 5) * WVEC(K, 2);
      TEMPC = PROD(K, 3) * WVEC(K, 4) + PROD(K, 4) * WVEC(K, 3);
      DEN(3) = DEN(3) - TEMPA - HALF * (TEMPB - TEMPC);
      DEN(7) = DEN(7) - HALF * (TEMPB + TEMPC);
      TEMPA = PROD(K, 1) * WVEC(K, 4) + PROD(K, 4) * WVEC(K, 1);
      DEN(4) = DEN(4) - TEMPA - PAR(2) + PAR(3);
      TEMPA = PROD(K, 1) * WVEC(K, 5) + PROD(K, 5) * WVEC(K, 1);
      TEMPB = PROD(K, 2) * WVEC(K, 3) + PROD(K, 3) * WVEC(K, 2);
      DEN(5) = DEN(5) - TEMPA - HALF * TEMPB;
      DEN(8) = DEN(8) - PAR(4) + PAR(5);
      TEMPA = PROD(K, 4) * WVEC(K, 5) + PROD(K, 5) * WVEC(K, 4);
      DEN(9) = DEN(9) - HALF * TEMPA;
    end
    %     Extend DEN so that it holds all the coefficients of DENOM.
    SUM = ZERO;
    for I = 1:5
      PAR(I) = HALF * PROD(KNEW, I)^2;
      SUM = SUM + PAR(I);
    end
    DENEX(1) = ALPHA * DEN(1) + PAR(1) + SUM;
    TEMPA = TWO * PROD(KNEW, 1) * PROD(KNEW, 2);
    TEMPB = PROD(KNEW, 2) * PROD(KNEW, 4);
    TEMPC = PROD(KNEW, 3) * PROD(KNEW, 5);
    DENEX(2) = ALPHA * DEN(2) + TEMPA + TEMPB + TEMPC;
    DENEX(6) = ALPHA * DEN(6) + TEMPB - TEMPC;
    TEMPA = TWO * PROD(KNEW, 1) * PROD(KNEW, 3);
    TEMPB = PROD(KNEW, 2) * PROD(KNEW, 5);
    TEMPC = PROD(KNEW, 3) * PROD(KNEW, 4);
    DENEX(3) = ALPHA * DEN(3) + TEMPA + TEMPB - TEMPC;
    DENEX(7) = ALPHA * DEN(7) + TEMPB + TEMPC;
    TEMPA = TWO * PROD(KNEW, 1) * PROD(KNEW, 4);
    DENEX(4) = ALPHA * DEN(4) + TEMPA + PAR(2) - PAR(3);
    TEMPA = TWO * PROD(KNEW, 1) * PROD(KNEW, 5);
    DENEX(5) = ALPHA * DEN(5) + TEMPA + PROD(KNEW, 2) * PROD(KNEW, 3);
    DENEX(8) = ALPHA * DEN(8) + PAR(4) - PAR(5);
    DENEX(9) = ALPHA * DEN(9) + PROD(KNEW, 4) * PROD(KNEW, 5);
    %     Seek the value of the angle that maximizes the modulus of DENOM.
    SUM = DENEX(1) + DENEX(2) + DENEX(4) + DENEX(6) + DENEX(8);
    DENOLD = SUM;
    DENMAX = SUM;
    ISAVE = 0;
    IU = 49;
    TEMP = TWOPI / (IU + 1);
    PAR(1) = ONE;
    for I = 1:IU
      ANGLE = (I) * TEMP;
      PAR(2) = cos(ANGLE);
      PAR(3) = sin(ANGLE);
      for J = 4:2:8
        PAR(J) = PAR(2) * PAR(J - 2) - PAR(3) * PAR(J - 1);
        PAR(J + 1) = PAR(2) * PAR(J - 1) + PAR(3) * PAR(J - 2);
      end
      SUMOLD = SUM;
      SUM = ZERO;
      for J = 1:9
        SUM = SUM + DENEX(J) * PAR(J);
      end
      if (abs(SUM) > abs(DENMAX))
        DENMAX = SUM;
        ISAVE = I;
        TEMPA = SUMOLD;
      elseif (I == ISAVE + 1)
        TEMPB = SUM;
      end
    end
    if (ISAVE == 0)
      TEMPA = SUM;
    end
    if (ISAVE == IU)
      TEMPB = DENOLD;
    end
    STEP = ZERO;
    if (TEMPA ~= TEMPB)
      TEMPA = TEMPA - DENMAX;
      TEMPB = TEMPB - DENMAX;
      STEP = HALF * (TEMPA - TEMPB) / (TEMPA + TEMPB);
    end
    ANGLE = TEMP * ((ISAVE) + STEP);
    %     Calculate the new parameters of the denominator, the new VLAG vector
    %     and the new D. Then test for convergence.
    PAR(2) = cos(ANGLE);
    PAR(3) = sin(ANGLE);
    for J = 4:2:8
      PAR(J) = PAR(2) * PAR(J - 2) - PAR(3) * PAR(J - 1);
      PAR(J + 1) = PAR(2) * PAR(J - 1) + PAR(3) * PAR(J - 2);
    end
    BETA = ZERO;
    DENMAX = ZERO;
    for J = 1:9
      BETA = BETA + DEN(J) * PAR(J);
      DENMAX = DENMAX + DENEX(J) * PAR(J);
    end
    for K = 1:NDIM
      VLAG(K) = ZERO;
      for J = 1:5
        VLAG(K) = VLAG(K) + PROD(K, J) * PAR(J);
      end
    end
    TAU = VLAG(KNEW);
    DD = ZERO;
    TEMPA = ZERO;
    TEMPB = ZERO;
    for I = 1:N
      D(I) = PAR(2) * D(I) + PAR(3) * S(I);
      W(I) = XOPT(I) + D(I);
      DD = DD + D(I)^2;
      TEMPA = TEMPA + D(I) * W(I);
      TEMPB = TEMPB + W(I) * W(I);
    end
    if (ITERC >= N)
      break
    end
    if (ITERC > 1)
      DENSAV = max(DENSAV, DENOLD);
    end
    if (abs(DENMAX) <= 1.1e0 * abs(DENSAV))
      break
    end
    DENSAV = DENMAX;
    %     Set S to half the gradient of the denominator with respect to D.
    %     Then branch for the next iteration.
    for I = 1:N
      TEMP = TEMPA * XOPT(I) + TEMPB * D(I) - VLAG(NPT + I);
      S(I) = TAU * BMAT(KNEW, I) + ALPHA * TEMP;
    end
    for K = 1:NPT
      SUM = ZERO;
      for J = 1:N
        SUM = SUM + XPT(K, J) * W(J);
      end
      TEMP = (TAU * W(N + K) - ALPHA * VLAG(K)) * SUM;
      for I = 1:N
        S(I) = S(I) + TEMP * XPT(K, I);
      end
    end
    SS = ZERO;
    DS = ZERO;
    for I = 1:N
      SS = SS + S(I)^2;
      DS = DS + D(I) * S(I);
    end
    SSDEN = DD * SS - DS * DS;
    if (SSDEN >= 1.0e-8 * DD * SS)
      continue
    else
      break
    end
  end
  %     Set the vector W before the return from the subroutine.
  for K = 1:NDIM
    W(K) = ZERO;
    for J = 1:5
      W(K) = W(K) + WVEC(K, J) * PAR(J);
    end
  end
  VLAG(KOPT) = VLAG(KOPT) + ONE;

end