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
function [ALPHA, D, HCOL, GC, GD, S, W] = ...
    biglag (N, NPT, XOPT, XPT, BMAT, ZMAT, IDZ, KNEW, ...
    DELTA, D, HCOL, GC, GD, S, W)
  % IMPLICIT REAL(8) (A-H,O-Z)
  % DIMENSION XOPT(*),XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),D(*), &
  %    HCOL(*),GC(*),GD(*),S(*),W(*)

  %     N is the number of variables.
  %     NPT is the number of interpolation equations.
  %     XOPT is the best interpolation point so far.
  %     XPT contains the coordinates of the current interpolation points.
  %     BMAT provides the last N columns of H.
  %     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
  %     NDIM is the first dimension of BMAT and has the value NPT+N.
  %     KNEW is the index of the interpolation point that is going to be moved.
  %     DELTA is the current trust region bound.
  %     D will be set to the step from XOPT to the new point.
  %     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
  %     HCOL, GC, GD, S and W will be used for working space.

  %     The step D is calculated in a way that attempts to maximize the modulus
  %     of LFUNC(XOPT+D), subject to the bound ||D|| <= DELTA, where LFUNC is
  %     the KNEW-th Lagrange function.

  %     Set some constants.

  HALF = 0.5e0;
  ONE = 1.0e0;
  ZERO = 0.0e0;
  TWOPI = 8.0e0 * atan(ONE);
  DELSQ = DELTA * DELTA;
  NPTM = NPT - N - 1;

  %     Set the first NPT components of HCOL to the leading elements of the
  %     KNEW-th column of H.

  ITERC = 0;
  for K = 1:NPT
    HCOL(K) = ZERO;
  end
  for J = 1:NPTM
    TEMP = ZMAT(KNEW, J);
    if (J < IDZ)
      TEMP = -TEMP;
    end
    for K = 1:NPT
      HCOL(K) = HCOL(K) + TEMP * ZMAT(K, J);
    end
  end
  ALPHA = HCOL(KNEW);

  %     Set the unscaled initial direction D. Form the gradient of LFUNC at
  %     XOPT, and multiply D by the second derivative matrix of LFUNC.

  DD = ZERO;
  for I = 1:N
    D(I) = XPT(KNEW, I) - XOPT(I);
    GC(I) = BMAT(KNEW, I);
    GD(I) = ZERO;
    DD = DD + D(I)^2;
  end
  for K = 1:NPT
    TEMP = ZERO;
    SUM = ZERO;
    for J = 1:N
      TEMP = TEMP + XPT(K, J) * XOPT(J);
      SUM = SUM + XPT(K, J) * D(J);
    end
    TEMP = HCOL(K) * TEMP;
    SUM = HCOL(K) * SUM;
    for I = 1:N
      GC(I) = GC(I) + TEMP * XPT(K, I);
      GD(I) = GD(I) + SUM * XPT(K, I);
    end
  end

  %     Scale D and GD, with a sign change if required. Set S to another
  %     vector in the initial two dimensional subspace.

  GG = ZERO;
  SP = ZERO;
  DHD = ZERO;
  for I = 1:N
    GG = GG + GC(I)^2;
    SP = SP + D(I) * GC(I);
    DHD = DHD + D(I) * GD(I);
  end
  SCALE = DELTA / sqrt(DD);
  if (SP * DHD < ZERO)
    SCALE = -SCALE;
  end
  TEMP = ZERO;
  if (SP * SP > 0.99e0 * DD * GG)
    TEMP = ONE;
  end
  TAU = SCALE * (abs(SP) + HALF * SCALE * abs(DHD));
  if (GG * DELSQ < 0.01e0 * TAU * TAU)
    TEMP = ONE;
  end
  for I = 1:N
    D(I) = SCALE * D(I);
    GD(I) = SCALE * GD(I);
    S(I) = GC(I) + TEMP * GD(I);
  end
  %     Begin the iteration by overwriting S with a vector that has the
  %     required length and direction, except that termination occurs if
  %     the given D and S are nearly parallel.
  while (1)
    ITERC = ITERC + 1;
    DD = ZERO;
    SP = ZERO;
    SS = ZERO;
    for I = 1:N
      DD = DD + D(I)^2;
      SP = SP + D(I) * S(I);
      SS = SS + S(I)^2;
    end
    TEMP = DD * SS - SP * SP;
    if (TEMP <= 1.0e-8 * DD * SS)
      break
    end
    DENOM = sqrt(TEMP);
    for I = 1:N
      S(I) = (DD * S(I) - SP * D(I)) / DENOM;
      W(I) = ZERO;
    end
    %     Calculate the coefficients of the objective function on the circle,
    %     beginning with the multiplication of S by the second derivative matrix.
    for K = 1:NPT
      SUM = ZERO;
      for J = 1:N
        SUM = SUM + XPT(K, J) * S(J);
      end
      SUM = HCOL(K) * SUM;
      for I = 1:N
        W(I) = W(I) + SUM * XPT(K, I);
      end
    end
    CF1 = ZERO;
    CF2 = ZERO;
    CF3 = ZERO;
    CF4 = ZERO;
    CF5 = ZERO;
    for I = 1:N
      CF1 = CF1 + S(I) * W(I);
      CF2 = CF2 + D(I) * GC(I);
      CF3 = CF3 + S(I) * GC(I);
      CF4 = CF4 + D(I) * GD(I);
      CF5 = CF5 + S(I) * GD(I);
    end
    CF1 = HALF * CF1;
    CF4 = HALF * CF4 - CF1;
    %     Seek the value of the angle that maximizes the modulus of TAU.
    TAUBEG = CF1 + CF2 + CF4;
    TAUMAX = TAUBEG;
    TAUOLD = TAUBEG;
    ISAVE = 0;
    IU = 49;
    TEMP = TWOPI / (IU + 1);
    for I = 1:IU
      ANGLE = (I) * TEMP;
      CTH = cos(ANGLE);
      STH = sin(ANGLE);
      TAU = CF1 + (CF2 + CF4 * CTH) * CTH + (CF3 + CF5 * CTH) * STH;
      if (abs(TAU) > abs(TAUMAX))
        TAUMAX = TAU;
        ISAVE = I;
        TEMPA = TAUOLD;
      elseif (I == ISAVE + 1)
        TEMPB = TAU;
      end
      TAUOLD = TAU;
    end
    if (ISAVE == 0)
      TEMPA = TAU;
    end
    if (ISAVE == IU)
      TEMPB = TAUBEG;
    end
    STEP = ZERO;
    if (TEMPA ~= TEMPB)
      TEMPA = TEMPA - TAUMAX;
      TEMPB = TEMPB - TAUMAX;
      STEP = HALF * (TEMPA - TEMPB) / (TEMPA + TEMPB);
    end
    ANGLE = TEMP * ((ISAVE) + STEP);
    %     Calculate the new D and GD. Then test for convergence.
    CTH = cos(ANGLE);
    STH = sin(ANGLE);
    TAU = CF1 + (CF2 + CF4 * CTH) * CTH + (CF3 + CF5 * CTH) * STH;
    for I = 1:N
      D(I) = CTH * D(I) + STH * S(I);
      GD(I) = CTH * GD(I) + STH * W(I);
      S(I) = GC(I) + GD(I);
    end
    if (abs(TAU) <= 1.1e0 * abs(TAUBEG))
      break
    elseif (ITERC < N)
      continue
    else
      break
    end
  end

end
