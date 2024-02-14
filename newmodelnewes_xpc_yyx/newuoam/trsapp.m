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
function [STEP, D, G, HD, HS, CRVMIN] = ...
    trsapp (N, NPT, XOPT, XPT, GQ, HQ, PQ, DELTA, STEP, ...
    D, G, HD, HS, CRVMIN)
  % IMPLICIT REAL(8) (A-H,O-Z)
  % DIMENSION XOPT(*),XPT(NPT,*),GQ(*),HQ(*),PQ(*),STEP(*), &
  %    D(*),G(*),HD(*),HS(*)

  %     N is the number of variables of a quadratic objective function, Q say.
  %     The arguments NPT, XOPT, XPT, GQ, HQ and PQ have their usual meanings,
  %       in order to define the current quadratic model Q.
  %     DELTA is the trust region radius, and has to be positive.
  %     STEP will be set to the calculated trial step.
  %     The arrays D, G, HD and HS will be used for working space.
  %     CRVMIN will be set to the least curvature of H along the conjugate
  %       directions that occur, except that it is set to zero if STEP goes
  %       all the way to the trust region boundary.

  %     The calculation of STEP begins with the truncated conjugate gradient
  %     method. If the boundary of the trust region is reached, then further
  %     changes to STEP may be made, each one being in the 2D space spanned
  %     by the current STEP and the corresponding gradient of Q. Thus STEP
  %     should provide a substantial reduction to Q within the trust region.

  %     Initialization, which includes setting HD to H times XOPT.

  HALF = 0.5e0;
  ZERO = 0.0e0;
  TWOPI = 8.0e0 * atan(1.0e0);
  DELSQ = DELTA * DELTA;
  ITERC = 0;
  ITERMAX = N;
  ITERSW = ITERMAX;

  for I = 1:N
    D(I) = XOPT(I);
  end

  %     The following instructions act as a subroutine for setting the vector
  %     HD to the vector D multiplied by the second derivative matrix of Q.
  %     They are called from three different places, which are distinguished
  %     by the value of ITERC.
  while (1)
    for I = 1:N
      HD(I) = ZERO;
    end

    for K = 1:NPT
      TEMP = ZERO;
      for J = 1:N
        TEMP = TEMP + XPT(K, J) * D(J);
      end
      TEMP = TEMP * PQ(K);
      for I = 1:N
        HD(I) = HD(I) + TEMP * XPT(K, I);
      end
    end
    IH = 0;
    for J = 1:N
      for I = 1:J
        IH = IH + 1;
        if (I < J)
          HD(J) = HD(J) + HQ(IH) * D(I);
        end
        HD(I) = HD(I) + HQ(IH) * D(J);
      end
    end

    if (ITERC == 0)
      %     Prepare for the first line search.
      QRED = ZERO;
      DD = ZERO;
      for I = 1:N
        STEP(I) = ZERO;
        HS(I) = ZERO;
        G(I) = GQ(I) + HD(I);
        D(I) = -G(I);
        DD = DD + D(I)^2;
      end
      CRVMIN = ZERO;
      if (DD == ZERO)
        return
      end
      DS = ZERO;
      SS = ZERO;
      GG = DD;
      GGBEG = GG;
      %     Calculate the step to the trust region boundary and the product HD.
      ITERC = ITERC + 1;
      TEMP = DELSQ - SS;
      BSTEP = TEMP / (DS + sqrt(DS * DS + DD * TEMP));
      continue
    end

    if (ITERC <= ITERSW)
      DHD = ZERO;
      for J = 1:N
        DHD = DHD + D(J) * HD(J);
      end
      %     Update CRVMIN and set the step-length ALPHA.
      ALPHA = BSTEP;
      if (DHD > ZERO)
        TEMP = DHD / DD;
        if (ITERC == 1)
          CRVMIN = TEMP;
        end
        CRVMIN = min(CRVMIN, TEMP);
        ALPHA = min(ALPHA, GG / DHD);
      end
      QADD = ALPHA * (GG - HALF * ALPHA * DHD);
      QRED = QRED + QADD;
      %     Update STEP and HS.
      GGSAV = GG;
      GG = ZERO;
      for I = 1:N
        STEP(I) = STEP(I) + ALPHA * D(I);
        HS(I) = HS(I) + ALPHA * HD(I);
        GG = GG + (G(I) + HS(I))^2;
      end
      %     Begin another conjugate direction iteration if required.
      if (ALPHA < BSTEP)
        if (QADD <= 0.01e0 * QRED)
          return;
        end
        if (GG <= 1.0e-4 * GGBEG)
          return;
        end
        if (ITERC == ITERMAX)
          return;
        end
        TEMP = GG / GGSAV;
        DD = ZERO;
        DS = ZERO;
        SS = ZERO;
        for I = 1:N
          D(I) = TEMP * D(I) - G(I) - HS(I);
          DD = DD + D(I)^2;
          DS = DS + D(I) * STEP(I);
          SS = SS + STEP(I)^2;
        end
        if (DS <= ZERO)
          return
        end
        if (SS < DELSQ)
          ITERC = ITERC + 1;
          TEMP = DELSQ - SS;
          BSTEP = TEMP / (DS + sqrt(DS * DS + DD * TEMP));
          continue
        end
      end
      CRVMIN = ZERO;
      ITERSW = ITERC;
      %     Test whether an alternative iteration is required.
      if (GG <= 1.0e-4 * GGBEG)
        return
      end
      SG = ZERO;
      SHS = ZERO;
      for I = 1:N
        SG = SG + STEP(I) * G(I);
        SHS = SHS + STEP(I) * HS(I);
      end
      SGK = SG + SHS;
      ANGTEST = SGK / sqrt(GG * DELSQ);
      if (ANGTEST <= -0.99e0)
        return
      end
      %     Begin the alternative iteration by calculating D and HD and some
      %     scalar products.
      ITERC = ITERC + 1;
      TEMP = sqrt(DELSQ * GG - SGK * SGK);
      TEMPA = DELSQ / TEMP;
      TEMPB = SGK / TEMP;
      for I = 1:N
        D(I) = TEMPA * (G(I) + HS(I)) - TEMPB * STEP(I);
      end
      continue
    end
    DG = ZERO;
    DHD = ZERO;
    DHS = ZERO;
    for I = 1:N
      DG = DG + D(I) * G(I);
      DHD = DHD + HD(I) * D(I);
      DHS = DHS + HD(I) * STEP(I);
    end
    %     Seek the value of the angle that minimizes Q.
    CF = HALF * (SHS - DHD);
    QBEG = SG + CF;
    QSAV = QBEG;
    QMIN = QBEG;
    ISAVE = 0;
    IU = 49;
    TEMP = TWOPI / (IU + 1);
    for I = 1:IU
      ANGLE = (I) * TEMP;
      CTH = cos(ANGLE);
      STH = sin(ANGLE);
      QNEW = (SG + CF * CTH) * CTH + (DG + DHS * CTH) * STH;
      if (QNEW < QMIN)
        QMIN = QNEW;
        ISAVE = I;
        TEMPA = QSAV;
      elseif (I == ISAVE + 1)
        TEMPB = QNEW;
      end
      QSAV = QNEW;
    end
    if (ISAVE == ZERO)
      TEMPA = QNEW;
    end
    if (ISAVE == IU)
      TEMPB = QBEG;
    end
    ANGLE = ZERO;
    if (TEMPA ~= TEMPB)
      TEMPA = TEMPA - QMIN;
      TEMPB = TEMPB - QMIN;
      ANGLE = HALF * (TEMPA - TEMPB) / (TEMPA + TEMPB);
    end
    ANGLE = TEMP * ((ISAVE) + ANGLE);
    %     Calculate the new STEP and HS. Then test for convergence.
    CTH = cos(ANGLE);
    STH = sin(ANGLE);
    REDUC = QBEG - (SG + CF * CTH) * CTH - (DG + DHS * CTH) * STH;
    GG = ZERO;
    for I = 1:N
      STEP(I) = CTH * STEP(I) + STH * D(I);
      HS(I) = CTH * HS(I) + STH * HD(I);
      GG = GG + (G(I) + HS(I))^2;
    end
    QRED = QRED + REDUC;
    RATIO = REDUC / QRED;
    if (ITERC < ITERMAX && RATIO > 0.01e0)
      if (GG <= 1.0e-4 * GGBEG)
        return
      end
      SG = ZERO;
      SHS = ZERO;
      for I = 1:N
        SG = SG + STEP(I) * G(I);
        SHS = SHS + STEP(I) * HS(I);
      end
      SGK = SG + SHS;
      ANGTEST = SGK / sqrt(GG * DELSQ);
      if (ANGTEST <= -0.99e0)
        return
      end
      %     Begin the alternative iteration by calculating D and HD and some
      %     scalar products.
      ITERC = ITERC + 1;
      TEMP = sqrt(DELSQ * GG - SGK * SGK);
      TEMPA = DELSQ / TEMP;
      TEMPB = SGK / TEMP;
      for I = 1:N
        D(I) = TEMPA * (G(I) + HS(I)) - TEMPB * STEP(I);
      end
    else
      return
    end
  end
end
