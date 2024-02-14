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
function [BMAT, ZMAT, IDZ, VLAG, W] = ...
    update (N, NPT, BMAT, ZMAT, IDZ, VLAG, BETA, KNEW, W)
  % IMPLICIT REAL(8) (A-H,O-Z)
  % DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)

  %     The arrays BMAT and ZMAT with IDZ are updated, in order to shift the
  %     interpolation point that has index KNEW. On entry, VLAG contains the
  %     components of the vector Theta*Wcheck+e_b of the updating formula
  %     (6.11), and BETA holds the value of the parameter that has this name.
  %     The vector W is used for working space.

  %     Set some constants.
  ONE = 1.0e0;
  ZERO = 0.0e0;
  NPTM = NPT - N - 1;
  %     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
  JL = 1;
  for J = 2:NPTM
    if (J == IDZ)
      JL = IDZ;
    elseif (ZMAT(KNEW, J) ~= ZERO)
      TEMP = sqrt(ZMAT(KNEW, JL)^2 + ZMAT(KNEW, J)^2);
      TEMPA = ZMAT(KNEW, JL) / TEMP;
      TEMPB = ZMAT(KNEW, J) / TEMP;
      for I = 1:NPT
        TEMP = TEMPA * ZMAT(I, JL) + TEMPB * ZMAT(I, J);
        ZMAT(I, J) = TEMPA * ZMAT(I, J) - TEMPB * ZMAT(I, JL);
        ZMAT(I, JL) = TEMP;
      end
      ZMAT(KNEW, J) = ZERO;
    end
  end

  %     Put the first NPT components of the KNEW-th column of HLAG into W,
  %     and calculate the parameters of the updating formula.

  TEMPA = ZMAT(KNEW, 1);
  if (IDZ >= 2)
    TEMPA = -TEMPA;
  end
  if (JL > 1)
    TEMPB = ZMAT(KNEW, JL);
  end
  for I = 1:NPT
    W(I) = TEMPA * ZMAT(I, 1);
    if (JL > 1)
      W(I) = W(I) + TEMPB * ZMAT(I, JL);
    end
  end
  ALPHA = W(KNEW);
  TAU = VLAG(KNEW);
  TAUSQ = TAU * TAU;
  DENOM = ALPHA * BETA + TAUSQ;
  VLAG(KNEW) = VLAG(KNEW) - ONE;
  %     Complete the updating of ZMAT when there is only one nonzero element
  %     in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to one,
  %     then the first column of ZMAT will be exchanged with another one later.
  IFLAG = 0;
  if (JL == 1)
    TEMP = sqrt(abs(DENOM));
    TEMPB = TEMPA / TEMP;
    TEMPA = TAU / TEMP;
    for I = 1:NPT
      ZMAT(I, 1) = TEMPA * ZMAT(I, 1) - TEMPB * VLAG(I);
    end
    if (IDZ == 1 && TEMP < ZERO)
      IDZ = 2;
    end
    if (IDZ >= 2 && TEMP >= ZERO)
      IFLAG = 1;
    end
  else
    %     Complete the updating of ZMAT in the alternative case.
    JA = 1;
    if (BETA >= ZERO)
      JA = JL;
    end
    JB = JL + 1 - JA;
    TEMP = ZMAT(KNEW, JB) / DENOM;
    TEMPA = TEMP * BETA;
    TEMPB = TEMP * TAU;
    TEMP = ZMAT(KNEW, JA);
    SCALA = ONE / sqrt(abs(BETA) * TEMP * TEMP + TAUSQ);
    SCALB = SCALA * sqrt(abs(DENOM));
    for I = 1:NPT
      ZMAT(I, JA) = SCALA * (TAU * ZMAT(I, JA) - TEMP * VLAG(I));
      ZMAT(I, JB) = SCALB * (ZMAT(I, JB) - TEMPA * W(I) - TEMPB * VLAG(I));
    end
    if (DENOM <= ZERO)
      if (BETA < ZERO)
        IDZ = IDZ + 1;
      else
        IFLAG = 1;
      end
    end
  end

  %     IDZ is reduced in the following case, and usually the first column
  %     of ZMAT is exchanged with a later one.

  if (IFLAG == 1)
    IDZ = IDZ - 1;
    for I = 1:NPT
      TEMP = ZMAT(I, 1);
      ZMAT(I, 1) = ZMAT(I, IDZ);
      ZMAT(I, IDZ) = TEMP;
    end
  end

  %     Finally, update the matrix BMAT.

  for J = 1:N
    JP = NPT + J;
    W(JP) = BMAT(KNEW, J);
    TEMPA = (ALPHA * VLAG(JP) - TAU * W(JP)) / DENOM;
    TEMPB = (-BETA * W(JP) - TAU * VLAG(JP)) / DENOM;
    for I = 1:JP
      BMAT(I, J) = BMAT(I, J) + TEMPA * VLAG(I) + TEMPB * W(I);
      if (I > NPT)
        BMAT(JP, I - NPT) = BMAT(I, J);
      end
    end
  end

end
