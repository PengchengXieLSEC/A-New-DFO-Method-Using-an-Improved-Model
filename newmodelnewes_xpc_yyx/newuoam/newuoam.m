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
function [X, NF, F] = newuoam(Func, N, NPT, X, RHOBEG, RHOEND, IPRINT, MAXFUN)
  % implicit none
  % integer, intent(in) :: N, NPT, IPRINT, MAXFUN
  % real(kind=8), intent(inout) :: X(:), W(*)
  % real(kind=8), intent(in) :: RHOBEG,RHOEND
  % integer, intent(inout) :: NF
  % interface
  %    subroutine CALFUN(i_x, o_f,i_k,i_problem)
  %       real(kind=8), dimension(:) :: i_x
  %       real(kind=8) :: o_f
  %       integer(kind=4) :: i_k
  %       character(len=15) :: i_problem
  %    end subroutine
  % end interface
  % local variables:
  % integer :: NP, NPTM, NDIM, IXB, IXO, IXN, IXP, IFV, IGQ, &
  %    IHQ, IPQ, IBMAT, IZMAT, ID, IVL, IW
  % character(len=15) :: PROBLEM

  %%
  %     This subroutine seeks the least value of a function of many variables,
  %     by a trust region method that forms quadratic models by interpolation.
  %     There can be some freedom in the interpolation conditions, which is
  %     taken up by minimizing the Frobenius norm of the change to the second
  %     derivative of the quadratic model, beginning with a zero matrix. The
  %     arguments of the subroutine are as follows.

  %     N must be set to the number of variables and must be at least two.
  %     NPT is the number of interpolation conditions. Its value must be in the
  %       interval [N+2,(N+1)(N+2)/2].
  %     Initial values of the variables must be set in X(1),X(2),...,X(N). They
  %       will be changed to the values that give the least calculated F.
  %     RHOBEG and RHOEND must be set to the initial and final values of a trust
  %       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
  %       RHOBEG should be about one tenth of the greatest expected change to a
  %       variable, and RHOEND should indicate the accuracy that is required in
  %       the final values of the variables.
  %     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
  %       amount of printing. Specifically, there is no output if IPRINT=0 and
  %       there is output only at the return if IPRINT=1. Otherwise, each new
  %       value of RHO is printed, with the best vector of variables so far and
  %       the corresponding value of the objective function. Further, each new
  %       value of F with its variables are output if IPRINT=3.
  %     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
  %     The array W will be used for working space. Its length must be at least
  %     (NPT+13)*(NPT+N)+3*N*(N+3)/2.

  %     SUBROUTINE CALFUN (X,F) must be provided by the user. It must set F to
  %     the value of the objective function for the variables X(1),X(2),...,X(N).

  %     Partition the working space array, so that different parts of it can be
  %     treated separately by the subroutine that performs the main calculation.
  NP = N + 1;

  if (NPT < N + 2 || NPT > ((N + 2) * NP) / 2)
    error('return from NEWUOA because NPT is not in' + ...
    ' the required interval')
    return
  end

  %     The above settings provide a partition of W for subroutine NEWUOB.
  %     The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
  %     W plus the space that is needed by the last array of NEWUOB.
  [NF, X, F] = newuob(Func, N, NPT, X, RHOBEG, RHOEND, IPRINT, MAXFUN);
end