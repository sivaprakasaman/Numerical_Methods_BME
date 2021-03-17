function [iseed,u16807d]=u16807d( iseed )
%..........................................................................
%     bruce schmeiser     january 1989.                                   .
%     a linear congruential pseudorandom number generator                 .
%       using constant 16807 and modulus (2**31)-1.                       .
%              iseed = iseed*16807 (mod 2^31 -1)                          .
%     for implementations that don't require double precision, see        .
%       s.k. park and k.w. miller, "random numbers generators: good       .
%         ones are hard to find," cacm, 31, 10 (October 1988), 1192-1201. .
%     in correct implementations, starting with a seed value of 1 will    .
%       result in a seed value of 1043618065 on call number 10001.        .
%..........................................................................
%     input:  iseed.   integer.                                           . 
%                        chosen from [1,2147483646] on the first call.    .
%                        thereafter, the value returned from the last call.
%     output: iseed.   integer.                                           .
%                        to be used in the next call.                     .
%     output: u16807d.  real.                                             .
%                        a pseudorandom number in (0,1).                  .
%..........................................................................
%iseed=0;     
u16807d=0;
while (u16807d<=0 || u16807d>=1)
    iseed = mod (iseed * 16807,2147483647);
    u16807d = iseed / 2147483648;
end
%--------------------------------------------------------------------------