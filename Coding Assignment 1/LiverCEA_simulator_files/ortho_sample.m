%Andrew Sivaprakasam
%Calculate orthogonal samples by brute forcing LHS. No doubt could be
%optimized, but fast enough.

function [dist] = ortho_sample(n,divs)
%Brute force calculation of orthogonal sample set for 2D sampling
%MUST MAKE SURE n is evenly divided by divs!!!!

if(mod(n,divs^2)~=0)
    error('MUST MAKE SURE n is evenly divided by divs*divs!!!!');
end
satisfied = 0;

while(~satisfied)
   
   dist = lhsdesign(n,2);
   
   bin_size = 1/divs;
   subspaces = divs^2;
   
   for i = 1:divs
       lb_y = (i-1)*bin_size;
       ub_y = i*bin_size;
        
       for j = 1:divs
              lb_x = (j-1)*bin_size;
              ub_x = j*bin_size;
              sum_space(i,j) = sum((lb_y<dist(:,1) & dist(:,1)<ub_y)&(lb_x<dist(:,2) & dist(:,2)<ub_x));
       end
   end 
   
   if(prod(prod(sum_space == sum_space(1,1))))
        satisfied = 1;
   end
end

end

