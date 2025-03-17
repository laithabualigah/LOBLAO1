function [c] = ceniti(k,s,D)
  c=zeros(1,size(D,2));
for k1=1:k
    coun1=1; c1=zeros(1,size(D,2));
 for i=1:size(s,2)
     if s(1,i)==k1
         c1(coun1,:)=D(i,:);
         coun1=coun1+1;
     end     
 end
 
c(k1,:)=mean(c1);
end
 
end

