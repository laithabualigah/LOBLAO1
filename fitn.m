function [dis] = fitn(k,s,D)
[c] = ceniti(k,s,D);
dis=0;

for i=1:k
    for ii=1:size(s,2)
        if s(1,ii)==i
           for y=1:size(D,2)
             dis=dis+(abs(D(ii,y)-c(i,y)))^2;
           end
        end  
    end
end
 dis=sqrt(dis);
end




% [c] = ceniti(k,s,D);
% dis=0;
% 
% for i=1:k
%     for ii=1:size(s,2)
%         if s(1,ii)==i
%            for y=1:size(D,2)
%              dis=dis+(abs(D(ii,y)-c(i,y)))^2;
%             %dis=dis+abs(D(ii,y)-c(i,y));
%            end
%         end  
%     end
% end
%  dis=sqrt(dis);%*(k+1);