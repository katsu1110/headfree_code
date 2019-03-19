function [newvec] = HiPaFi(oldvec, take, varargin)
% high-pass filter by sliding window averaging

if nargin < 2
    take = 10;
end
newvec = oldvec;

%     for i = 1:length(oldvec)
%         if i >= take+1 && i <= length(oldvec) - take;
%             vec = oldvec(i-take:i+take);    
%         elseif i<take+1
%             vec = oldvec(1:2*take+1);
%         elseif i>length(oldvec)-take
%             vec = oldvec(end-2*take:end);
%         else newvec(i) = oldvec(i);
%         end
%         vec(find(vec==oldvec(i))) = [];
%         newvec(i) = (sum(vec))/2*take; 
%     end

%     for i = 1:length(oldvec)
%         if i >= take+1 && i <= length(oldvec) - take;
%             newvec(i) = (sum([oldvec(i-take:i+take)])-oldvec(i))/(2*take);    
%         else newvec(i) = oldvec(i);
%         end
%     end
    
    for i = 1:length(oldvec)
        if i >= take+1 && i <= length(oldvec) - take
            newvec(i) = (sum(oldvec(i-take:i+take))-oldvec(i))/(2*take);    
        elseif i<take+1
            newvec(i) = (sum(oldvec(1:2*take+1))-oldvec(i))/(2*take);    
        elseif i>length(oldvec)-take
            newvec(i) = (sum(oldvec(end-2*take:end))-oldvec(i))/(2*take);    
        else newvec(i) = oldvec(i);
        end
    end