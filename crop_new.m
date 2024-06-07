function res = crop_new(x,s)
%  res = crop(x,sx,sy)
%  crops a 2D matrix around its center.
%
%
%  res = crop(x,sx,sy,sz,st)
%  crops a 4D matrix around its center
%
%  
%  res = crop(x,[sx,sy,sz,st])
%  same as the previous example
%
%
%
%
% % (c) Michael Lustig 2007
% 
% if nargin < 2
% 	error('must have a target size')
% end
% 
% if nargin == 2
% 	s = sx;
% end
% 
% if nargin == 3
%     s = [sx,sy];
% end
% 
% if nargin == 4
%     s = [sx,sy,sz];
% end
% 
% if nargin == 5
%     s = [sx,sy,sz,st];
% end
% 

m=size(x);
if ndims(x)==2
    sx = s(1); sy = s(2);
    res = x(floor((size(x,1)-sx)/2)+1:1:floor((size(x,1)-sx)/2)+sx,floor((size(x,2)-sy)/2)+1:1:floor((size(x,2)-sy)/2)+sy);
end

if ndims(x)==3
    sx = s(1); sy = s(2); sz = s(3);
    res = x(floor((size(x,1)-sx)/2)+1:1:floor((size(x,1)-sx)/2)+sx,floor((size(x,2)-sy)/2)+1:1:floor((size(x,2)-sy)/2)+sy,floor((size(x,3)-sz)/2)+1:1:floor((size(x,3)-sz)/2)+sz);
end

if ndims(x)==4
    sx = s(1); sy = s(2); sz = s(3); st = s(4);
    res = x(floor((size(x,1)-sx)/2)+1:1:floor((size(x,1)-sx)/2)+sx,floor((size(x,2)-sy)/2)+1:1:floor((size(x,2)-sy)/2)+sy,floor((size(x,3)-sz)/2)+1:1:floor((size(x,3)-sz)/2)+sz,floor((size(x,4)-st)/2)+1:1:floor((size(x,4)-st)/2)+st);
end
    % m = size(x);
    % if length(s) < length(m)
	%     s = [s, ones(1,length(m)-length(s))];
    % end
    % 
    % % if sum(m==s)==length(m)
    % if sum(m==[s(1),s(2)])==length(m)
	% res = x;
	% return;
    % end
    % 
    % 
    % for n=1:length(s)
	%     idx{n} = floor(m(n)/2)+1+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2);
    % end
    % 
    % % this is a dirty ugly trick
    % cmd = 'res = x(idx{1}';
    % for n=2:length(s)
    % 	cmd = sprintf('%s,idx{%d}',cmd,n);
    % end
    % cmd = sprintf('%s);',cmd);
    % eval(cmd);
