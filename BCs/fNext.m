function f = fNext(x,y,n)

% EX 1
grad = [2*(x-0.5),2*(y-0.5)];
f = dot(grad,n);

% EX 2
%n = -[2*(x-0.5),2*(y-0.5)];
%n = n/norm(n);

% f = dot([2,1],n);

end