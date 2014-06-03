function x = fourpts(n, dom)

if ( nargin == 1 )
    dom = [-1,1];
end

x = linspace(dom(1), dom(end), n+1).';
x(end) = [];

end