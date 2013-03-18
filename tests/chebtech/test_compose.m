% Test file for chebtech/compose.

function pass = test_compose(pref)

if ( nargin < 1 )
    pref = chebtech.pref;
end

for ( n = 1:2 )
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    tol = 10*pref.chebtech.eps;

    % Compose a scalar-valued CHEBTECH object with sin(x):
    f = testclass.make(@(x) x);
    g = compose(f, @sin, pref);
    h = testclass.make(@sin);
    pass(n, 1) = norm(h.values - g.values, inf) < tol;

    % Compose a multi-valued CHEBTECH object with sin(x):
    f = testclass.make(@(x) [x x]);
    g = compose(f, @sin, pref);
    h = testclass.make(@(x) [sin(x) sin(x)]);
    pass(n, 2) = norm(h.values - g.values, inf) < tol;

    % Compose a multi-valued CHEBTECH object with sin(x):
    f = testclass.make(@(x) [x x.^2]);
    g = compose(f, @sin, pref);
    x = g.points();
    pass(n, 3) = norm(sin([x, x.^2])- g.values, inf) < tol;
    
    % Compose a multivalued CHEBTECH object with sin(x):
    f = testclass.make(@(x) [x x x.^2]);
    g = compose(f, @sin, pref);
    x = g.points();
    pass(n, 4) = norm(sin([x x x.^2]) - g.values, inf) < tol;
    
    % Compose 2 CHEBTECH objects with a binary function:
    f1 = testclass.make(@(x) sin(x));
    f2 = testclass.make(@(x) cos(x));
    g = compose(f1, @plus, f2, pref);
    h = testclass.make(@(x) sin(x) + cos(x));
    pass(n, 5) = norm(h.values - g.values, inf) < tol;
    
    % Compose 2 multivalued CHEBTECH objects with a binary function:
    f1 = testclass.make(@(x) [sin(x) cos(x)]);
    f2 = testclass.make(@(x) [cos(x) exp(x)]);
    g = compose(f1, @times, f2, pref);
    h = testclass.make(@(x) [sin(x).*cos(x) cos(x).*exp(x)]);
    pass(n, 6) = norm(h.values - g.values, inf) < tol;
    
    % Compose f(g), when f and g are CHEBTECH objects:
    f = testclass.make(@(x) x.^2);
    g = testclass.make(@(x) sin(x));
    h = compose(f, g);
    x = testclass.chebpts(length(h));
    pass(n, 7) = norm(h.values - sin(x.^2), inf) < tol;
    
    % Compose f(g), when f and g are CHEBTECH objects and g is multivalued:
    f = testclass.make(@(x) x.^2);
    g = testclass.make(@(x) [sin(x) cos(x)]);
    h = compose(f, g);
    x = testclass.chebpts(length(h));
    pass(n, 8) = norm(h.values - [sin(x.^2) cos(x.^2)], inf) < tol;
    
    % Compose f(g), when f and g are CHEBTECH objects and f is multivalued:
    f = testclass.make(@(x) [x x.^2]);
    g = testclass.make(@(x) sin(x));
    h = compose(f, g);
    x = testclass.chebpts(length(h));
    pass(n, 9) = norm(h.values - [sin(x) sin(x.^2)], inf) < tol;
    
    % We cannot expect to compose two multivalued CHEBTECH objects f(g):
    try 
        f = testclass.make(@(x) [x x.^2]);
        g = testclass.make(@(x) [sin(x), cos(x)]);
        compose(f, g);
        pass(n, 10) = false;
    catch ME 
        pass(n, 10) = strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:compose:multival');
    end
    
    % We cannot expect to compose two multivalued CHEBTECH objects in this way:
    try 
        f = testclass.make(@(x) [x x.^2]);
        g = testclass.make(@(x) sin(x));
        compose(f, @plus, g)
        pass(n, 11) = false;
    catch ME
        pass(n, 11) = strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:compose:dim');
    end
    
    try
        f = testclass.make(@(x) sin(x));
        g = testclass.make(@(x) 100*cos(x));
        compose(g,f);
        pass(n, 12) = false;
    catch ME
        pass(n, 12) = strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:compose:range');
    end
end

end