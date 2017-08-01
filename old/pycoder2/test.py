import splitCoder as split
import numpy as np
from   numpy.linalg import inv

import compiler

x = split.Var('x', 3)
y = split.Const('y', np.array([[4],[5],[6]]))

# vars = [x,y]
# shapes = [var.shape for var in vars]

# f = lambda i: reduce(lambda x,y: max(x[i],y[i]), shapes)
# print (f(0), f(1))

# print shapes

# shapes = [(2,1),(1,4),(10,-3)]
# print reduce(lambda x, y: (max(x[0],y[0]),max(x[1],y[1])), shapes)


split.gen.sum((x, -1.0),y,5.6)

# str = 'invK * (rho * (L.T * sum([-l, y, -lam]) - f))'
# expr = compiler.parse(str,'eval')

# print expr
