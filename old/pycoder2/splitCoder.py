from genData import dat, args, constants, settings
import numpy as np

class Gen:
  def sum(self, *args):
    """Weighted sum of vectors

    arg[i] = Var or Const or (Var, alpha) or (Const, alpha)

    If given with an alpha, then it's a weighted sum
    """
    # Extract arguments and weights
    f = lambda x: x[0] if isinstance(x, tuple) else x
    vars = [f(arg) for arg in args]

    f = lambda x: x[1] if isinstance(x, tuple) else 1.0
    weights = [f(arg) for arg in args]

    # Get the size of the sum
    shapes = [var.shape for var in vars]
    f = lambda i: reduce(lambda x,y: max(x[i],y[i]), shapes)
    shape = (f(0), f(1))

    if shape[1] > 1:
      raise Exception('Currently can only sum over vectors')

    # Generate the c-code
    cstr  = 'for(i=0; i<%i; i++)'
    cstr += 


gen = Gen()

class Var:
  def __init__(self, name, rows=1, cols=1, dtype=float):
    self.name  = name
    self.rows  = rows
    self.cols  = cols
    self.dtype = dtype
    self.shape = (rows,cols)


class Const:
  def __init__(self, name, x=None):
    if x is None:
      x = dat[name]
    self.name = name
    self.x    = x
    self.shape = x.shape


class Arg:
  def __init__(self, name, rows=1):
    self.name = name
    self.rows = rows

# Stores list of commands to be generated
class Prog:
  def __init__(self):
    self.code = []
    self.vars = []
    self.const = []

  def forloop(self, ):
    pass


# Helper functions
def concat(M):
  K = np.array([])
  for row in M:
    if K.size == 0:
      K = np.hstack(row)
    else:
      K = np.vstack([K,np.hstack(row)])
  return K

