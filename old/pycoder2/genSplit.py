import datetime
import sys
import numpy as np
import inspect
import logging as lg

try:
    import cog
    use_cog = True
except ImportError:
    import sys
    use_cog = False

lg.basicConfig(format='%(levelname)-8s[%(funcName)-25s] %(message)s', level=lg.DEBUG)

# Load problem generation data
from genData import dat, args, constants, settings

########################################################
# Write out a dynamic header
########################################################
def gen_header():
  lg.info('Generating header')
  prtl("* Automatically generated on %s" % datetime.datetime.now())

########################################################
# Writes data to a static c-vector ROWWISE
########################################################
def load_data(var, a=None, useDefine=False, defName=None):
  lg.info('Loading data ' + var)
  if a is None:
    a = dat[var]
  if defName is not None:
    var = defName # Use a different name from var to define this variable
  if isinstance(a, np.ndarray):
    if a.size == 0:
      # The variable is empty, so it will never be used anyway
      prt('static const double ' + var + '[1] = {0.0') 
    else:
      prt('static const double ' + var + '[%i] = {' % (a.shape[0]*a.shape[1]))
      sep = ''
      for val in a:
        for x in val:
          prt(sep + '%.10g' % x)
          sep = ','
    prtl('};')
  elif isinstance(a, str):
    prtl('static const char %s[%i] = "%s";' % (var, len(a), a))
  elif isinstance(a, int):
    if useDefine:
      prtl('#define %s %g' % (var, a))
    else:
      prtl('static const int %s = %g;' % (var, a))
  elif isinstance(a, float):
    if useDefine:
      prtl('#define %s %g' % (var, a))
    else:
      prtl('static const double %s = %g;' % (var, a))
  else:
    raise NameError('Unknown type')

########################################################
# Loads all variables in 
########################################################
def load_constants():
  for key,val in constants.iteritems():
    load_data(key, a=val, useDefine=True)

########################################################
# Generate the calling structure of the solver function
########################################################
def gen_function_prototype():
  lg.info('Generating function prototype')
  prt("void %s(Sol *sol" % settings['functionName'])
  for key,val in args.iteritems():
    prt(", const double " + key + "[%i]" % val)
  prt(")")

########################################################
# Create a #define statement for the size of the given 
# var (which must be in dat)
########################################################
def define_size(var):
  lg.info('Defining variable size of ' + var)
  n,m = dat[var].shape
  prtl("#define %s_rows %i" % (var, n))
  prtl("#define %s_cols %i\n" % (var, m))

########################################################
# Generate code to multiply a constant matrix matName 
# from the dictionary with a vector vecName
########################################################
def matvec_product(matName, vecRHSName, vecLHSName, transpose=False):
  lg.info('Producing sparse code for ' + vecLHSName + '=' + matName + '*' + vecRHSName)
  M = dat[matName] if transpose==False else dat[matName].T
  density = float(np.count_nonzero(M)) / (M.shape[0]*M.shape[1])
  if density > 0.5:
    lg.warning('Matrix has %.2f%% non-zeros => Write dense multiplication code' % (density*100))
  for i,r in enumerate(M): # Iterate over rows of the matrix
    prt('%s[%i] = ' % (vecLHSName, i))
    nz = r.nonzero()
    if nz[0].size == 0:
      prt('0.0')
    else:
      prtNice = lambda x: '' if x == 1.0 else '%.10g*' % x
      [prt('+%s%s[%i]' % (prtNice(M[i,j]), vecRHSName, j)) for j in nz[0]]
    prtl(';')

########################################################
# Solve a system of linear equations x = A\b
# For now just inverts the matrix x = inv(A)*b
########################################################
# def mldivide(A,b,x): # A = np.array, b = str, x = str
#   inv(A)

########################################################
# Generate the proximal operators
#  gen_prox(x, xprox, rho)
########################################################
def gen_prox(x, xprox, rho):
  proxes = dat['prox']
  if isinstance(proxes, tuple):
    for prox in proxes:
      _gen_prox(prox, x, xprox, rho)
  else:
    _gen_prox(proxes, x, xprox, rho)

def _gen_prox(prox, x, xprox, rho):
  # Grab the function call structure
  proxMap = {
  'box':                 'notImplemented',
  'ellipse':             'notImplemented',
  'ellipseConj':         'notImplemented',
  'secondOrderCone':     'proj_secondOrderCone(%(xprox)s, %(x)s, %(len)s);',
  'secondOrderConeConj': 'notImplemented',
  'nonPositive':         'proj_negative(%(xprox)s, %(x)s, %(len)s);',
  'nonNegative':         'proj_positive(%(xprox)s, %(x)s, %(len)s);',
  'normBall':            'proj_normball_%(normtype)s(%(xprox)s, %(x)s, %(c)g, %(len)s);',
  'normProx':            'prox_norm_%(normtype)s(%(xprox)s, %(x)s, %(rho)s, %(len)s);'
  }
  funcCall = proxMap[prox['type']]

  # Define the arguments to the function
  argMap = {
    'len':      '%i' % prox['len'],            # Length of this prox vector
    'rho':      rho,
    'xprox':    '%s+%i' % (xprox,prox['ind']), # Add the offset into the main prox vector
    'x':        '%s+%i' % (x,prox['ind']),
  }
  if isinstance(prox['dat'], dict):
    dat = prox['dat']
    if 'p' in dat:
      argMap['normtype'] = {
        "Inf": 'inf',
        "1": 'one',
        "2": 'two'
      }[dat['p']]
    if 'c' in dat:
      argMap['c'] = dat['c']

  # Merge function call with its arguments
  funcCall = funcCall % argMap

  # Produce output
  prtl(funcCall)
  lg.info("Generating prox of type " + prox['type'] + ": " + funcCall)


########################################################
# Debug functions
########################################################

# Print without a newline
def prt(s):
  if use_cog:
    cog.out(s)
  else:
    sys.stdout.write(s)

# Print and add a newline
def prtl(s):
  if use_cog:
    cog.outl(s)
  else:
    print s

