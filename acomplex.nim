#
#
#            Nim's Runtime Library
#        (c) Copyright 2010 Andreas Rumpf
#
#    See the file "copying.txt", included in this
#    distribution, for details about the copyright.
#

## This module implements complex numbers using arrays.
{.push checks:off, line_dir:off, stack_trace:off, debugger:off.}
# the user does not want to trace a part
# of the standard library!

import
  math

const
  EPS = 1.0e-7  ## Epsilon used for float comparisons.

type
  complex*[T] = array[2, T]
  complex32* = complex[float32]
  complex64* = complex[float64]
    ## a complex number, consisting of a real and an imaginary part

proc tocomplex*[T](x: T): complex[T] =
  ## Convert ``x`` to a complex number.
  result[0] = x
  result[1] = 0

proc `==` *(x, y: complex): bool =
  ## Compare two complex numbers `x` and `y` for equality.
  result = x[0] == y[0] and x[1] == y[1]

proc `=~` *(x, y: complex): bool =
  ## Compare two complex numbers `x` and `y` approximately.
  result = abs(x[0]-y[0])<EPS and abs(x[1]-y[1])<EPS

proc `+` *(x, y: complex): complex =
  ## Add two complex numbers.
  result[0] = x[0] + y[0]
  result[1] = x[1] + y[1]

proc `+` *[T](x: complex[T], y: T): complex[T] =
  ## Add complex `x` to `y`.
  result[0] = x[0] + y
  result[1] = x[1]

proc `+` *[T](x: T, y: complex[T]): complex[T] =
  ## Add `x` to complex `y`.
  result[0] = x + y[0]
  result[1] = y[1]

proc `-` *(z: complex): complex =
  ## Unary minus for complex numbers.
  result[0] = -z[0]
  result[1] = -z[1]

proc `-` *(x, y: complex): complex =
  ## Subtract two complex numbers.
  result[0] = x[0] - y[0]
  result[1] = x[1] - y[1]

proc `-` *[T](x: complex[T], y: T): complex[T] =
  ## Subtracts `y` from complex `x`.
  result = x + (-y)

proc `-` *[T](x: T, y: complex[T]): complex[T] =
  ## Subtracts complex `y` from `x`.
  result = x + (-y)

proc `/` *(x, y: complex): complex =
  ## Divide `x` by `y`.
  var
    r, den: float
  if abs(y[0]) < abs(y[1]):
    r = y[0] / y[1]
    den = y[1] + r * y[0]
    result[0] = (x[0] * r + x[1]) / den
    result[1] = (x[1] * r - x[0]) / den
  else:
    r = y[1] / y[0]
    den = y[0] + r * y[1]
    result[0] = (x[0] + r * x[1]) / den
    result[1] = (x[1] - r * x[0]) / den

proc `/` *[T](x : complex[T], y: T): complex[T] =
  ## Divide complex `x` by `y`.
  result[0] = x[0]/y
  result[1] = x[1]/y

proc `/` *[T](x : T, y: complex[T]): complex[T] =
  ## Divide `x` by complex `y`.
  var num : complex[T] = [x, 0.T]
  result = num/y

proc `*` *(x, y: complex): complex =
  ## Multiply `x` with `y`.
  result[0] = x[0] * y[0] - x[1] * y[1]
  result[1] = x[1] * y[0] + x[0] * y[1]

proc `*` *[T](x: T, y: complex[T]): complex[T] =
  ## Multiply `x` with complex `y`.
  result[0] = x * y[0]
  result[1] = x * y[1]

proc `*` *[T](x: complex[T], y: T): complex[T] =
  ## Multiply complex `x` with `y`.
  result[0] = x[0] * y
  result[1] = x[1] * y

proc `+=` *(x: var complex, y: complex) =
  ## Add `y` to `x`.
  x[0] += y[0]
  x[1] += y[1]

proc `+=` *[T](x: var complex[T], y: T) =
  ## Add `y` to the complex number `x`.
  x[0] += y

proc `-=` *(x: var complex, y: complex) =
  ## Subtract `y` from `x`.
  x[0] -= y[0]
  x[1] -= y[1]

proc `-=` *[T](x: var complex[T], y: T) =
  ## Subtract `y` from the complex number `x`.
  x[0] -= y

proc `*=` *(x: var complex, y: complex) =
  ## Multiply `y` to `x`.
  let im = x[1] * y[0] + x[0] * y[1]
  x[0] = x[0] * y[0] - x[1] * y[1]
  x[1] = im

proc `*=` *[T](x: var complex, y: T) =
  ## Multiply `y` to the complex number `x`.
  x[0] *= y
  x[1] *= y

proc `/=` *(x: var complex, y: complex) =
  ## Divide `x` by `y` in place.
  x = x / y

proc `/=` *[T](x : var complex[T], y: T) =
  ## Divide complex `x` by `y` in place.
  x[0] /= y
  x[1] /= y

proc abs*[T](z: complex[T]): T =
  ## Return the distance from (0,0) to `z`.
  # optimized by checking special cases (sqrt is expensive)
  var x, y, temp: T
  x = abs(z[0])
  y = abs(z[1])
  if x == 0.T:
    result = y
  elif y == 0.T:
    result = x
  elif x > y:
    temp = y / x
    result = x * sqrt(1 + temp * temp)
  else:
    temp = x / y
    result = y * sqrt(1 + temp * temp)

proc conjugate*(z: complex): complex =
  ## Conjugate of complex number `z`.
  result[0] = z[0]
  result[1] = -z[1]

proc sqrt*[T](z: complex[T]): complex[T] =
  ## Square root for a complex number `z`.
  var x, y, w, r: T
  if z[0] == 0.T and z[1] == 0.T:
    result = z
  else:
    x = abs(z[0])
    y = abs(z[1])
    if x >= y:
      r = y / x
      w = sqrt(x) * sqrt(0.5 * (1.0 + sqrt(1.0 + r * r)))
    else:
      r = x / y
      w = sqrt(y) * sqrt(0.5 * (r + sqrt(1.0 + r * r)))
    if z[0] >= 0.T:
      result[0] = w
      result[1] = z[1] / (w * 2.0)
    else:
      if z[1] >= 0.T:
        result[1] = w
      else:
        result[1] = -w
      result[0] = z[1] / (result[1] + result[1])

proc exp*(z: complex): complex =
  ## e raised to the power `z`.
  var rho   = exp(z[0])
  var theta = z[1]
  result[0] = rho*cos(theta)
  result[1] = rho*sin(theta)

proc ln*(z: complex): complex =
  ## Returns the natural log of `z`.
  result[0] = ln(abs(z))
  result[1] = arctan2(z[1],z[0])

proc log10*[T](z: complex[T]): complex[T] =
  ## Returns the log base 10 of `z`.
  result = ln(z)/ln(10.T)

proc log2*[T](z: complex[T]): complex[T] =
  ## Returns the log base 2 of `z`.
  result = ln(z)/ln(2.T)

proc pow*[T](x, y: complex[T]): complex[T] =
  ## `x` raised to the power `y`.
  if x[0] == 0.T  and  x[1] == 0.T:
    if y[0] == 0.T  and  y[1] == 0.T:
      result[0] = 1.T
      result[1] = 0.T
    else:
      result[0] = 0.T
      result[1] = 0.T
  elif y[0] == 1.T  and  y[1] == 0.T:
    result = x
  elif y[0] == -1.T  and  y[1] == 0.T:
    result = 1.T/x
  else:
    var rho   = sqrt(x[0]*x[0] + x[1]*x[1])
    var theta = arctan2(x[1],x[0])
    var s     = pow(rho,y[0]) * exp(-y[1]*theta)
    var r     = y[0]*theta + y[1]*ln(rho)
    result[0] = s*cos(r)
    result[1] = s*sin(r)

proc sin*(z: complex): complex =
  ## Returns the sine of `z`.
  result[0] = sin(z[0])*cosh(z[1])
  result[1] = cos(z[0])*sinh(z[1])

proc arcsin*[T](z: complex[T]): complex[T] =
  ## Returns the inverse sine of `z`.
  let i = [0.T,1.T]
  result = -i*ln(i*z + sqrt(1.T-z*z))

proc cos*(z: complex): complex =
  ## Returns the cosine of `z`.
  result[0] = cos(z[0])*cosh(z[1])
  result[1] = -sin(z[0])*sinh(z[1])

proc arccos*[T](z: complex[T]): complex[T] =
  ## Returns the inverse cosine of `z`.
  let i = [0.T,1.T]
  result = -i*ln(z + sqrt(z*z-1.T))

proc tan*(z: complex): complex =
  ## Returns the tangent of `z`.
  result = sin(z)/cos(z)

proc arctan*[T](z: complex[T]): complex[T] =
  ## Returns the inverse tangent of `z`.
  let i = [0.T,1.T]
  result = T(0.5)*i*(ln(1.T-i*z)-ln(1.T+i*z))

proc cot*(z: complex): complex =
  ## Returns the cotangent of `z`.
  result = cos(z)/sin(z)

proc arccot*[T](z: complex[T]): complex[T] =
  ## Returns the inverse cotangent of `z`.
  let i = (0.T,1.T)
  result = T(0.5)*i*(ln(1.T-i/z)-ln(1.T+i/z))

proc sec*[T](z: complex[T]): complex[T] =
  ## Returns the secant of `z`.
  result = 1.T/cos(z)

proc arcsec*[T](z: complex[T]): complex[T] =
  ## Returns the inverse secant of `z`.
  let i = (0.T,1.T)
  result = -i*ln(i*sqrt(1.T-1.T/(z*z))+1.T/z)

proc csc*[T](z: complex[T]): complex[T] =
  ## Returns the cosecant of `z`.
  result = 1.T/sin(z)

proc arccsc*[T](z: complex[T]): complex[T] =
  ## Returns the inverse cosecant of `z`.
  let i = (0.T,1.T)
  result = -i*ln(sqrt(1.T-1.T/(z*z))+i/z)

proc sinh*[T](z: complex[T]): complex[T] =
  ## Returns the hyperbolic sine of `z`.
  result = T(0.5)*(exp(z)-exp(-z))

proc arcsinh*[T](z: complex[T]): complex[T] =
  ## Returns the inverse hyperbolic sine of `z`.
  result = ln(z+sqrt(z*z+1.T))

proc cosh*[T](z: complex[T]): complex[T] =
  ## Returns the hyperbolic cosine of `z`.
  result = T(0.5)*(exp(z)+exp(-z))

proc arccosh*[T](z: complex[T]): complex[T] =
  ## Returns the inverse hyperbolic cosine of `z`.
  result = ln(z+sqrt(z*z-1.T))

proc tanh*(z: complex): complex =
  ## Returns the hyperbolic tangent of `z`.
  result = sinh(z)/cosh(z)

proc arctanh*[T](z: complex[T]): complex[T] =
  ## Returns the inverse hyperbolic tangent of `z`.
  result = T(0.5)*(ln((1.T+z)/(1.T-z)))

proc sech*[T](z: complex[T]): complex[T] =
  ## Returns the hyperbolic secant of `z`.
  result = 2.T/(exp(z)+exp(-z))

proc arcsech*[T](z: complex[T]): complex[T] =
  ## Returns the inverse hyperbolic secant of `z`.
  result = ln(1.T/z+sqrt(1.T/z+1.T)*sqrt(1.T/z-1.T))

proc csch*[T](z: complex[T]): complex[T] =
  ## Returns the hyperbolic cosecant of `z`.
  result = 2.T/(exp(z)-exp(-z))

proc arccsch*[T](z: complex[T]): complex[T] =
  ## Returns the inverse hyperbolic cosecant of `z`.
  result = ln(1.T/z+sqrt(1.T/(z*z)+1.T))

proc coth*(z: complex): complex =
  ## Returns the hyperbolic cotangent of `z`.
  result = cosh(z)/sinh(z)

proc arccoth*[T](z: complex[T]): complex[T] =
  ## Returns the inverse hyperbolic cotangent of `z`.
  result = T(0.5)*(ln(1.T+1.T/z)-ln(1.T-1.T/z))

proc phase*[T](z: complex[T]): T =
  ## Returns the phase of `z`.
  arctan2(z[1], z[0])

proc polar*[T](z: complex[T]): tuple[r, phi: T] =
  ## Returns `z` in polar coordinates.
  result.r = abs(z)
  result.phi = phase(z)

proc rect*[T](r: T, phi: T): complex[T] =
  ## Returns the complex number with polar coordinates `r` and `phi`.
  result[0] = r * cos(phi)
  result[1] = r * sin(phi)

proc `$`*(z: complex): string =
  ## Returns `z`'s string representation as ``"[re, im]"``.
  result = "[" & $z[0] & ", " & $z[1] & "]"

{.pop.}


when isMainModule:
  let
    z = [0.0, 0.0]
    oo = [1.0,1.0]
    a = [1.0, 2.0]
    b = [-1.0, -2.0]
    m1 = [-1.0, 0.0]
    i = [0.0,1.0]
    one = [1.0,0.0]
    tt = [10.0, 20.0]
    ipi = [0.0, -PI]
    t = polar(a)

  doAssert( a == a )
  doAssert( (a-a) == z )
  doAssert( (a+b) == z )
  doAssert( (a/b) == m1 )
  doAssert( (1.0/a) == [0.2, -0.4] )
  doAssert( (a*b) == [3.0, -4.0] )
  doAssert( 10.0*a == tt )
  doAssert( a*10.0 == tt )
  doAssert( tt/10.0 == a )
  doAssert( oo+(-1.0) == i )
  doAssert( (-1.0)+oo == i )
  doAssert( abs(oo) == sqrt(2.0) )
  doAssert( conjugate(a) == [1.0, -2.0] )
  doAssert( sqrt(m1) == i )
  doAssert( exp(ipi) =~ m1 )

  doAssert( pow(a,b) =~ [-3.72999124927876, -1.68815826725068] )
  doAssert( pow(z,a) =~ [0.0, 0.0] )
  doAssert( pow(z,z) =~ [1.0, 0.0] )
  doAssert( pow(a,one) =~ a )
  doAssert( pow(a,m1) =~ [0.2, -0.4] )

  doAssert( ln(a) =~ [0.804718956217050, 1.107148717794090] )
  doAssert( log10(a) =~ [0.349485002168009, 0.480828578784234] )
  doAssert( log2(a) =~ [1.16096404744368, 1.59727796468811] )

  doAssert( sin(a) =~ [3.16577851321617, 1.95960104142161] )
  doAssert( cos(a) =~ [2.03272300701967, -3.05189779915180] )
  doAssert( tan(a) =~ [0.0338128260798967, 1.0147936161466335] )
  doAssert( cot(a) =~ 1.0/tan(a) )
  doAssert( sec(a) =~ 1.0/cos(a) )
  doAssert( csc(a) =~ 1.0/sin(a) )
  doAssert( arcsin(a) =~ [0.427078586392476, 1.528570919480998] )
  doAssert( arccos(a) =~ [1.14371774040242, -1.52857091948100] )
  doAssert( arctan(a) =~ [1.338972522294494, 0.402359478108525] )

  doAssert( cosh(a) =~ [-0.642148124715520, 1.068607421382778] )
  doAssert( sinh(a) =~ [-0.489056259041294, 1.403119250622040] )
  doAssert( tanh(a) =~ [1.1667362572409199,-0.243458201185725] )
  doAssert( sech(a) =~ 1.0/cosh(a) )
  doAssert( csch(a) =~ 1.0/sinh(a) )
  doAssert( coth(a) =~ 1.0/tanh(a) )
  doAssert( arccosh(a) =~ [1.528570919480998, 1.14371774040242] )
  doAssert( arcsinh(a) =~ [1.469351744368185, 1.06344002357775] )
  doAssert( arctanh(a) =~ [0.173286795139986, 1.17809724509617] )
  doAssert( arcsech(a) =~ arccosh(1.0/a) )
  doAssert( arccsch(a) =~ arcsinh(1.0/a) )
  doAssert( arccoth(a) =~ arctanh(1.0/a) )

  doAssert( phase(a) == 1.1071487177940904 )
  doAssert( rect(t.r, t.phi) =~ a )
  doAssert( rect(1.0, 2.0) =~ [-0.4161468365471424, 0.9092974268256817] )