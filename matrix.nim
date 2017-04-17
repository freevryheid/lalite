import
  math,
  strutils,
  acomplex

type
  matrixshape* = tuple[rows: int, cols: int]
  matrix*[T] = object
    data*: seq[T]
    shape*: matrixshape
    order*: bool  # row [true] or col [false] major aligned
  MatrixError* = object of Exception

proc `[]=`*[T](m: var matrix[T], i, j: int, val: T) {. inline .} =
  if m.order:
    m.data[j * m.shape.rows + i] = val
  else:
    m.data[i * m.shape.cols + j] = val

proc `[]`*[T](m: matrix[T], i, j: int): T {. inline .} =
  if m.order:
    m.data[j * m.shape.rows + i]
  else:
    m.data[i * m.shape.cols + j]

proc toMatrix*[T](m: seq[T]; shape: matrixshape = (0,0)): matrix[T] =
  var mlen = m.len
  var t: float
  var rows, cols: int
  var k = 0
  if shape.rows < 0 or shape.cols < 0:
    raise newException(MatrixError, "cannot instantiate matrix - invalid shape dimensions")
  elif shape.rows == 0 and shape.cols == 0:
    var t = sqrt(mlen.float)
    if t.ceil != t:
      raise newException(MatrixError, "cannot instantiate square matrix - provide shape dimension")
    else:
      rows = t.int
      cols = t.int
  elif shape.rows == 0 and shape.cols > 0:
    t = mlen/shape.cols
    if t.ceil != t:
      raise newException(MatrixError, "cannot instantiate matrix - invalid shape dimensions")
    else:
      rows = t.int
      cols = shape.cols
  elif shape.rows > 0 and shape.cols == 0:
    t = mlen/shape.rows
    if t.ceil != t:
      raise newException(MatrixError, "cannot instantiate matrix - invalid shape dimensions")
    else:
      rows = shape.rows
      cols = t.int
  elif rows+cols != m.len:
    raise newException(MatrixError, "cannot instantiate matrix - invalid shape dimensions")
  else:
    rows = shape.rows
    cols = shape.cols

  result.data = @[]
  result.shape = (rows: rows, cols: cols)
  result.order = true
  for i in 0..<rows:
    for j in 0..<cols:
      result.data.add(m[k])
      k.inc

proc transpose*[T](m: matrix[T]): matrix[T] {.procvar.} =
  result.data = @[]
  var k: int
  for i in 0..<m.shape.cols:
    for j in 0..<m.shape.rows:
      k = j*m.shape.cols+i
      result.data.add(m.data[k])
  result.shape.rows = m.shape.cols
  result.shape.cols = m.shape.rows
  result.order = not m.order

proc eye*[T](dim: int, d: T): matrix[T] =
  result.data = @[]
  result.shape.rows = dim
  result.shape.cols = dim
  result.order = true
  for i in 0..<result.shape.rows:
    for j in 0..<result.shape.cols:
      if i == j:
        result.data.add(d)
      else:
        result.data.add(0)

proc fill*[T](shape: matrixshape, d: T): matrix[T] =
  result.data = @[]
  result.shape.rows = shape.rows
  result.shape.cols = shape.cols
  result.order = true
  for i in 0..<result.shape.rows:
    for j in 0..<result.shape.cols:
        result.data.add(d)

#HACK: strutils doesn't handle complex numbers
proc print*[T](m: matrix[T], decimals=2) =
  var k = 0
  var maxl = -1
  when T is float32 or T is float64:
    for d in m.data:
      k = ($(d.splitDecimal[0].int)).len
      if k > maxl: maxl = k
    k = 0
    for i in 0..<m.shape.rows:
      write stdout $"\n"
      for j in 0..<m.shape.cols:
        write stdout $formatFloat(m.data[k], ffDecimal, decimals).align(maxl+decimals+3)
        k.inc
    echo "\n"
  else:  # complex
    for d in m.data:
      var re = d[0]
      k = ($(re.splitDecimal[0].int)).len
      if k > maxl: maxl = k
    for d in m.data:
      var im = d[1]
      k = ($(im.splitDecimal[0].int)).len
      if k > maxl: maxl = k
    k = 0
    for i in 0..<m.shape.rows:
      write stdout $"\n"
      for j in 0..<m.shape.cols:
        write stdout $"["&formatFloat(m.data[k][0], ffDecimal, decimals).align(maxl+decimals+3), $formatFloat(m.data[k][1], ffDecimal, decimals).align(maxl+decimals+3)&']'
        k.inc
    echo "\n"
