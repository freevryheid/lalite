import acomplex
import matrix
import lapack
import blas

let
  ι = [0.0, 1.0]

proc solve*[T](
  # solution to linear systems of equation with real/complex coefficients
  a: matrix[T];
  b: matrix[T]): matrix[T] {.inline.} =
  var
    ma = a
    mb = b
  if ma.order:
    ma = ma.transpose
  if mb.order:
    mb = mb.transpose
  var
    n = ma.shape.rows
    nrhs = mb.shape.rows
    lda = ma.shape.rows
    ipiv: seq[int] = @[]
    ldb = mb.shape.cols
    info: int
  for i in 0..<n:
    ipiv.add(0)
  gesv(n.addr, nrhs.addr, ma.data[0].addr, lda.addr, ipiv[0].addr, mb.data[0].addr, ldb.addr, info.addr)
  if info != 0:
    nfo(info, "gesv")
  else:
    result.order = false
    result = mb.transpose

proc inv*[T](
  # matrix inversion for real/complex square matrices
  a: matrix[T]): matrix[T] {.inline.} =
  var
    ma = a
  if ma.order:
    ma = ma.transpose
  var
    m = ma.shape.rows
    n = ma.shape.cols
    lda = n
    ipiv: seq[int] = @[]
    info: int
  for i in 0..<n:
    ipiv.add(0)
  var
    ispec = 1
    name = "DGETRI"
    opts = "UN"
    n1 = n
    n2 = -1
    n3 = -1
    n4 = -1
    nb = ilaenv(addr ispec, addr name, addr opts, addr n1, addr n2, addr n3, addr n4)  # TODO: check UN param
    lwork = n*nb
    work: seq[T] = @[]
  for i in 0..<lwork:
    work.add(0)
  if nb < 1: nb = max(1, n)
  dgetrf(m.addr, n.addr, ma.data[0].addr, lda.addr, ipiv[0].addr, info.addr)
  if info != 0:
    nfo(info, "dgetrf")
  else:
    dgetri(n.addr, ma.data[0].addr, lda.addr, ipiv[0].addr, work[0].addr, lwork.addr, info.addr)
    if info != 0:
      nfo(info, "dgetri")
    else:
      result.order = false
      result = ma.transpose

proc eigvals*[T](
  # eigenvalues for general matrices
  a: matrix[T]): matrix[complex[T]] {.inline.} =
  var
    ma = a
  if ma.order:
    ma = ma.transpose
  var
    n = ma.shape.rows
    lda = ma.shape.cols
    ldvl = n
    ldvr = n
    lwork = 8*n  #TODO: can this size be optimized? query first?
    vl = ma.data
    vr = ma.data
    work: seq[T] = @[]
    wr: seq[T] = @[]
    wi: seq[T] = @[]
    info: int
    jobvl = 'N'
    jobvr = 'N'
  for i in 0..<lwork:
    work.add(0)
  for i in 0..<n:
    wr.add(0)
    wi.add(0)
  geev(jobvl.addr, jobvr.addr, n.addr, ma.data[0].addr, lda.addr, wr[0].addr, wi[0].addr, vl[0].addr, ldvl.addr, vr[0].addr, ldvr.addr, work[0].addr, lwork.addr, info.addr)
  if info != 0:
    nfo(info, "geev")
  else:
    result.data = @[]
    for i in 0..<n:
      result.data.add(wr[i] + ι*wi[i])
    result.shape.rows = n
    result.shape.cols = 1
    result.order = false
    result = result.transpose

proc eig*[T](
  # eigen-value/-vector problem for general matrices
  a: matrix[T]): tuple[val: matrix[complex[T]], vec: matrix[complex[T]]] {.inline.} =
  # type
  #   retType = tuple[val: matrix[complex[T]], vec: matrix[complex[T]]]
  var
    ma = a
    # ret: retType
  if ma.order:
    ma = ma.transpose
  var
    n = ma.shape.rows
    lda = ma.shape.cols
    ldvl = n
    ldvr = n
    lwork = 8*n  #TODO: can this size be optimized? query first?
    vl = ma.data
    vr = ma.data
    work: seq[T] = @[]
    wr: seq[T] = @[]
    wi: seq[T] = @[]
    w: seq[complex[T]] = @[]  # used for eigenvector result
    info: int
    jobvl = 'N'
    jobvr = 'V'
  for i in 0..<lwork:
    work.add(0)
  for i in 0..<n:
    wr.add(0)
    wi.add(0)
  geev(jobvl.addr, jobvr.addr, n.addr, ma.data[0].addr, lda.addr, wr[0].addr, wi[0].addr, vl[0].addr, ldvl.addr, vr[0].addr, ldvr.addr, work[0].addr, lwork.addr, info.addr)
  if info != 0:
    nfo(info, "geev")
  else:
    # eigenvalues
    result.val.data = @[]
    for i in 0..<n:
      result.val.data.add(wr[i] + ι*wi[i])
    result.val.shape.rows = n
    result.val.shape.cols = 1
    result.val.order = false
    result.val = result.val.transpose
    # right eigenvectors
    var v = vr.toMatrix
    for j in 0..<n:
      if wi[j] > 0.0:  # first of two conjugate eigenvalues
        for i in 0..<n:
          w.add(v[i,j] + ι*v[i, j+1])
      elif wi[j] < 0.0:  # second of two conjugate eigenvalues
        for i in 0..<n:
          w.add(v[i,j-1] - ι*v[i, j])
      else:
        for i in 0..<n:
          w.add([v[i,j], 0.0])
    result.vec = w.toMatrix
    result.vec.order = false
    result.vec = result.vec.transpose

proc matmul*[T](a,b: matrix[T]): matrix[T] =
  var
    ma = a
    mb = b
    mc = fill((a.shape.rows,b.shape.cols), 0.0)
  if ma.order:
    ma = ma.transpose
  if mb.order:
    mb = mb.transpose
  var
    transa = 'N'
    transb = 'N'
    m = a.shape.rows
    n = b.shape.cols
    k = a.shape.cols
    alpha = 1.0
    lda = a.shape.rows
    ldb = k
    beta = 0.0
    ldc = m
  gemm(transa.addr, transb.addr, m.addr, n.addr, k.addr, alpha.addr, ma.data[0].addr, lda.addr, mb.data[0].addr, ldb.addr, beta.addr, mc.data[0].addr, ldc.addr)
  mc.order = not mc.order
  return mc.transpose

proc abs*[T](a: matrix[T]): T =
  result = 0
  for i in a.shape.rows:
    for j in a.shape.cols:
      result += a[i,j][0]*a[i,j][0] + a[i,j][1]*a[i,j][1]
  result = sqrt(result)

when isMainModule:


  let a = @[-1.01, 0.86,-4.60, 3.31,-4.81,
             3.98, 0.53,-7.04, 5.29, 3.55,
             3.30, 8.26,-3.89, 8.20,-1.51,
             4.43, 4.96,-7.66,-7.33, 6.18,
             7.31,-6.43,-6.16, 2.47, 5.58]
  let ma = a.toMatrix
  let z = ma.eig

  "\nmatrix A:".echo
  ma.print
  "eigenvectors of matrix a:".echo
  z.vec.print

