when defined(windows):
  const lapackSuffix = ".dll"
elif defined(macosx):
  const lapackSuffix = ".dylib"
else:
  const lapackSuffix = ".so"

const lapackPrefix = "liblapack"
const lapackName = lapackPrefix & lapackSuffix

import acomplex

proc gesv*(
  n,nrhs: ptr int,
  a: ptr float32,
  lda,ipiv: ptr int,
  b: ptr float32,
  ldb,info: ptr int,
) {.cdecl, importc: "sgesv_", dynlib: lapackName.}

proc gesv*(
  n,nrhs: ptr int,
  a: ptr float64,
  lda,ipiv: ptr int,
  b: ptr float64,
  ldb,info: ptr int,
) {.cdecl, importc: "dgesv_", dynlib: lapackName.}

proc gesv*(
  n,nrhs: ptr int,
  a: ptr complex32,
  lda,ipiv: ptr int,
  b: ptr complex32,
  ldb,info: ptr int,
) {.cdecl, importc: "cgesv_", dynlib: lapackName.}

proc gesv*(
  n,nrhs: ptr int,
  a: ptr complex,
  lda,ipiv: ptr int,
  b: ptr complex,
  ldb,info: ptr int,
) {.cdecl, importc: "zgesv_", dynlib: lapackName.}

proc dgetrf*(
  m,n: ptr int,
  a: ptr float64,
  lda,ipiv, info: ptr int,
) {.cdecl, importc: "dgetrf_", dynlib: lapackName.}

proc dgetri*(
  n: ptr int,
  a: ptr float64,
  lda,ipiv: ptr int,
  work: ptr float64,
  lwork,info: ptr int,
) {.cdecl, importc: "dgetri_", dynlib: lapackName.}

proc geev*(
  jobvl, jobvr: ptr char,
  n: ptr int,
  a: ptr float64,
  lda: ptr int,
  wr,wi,vl: ptr float64,
  ldvl: ptr int,
  vr: ptr float64,
  ldvr: ptr int,
  work: ptr float64,
  lwork,info: ptr int,
) {.cdecl, importc: "dgeev_", dynlib: lapackName.}

proc ilaenv*(
  ispec: ptr int,
  name,opts: ptr string,
  n1,n2,n3,n4: ptr int
): int {.cdecl, importc: "ilaenv_", dynlib: lapackName.}

proc nfo*(info: int, f: string) =
  discard#for now