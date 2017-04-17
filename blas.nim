when defined(windows):
  const blasSuffix = ".dll"
elif defined(macosx):
  const blasSuffix = ".dylib"
else:
  const blasSuffix = ".so"

const blasPrefix = "libblas"
const blasName = blasPrefix & blasSuffix

proc gemm*(
  transa, transb: ptr char,
  m,n,k: ptr int,
  alpha,a: ptr float64,
  lda: ptr int,
  b: ptr float64,
  ldb: ptr int,
  beta,c: ptr float64,
  ldc: ptr int,
) {.cdecl, importc: "dgemm_", dynlib: blasName.}
