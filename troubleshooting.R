# Clean environment
rm(list = ls())
gc()

# load in library
library(nimble)
Sys.getenv()
path <- Sys.getenv('PATH')
newPath <- paste("C:\\rtools43\\usr\\bin\\;C:\\rtools43\\x86_64-w64-mingw32.static.posix\\bin\\;",
                 path, sep = "")
Sys.setenv(PATH = newPath) 

path <- Sys.getenv("PATH")
newPath <- paste("C:\\rtools43\\usr\\bin\\;C:\\rtools43\\mingw_64\\bin\\;", path, sep = "")
Sys.setenv(PATH = newPath)
Sys.setenv(BINPREF="C:\\rtools43\\x86_64-w64-mingw32.static.posix\\bin\\")

shortPathName("C:/Program Files/")
Sys.setenv(R_RTOOLS43_PATH = "C:\\rtools43\\x86_64-w64-mingw32.static.posix\\bin\\")

Sys.setenv("_files")
pkgbuild::find_rtools(debug = TRUE)

foo <- nimbleFunction( run = function(x = double(1)) {return(sum(x)); returnType(double())})
cfoo <- compileNimble(foo, showCompilerOutput = TRUE)
cfoo(1:10)


install.packages("Rcpp")
library("Rcpp") 
evalCpp("2 + 2")
Sys.setenv(BINPREF="")
Sys.setenv(RTOOLS40_HOME ="")
Sys.getenv('TEMP')
tempdir()
library(nimble,lib.loc='C:/Users/veron/AppData/Local/R/win-library')
pkgbuild::has_compiler(debug = TRUE)

RTOOLS_ROOT="C:/PROGRA~1/R/R-43~1.2/"
PATH="$rtools43/bin;${PATH}"
BINPREF="C:/rtools43/mingw64/bin"

shortPathName("C:/rtools43/mingw64/bin")

RTOOLS_ROOT="C:/PROGRA~1/R/Rtools-3.5"
PATH="${RTOOLS_ROOT}/bin;${PATH}"
BINPREF="${RTOOLS_ROOT}/mingw_$(WIN)/bin/"

# THE SOLUTION ----
# change temp file location
write("TMP = 'C:/Users/veron/AppData/Local/Temp'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))
