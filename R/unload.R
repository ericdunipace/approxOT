.onUnload <- function (libpath) {
  library.dynam.unload("approxOT", libpath)
}