.onUnload <- function (libpath) {
  library.dynam.unload("limbs", libpath)
}