".First.lib" <-
function(lib,pkg)
{
   cat("load changeLOS: ", lib, "...\n")
   library.dynam("changeLOS",pkg,lib)
}
".Last.lib" <-
function(lib)
{
   cat("unload changeLOS: ", lib, "...\n")
   library.dynam.unload("changeLOS",lib)
}
