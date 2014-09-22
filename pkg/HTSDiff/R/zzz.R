.onAttach <- function(libname, pkgname)
{
	if(.Platform$OS.type=="windows" && 
		.Platform$GUI=="Rgui") {

		winMenuAddItem("Vignettes", "HTSDiff",
		"shell.exec(system.file(\"doc\",\"HTSDiff.pdf\", package=\"HTSDiff\"))")

	}
}