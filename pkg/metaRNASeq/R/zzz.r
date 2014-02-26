.onAttach <- function(libname, pkgname)
{
	if(.Platform$OS.type=="windows" && 
		.Platform$GUI=="Rgui") {

		winMenuAddItem("Vignettes", "metaRNASeq",
		"shell.exec(system.file(\"vignettes\",\"metaRNASeq.pdf\", package=\"metaRNASeq\"))")

	}
}