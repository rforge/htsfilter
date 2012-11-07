## Register vignette with Windows GUI
.onAttach <- function(libname, pkgname) {
	if (.Platform$OS.type == "windows" && 
		require("Biobase") && interaction() && 
		.Platform$GUI == "Rgui") {
		
		addVigs2WinMenu("HTSFilter")
	}
}