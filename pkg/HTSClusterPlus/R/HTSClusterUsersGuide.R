#' View HTSCluster User's Guide
#' 
#' Finds the location of the HTSCluster User's Guide and optionally opens it.
#' 
#' The function \code{vignette("HTSCluster")} will find the short HTSCluster
#' Vignette which describes how to obtain the HTSCluster User's Guide. The
#' User's Guide is not itself a true vignette because it is not automatically
#' generated using \code{\link{Sweave}} during the package build process. This
#' means that it cannot be found using \code{vignette}, hence the need for this
#' special function.
#' 
#' If the operating system is other than Windows, then the PDF viewer used is
#' that given by \code{Sys.getenv("R_PDFVIEWER")}. The PDF viewer can be
#' changed using \code{Sys.putenv(R_PDFVIEWER=)}.
#' 
#' Note that this function was adapted from that defined by Gordon Smyth in the
#' edgeR package.
#' 
#' @param view logical, should the document be opened using the default PDF
#' document reader?
#' @return Character string giving the file location. If \code{view=TRUE}, the
#' PDF document reader is started and the User's Guide is opened, as a side
#' effect.
#' @author Gordon Smyth
#' @seealso \code{\link{system}}
#' @keywords documentation
#' @examples
#' 
#' # To get the location:
#' HTSClusterUsersGuide(view=FALSE)
#' # To open in pdf viewer:
#' \dontrun{HTSClusterUsersGuide()}
#' 
#' @export HTSClusterUsersGuide
HTSClusterUsersGuide <- function(view=TRUE)
#	Adapted from edgeR User's Guide
{
	f <- system.file("doc","HTSClusterUsersGuide.pdf",package="HTSCluster")
	if(view) {
		if(.Platform$OS.type == "windows") 
			shell.exec(f)
		else
			system(paste(Sys.getenv("R_PDFVIEWER"),f,"&"))
	}
	return(f)
}
