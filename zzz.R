.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "WARNING: '", pkgname, "' is no longer maintained.\n",
    "The package has been moved and expanded. Please install the new version:\n",
    "  - CRAN:   install.packages(\"multiDEGGs\")\n",
    "  - GitHub: https://github.com/elisabettasciacca/multiDEGGs\n",
    "This version will not receive further updates."
  )
}
