
cppLiteral <- function(text) {
    template <- quote(cppLiteral(text))
    template[[2]] <- text
    template
}

writeCode <- function(x, ...) writeLines(unlist(x), ...)

generateAll <- function(x, ...) lapply(x, function(y) y$generate(...))

codeSubstitute <- function(code, subList) {
    eval(substitute(substitute(code, subList), list(code = code)))
}
## This allows definition of a template and substitution into it with lighter syntax
## examples
## template <- quote({RNAME <- foo(CNAME)})
## codeSubstitute(template, list(RNAME = 'Ra', CNAME = 'ca'))
## l1 <- codeSubstitute(template, list(RNAME = 'Ra', CNAME = as.name('ca')))

putCodeLinesInBrackets <- function(codeLines) {
    as.call(c(as.name('{'), codeLines))
}

# This is the location of the RcppUtils.cpp, etc. files.
IncludeCodeDir = character()
# NimbleCodeDir = system.file("CppCode", package = "nimble")

nimbleIncludeFile =
function(file, path = IncludeCodeDir)
{
  if(length(path)) 
     sprintf('"%s"', normalizePath(sprintf("%s/%s", sub("/$", "", path), file)))
  else
     sprintf("<nimble/%s>", file)
}

makeDefaultDirName <- function()  normalizePath(file.path(tempdir(), 'nimble_generatedCode'), winslash = "\\", mustWork=FALSE)
