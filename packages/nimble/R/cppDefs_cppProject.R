## Classes for a cppProject
cppCodeFileClass <- setRefClass('cppCodeFileClass',
                             fields = list(
                                 filename = 'ANY',	#character
                                 includes = 'ANY',	#character
                                 usings = 'ANY',	#character
                                 cppDefs = 'ANY'	#list()
                                 ),
                             methods = list(
                             	initialize = function(...){filename <<- character(); includes <<- character(); usings <<- character(); cppDefs <<- list(); callSuper(...)},

                                 writeIncludes = function(con = stdout()) {
                                     writeLines(c('#ifndef R_NO_REMAP', '#define R_NO_REMAP', '#endif'), con)
                                     if(length(includes) > 0) writeLines(paste0('#include ', includes), con)
                                 },
                                 writeUsings = function(con = stdout()) {
                                     if(length(usings) > 0) writeLines(paste0('using ', usings,';'), con)
                                 },
                                 writeDecs = function(con = stdout()) {
                                     lapply(cppDefs, function(x) {writeLines("", con); writeCode(x$generate(declaration = TRUE), con)})
                                 },
                                 writeDefs = function(con = stdout()) {
                                     lapply(cppDefs, function(x) {writeLines("", con); writeCode(x$generate(declaration = FALSE), con)})
                                 }
                                 )
                             )

cppHfileClass <- setRefClass('cppHfileClass',
                             contains = 'cppCodeFileClass',
                             fields = list(
                                 ifndefName = 'ANY' #'character'
                                 ),
                             methods = list(
                             	initialize = function(...){ifndefName <<- character(); callSuper(...)},
                                 writeFile = function(con = filename, dir = character()) {
                                     if(is.character(con)) {
                                         con <- normalizePath(file.path(dir, paste0(con, '.h')), winslash = "\\", mustWork=FALSE)
                                         zz <- file(con, open = 'w')
                                         on.exit(close(zz))
                                     } else
                                         zz <- con
                                     writeIfndef(zz)
                                     writeIncludes(zz)
                                     writeUsings(zz)
                                     writeDecs(zz)
                                     closeIfndef(zz)
                                     invisible(NULL)
                                 },
                                 writeIfndef = function(con = stdout()) {
                                     writeLines(paste0('#ifndef ',ifndefName,
                                                       '\n#define ',ifndefName),
                                                con)
                                 },
                                 closeIfndef = function(con = stdout()) {
                                     writeLines('#endif', con)
                                 }
                                 )
                             )

cppCPPfileClass <- setRefClass('cppCPPfileClass',
                               contains = 'cppCodeFileClass',
                               fields = list(
                                   ifndefName = 'ANY'
                                   ),
                               methods = list(
                               		initialize = function(...){ifndefName <<- character(); callSuper(...)},
                                   writeFile = function(con = filename, dir = character()) {
                                       if(is.character(con)) {
                                           con <- normalizePath(file.path(dir, paste0(con,'.cpp')), winslash = "\\", mustWork=FALSE)
                                           zz <- file(con, open = 'w')
                                           on.exit(close(zz))
                                       } else
                                           zz <- con
                                       writeIfndef(zz)
                                       writeIncludes(zz)
                                       writeUsings(zz)
                                       writeDefs(zz)
                                       closeIfndef(zz)
                                       invisible(NULL)
                                   },
                                   writeIfndef = function(con = stdout()) {
                                       if(length(ifndefName) > 0)
                                           writeLines(paste0('#ifndef ',ifndefName,
                                                             '\n#define ',ifndefName),
                                                      con)
                                 },
                                   closeIfndef = function(con = stdout()) {
                                       if(length(ifndefName) > 0) writeLines('#endif', con)
                                   }
                                   ))

cppProjectClass <- setRefClass('cppProjectClass',
                               fields = list(
                                   dirName = 'ANY', #'character',
                                   cppDefs = 'ANY', #'list',
                                   dll = 'ANY', #"DLLInfoOrNULL",
                                   outputSOfile = 'ANY',# "character"
                                   Oincludes = 'ANY'
                                   ),
                               methods = list(
                               		initialize = function(...){dirName <<- character(); cppDefs <<- list(); dll <<- NULL; outputSOfile <<- character();callSuper(...)},
                                   addFunction = function(funDef, name, filename) {
                                       if(missing(name)) name <- funDef$name
                                       cppDefs[[name]] <<- funDef
                                       if(!missing(filename)) {
                                           filename <- Rname2CppName(filename); funDef$filename <- filename
                                       } else {
                                           if(length(funDef$filename)==0) {
                                               filename <- Rname2CppName(name); funDef$filename <- filename
                                           }
                                       }
                                   },
                                   addClass = function(classDef, name, filename, includeNeededTypeDefs = TRUE) {
                                       if(missing(name)) name <- classDef$name
                                       if(!missing(filename)) {
                                           filename <- Rname2CppName(filename)
                                       } else {
                                           if(length(classDef$filename)==0) {
                                               filename <- Rname2CppName(name)
                                           } else {
                                               stop("Error in addClass: Can't determine filename")
                                           }
                                       }
                                       if(includeNeededTypeDefs) {
                                           for(iNTD in classDef$neededTypeDefs) addClass(iNTD, filename = filename)
                                       }
                                       classDef$filename <- filename
                                       cppDefs[[name]] <<- classDef
                                   },
                                   writeFiles = function(filename, con = filename) {
                                       filename <- Rname2CppName(filename)

                                       whichDefs <- which(unlist(lapply(cppDefs, `[[`, 'filename')) == filename)

                                       defs <- cppDefs[whichDefs]
                                       Hincludes <- unlist(lapply(defs, function(x) lapply(x$getHincludes(), function(xx) if(is.character(xx)) xx else paste0('\"',xx$filename,'.h\"'))))
                                       Hincludes <- unique(Hincludes)
                                       CPPincludes <- unlist(lapply(defs, function(x) lapply(x$getCPPincludes(), function(xx) if(is.character(xx)) xx else paste0('\"', xx$filename,'.cpp\"'))))
                                       CPPincludes <- unique(CPPincludes)
                                       selfCPP <- if(is.character(con)) paste0('"', con, '.cpp"') else '"[FILENAME].cpp"'
                                       CPPincludes <- CPPincludes[ CPPincludes != selfCPP ]

                                      ## TODO Simplify Eigen include logic now that Nimble defines `R_NO_REMAP`.
                                       ## similar for cppad
                                       iEigenInclude <- grep("EigenTypedefs", CPPincludes)
                                       if(length(iEigenInclude) > 0) {
                                           CPPincludes <- c(CPPincludes[iEigenInclude], CPPincludes[-iEigenInclude])
                                       }
                                       iCppInclude <- grep("cppad", CPPincludes)
                                       if(length(iCppInclude) > 0) {
                                           CPPincludes <- c(CPPincludes[iCppInclude], CPPincludes[-iCppInclude])
                                       }

                                       ## at this point strip out CPPincludes other than EigenTypedefs that have .cpp and gsub .cpp to .o
                                       boolConvertCppIncludeToOinclude <- grepl("\\.cpp", CPPincludes)
                                       Oincludes <<- gsub("\\.cpp", ".o", CPPincludes[boolConvertCppIncludeToOinclude])
                                       CPPincludes <- CPPincludes[!boolConvertCppIncludeToOinclude]

                                       CPPusings <- unlist(lapply(defs, function(x) x$getCPPusings()))
                                       CPPusings <- unique(CPPusings)

                                       ifndefName <- if(is.character(con)) toupper(paste0('__', con)) else '"__IFNDEFNAME"'
                                       cppPieces <- do.call('c', lapply(defs, function(x) x$getDefs()))

                                       hFile <- cppHfileClass(filename = filename,
                                                              includes = Hincludes,
                                                              cppDefs = cppPieces,
                                                              ifndefName = ifndefName)
                                       selfInclude <- if(is.character(con)) paste0('"', con, '.h', '"') else '"[FILENAME].h"'
                                       CPPincludes <- c(CPPincludes, selfInclude) ## selfInclude has to come last because Rinternals.h makes a name conflict with Eigen (this may be moot, 7/17)
                                       
                                       cppIfndefName <- paste0(ifndefName,'_CPP')
                                       cppFile <- cppCPPfileClass(filename = filename,
                                                                  includes = CPPincludes,
                                                                  usings = CPPusings,
                                                                  cppDefs= cppPieces,
                                                                  ifndefName = cppIfndefName
                                                                  )

                                       if(is.character(con)) createDir()
                                       hFile$writeFile(con = con, dir = dirName)
                                       cppFile$writeFile(con = con, dir = dirName)
                                   },
                                   writeDynamicRegistrationsDotCpp = function(dynamicRegistrationsCppName, dllName) {
                                       ## this writes dynamicRegistrations.cpp to include <nimble/dynamicRegistrations.h> and
                                       ## then add an R_init function, which must have the name R_init_[outputSOfile]
                                       ##
                                       ## At the moment this step is called immediately before doing R CMD SHLIB
                                       ## And it is only at that stage that the outputSOfile is known, since the SOname includes a time stamp
                                       ## for uniqueness that is created right before R CMD SHLIB.
                                       ## However this content for this file could be integrated into the cppDefs
                                       ## and then automatically written as part of writeFiles.  Doing so would require
                                       ## that the SOname be generated earlier, which would probably be fine.
                                       contentLines <- c(
                                            "#include <nimble/dynamicRegistrations.h>",
                                           "",
                                           "extern \"C\"",
                                           paste0("void R_init_", dllName, "(DllInfo *dll) {"),
                                           "  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);",
                                           "}")
                                       writeLines(contentLines, con = dynamicRegistrationsCppName)
                                   },
                                   compileStaticCode = function(dllName, cppName, showCompilerOutput) {
                                       ssDllName <- normalizePath(file.path(dirName, paste0(dllName, .Platform$dynlib.ext)), winslash = "\\", mustWork=FALSE)
                                       ssdSHLIBcmd <- normalizePath(file.path(R.home('bin'), 'R'),
                                                                          winslash = "\\", mustWork=FALSE)
                                       ssdSHLIBargs <- paste('CMD SHLIB',
                                                             ifelse(nimbleOptions()$precleanCompilation, '--preclean', ''),
                                                             cppName, '-o', basename(ssDllName))

                                       logFile <- paste0(dllName, ".log")
                                       errorFile <- paste0(dllName, ".err")
                                       status = system2(ssdSHLIBcmd, ssdSHLIBargs, stdout = logFile, stderr = errorFile)

                                       if(status == 0 && showCompilerOutput) {
                                           output <- c(readLines(logFile), readLines(errorFile))
                                           cat(output, sep = '\n')
                                       }
                                       if(status != 0) {
                                           if(showCompilerOutput) {
                                               warnLength <- options()$warning.length
                                               options(warning.length = 8000)
                                               on.exit(options(warning.length = warnLength))
                                               
                                               output <- c(readLines(logFile), readLines(errorFile))
                                               stop(structure(simpleError(paste0("Failed to create the shared library.\n", paste0(output, collapse = "\n"))),
                                                              class = c("SHLIBCreationError", "ShellError", "simpleError", "error", "condition")))

                                           } else {
                                               nimbleUserNamespace$errorFile <- normalizePath(file.path(dirName, errorFile), winslash = "\\", mustWork=FALSE)
                                               stop(structure(simpleError(paste0("Failed to create the shared library. Run 'printErrors()' to see the compilation errors.\n")),
                                                              class = c("SHLIBCreationError", "ShellError", "simpleError", "error", "condition")))
                                           }
                                       }
                                       return(dyn.load(basename(ssDllName), local = TRUE))
                                   },
                                   compileDynamicRegistrations = function(showCompilerOutput = nimbleOptions('showCompilerOutput')) {
                                       timeStamp <- format(Sys.time(), "%m_%d_%H_%M_%S")
                                       dllName <- paste0("dynamicRegistrations_", timeStamp)
                                       cppName <- paste0(dllName, ".cpp")
                                       writeDynamicRegistrationsDotCpp(cppName, dllName)
                                       nimbleUserNamespace$sessionSpecificDll <- compileStaticCode(dllName, cppName, showCompilerOutput)
                                   },
                                 compile_nimbleCppADbaseClass = function(showCompilerOutput = nimbleOptions('showCompilerOutput')) {
                                   timeStamp <- format(Sys.time(), "%m_%d_%H_%M_%S")
                                   dllName <- paste0("nimbleCppADbaseClass_", timeStamp)
                                   cppName <- "nimbleCppADbaseClass.cpp"
                                   origFile <-  system.file(file.path("include","nimble", cppName), package = "nimble")
                                   if(nchar(origFile) == 0) {
                                     warning("Could not compile nimbleCppADbaseClass. Subsequent steps will likely generate errors.")
                                     return(NULL)
                                   }
                                   file.copy(origFile, cppName)
                                   # We don't need this DLL, but we might as well hold onto it since it exists.
                                   nimbleUserNamespace$nimbleCppADbaseClassDll <- compileStaticCode(dllName, cppName, showCompilerOutput)
                                 },
                                   compileFile = function(names, showCompilerOutput = nimbleOptions('showCompilerOutput'),
                                                          .useLib = UseLibraryMakevars) {
                                       names <- Rname2CppName(names)
                                       isWindows = (.Platform$OS.type == "windows")

                                       includes <- character()
                                       timeStamp <- format(Sys.time(), "%m_%d_%H_%M_%S")

                                       mainfiles <- paste(basename(
                                           normalizePath(file.path(dirName, paste0(names,'.cpp')), winslash = "\\", mustWork=FALSE)
                                       ),
                                       collapse = ' ')

				       if(!file.exists(normalizePath(file.path(dirName, sprintf("Makevars%s", if(isWindows) ".win" else "")), winslash = "\\", mustWork=FALSE)) && NeedMakevarsFile) # should reverse the order here in the long term.
				           createMakevars(.useLib = .useLib, dir = dirName)

                                       dllName <- paste0(names[1], "_", timeStamp)
                                                                             
                                       outputSOfile <<- normalizePath(file.path(dirName, paste0(dllName, .Platform$dynlib.ext)), winslash = "\\", mustWork=FALSE)

                                       if(!inherits(Oincludes, 'uninitializedField')) { ## will only be uninitialized if writeFiles was skipped due to specialHandling (developer backdoor)
                                           includes <- c(includes, Oincludes) ## normal operation will have Oincludes.
                                       }
                                       SHLIBcmd <- normalizePath(file.path(R.home('bin'), 'R'), winslash = "\\", mustWork=FALSE)
                                       SHLIBargs <- paste('CMD SHLIB',
                                                          ifelse(nimbleOptions()$precleanCompilation, '--preclean', ''),
                                                          paste(c(mainfiles, includes), collapse = ' '), '-o', basename(outputSOfile))

                                       cur = getwd()
                                       setwd(dirName)
                                       on.exit(setwd(cur))

                                       if(is.null(nimbleUserNamespace$sessionSpecificDll)) {
                                           compileDynamicRegistrations(showCompilerOutput = showCompilerOutput)
                                       }
                                       if(isTRUE(nimbleOptions("enableDerivs"))) {
                                         if(any(grepl("^nimbleCppADbaseClass.o$", Oincludes))) {
                                           if(is.null(nimbleUserNamespace$nimbleCppADbaseClassDll)) {
                                             compile_nimbleCppADbaseClass(showCompilerOutput = showCompilerOutput)
                                           }
                                         }
                                       }
                                       origSHLIBcmd <- SHLIBcmd
                                       if(isTRUE(nimbleOptions('stopCompilationBeforeLinking'))) {## used only for testing, when we want to go quickly and skip linking and bail out
                                           ## get the dry run commands, run only those that contain -c for compile-only (don't link)
                                           ## this has only been tested with single .cpp files, not multiple .cpp files
                                           stop("Option 'stopCompilationBeforeLinking' has been disabled.")
                                           dryRunCmd <- paste0(SHLIBcmd, " -n")
                                           dryRunResult <- system(dryRunCmd, intern = TRUE)
                                           compileOnlyLines <- dryRunResult[ grepl("-c", dryRunResult) ]
                                           SHLIBcmd <- paste0(compileOnlyLines, collapse =  ";" )
                                       }

                                       if(isTRUE(nimbleOptions('forceO1'))) { ## replace -On flags with -O1 to reduce compiler time due to higher optimization levels 
                                           ## If forceO1 is TRUE and we did not already strip out -c flags, do so now
                                           stop("Option 'forceO1' has been disabled.")
                                           if(!isTRUE(nimbleOptions('stopCompilationBeforeLinking'))) {
                                               dryRunCmd <- paste0(SHLIBcmd, " -n")
                                               dryRunResult <- system(dryRunCmd, intern = TRUE)
                                               compileOnlyLines <- dryRunResult[ grepl("-c", dryRunResult) ]
                                               SHLIBcmd <- paste0(compileOnlyLines, collapse =  ";" )
                                           }
                                           SHLIBcmd <- gsub("-O[1-9]", "-O1", SHLIBcmd)
                                           SHLIBcmd <- paste0(SHLIBcmd, "; ", origSHLIBcmd)
                                       }
 
                                       
                                       logFile <- paste0(names[1], "_", format(Sys.time(), "%m_%d_%H_%M_%S"), ".log")
                                       errorFile <- paste0(names[1], "_", format(Sys.time(), "%m_%d_%H_%M_%S"), ".err")

                                       if(nimbleOptions('pauseAfterWritingFiles')) browser()
                                       ## We formerly used ignore.stdout = !showCompilerOutput, ignore.stderr = !showCompilerOutput
                                       ## but when ignore.stdout and ignore.stderr are TRUE nothing gets printed to stdout and stderr so
                                       ## .log and .err files are empty.
                                       status = system2(SHLIBcmd, SHLIBargs, stdout = logFile, stderr = errorFile)
                                       if(status == 0 && showCompilerOutput) {
                                           output <- c(readLines(logFile), readLines(errorFile))
                                           cat(output, sep = '\n')
                                       }
                                       if(status != 0) {
                                           if(showCompilerOutput) {
                                               warnLength <- options()$warning.length
                                               options(warning.length = 8000)
                                               on.exit(options(warning.length = warnLength))

                                               output <- c(readLines(logFile), readLines(errorFile))
                                               stop(structure(simpleError(paste0("Failed to create the shared library.\n", paste0(output, collapse = "\n"))),
                                                              class = c("SHLIBCreationError", "ShellError", "simpleError", "error", "condition")))

                                           } else {
                                               nimbleUserNamespace$errorFile <- normalizePath(file.path(dirName, errorFile), winslash = "\\", mustWork=FALSE)
                                               stop(structure(simpleError(paste0("Failed to create the shared library. Run 'printErrors()' to see the compilation errors.\n")),
                                                              class = c("SHLIBCreationError", "ShellError", "simpleError", "error", "condition")))
                                           }
                                       }
                                       if(isTRUE(nimbleOptions()$stopCompilationBeforeLinking)) stop("safely stopping before linking", call.=FALSE)
                                   },
                                   loadSO = function(name) {
                                       dll <<- dyn.load(getSOName(), local = TRUE)
                                   },
                                   unloadSO = function(check = TRUE, force = FALSE) { ## The book-keeping on different names isn't quite connected to here yet.  Instead we just unload dll.
				       if(!is.null(dll)) {
                                           objectNames <- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$RNimble_Ptr_CheckAndRunAllDllFinalizers, dll[['handle']], force))
                                           if(length(objectNames) > 0 & check) {
                                               warning(paste0("A DLL to be unloaded has non-zero (", paste(objectNames, collapse = ", "), ") objects that need to be finalized first. ", if(force) "It's objects were cleared and it was unloaded anyway." else "It was not unloaded." ))
                                               browser()
                                               if(!force) return(NULL)
                                           }
                                           status = dyn.unload(dll[["path"]])
                                           dll <<- NULL
                                           status
                                       } else
                                           FALSE
                                   },
                                   getSOName = function() {
                                           return(outputSOfile)
                                       },
                                   createDir = function() {
                                       if(!file.exists(dirName)) dir.create(dirName)
                                   })
                               )

