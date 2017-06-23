## C++ definition for a namedObjects class
## This is a permanaent C++ base class to provide one interface
## to get external pointers to class member data by character names

cppNamedObjectsClass <- setRefClass('cppNamedObjectsClass',
                               contains = 'cppClassDef',
                               fields = list(Rnames2CppNames = 'ANY'),
                               methods = list(
                                   initialize = function(...) {
                                   		Rnames2CppNames <<- list()
                                       inheritance <<- c(inheritance, list('NamedObjects'))
                                       Hincludes <<- c(Hincludes, nimbleIncludeFile("NamedObjects.h"))
                                       callSuper(...)
                                   },
                                   namedObjectsConstructorCodeBlock = function() {
                                       template <- quote( {namedObjects[RNAME] = cppReference(CNAME)} )[[2]]
                                       codeLines <- vector('list', length = length(Rnames2CppNames))
                                       for(i in seq_along(Rnames2CppNames))
                                           codeLines[[i]] <- codeSubstitute(template, list(RNAME = names(Rnames2CppNames)[i],
                                                                                           CNAME = as.name(Rnames2CppNames[[i]])))
                                       cppCodeBlock(code = putCodeLinesInBrackets(codeLines), skipBrackets = TRUE)
                                   }
                               )
                               )
