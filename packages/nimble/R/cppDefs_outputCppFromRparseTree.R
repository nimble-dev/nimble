## We have two systems for outputting C++ code
## This is the older system.
## It uses R parse trees that have appropriate contents (not necessarily evaluatable in R)

## The only two parameters we need for any function here is the modified code and Indentation and keeps increasing or decreasing accordingly

# NOTE: Jagadish had some code embedded in the roxygen help in the description fields that was messing up roxygen processing, so I moved those lines out of the roxygen documentation - CJP 6/17/14

makeSpaces <- function(num) paste(rep(' ', num), collapse = '')

pasteSemicolon <- function(x, indent = '') {
    if(is.numeric(indent)) indent <- makeSpaces(indent)
    if(is.character(x)) return(paste0(indent, x, ';'))
    if(is.list(x)) return(lapply(x, function(ix) paste0(indent, ix, ';')))
    if(is.null(x)) return(character())
    stop(paste0('Error, pasteSemicolon called for object of class ', class(x), '. Must be character or list.')) 
}

#Special Lists
cppKeywordsThatFillSemicolons <- c("cppCommentBlock","scopeBlock", "{", "cppFor", "if", "cppWhile", "cppLiteral")

CppBinaryOperators2 <- list('+' = '+', '-' = '-', '/' = '/', '*' = '*',
                             '%+=%' = '+=', '%-=%' = '-=', '%/=%' = '/=', '%*=%' = '*=',
                             '=' = '=',  '==' = '==',  '!=' = '!=',  '<=' = '<=',
                             '>=' = '>=',  '>' = '>',  '<' = '<', '%%' = '%',
                             '&&' = '&&', '||' = '||', '&' = '&', '|' = '|',
                           '%.%' = '.', '%->%' = '->') ## used to have a separation function "withoutSpace" for these

ModifiedRmmParseKeywords2 <- c( structure( as.list(rep('outputCppBinaryOperation2', length(CppBinaryOperators2))), names = names(CppBinaryOperators2)),
                               list(
                                   '(' = 'outputCppOpenParentheses2',
                                   'cppOpenParentheses' = 'outputCppOpenParentheses2', 'if' = 'outputCppIfElse2', '{' = 'outputCppOpenBracketInvisible',
                                   scopeBlock = 'outputCppOpenBracket2',
                                   ##'[[' = 'outputCppArrayIndex2',
                                   '[' = 'outputCppArrayIndex2', 'cppFor' = 'outputCppFor2', 'cppDereference' = 'outputCppDereference2',
                                   'cppReference' = 'outputCppReference2', 'cppCommentBlock' = 'outputCppCommentBlock2', 'cppCout' = 'outputCppCout2',
                                   'c' = 'outputCppTempArray2', 'cppExit' = 'outputCppExit2', 'cppWhile' = 'outputCppWhile2', 'cppScopeResolution' = 'outputCppScopeResolution2',
                                   'EigenOperation' = 'outputEigenOperation2', 'TypeCast' = 'outputTypeCast2',
                                   'cppLiteral' = 'outputCppLiteral'
                                   )
                               )
## '[[' should not be there.  Looks like '[' was for Nim Array index only.
## 'c' is still there for array initialization I guess
##

outputCppLiteral <- function(modifiedRmmCode, indent = '') {
    contents <- eval(modifiedRmmCode[[2]])
    if(length(contents) > 1) as.list(contents) else contents
}

# e.g. outputEigenOperation(quote(TypeCast("double", A)))


# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is
# a TypeCast
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
# @export
outputTypeCast2 <- function(ModifiedRmmCode, Indentation = '\t') {
  paste0(ModifiedRmmCode[[2]], '(', outputCppParseTree2(ModifiedRmmCode[[3]], Indentation), ')')
}

# e.g. outputEigenOperation(quote(EigenOperation("array", A)))
# e.g. outputEigenOperation(quote(EigenOperation("matrix", A)))

# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is
# an Eigen Operation
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
# @export
outputEigenOperation2 <- function(ModifiedRmmCode, Indentation = '\t') {
  ## to avoid extra brackets like (A.matrix()).transpose() we just want A.matrix().transpose()
  if(is.call(ModifiedRmmCode[[3]]) && !(deparse(ModifiedRmmCode[[3]][[1]]) %in% c('(', 'EigenOperation', 'TypeCast'))) 
    withinBrackets <- TRUE
   else withinBrackets <- FALSE
  paste0(if(withinBrackets) '(',  outputCppParseTree2(ModifiedRmmCode[[3]], Indentation), if(withinBrackets) ')',
         '.', ModifiedRmmCode[[2]], '(',
         if(length(ModifiedRmmCode) > 3) {
           Parameters <- list() ## handling all parameters within the function eg A.pow(3), so here 3 is the parameter, all we do is separate with a comma
           for(i in 4:length(ModifiedRmmCode))
            Parameters[length(Parameters) + 1] <- outputCppParseTree2(ModifiedRmmCode[[i]], Indentation)
           paste0(Parameters, collapse = ', ')
           },
         ')'
  )
}

# e.g. outputCppExit(quote(cppExit('Error.')))


# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is cppExit
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
# @export
outputCppExit2 <- function(ModifiedRmmCode, Indentation = '\t') {
  paste0('RBREAK("', ModifiedRmmCode[[2]], '")')   
}

# e.g. outputCppCout(quote(cppCout(a)))

# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is cppOut
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
# @export
outputCppCout2 <- function(ModifiedRmmCode, Indentation = '') {
    paste0(Indentation, 'cout<<', outputCppParseTree2(ModifiedRmmCode[[2]], ''), '<<"\\n"')
}

# e.g.  outputCppCommentBlock(quote(comment('...of the world (Still singing eh?)')))

# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is a CppCommentBlock
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
# @export
outputCppCommentBlock2 <- function(ModifiedRmmCode, Indentation = '') {
  paste0('/*\n', Indentation, if(is.character(ModifiedRmmCode[[2]])) ModifiedRmmCode[[2]] else deparse(ModifiedRmmCode[[2]]),'\n', '*/')
}

# e.g. outputCppArrayIndex(quote(mv[i - 1]))

# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is CppArrayIndex
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
outputCppArrayIndex2 <- function(ModifiedRmmCode, Indentation = '') { ## this has only one element inside the []
    partThree <- if(length(ModifiedRmmCode) == 3) outputCppParseTree2(ModifiedRmmCode[[3]],'') else NULL
    paste0(Indentation, outputCppParseTree2(ModifiedRmmCode[[2]], ''), '[', partThree, ']')
}


# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is a CppOpenParentheses
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
# @export
outputCppOpenParentheses2 <- function(ModifiedRmmCode, Indentation = '') {
    BlockLength <- length(ModifiedRmmCode)
    outputCppCode <- character(BlockLength - 1)
    if(BlockLength > 1)
        for(i in 2:BlockLength)
            outputCppCode[i-1] <- outputCppParseTree2(ModifiedRmmCode[[i]], '')
    return(paste0(Indentation,'(', paste0(outputCppCode, collapse = ', '), ')'))
}


# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is a concatenate i.e. c()
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
# @export
outputCppTempArray2 <- function(ModifiedRmmCode, Indentation = '') {
  ## c(1,2,3) becomes {1,2,3}
    BlockLength <- length(ModifiedRmmCode)
    outputCppCode <- character(BlockLength-1)
    if(BlockLength > 1)
        for(i in 2:BlockLength)
            outputCppCode[i-1] <- outputCppParseTree2(ModifiedRmmCode[[i]], '')
    return(paste0(Indentation, '{', paste0(outputCppCode, collapse = ', '), '}'))
}


# e.g. outputCppReference(quote(cppReference(x)))

# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is a CppReference
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
# @export
outputCppReference2 <- function(ModifiedRmmCode, Indentation = '') {
##  if(is.recursive(ModifiedRmmCode[[2]])) {
    if(!is.name(ModifiedRmmCode[[2]])) {
      return(paste0(Indentation, '&(', outputCppParseTree2(ModifiedRmmCode[[2]], ''), ')'))
  }
  else
      return(paste0(Indentation, '&', outputCppParseTree2(ModifiedRmmCode[[2]], '')))  
}

# e.g. outputCppDeReference(quote(cppDeReference(x)))


# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is a CppDeReference
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
# @export
outputCppDereference2 <- function(ModifiedRmmCode, Indentation = '') {
    if(!is.name(ModifiedRmmCode[[2]]))
        return(paste0(Indentation, '(*(', outputCppParseTree2(ModifiedRmmCode[[2]], ''), '))'))
    else
        return(paste0(Indentation, '(*', outputCppParseTree2(ModifiedRmmCode[[2]], ''), ')'))
}

# e.g. outputGeneralCppFunction(quote(some_random_function(x)))

# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is a GeneralCppFunction
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
# @export
outputGeneralCppFunction2 <- function(ModifiedRmmCode, indent = '') {
  BlockLength <- length(ModifiedRmmCode)
  outputCppCode <- character(BlockLength-1)
  if(BlockLength > 1)
    for(i in 2:BlockLength)
      outputCppCode[i-1] <- outputCppParseTree2(ModifiedRmmCode[[i]], '')

  functionName <- if(is.name(ModifiedRmmCode[[1]])) genName2(deparse(ModifiedRmmCode[[1]]))
  else outputCppParseTree2(ModifiedRmmCode[[1]], '')
  
  return(paste0(paste0(indent, functionName, '('),
                paste0(outputCppCode, collapse = ', '),
                ')'))
}

genName2 <- function(x) x

# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is a CppOpenBracket
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
# @export
outputCppOpenBracketInvisible  <- function(ModifiedRmmCode, Indentation = '') {
  ## curly bracket, all elements in the parse tree are new lines
  BlockLength <- length(ModifiedRmmCode)
  if(BlockLength > 1)
      outputCppCode <- vector('list', length = BlockLength-1)
  else
      return(list())
  for(i in 2:BlockLength) {
      outputCppCode[[i-1]] <- outputCppParseTree2(ModifiedRmmCode[[i]], Indentation)
      if(!((is.call(ModifiedRmmCode[[i]]) &&
            as.character(ModifiedRmmCode[[i]][[1]]) %in% cppKeywordsThatFillSemicolons) ||
           inherits(ModifiedRmmCode[[i]], 'cppCodeBlock')))  ## This prevents adding semicolons to blocks that should have added their own.
          outputCppCode[[i-1]] <- pasteSemicolon(outputCppCode[[i-1]])
  }
  return(outputCppCode)
}

outputCppOpenBracket2  <- function(ModifiedRmmCode, Indentation = '') {
  ## curly bracket, all elements in the parse tree are new lines
  BlockLength <- length(ModifiedRmmCode)
  outputCppCode <- vector('list', length = BlockLength+1) 
  outputCppCode[[1]] <- paste0(Indentation, '{')
  if(BlockLength > 1)
    for(i in 2:BlockLength) {
        outputCppCode[[i]] <- outputCppParseTree2(ModifiedRmmCode[[i]], paste0(Indentation,' '))
      if(!(is.call(ModifiedRmmCode[[i]]) &&
           as.character(ModifiedRmmCode[[i]][[1]]) %in% cppKeywordsThatFillSemicolons))
          outputCppCode[[i]] <- pasteSemicolon(outputCppCode[[i]])
    }
  outputCppCode[[BlockLength+1]] <- paste0(Indentation, '}')
  return(outputCppCode)
}

# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is CppBinaryOperation
# e.g. outputCppBinaryOperation(quote(a + b))
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
# @export
outputCppBinaryOperation2 <- function(ModifiedRmmCode, Indentation = '') {
    if(length(ModifiedRmmCode) != 3) stop(paste0('Error generating C++ output for expression \"', deparse(ModifiedRmmCode), '\" which was expected to have length 3.'))
    return(paste0(Indentation, outputCppParseTree2(ModifiedRmmCode[[2]], ''), CppBinaryOperators2[[ deparse(ModifiedRmmCode[[1]]) ]],
        outputCppParseTree2(ModifiedRmmCode[[3]], '')))
}


# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is the
# Scope Resolution Operator
# e.g. outputCppScopeResolution(quote(cppScopeResolution(Pump,theta)))
# e.g. outputCppScopeResolution(quote(cppScopeResolution(x)))
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
# @export
outputCppScopeResolution2 <- function(ModifiedRmmCode, Indentation = '\t') {
  ## this handles A::B and also global objects ::Z
  if(length(ModifiedRmmCode) == 2)
    return(paste0('::', ModifiedRmmCode[2], collapse = ''))
  return(paste0(ModifiedRmmCode[-1], collapse = '::'))
}


# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is a CppFor
# e.g. outputCppFor(quote(cppFor(j, 1, m, {sum = sum + j})))
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
# @export
outputCppFor2 <- function(ModifiedRmmCode, Indentation = '') {
  loopVar <- outputCppParseTree2(ModifiedRmmCode[[2]], '')
  Begin <-  outputCppParseTree2(ModifiedRmmCode[[3]], '')
  End <-  outputCppParseTree2(ModifiedRmmCode[[4]], '')
  cppForBody <- outputCppParseTree2(ModifiedRmmCode[[5]], paste0(Indentation,' '))
  ans <- list(paste0(Indentation, 'for(', loopVar, ' = ', Begin, '; ', loopVar, ' <= ', End, '; ', loopVar, '++) {'),
              cppForBody,
              paste0(Indentation, '}'))
  ans
}


# e.g. outputCppFor(quote(cppWhile(a < 5, {a %+=% 2})))

# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is a CppWhile
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}
# @export
outputCppWhile2 <- function(ModifiedRmmCode, Indentation = '\t') {
  condition <- outputCppParseTree2(ModifiedRmmCode[[2]], Indentation)
  cppWhileBody <- outputCppParseTree2(ModifiedRmmCode[[3]], Indentation)
  return(paste0('while(', condition, ')', cppWhileBody))
}

# e.g. outputCppIfElse(quote(if(a == 5, {a %+=% 2})))
# e.g. outputCppIfElse(quote(if(a == 5, {a %+=% 2}, {a %-=% 2})))

# Outputs the Low Level Nimble Parse Tree after modifying it to C++ format
# 
# To print the modified the Rmm Parse Code we call the Function "outputCppParseTree" which in turn will
# call an appropriate function to print the Code. This function is called when the current element is an 'if'
#
# @param ModifiedRmmCode The Modified Low Level Nimble Parse Tree 
# @param Indentation for spacing in the Cpp code to make it readable (default: One Tab Space)
# @return returns C++ Output Code
# @author Jagadish Babu
# @keywords OutputCppParseTree
# @seealso \code{\link{outputCppParseTree}}, \code{\link{outputCppIf}}
# @export
outputCppIfElse2  <- function(ModifiedRmmCode, Indentation = '') {
    codeLength <- length(ModifiedRmmCode)
    if(!(codeLength == 3 | codeLength == 4)) stop(paste0('Error outputting \"', deparse(ModifiedRmmCode), '\". Length of if-then-else code must be 3 or 4'))
    outputCppBlock <- vector('list', if(codeLength == 3) 3 else 5)
    outputCppBlock[[1]] <- paste0(Indentation, 'if(', outputCppParseTree2(ModifiedRmmCode[[2]]), ') {')
    outputCppBlock[[2]] <- outputCppParseTree2(ModifiedRmmCode[[3]], paste0(' ', Indentation))
    if(length(ModifiedRmmCode)==3) {
        outputCppBlock[[3]] <- paste0(Indentation, '}')
    } else {
        outputCppBlock[[3]] <- paste0(Indentation, '} else {')
        outputCppBlock[[4]] <- outputCppParseTree2(ModifiedRmmCode[[4]], paste0(' ',Indentation))
        outputCppBlock[[5]] <- paste0(Indentation, '}')
    }
    outputCppBlock
}


## new version of outputCppParseTree
outputCppParseTree2 <- function(code, indent = '') {
    if(is.call(code)) {
        if(deparse(code[[1]]) %in% names(ModifiedRmmParseKeywords2))
            return(eval(call(ModifiedRmmParseKeywords2[[ deparse(code[[1]]) ]], quote(code), quote(indent))))
        else
            return(outputGeneralCppFunction2(code, indent))
    }
    
    dpc <- deparse(code)
    if(is.name(code)) {
        return(dpc)
    }

    if(inherits(code, 'cppCodeBlock')) return(code$generate())
  ##If it is not any of these cases then just return the ModifiedRmmCode back as text 
    return(dpc)

}
