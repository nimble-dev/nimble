makeEnvName <- function(name) paste0('.env_',name) ## already in modelBUGS.R for old system
makeNameName <- function(name) paste0('.name_',name)
makeRowName <- function(name) paste0('.row_', name)

makeLogProbName <- function(name) paste0('logProb_', name)
makeLogProbEnvName <- function(name) makeEnvName(makeLogProbName(name))

makeVecName <- function(name){
	if(length(name) == 0)	return(name)		#This was causing a problem when we have blank modelValues (i.e. MCMCconf if we have no thinning monitors)
	paste0(name,'_Vec')
}
removeLogProbName <- function(name)  sub('^logProb_', '', name)

is.indexed <- function(node) {
    if(is.character(node)) return(grepl('\\[', node))
    else return('[' %in% all.names(node))
}

is.vectorized <- function(node) {
    if(is.character(node)) return(grepl('\\:', node))
    else return(':' %in% all.names(node))
}


is.blank <- function(arg) {
    if(is.null(arg)) return(FALSE)
    return(identical(arg, quote(x[])[[3]]))
}

removeIndexing <- function(nodes) {
    return(gsub('\\[.*', '', nodes))
}

## this utility function is used in the setup code of conjugate sampler functions.
## the 'd' dimension variables used to come from nodeInfo object,
## but now we calculate it (from the targetNode name) using this function.
determineNodeIndexSizes <- function(node) {
    if(!is.indexed(node)) return(numeric(0))
    indString <- gsub('.*\\[', '', gsub('\\]', '', node))
    indStringList <- as.list(strsplit(indString, ', ')[[1]])
    indExprList <- lapply(indStringList, function(ind) parse(text=ind)[[1]])
    sizes <- unlist(lapply(indExprList, function(ind) max(eval(ind)) - min(eval(ind)) + 1))
    return(sizes)
}

## This should add model$ in front of any names that are not already part of a '$' expression
addModelDollarSign <- function(expr, exceptionNames = character(0)) {
    if(is.numeric(expr)) return(expr)
    if(is(expr, 'srcref')) return(expr)
    if(is.name(expr)) {
        if(as.character(expr) %in% exceptionNames)    return(expr)
        proto <- quote(model$a)
        proto[[3]] <- expr
        return(proto)
    }
    if(is.call(expr)) {
        if(expr[[1]] == '$'){
          if(length(expr[[2]] > 0)){
            expr[2] <- lapply(expr[2], function(listElement) addModelDollarSign(listElement, exceptionNames))
            return(expr)
          } 
          else{
            return(expr)
          }
        }
        if(expr[[1]] == 'returnType')
            return(expr)
        if(length(expr) > 1) {
            expr[2:length(expr)] <- lapply(expr[-1], function(listElement) addModelDollarSign(listElement, exceptionNames))
            return(expr)
        }
    }
    return(expr)
}

# Determine if a piece of code contains a '['
hasBracket <- function(code) {
    if(length(code) < 2) return(FALSE)
    if(code[[1]] == '[') return(TRUE)
    FALSE
}


evalBracketArgs <- function(code, constantEnv) {
    if(hasBracket(code)) {
        for(i in 3:length(code)) {
            if(!is.vectorized(code[[i]])) {
                ## no vectorization; just evaluate it
                code[[i]] <- as.numeric(eval(code[[i]], constantEnv))
            } else {
                ## vectorized index. evaluate it, then set it to the resulting expression:  MIN:MAX
                indicies <- as.numeric(eval(code[[i]], constantEnv))
                code[[i]] <- substitute(MIN:MAX, list(MIN=min(indicies), MAX=max(indicies)))
            }
        }
    }
    return(code)
}

evalBracketArgsKnownBracket <- function(code, constantEnv, isVectorized) {
    ## Same as above but when hasBracket(code) is known to be true
    for(i in 3:length(code)) {
        if(!isVectorized[i]) {
            ## no vectorization; just evaluate it
            code[[i]] <- as.numeric(eval(code[[i]], constantEnv))
        } else {
                ## vectorized index. evaluate it, then set it to the resulting expression:  MIN:MAX
            indicies <- as.numeric(eval(code[[i]], constantEnv))
            code[[i]] <- substitute(MIN:MAX, list(MIN=min(indicies), MAX=max(indicies)))
        }
    }
    code
}



exprAsListDrop2 <- function(expr) {
    lem2 <- length(expr) - 2
    if(lem2 < 1) return(NULL)
    ans <- switch(lem2,
                  list(expr[[3]]),
                  list(expr[[3]], expr[[4]]),
                  list(expr[[3]], expr[[4]], expr[[5]]),
                  list(expr[[3]], expr[[4]], expr[[5]], expr[[6]])
                  )
    if(is.null(ans)) as.list(expr[-c(1,2)]) else ans
}

getCallText <- function(code)      as.character(code[[1]])

parseTreeSubstitute <- function(pt, pattern, replacement) {
    if(identical(pt, pattern))    return(replacement)
    if(class(pt) == 'name')       return(pt)
    if(is.numeric(pt))            return(pt)
    if(is.logical(pt))            return(as.numeric(pt)) # for now turn logicals into numerics until we propagate logicals to C++
    for(i in seq_along(pt))    { pt[[i]] <- parseTreeSubstitute(pt[[i]], pattern, replacement) }
    return(pt)
}

# no longer documented in Rd
# Generates a valid C++ name from an R Name
# 
# replaces [ ( $ and a few other symbols with underscores, and removes ] ) and spaces in a string
# 
# @param rName A String
# @return returns a string representing the modified rName
# @author Jagadish Babu
# @keywords Name
# @seealso \code{\link{character}} 
# @export
# @examples
#  genName('theta[1]')
Rname2CppName <- function(rName, colonsOK = TRUE) {
    ## This will serve to replace and combine our former Rname2CppName and nameMashupFromExpr
    ## which were largely redundant
    if (!is.character(rName))
        rName <- deparse(rName)
    if( colonsOK)
        rName <- gsub(':', 'to', rName)  # replace colons with 'to'
    else
        if(grepl(':', rName))
            stop(paste0('trying to do name mashup on expression with colon (\':\') from ', rName))
    rName <- gsub(' ', '', rName)
    rName <- gsub('\\.', 'p', rName) 
    rName <- gsub("\"", "", rName)
    rName <- gsub(',', '_', rName)   
    rName <- gsub("`", "" , rName)
    rName <- gsub('\\[', '_oB', rName)
    rName <- gsub('\\]', '_cB', rName)
    rName <- gsub('\\(', '_oP', rName)
    rName <- gsub('\\)', '_cP', rName)
    rName <- gsub("\\$", "_" , rName)
    rName <- gsub(">=", "_gte_", rName)
    rName <- gsub("<=", "_lte_", rName)
    rName <- gsub("<=", "_eq_", rName)
    rName <- gsub("!=", "_neq_", rName)
    rName <- gsub(">", "_gt_", rName)
    rName <- gsub("<", "_lt_", rName)
    rName <- gsub("!", "_not_", rName)
    rName <- gsub("\\|\\|", "_or2_", rName)
    rName <- gsub("&&", "_and2_", rName)
    rName <- gsub("\\|", "_or_", rName)
    rName <- gsub("&", "_and_", rName)
    rName <- gsub("%%", "_mod_", rName)
    rName <- gsub("%\\*%", "_matmult_", rName)
    rName <- gsub("=", "_" , rName)
    rName <- gsub("\\(", "_" , rName)
    rName <- gsub("\\+", "_plus_" , rName)
    rName <- gsub("-", "_minus_" , rName)
    rName <- gsub("\\*", "_times_" , rName)
    rName <- gsub("/", "_over_" , rName)
    rName <- gsub('\\^', '_tothe_', rName)
    rName <- gsub('^_+', '', rName) # remove leading underscores.  can arise from (a+b), for example
    rName <- gsub('^([[:digit:]])', 'd\\1', rName)    # if begins with a digit, add 'd' in front
    rName
    
}


vectorIndex_2_flat <- function(index, strides){
	return(sum((index-1) * strides) + 1)
}

extractFlatIndices_wVarInfo <- function(nodeNames, varInfo){
	varName <- varInfo[['varName']]
	numNodes <- length(nodeNames)
	firstDropNumber <- rep(nchar(varName) + 2, numNodes) 
	lastDropNumber <- nchar(nodeNames) - 1
	charIndices_wCommas <- substring(nodeNames, firstDropNumber, lastDropNumber)
	charIndices_nCommas <- strsplit(charIndices_wCommas, ',') 
	numIndices <- lapply(charIndices_nCommas, as.numeric)
	if(varInfo$nDim == 1)
		return(unlist(numIndices) )
	strides <- rep(1, varInfo$nDim)
	for(i in 2:(varInfo$nDim) ){
		strides[i] = strides[i-1] * varInfo$maxs[i-1]
		}
	flatIndices <- unlist(lapply(numIndices, vectorIndex_2_flat, strides = strides) )
	return(flatIndices)
}

expandIndexSet4sapply <- function(nextIndices, prevIndices, stride)
	return(prevIndices + (nextIndices-1) * stride)

combineIndices2Flat <- function(prevIndices, nextIndices, stride)
	as.numeric(sapply(nextIndices, expandIndexSet4sapply, prevIndices, stride) ) 
	
character2index <- function(thisChar){
	splitValues <- strsplit(thisChar, split = ":")[[1]]
		if(length(splitValues) == 1)
		return(as.numeric(splitValues) ) 
	if(length(splitValues)){
		begin = as.numeric(splitValues[1])
		end = as.numeric(splitValues[2])
		return(begin:end)
	}
	else
		stop("Error: too many :'s in index")
}

# extracts dimension from character vec of form such as c("double(0)", "integer(1)")
## getDimFromType <- function(text) {
##     sapply(parse(text = text), '[[', 2)
## }

