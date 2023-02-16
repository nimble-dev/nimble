
#' EXPERIMENTAL: Turn a function into a model macro builder
#' A model macro expands one line of code in a nimbleModel into one or
#' more new lines.  This supports compact programming by defining
#' re-usable modules.  \code{model_macro_builder} takes as input a
#' function that constructs new lines of model code from the original
#' line of code.  It returns a function suitable for internal use by
#' \code{nimbleModel} that arranges arguments for input function.  Macros
#' are an experimental feature and are available only after setting
#' \code{nimbleOptions(enableModelMacros = TRUE)}.
#'
#' @param fun A function written to construct new lines of model code.
#'
#' @param use3pieces (TRUE or FALSE) Should the arguments from the input
#' line be split into pieces for the LHS (left-hand side), RHS
#' (right-hand side, possibly further split depending on
#' \code{unpackArgs}), and \code{stoch} (TRUE if the line uses a
#' \code{~}, FALSE otherwise)?  See details and examples.
#'
#' @param unpackArgs (TRUE or FALSE) Should arguments be passed as a list
#' (FALSE) or as separate arguments (TRUE)?  See details and examples.
#'
#' @details The arguments \code{use3pieces} and \code{unpackArgs}
#' indicate how \code{fun} expects to have arguments arranged from an
#' input line of code (processed by \code{nimbleModel}).
#'
#' Consider the defaults \code{use3pieces = TRUE} and \code{unpackArgs =
#' TRUE}, for a macro called \code{macro1}.  In this case, the line of
#' model code \code{x ~ macro1(arg1 = z[1:10], arg2 = "hello")} will be
#' passed to \code{fun} as \code{fun(stoch = TRUE, LHS = x, arg1 =
#' z[1:10], arg2 = "hello")}.
#'
#' If \code{use3pieces = TRUE} but \code{unpackArgs = FALSE}, then the
#' RHS will be passed as is, without unpacking its arguments into
#' separate arguments to \code{fun}.  In this case, \code{x ~ macro1(arg1
#' = z[1:10], arg2 = "hello")} will be passed to \code{fun} as
#' \code{fun(stoch = TRUE, LHS = x, RHS = macro1(arg1 = z[1:10], arg2 =
#' "hello"))}.
#'
#' If \code{use3pieces = FALSE} and \code{unpackArgs = FALSE}, the entire
#' line of code is passed as a single object.  In this case, \code{x ~
#' macro1(arg1 = z[1:10], arg2 = "hello")} will be passed to \code{fun}
#' as \code{fun(x ~ macro1(arg1 = z[1:10], arg2 = "hello"))}.  It is also
#' possible in this case to pass a macro without using a \code{~} or
#' \code{<-}.  For example, the line \code{macro1(arg1 = z[1:10], arg2 =
#' "hello")} will be passed to \code{fun} as \code{fun(macro1(arg1 =
#' z[1:10], arg2 = "hello"))}.
#'
#' If \code{use3pieces = FALSE} and \code{unpackArgs = TRUE}, it
#' won't make sense to anticipate a declaration using \code{~} or \code{<-}.  Ins#' tead, arguments from an arbitrary call will be passed as separate arguments.  #' For example, the line \code{macro1(arg1 = z[1:10], arg2 = "hello")} will be pa#' ssed to \code{fun} as \code{fun(arg1 = z[1:10], arg2 = "hello")}.
#'
#' It is extremely useful to be familiar with processing R code as an
#' object to write \code{fun} correctly.  Functions such as
#' \code{\link{substitute}} and \code{\link{as.name}}
#' (e.g. \code{as.name('~')}), \code{\link{quote}}, \code{\link{parse}}
#' and \code{\link{deparse}} are particularly handy.
#'
#' Multiple lines of new code should be contained in \code{ {} }. Extra
#' curly braces are not a problem. See example 2.
#'
#' Macro expansion is done recursively: One macro can return code that
#' invokes another macro.
#' 
#' @return A list with a named element \code{code} that contains the
#' replacement code.
#'
#' @export
#' 
#' @examples
#' nimbleOptions(enableModelMacros = TRUE)
#' nimbleOptions(verbose = FALSE)
#' 
#' ## Example 1: Say one is tired of writing "for" loops.
#' ## This macro will generate a "for" loop with dnorm declarations
#' all_dnorm <- model_macro_builder(
#'     function(stoch, LHS, RHSvar, start, end, sd = 1) {
#'         newCode <- substitute(
#'             for(i in START:END) {
#'                 LHS[i] ~ dnorm(RHSvar[i], SD)
#'             },
#'             list(START = start,
#'                  END = end,
#'                  LHS = LHS,
#'                  RHSvar = RHSvar,
#'                  SD = sd))
#'         list(code = newCode)
#'     },
#'     use3pieces = TRUE,
#'     unpackArgs = TRUE 
#' )
#' 
#' model1 <- nimbleModel(
#'     nimbleCode(
#'     {
#'         ## Create a "for" loop of dnorm declarations by invoking the macro
#'         x ~ all_dnorm(mu, start = 1, end = 10)
#'     }
#'     ))
#' 
#' ## show code from expansion of macro
#' model1$getCode()
#' ## The result should be:
#' ## {
#' ##     for (i in 1:10) {
#' ##         x[i] ~ dnorm(mu[i], 1)
#' ##     }
#' ## }
#' 
#' ## Example 2: Say one is tired of writing priors.
#' ## This macro will generate a set of priors in one statement
#' flat_normal_priors <- model_macro_builder(
#'     function(...) {
#'         allVars <- list(...)
#'         priorDeclarations <- lapply(allVars,
#'                                     function(x)
#'                                         substitute(VAR ~ dnorm(0, sd = 1000),
#'                                                    list(VAR = x)))
#'         newCode <- quote({})
#'         newCode[2:(length(allVars)+1)] <- priorDeclarations
#'         list(code = newCode)
#'     },
#'     use3pieces = FALSE,
#'     unpackArgs = TRUE
#' )
#' 
#' model2 <- nimbleModel(
#'     nimbleCode(
#'     {
#'         flat_normal_priors(mu, beta, gamma)
#'     }
#'     ))
#' 
#' ## show code from expansion of macro
#' model2$getCode()
#' ## The result should be:
#' ## {
#' ##     {
#' ##         mu ~ dnorm(0, sd = 1000)
#' ##         beta ~ dnorm(0, sd = 1000)
#' ##         gamma ~ dnorm(0, sd = 1000)
#' ##     }
#' ## }
#' ## Extra curly braces do not matter.
model_macro_builder <- function(fun,
                                use3pieces = TRUE,
                                unpackArgs = TRUE ) {
    if(use3pieces) {
        wrapper <- function(code, .constants, .env) {
            args <- as.list(code)
            args[[1]] <- args[[1]] == '~'
            names(args) <- c('stoch','LHS','RHS')
            if(unpackArgs) {
                RHSargs <- as.list(args[[3]])[-1]
                args <- c(args[1:2], RHSargs)
            }
            ## Since code contains unevaluated expressions, they need
            ## to be wrapped in another layer of quote() that do.call
            ## will strip.  Otherwise there will be an attempted
            ## evaluation of variable names.
            args <-  lapply(args,
                            function(x)
                                substitute(quote(X), list(X = x)))
            do.call(fun, c(args, list(.constants = .constants, .env = .env)))
        }
    } else {
        if(unpackArgs) {
            wrapper <- function(code, .constants, .env) {
                args <- as.list(code)[-1]
                args <-  lapply(args,
                            function(x)
                                substitute(quote(X), list(X = x)))
                do.call(fun, c(args, list(.constants = .constants, .env = .env)))
            }
        } else {
            wrapper <- fun
        }
    }
    ans <- structure(list(process = wrapper), class = "model_macro")
    ans
}

