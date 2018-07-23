#' Latex to R markdown converter
#' 
#' @description Does 90 percent of the job, the rest is manual edits.
#' 
#' @param input Character Filename of input file
#' @param output character Filename of output file
#' 
#' @details Doesn't insert figures from external sources.
#' 
#' @author Loic Dutrieux, modified by Chris Paciorek
#' 
#' @import stringr
#'
#'
library(stringr)

## check html output; title page for html output

## check new version of eigen API and if tables for eigen are accurate

## % mcmc-litters chunk issue in chapter_MCMC

# for \n included with . could do regex("<<.*>>", dotall=TRUE)

latex2Rmd <- function(input, save = TRUE) {
    inputFile <- paste0(input, ".Rnw")
    text <- readChar(inputFile, file.info(inputFile)$size)

    ## Order these in likelihood of nestedness
    
#   text <- str_replace_all(string = text, pattern = ",results='markup;\\$=\\\\\\$'", '') # can't get it to deal with \$ correctly
    text <- str_replace_all(string = text, pattern = ",results='markup;.*?'", '')
                                       
    text <- str_replace_all(string = text, pattern = "^\\s*%(.*?)\\n", "<!--- \\1 -->\n")
    text <- str_replace_all(string = text, pattern = "\\n\\s*%(.*?)\\n", "\n<!--- \\1 -->\n")
    text <- str_replace_all(string = text, pattern = '<<([\\s\\S]*?)>>=([\\s\\S]*?)@', replacement = '```{r, \\1}\\2```')
    text <- str_replace_all(string = text, pattern = '\\\\href\\{([\\s\\S]*?)\\}\\{([\\s\\S]*?)\\}', replacement = '[\\2](\\1)')
    text <- str_replace_all(string = text, pattern = '\\\\url\\{([\\s\\S]*?)\\}', replacement = '[\\1]')
    text <- str_replace_all(string = text, pattern = '\\\\(texttt|code|file)\\{(.*?)\\}', replacement = '`\\2`')
    text <- str_replace_all(string = text, pattern = '\\\\(textit|emph|nm)\\{(.*?)\\}', replacement = '*\\2*')
   
    text <- str_replace_all(string = text, pattern = '\\\\chapter\\*?\\{(.*?)\\}', replacement = '# \\1') 
    text <- str_replace_all(string = text, pattern = '\\\\section\\*?\\{(.*?)\\}', replacement = '## \\1')
    text <- str_replace_all(string = text, pattern = '\\\\subsection\\*?\\{(.*?)\\}', replacement = '### \\1')
    text <- str_replace_all(string = text, pattern = '\\\\subsubsection\\*?\\{(.*?)\\}', replacement = '#### \\1')
    text <- str_replace_all(string = text, pattern = '\\\\paragraph\\*?\\{(.*?)\\}', replacement = '##### \\1')

    text <- str_replace_all(string = text, pattern = '\\n\\\\label\\*?\\{(.*?)\\}', replacement = ' \\{#\\1\\}')
    text <- str_replace_all(string = text, pattern = '\\\\ref\\*?\\{(.*?)\\}', replacement = '\\\\@ref\\(\\1\\)')
    
    text <- str_replace_all(string = text, pattern = '\\\\cd\\*?\\{([\\s\\S]*?)\\}', replacement = '`\\1`') 
    text <- str_replace_all(string = text, pattern = '\\\\begin\\{(.*?)\\}', replacement = '')
    text <- str_replace_all(string = text, pattern = '\\\\end\\{(.*?)\\}', replacement = '')
    text <- str_replace_all(string = text, pattern = '\\\\item', replacement = '  1.')

    text <- str_replace_all(string = text, pattern = '\\\\verb\\|(.*?)\\|', replacement = '`\\1`')

    text <- str_replace_all(string = text, pattern = '\\\r', replacement = '')
    text <- str_replace_all(string = text, pattern = '\\\\ver', replacement = '`r version`')

    text <- str_replace_all(string = text, pattern = '\\\\footnote\\{([\\s\\S]*?)\\}', replacement = '^[\\1]')
    text <- str_replace_all(string = text, pattern = '\\\\Sexpr\\{(.*?)\\}', replacement = '')
    text <- str_replace_all(string = text, pattern = "``(.*?)''", "'\\1'")
    text <- str_replace_all(string = text, pattern = '\\\\_', '_')
    text <- str_replace_all(string = text, pattern = '\\\\\\$', '$')
    text <- str_replace_all(string = text, pattern = '\\\\\\\\', '\n')
    text <- str_replace_all(string = text, pattern = '\\.tex', '\\.md')
    text <- str_replace_all(string = text, pattern = '\\$times\\$', ' $ \times $ ')
    text <- str_replace_all(string = text, pattern = '\\\\cite[pt]\\*?\\{(.*?)\\}', replacement = '[@\\1]')

    text <- str_replace_all(string = text, pattern = '\\\\input\\{(.*?)\\}', "```{r, child = '\\1'}\n```")

    if(save)  writeLines(text, paste0(input, ".Rmd"), sep = '')
}
