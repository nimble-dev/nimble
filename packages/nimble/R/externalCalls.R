## system to allow compiled nimbleFunctions to call semi-arbitrary external compiled (or compilable) code

externalCallClass <- setRefClass(
    fields = list(
        hFile = 'ANY',
        cFile = 'ANY',
        oFile = 'ANY',
        returnType = 'ANY',
        argTypes = 'ANY'
    )
)

