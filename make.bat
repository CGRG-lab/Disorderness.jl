@ECHO OFF

if "%1" == "" goto help

if "%1" == "help" (
	:help
	echo.Please use `make ^<target^>` where ^<target^> is one of
	echo.  instantiate     Download all the packages declared in that manifest.toml
	goto end
)

if "%1" == "instantiate" (
    julia --project=. --color=yes -E "using Pkg; Pkg.instantiate();"
    goto end
)

:end