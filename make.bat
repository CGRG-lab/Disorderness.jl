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

if "%1" == "pluto" (
	if "%2" == "sde" (
		@REM Before Pluto releases a new version, we should reload `Pluto.open_in_default_browser`.
		@REM Original issue: https://github.com/fonsp/Pluto.jl/issues/936
		start julia --project=. -e "using Pluto; Pluto.open_in_default_browser(url::AbstractString)::Bool = (try (Base.run(`powershell.exe Start \"'$url'\"`); true) catch ex false end); Pluto.run(Pluto.Configuration.Options(server=Pluto.Configuration.ServerOptions(notebook=\"Pluto_SDE.jl\")))"
		exit
		goto end
	)
)

:end