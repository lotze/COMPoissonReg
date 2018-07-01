equitest <- function(object, ...)
{
	UseMethod("equitest")
}

leverage <- function(object, ...)
{
	UseMethod("leverage")
}

nu <- function(object, ...)
{
	UseMethod("nu")
}

sdev <- function (object, ...)
{
	UseMethod("sdev")
}

parametric_bootstrap <- function(object, reps = 1000, report.period = reps + 1, ...)
{
	UseMethod("parametric_bootstrap")
}
