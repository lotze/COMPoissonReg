PROJNAME := COMPoissonReg
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_path := $(dir $(mkfile_path))

## See http://kbroman.org/pkg_primer/pages/docs.html
## for help on generating Roxygen2 documentation

pdf:
	R CMD Rd2pdf --force --no-preview --output=$(PROJNAME).pdf .

#manual:
#	R -e 'library(devtools); document(roclets = c("collate", "rd"))'

clean:
	rm -rf man $(PROJNAME).pdf src/$(PROJNAME).so
	rm -rf autom4te.cache config.log config.status src/Makevars
	find $(current_path)/src -type f -name '*.o' -exec rm {} \;

