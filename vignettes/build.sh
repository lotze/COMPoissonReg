#!/bin/bash

# This bash script compiles any vignettes present in the `vignettes` folder of
# the project, relative to the current working directory.
#
# You can also run the lines below within R from any platform to get the same
# effect.


read -r -d '' SCRIPT << EOF
require("R.rsp")
devtools::document(roclets = "vignette")
EOF

Rscript -e "${SCRIPT}"

