#!/bin/bash
#
FC=gfortran
DEFINES='-D_BYTESWAPIO -DWORDLENGTH=4 -DNML_EXTENDED_F77'
CPP='cpp  -traditional -P'
NOOPTFLAGS='-O0'
EXTENDED_SRC_FLAG='-ffixed-line-length-132'
GET_FC_VERSION="--version"

#  For IEEE, use the "-ffloat-store" option
if test "x$IEEE" = x ; then
    FFLAGS='-Wimplicit -Wunused -Wuninitialized'
    FOPTIM='-O3 -funroll-loops'
else
    FFLAGS='-Wimplicit -Wunused -ffloat-store'
    FOPTIM='-O0'
fi



