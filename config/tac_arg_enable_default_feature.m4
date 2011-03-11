dnl @synopsis TAC_ARG_ENABLE_DEFAULT_FEATURE(FEATUREE_NAME, FEATURE_DESCRIPTION, HAVE_NAME)
dnl
dnl Test for --enable-${FEATUER_NAME} and set to yes
dnl (unless --disable-default-packages is specified) if not specified.
dnl Also calls AC_DEFINE to define HAVE_${HAVE_NAME} if value is not equal to "no"
dnl 
dnl Use this macro to help define whether or not packages that are built by 
dnl default (unless --disable-default-packages is specified) should be built.
dnl  For example:
dnl
dnl TAC_ARG_ENABLE_DEFAULT_FEATURE(epetra, [Configure and build epetra], EPETRA)
dnl 
dnl will test for --enable-epetra when configure is run.  If it is defined 
dnl and not set to "no" or not defined and --disable-default-packages is not 
dnl specified then HAVE_EPETRA will be defined, if --enable-epetra is not 
dnl defined to be "yes" and --disable-default-packages is specified, 
dnl HAVE_EPETRA will not be defined.
dnl
dnl *NOTE: epetra, aztecoo, komplex, ifpack, and other software found in
dnl subdirectories of Trilinos/packages are "packages" in their own right.
dnl However, these packages are also "features" of the larger package
dnl "Trilinos".  Therefore, when configuring from the Trilinos directory,
dnl it is appropriate to refer to these software packages as "features".
dnl
dnl This file was based on tac_arg_enable_feature.m4
dnl @author James Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_ENABLE_DEFAULT_FEATURE],
[
AC_ARG_ENABLE([$1],
AC_HELP_STRING([--enable-$1],[$2 (default is yes unless the --disable-default-packages option is used]),
ac_cv_use_$1=$enableval
ac_cv_use_$1_explicit=$enableval,
ac_cv_use_$1=$ac_cv_use_default_packages
ac_cv_use_$1_explicit=no)

AC_MSG_CHECKING(whether to use [$1])

if test "X$ac_cv_use_$1" != "Xno"; then
  AC_MSG_RESULT(yes)  
  AC_DEFINE([HAVE_$3],,[Define if want to build $1])
else
  AC_MSG_RESULT(no)
fi
])

