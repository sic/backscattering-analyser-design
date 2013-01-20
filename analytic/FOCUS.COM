$ write sys$output "Compiling..."
$ f90 focus
$ write sys$output "Linking..."
$ link90 focus
