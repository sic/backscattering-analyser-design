$ write sys$output "Compiling..."
$ f90 hybrid
$ write sys$output "Linking..."
$ link90 hybrid
