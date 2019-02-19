#!/bin/bash

# This script just runs some of the driver programs to verify the build.
./drivers/print_gauss 10
./drivers/print_gauss_lobatto 10
./drivers/jacobi_rule 10
./drivers/conical_product_2D 3
./drivers/conical_product_3D 3
./drivers/gm_rule -s 2
./drivers/dubiner_verify 5
./drivers/lu --matrix-size 10
./drivers/pyramid_conical_product 4
