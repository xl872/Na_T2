#!/bin/bash
# Scan 613 (1H DtiStandard, 64x64x80 @ 0.10 mm, 7 b-values) -> AFNI, same pipeline as scan 599:
#   to3d -> (axis reorder) -> 3dMean (allave_t) -> 3dTstat -mean (allave)
#
# Orientation, from 613/method + pdata/1/visu_pars (identical prescription to 599):
#   PVM_SPackArrSliceOrient = sagittal, PVM_SPackArrReadOrient = H_F, RECO_transposition = 1
#   VisuCoreOrientation rows (subject coords, +x=L +y=P +z=S):
#       2dseq dim1 (64, 6.4 mm) = ( 0     1  0    )  -> A->P   (phase1)
#       2dseq dim2 (64, 6.4 mm) = (-0.055 0 -0.998)  -> S->I   (read, H_F)
#       2dseq dim3 (80, 8.0 mm) = (-0.998 0  0.055)  -> L->R   (slice pack)
#   599 was converted with generic labels -xFOV R-L -yFOV S-I -zFOV A-P (dataset RSA): there the tube
#   happened to lie along dim3, so the tube axis got the A-P label and AFNI coronal = round tube view.
#   In 613 the sample was placed 90 deg rotated: the tube lies along dim1 (true A-P).  Using the same
#   generic labels would put the tube on R-L and the round view in sagittal (that is what an earlier
#   conversion of scan 613 shows).  So here the true labels are used in to3d (dataset ASL), then the axes are
#   reordered to RSA with 3daxialize (pure permutation/flip, no interpolation) so the result has the
#   same layout as 599:  x = R-L (80), y = S-I (64), z = A-P (64, tube axis)  ->  coronal = round view,
#   and in MATLAB (BrikLoad) the round view is again data(:,:,slice).
HERE=$(cd "$(dirname "$0")" && pwd)   # release layout: paths relative to this script
RAW=$HERE/bruker/20250801_201306_Mouse_23Naplus1_08012025_1_81
OUT=$HERE/DWI_scan613
mkdir -p "$OUT"
NAME=20250801_201306_Mouse_23Naplus1_08012025_1_81.613.dif
cd $OUT

# 1) 2dseq -> AFNI.  64x64 images, 80 slices x 7 b-values = 560 images, 16-bit little-endian ("3D:" = short)
to3d -time:zt 80 7 500 zero -omri -prefix ${NAME}_asl -session $OUT -view orig \
     -xFOV 3.20A-P -yFOV 3.20S-I -zFOV 4.00L-R -skip_outliers \
     3D:0:0:64:64:560:$RAW/613/pdata/1/2dseq

# 2) reorder axes to RSA (same as the 599 dataset)
3daxialize -orient RSA -prefix $NAME ${NAME}_asl+orig
rm -f ${NAME}_asl+orig.HEAD ${NAME}_asl+orig.BRIK

# 3) same averaging steps as 599
3dMean  -prefix allave_t ${NAME}+orig
3dTstat -mean -prefix allave allave_t+orig
