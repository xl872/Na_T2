#!/bin/bash
# Scan 627 -> AFNI, same pipeline as DWI_scan613:
#   to3d -> 3daxialize (axis order/flip) -> 3dMean (allave_t) -> 3dTstat -mean (allave)
#
# Orientation, from 627/method + pdata/1/visu_pars:
#   PVM_SPackArrSliceOrient = coronal (rotated about the bore axis), PVM_SPackArrReadOrient = H_F, RECO_transposition = 1,
#   VisuCoreDiskSliceOrder = disk_reverse_slice_order
#   VisuCoreOrientation rows (subject coords, +x=L +y=P +z=S):
#       2dseq dim1  ~ towards L   (phase1)
#       2dseq dim2  ~ S->I        (read, H_F)
#       2dseq dim3  ~ A->P geometrically, but stored in REVERSE order on disk, so the
#                                  slice index k on disk runs P->A
#   The tube lies along dim3 (checked in the raw 2dseq), so with dim3 on the A-P axis the coronal
#   view is the round tube view, as in scan 599 and DWI_scan613.
#   to3d gets the true labels of the axes AS STORED (R-L, S-I, P-A -> dataset RSP), then 3daxialize
#   flips k to give RSA (x = R-L, y = S-I, z = A-P), so the physical handedness matches 613New.
#   Result: AFNI (i,j,k) = 2dseq (d1, d2, 80-d3_disk); AFNI k equals the geometric slice index.
HERE=$(cd "$(dirname "$0")" && pwd)   # release layout: paths relative to this script
RAW=$HERE/bruker/20250801_201306_Mouse_23Naplus1_08012025_1_81
OUT=$HERE/DWI_scan627
mkdir -p "$OUT"
NAME=20250801_201306_Mouse_23Naplus1_08012025_1_81.627.dif
cd $OUT

# 1) 2dseq -> AFNI, 16-bit little-endian ("3D:" = short)
to3d -time:zt 81 8 500 zero -omri -prefix ${NAME}_rsp -session $OUT -view orig \
     -xFOV 3.00R-L -yFOV 3.00S-I -zFOV 4.05P-A -skip_outliers \
     3D:0:0:60:60:648:$RAW/627/pdata/1/2dseq

# 2) flip the slice axis to RSA (same layout as DWI_scan613)
3daxialize -orient RSA -prefix $NAME ${NAME}_rsp+orig
rm -f ${NAME}_rsp+orig.HEAD ${NAME}_rsp+orig.BRIK

# 3) same averaging steps as 599 / 613New
3dMean  -prefix allave_t ${NAME}+orig
3dTstat -mean -prefix allave allave_t+orig
