SHELL=/bin/sh
set -x

##################################################################
# wafs using module compile standard
# 06/12/2018 yali.ma@noaa.gov:    Create module load version
##################################################################

module purge
set -x
mac=$(hostname | cut -c1-1)
mac2=$(hostname | cut -c1-2)
if [ $mac = v -o $mac = m  ] ; then            # For Dell
 machine=dell
 . $MODULESHOME/init/bash                 
elif [ $mac2 = hf ] ; then                        # For Hera
 machine=hera
 . /etc/profile
 . /etc/profile.d/modules.sh
elif [ $mac = O ] ; then           # For Orion
 machine=orion
 . /etc/profile
fi

moduledir=`dirname $(readlink -f ../modulefiles/wafs)`
module use ${moduledir}
module load wafs/wafs_v7.0.0-${machine}
module list

 curdir=`pwd`
 export INC="${G2_INC4}"
 export FC=ifort

# track="-O3 -g -traceback -ftrapuv -check all -fp-stack-check "
# track="-O2 -g -traceback"

 export FFLAGSawc="-FR -I ${G2_INC4} -I ${IP_INC4} -O2 -convert big_endian -check bounds -traceback -g -assume noold_ldout_format"
 export FFLAGSblnd="-O -I ${G2_INC4} -g -check bounds -traceback"
 export FFLAGST="-O -FR -I ${G2_INC4} -g -check bounds -traceback"
 export FFLAGSgcip="-FR -I ${G2_INC4} -I ${IP_INC4} -g -O3"
# export FFLAGSgcip="-FR -I ${G2_INC4} -I ${IP_INC4} ${track}"

 export FFLAGScnv="-O3 -g -I ${G2_INC4}"
 export FFLAGSmkwfs="-O3 -g -r8 -i8"

if [ ! -d "../exec" ] ; then
  mkdir -p ../exec
fi

for dir in wafs_awc_wafavn.fd wafs_gcip.fd wafs_blending.fd wafs_makewafs.fd wafs_cnvgrib2.fd wafs_setmissing.fd wafs_grib2_0p25.fd wafs_blending_0p25.fd ; do
#for dir in wafs_blending.fd ; do 
#for dir in wafs_blending_0p25.fd ; do 
#for dir in wafs_gcip.fd ; do
#for dir in wafs_makewafs.fd ; do
#for dir in wafs_awc_wafavn.fd ; do
#for dir in wafs_grib2_0p25.fd wafs_blending_0p25.fd ; do
 if [ $dir = "wafs_makewafs.fd" ] ; then
     export LIBS="${W3NCO_LIB8}  ${W3EMC_LIB8} ${BACIO_LIB8}"
 else
     export LIBS="${G2_LIB4} ${W3NCO_LIB4} ${BACIO_LIB4} ${IP_LIB4} ${SP_LIB4} ${JASPER_LIBRARIES}/libjasper.a ${PNG_ROOT}/lib64/libpng.a ${ZLIB_LIBRARIES}/libz.a  ${BUFR_LIB4}"
 fi
 cd ${curdir}/$dir
 make clean
 make
 make install
 make clean
done
