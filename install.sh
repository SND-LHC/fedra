#!/bin/bash
# Installation script for the fedra-package   #
# 14.11.03                                    #

if [[ -z "${FEDRA_INSTALL_DIR}" ]]; then
 installdir=`pwd`
else
 installdir="${FEDRA_INSTALL_DIR}"
 scriptdir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
 cp -r ${scriptdir}/* ${installdir}
fi
if [[ ${installdir##*/} == "src" ]] ; then
 installdir=`cd ..; pwd`
fi
export installdir=$installdir

# making directories, if not already existing
if [[ ! -a $installdir/bin ]] ; then
  mkdir $installdir/bin
fi
if [[ ! -a $installdir/include ]] ; then
  mkdir $installdir/include
fi
if [[ ! -a $installdir/lib ]] ; then
  mkdir $installdir/lib
fi

# Define links from include-dir to src-files; former "setlinks.sh"
PROJECT_INC=$installdir/include
PROJECT_SRC=$installdir/src

LIBS="libEdb libEbase libERootTools libDataConversion libEmath libEGA libEphys libEdr libEIO libEdd libEMC libAlignment libScan libShower libEOracle libEmr libEDA libMosaic libShowRec"

for lib in ${LIBS} ; do
    echo "set links for ${lib} ..."
    listLib=`ls $installdir/src/${lib}/*.h`
    for hname in ${listLib} ; do
	ln -fs $PROJECT_SRC/${lib}/${hname##*/} $PROJECT_INC/${hname##*/}
    done
done
echo "set links for libVt++ ..."
ln -fs $PROJECT_SRC/libVt++/vt++/include $PROJECT_INC/vt++
ln -fs $PROJECT_SRC/libVt++/smatrix/include $PROJECT_INC/smatrix
echo "set dataIO link for libDataConversion ..."
ln -fs $PROJECT_SRC/libDataConversion/dataIO $PROJECT_INC/dataIO

chmod +x $PROJECT_SRC/libVt++/smatrix/CreateBinaryOp.sh
chmod +x $PROJECT_SRC/libVt++/smatrix/CreateUnaryOp.sh
chmod +x $PROJECT_SRC/libVt++/version.sh

# create ProjectDef.mk in ./config-directory
ProjectDef="$installdir/src/config/ProjectDef.mk"

echo "PROJECT_ROOT=$installdir"  > $ProjectDef
echo >> $ProjectDef
echo 'BIN_DIR = $(PROJECT_ROOT)/bin' >> $ProjectDef
echo 'LIB_DIR = $(PROJECT_ROOT)/lib' >> $ProjectDef
echo 'INC_DIR = $(PROJECT_ROOT)/include' >> $ProjectDef
echo >> $ProjectDef
#echo 'PROJECT_LIBS = -L$(LIB_DIR)' >> $ProjectDef

cat $ProjectDef
# exit 0

# create RootDef.mk and TargetsDef.mk links
configdir="$installdir/src/config"
ln -fs $configdir/RootDef.linux.mk $configdir/RootDef.mk 
ln -fs $configdir/TargetsDef.linux.mk $configdir/TargetsDef.mk 

ln -fs $installdir/src/libScan/Makefile.linux $installdir/src/libScan/Makefile 

# create setup_new.sh, depending on default-shell used
usedshell=`echo $SHELL`
if [[ (${usedshell##*/} == "tcsh") ||  (${usedshell##*/} == "csh") ]] ; then 
 echo 'setenv  FEDRA_ROOT '"$installdir" > $installdir/setup_new.sh
 echo 'setenv  LD_LIBRARY_PATH $FEDRA_ROOT/lib:${LD_LIBRARY_PATH}' >> $installdir/setup_new.sh
 echo 'setenv  PATH ${PATH}:$FEDRA_ROOT/bin' >> $installdir/setup_new.sh
 echo 'setenv  PYTHONPATH ${PYTHONPATH}:$FEDRA_ROOT/python' >> $installdir/setup_new.sh
fi

if [[ (${usedshell##*/} == "bash") ||  ${usedshell##*/} == "sh" || ${usedshell##*/} == "ksh" ]] ; then
 echo 'export FEDRA_ROOT='"$installdir" > $installdir/setup_new.sh
 echo 'export LD_LIBRARY_PATH=$FEDRA_ROOT/lib:${LD_LIBRARY_PATH}' >> $installdir/setup_new.sh
 echo 'export PATH=$FEDRA_ROOT/bin:${PATH}' >> $installdir/setup_new.sh
 echo 'export PYTHONPATH=${PYTHONPATH}:$FEDRA_ROOT/python' >> $installdir/setup_new.sh
fi
export FEDRA_ROOT="$installdir"
export LD_LIBRARY_PATH=$FEDRA_ROOT/lib:${LD_LIBRARY_PATH}
export PATH=${PATH}:$FEDRA_ROOT/bin

# compilation of libraries
#echo "Do you want to compile the libraries: [y/n]"
#read yesno
#if [[ $yesno == 'y' ]] ; then
  chmod u+x $installdir/src/makeall.sh
  cd $installdir/src
  $installdir/src/makeall.sh
#fi

# install fedra macros
cp -rf $PROJECT_SRC/appl/macros $installdir/macros
mv $PROJECT_SRC/*/*.pcm $installdir/lib
mv $PROJECT_SRC/*/*/*.pcm $installdir/lib
cd $installdir

# install config files

#test -d $FEDRA_ROOT/config || mkdir $FEDRA_ROOT/config
#cp  $PROJECT_SRC/appl/config/*cfg $FEDRA_ROOT/config
#cd $installdir


