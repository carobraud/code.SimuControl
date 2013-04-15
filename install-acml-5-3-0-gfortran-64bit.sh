#!/bin/sh
#
# Installation script for ACML version 5.3.0

toolname=`basename ${0}`
# Check for any switch arguments
doneswitches=0
badswitch=0
acceptlicense=0
installdir=""
while test "${doneswitches}" = "0"
do
  case "${1-}" in
    -accept*)
       acceptlicense=1
       shift;;
    -installdir=*)
       installdir=`echo ${1} | sed -e 's%.*=%%'`
       shift;;
    -*)
       echo "${toolname}: unrecognised switch: ${1}"
       badswitch=1
       exit
       shift;;
     *)
       doneswitches=1;;
  esac
done

showLicense()
{
if [ "${acceptlicense}" = "1" ]; then
  echo "Warning: by installing this software you have accepted"
  echo "the license agreement in ACML-EULA.txt"
  reply="accept"
else
  more ACML-EULA.txt

  reply=""
  while [ "${reply}" != "accept" -a "${reply}" != "decline" ]; do
    echo -e "[accept/decline]? : \c"
    read reply
    reply=`echo ${reply} | tr [:upper:] [:lower:]`
  done
fi
}

get_yes_no()
{
reply=""
while [ "$reply" != "y" -a "$reply" != "n" ]; do
  echo -e "$1 ? [y/n] : \c"
  read reply
  reply=`echo ${reply} | tr [:upper:] [:lower:]`
done
}

echo -e "                   ACML Installation  "
echo -e "                   =================  "
echo -e ""
echo -e "This script will install ACML version 5.3.0"

showLicense

if [ "${reply}" != "accept" ]; then
  echo "Installation declined. ACML not installed."
  exit
fi

echo -e ""
echo -e "Where do you want to install ACML?  Press return to use"
echo -e "the default location (/opt/acml5.3.0), or enter an alternative path."
echo -e "The directory will be created if it does not already exist."
if [ "${installdir}" != "" ]; then
  INSTALLDIR=${installdir}
else
  INSTALLDIR=""
  while [ "${INSTALLDIR}" = "" ]; do
    echo -e "> \c"
    read ans
    if [ $ans ]
    then
        case $ans in
        *) INSTALLDIR=$ans ;;
        esac
    else
        INSTALLDIR=/opt/acml5.3.0
    fi
  done
fi

# Replace any ~ by ${HOME} otherwise you end up with a
# subdirectory named ~ (dangerous if you then try to remove it!)
INSTALLDIR=`echo ${INSTALLDIR} | sed -e "s%~%${HOME}%g"`

echo -e ""
echo -e "Installing to : ${INSTALLDIR}"
echo -e ""

if [ ! -d ${INSTALLDIR} ]
then
  mkdir -p ${INSTALLDIR}
  if [ $? -ne 0 ]
  then
    echo -e "***** Cannot create installation directory, installation failed *****"
    exit
  fi
fi

# Extract everything from the compressed tar file
fromdir=`pwd`
cd ${INSTALLDIR}
fullpathinstalldir=`pwd`
tarfile=${fromdir}/contents-acml-5-3-0-gfortran-64bit.tgz
gzip -dc ${tarfile} | tar oxvf -

# Edit any example GNUmakefiles that just got installed so that they
# point at the actual installation directory. We find the names of
# the GNUmakefiles from the compressed tar file so that we don't
# edit any GNUmakefiles that happen already to have been installed.
gnumakefiles=`gzip -dc ${tarfile} | tar tvf - | grep examples | grep -i makefile`
for gnumakefile in ${gnumakefiles}
do
  fpgnumakefile=${fullpathinstalldir}/${gnumakefile}
  if [ -f "${fpgnumakefile}" ]; then
    echo ""
    echo "Editing the example program GNUmakefile:"
    echo "  ${fpgnumakefile}"
    echo "to point at the ACML installation directory ..."
    sed -e "s%ACMLDIR \(:*\)= /opt/acml5.3.0/%ACMLDIR \1= ${fullpathinstalldir}/%" < \
         ${fpgnumakefile} > ${fpgnumakefile}.tmp
    if [ -f "${fpgnumakefile}.tmp" ]; then
      mv ${fpgnumakefile}.tmp ${fpgnumakefile}
    fi
  fi
done

echo -e ""
echo -e "====== ACML installation complete ======"
