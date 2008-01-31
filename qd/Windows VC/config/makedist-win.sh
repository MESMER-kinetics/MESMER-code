#!/bin/bash
PACKAGE_NAME=qd
MAJOR_VERSION=2
MINOR_VERSION=3
PATCH_LEVEL=$1
if [[ -z "$PATCH_LEVEL" ]]; then
  echo "Usage: makedist-win.sh patch-level"
  exit
fi

VERSION=$MAJOR_VERSION.$MINOR_VERSION.$PATCH_LEVEL
echo "Creating $PACKAGE_NAME-$VERSION distribution..."

DIR=/var/tmp/$PACKAGE_NAME-$VERSION
ORIG_DIR=`pwd`

rm -rf $DIR &&
mkdir -p $DIR &&
cp -a . $DIR &&
cd $DIR &&
find . -name 'Makefile.in' -exec rm '{}' ';' &&
rm -f configure &&
rm -rf Debug Release &&
rm -f *.ncb *.suo *.user &&
mv config.h config.h.old &&
sed "s/devel/$PATCH_LEVEL/" config.h.old >config.h &&
rm -f config.h.old &&
git log --no-merges >ChangeLog &&
find . -name '.gitignore' -exec rm '{}' ';' &&
rm -rf .git &&
cd .. &&
zip -r $PACKAGE_NAME-$VERSION.zip $PACKAGE_NAME-$VERSION &&
mv -f $PACKAGE_NAME-$VERSION.zip $ORIG_DIR &&
rm -rf $DIR
