#!/bin/bash
set -e
if [ $# -ne 2 ] ;
then
  echo "Usage: $0 <new_version> <from_branch>"
  echo ""
  echo "Merge <from_branch> into stable to make <new_version>"
  echo "Last tag:" "stable_r"`git tag | cut -dr -f2 | sort -n | tail -1`
  exit 1
fi
cd $GACODE_ROOT
git checkout $2
git pull --rebase=false
git checkout stable
git pull
git merge $2
git submodule update
python shared/bin/gacode_regression.py -clean || git reset --hard origin/stable; exit
git tag -a $1
git push
git push --tags
git checkout master
git pull --rebase=false
git merge stable
git push
