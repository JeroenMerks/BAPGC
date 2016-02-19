#!/usr/bin/env bash

# Automatisch herkennen welke benodigde software-paketten van het systeem missen en indien nodig installeren.
# Pas op, het kan voorkomen dat je het wachtwoord van je gebruikersaccount dient op te geven!

if [ "$(dpkg -l | grep python-pip)" == "" ]; then
  echo "Installing python-pip..."
  sudo apt-get --assume-yes install python-pip
fi

if [ "$(dpkg -l | grep python-dev)" == "" ]; then
  echo "Installing python-dev..."
  sudo apt-get --assume-yes install python-dev
fi

if [ "$(dpkg -l | grep libfreetype6-dev)" == "" ]; then
  echo "Installing build-essential..."
  sudo apt-get --assume-yes install libfreetype6-dev
fi

if [ "$(dpkg -l | grep libxft-dev)" == "" ]; then
  echo "Installing libxft-dev..."
  sudo apt-get --assume-yes install libxft-dev
fi

if [ "$(dpkg -l | grep build-essential)" == "" ]; then
  echo "Installing build-essential..."
  sudo apt-get --assume-yes install build-essential
fi

if [ "$(dpkg -l | grep clustalw)" == "" ]; then
  echo "Installing clustalw..."
  sudo apt-get --assume-yes install clustalw
fi

#python -c 'import scipy'
#if [ "$(echo $?)" == 1 ]; then
#  echo "Installing Python module scipy..."
#  sudo pip install scipy -q
#fi

python -c 'import numpy'
if [ "$(echo $?)" == 1 ]; then
  echo "Installing Python module numpy..."
  sudo pip install numpy -q
fi

python -c 'import pandas'
if [ "$(echo $?)" == 1 ]; then
  echo "Installing Python module pandas..."
  sudo pip install pandas -q
fi

python -c 'import mygene'
if [ "$(echo $?)" == 1 ]; then
  echo "Installing Python module mygene..."
  sudo pip install mygene -q
fi

python -c 'import pylab'
if [ "$(echo $?)" == 1 ]; then
  echo "Installing Python module pylab..."
  sudo pip install pylab -q
fi

python -c 'import Bio'
if [ "$(echo $?)" == 1 ]; then
  echo "Installing Python module biopython..."
  sudo pip install biopython -q
fi

python -c 'import ete3'
if [ "$(echo $?)" == 1 ]; then
  echo "Installing Python module ete3..."
  sudo pip install --pre ete3 -q
fi

python -c 'import matplotlib'
if [ "$(echo $?)" == 1 ]; then
  echo "Installing Python module matplotlib..."
  sudo pip install matplotlib -q
fi

