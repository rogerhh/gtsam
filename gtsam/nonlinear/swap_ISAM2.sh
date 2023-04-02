#!/bin/bash

if [ $1 == "old" ] 
then
    rm ISAM2.cpp
    ln -s ISAM2-old.bak ISAM2.cpp
elif [ $1 == "new" ]
then
    rm ISAM2.cpp
    ln -s ISAM2-new.bak ISAM2.cpp
fi

touch ISAM2.cpp
