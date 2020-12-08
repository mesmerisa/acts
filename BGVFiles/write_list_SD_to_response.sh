#!/bin/bash

# variable to save strings
listSD=""

if [ $1 = "MALTA" ]; then
    # read strings from file and add | between them
    while read F  ; do
            listSD+="|"$F
    done < "MALTA_SD_list_"$2".dat"

    # print the list to screen
    echo $listSD

    # remove first character (which is a |)
    listSD="${listSD:1}"

    # replace strings in response files accordingly:

    sed "s/XXXLISTSDXXX/${listSD}/g" tgeo-bgv-MALTA_bare.response > tgeo-bgv-MALTA_$2.response
    sed -i "s/XXXMALTAXXX/${2}/g" tgeo-bgv-MALTA_$2.response

    sed "s/XXXLISTSDXXX/${listSD}/g" tgeo-bgv-matmap-MALTA_bare.response > tgeo-bgv-matmap-MALTA_$2.response
    sed -i "s/XXXMALTAXXX/${2}/g" tgeo-bgv-matmap-MALTA_$2.response

    sed "s/XXXLISTSDXXX/${listSD}/g" tgeo-bgv-valmatmap-MALTA_bare.response > tgeo-bgv-valmatmap-MALTA_$2.response
    sed -i "s/XXXMALTAXXX/${2}/g" tgeo-bgv-valmatmap-MALTA_$2.response
fi

if [ $1 = "SIDU" ]; then
    echo "SIDU!"
    # read strings from file and add | between them
    while read F  ; do
            listSD+="|"$F
    done < "SIDU_SD_list_"$2".dat"

    # print the list to screen
    echo $listSD

    # remove first character (which is a |)
    listSD="${listSD:1}"

    # replace strings in response files accordingly:

    sed "s/XXXLISTSDXXX/${listSD}/g" tgeo-bgv-SIDU_bare.response > tgeo-bgv-SIDU_$2.response
    sed -i "s/XXXSIDUXXX/${2}/g" tgeo-bgv-SIDU_$2.response

    sed "s/XXXLISTSDXXX/${listSD}/g" tgeo-bgv-matmap-SIDU_bare.response > tgeo-bgv-matmap-SIDU_$2.response
    sed -i "s/XXXSIDUXXX/${2}/g" tgeo-bgv-matmap-SIDU_$2.response

    sed "s/XXXLISTSDXXX/${listSD}/g" tgeo-bgv-valmatmap-SIDU_bare.response > tgeo-bgv-valmatmap-SIDU_$2.response
    sed -i "s/XXXSIDUXXX/${2}/g" tgeo-bgv-valmatmap-SIDU_$2.response
    
    sed "s/XXXLISTSDXXX/${listSD}/g" tgeo-bgv-truthtrack-SIDU_bare.response > tgeo-bgv-truthtrack-SIDU_$2.response
    sed -i "s/XXXSIDUXXX/${2}/g" tgeo-bgv-truthtrack-SIDU_$2.response
        
fi


