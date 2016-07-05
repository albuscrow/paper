if [ "$1" == "4" ]
then
    make 4 && foxitreader paper.pdf
else
    make quick && foxitreader paper.pdf
fi
