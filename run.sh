readder=evince
#readder=foxitreader
if [ "$1" == "4" ]
then
    make 4 && $readder paper.pdf
else
    make quick && $readder paper.pdf
fi
