#!/bin/bash



## set up environment

now=$(date +"%m-%d-%Y_%H:%M:%S")

. /etc/profile.d/modules.sh
module load java/1.7.0
export PATH=/data/.bds:$PATH

## submit job 

bds -c /group/rbao/ExScaliburGMD/example/myProject/NA12878trio.bds.cfg -reportHtml -reportYaml -y 0  -l  -s sge Run_ExScaliburGMD.bds -aligners bwaaln bwamem novoalign -callers freebayes gatkhc gatkug ivc mpileup platypus -projdir /group/rbao/ExScaliburGMD/example/myProject -project NA12878trio -samples NA12878 NA12891 NA12892 > Run_ExScaliburGMD.NA12878trio.$now.log.out 2> Run_ExScaliburGMD.NA12878trio.$now.log.err
