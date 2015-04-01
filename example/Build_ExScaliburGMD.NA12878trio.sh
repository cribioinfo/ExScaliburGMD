

## build pipeline scripts 

now=$(date +"%m-%d-%Y_%H:%M:%S")

module load exscaliburgmd/0.5.0

## script 
# BuildExScaliburGMDExe=/group/rbao/ExScaliburGMD/Build_ExScaliburGMD.pl
BuildExScaliburGMDExe=../Build_ExScaliburGMD.pl

## project info 
project=NA12878trio
metadata=NA12878trio.metadata.txt
config=NA12878trio.pipeline.yaml
projdir=$PWD/myProject

## options
threads=2
mapq=30
force="--force"
split="--split"
tree="--tree tree"

## bds 
platform="--cluster"
scheduler="--sge"
bdscfg=""
retry=0
log="--log"

## command 
echo "START" `date`
$BuildExScaliburGMDExe \
	--project $project \
	--metadata $metadata \
	--config $config \
	--projdir $projdir \
	--threads $threads \
	--mapq $mapq \
	--retry $retry \
	$force \
	$split \
	$tree \
	$platform \
	$scheduler \
	$bdscfg \
	$log \
	>& Build_ExScaliburGMD.$project.$now.log
echo "END" `date`
