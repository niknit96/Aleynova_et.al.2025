DIR="$(dirname "${BASH_SOURCE[0]}")"
DIR="$(realpath "${DIR}")"


# Download the SILVA 138 pre-trained classifier for 16s sequences (99% OTUs from V4 region of sequences) 
# from https://docs.qiime2.org/2023.5/data-resources/ (accessed on 19 July 2024)
ls $DIR/16s/ | ( grep silva-138-99-515-806-nb-classifier.qza > /dev/null && echo "SILVA 138 pre-trained classifier already downloaded" ) \
|| ( echo "Download the SILVA 138 pre-trained classifier" && wget -P $DIR/16s/ https://data.qiime2.org/2023.5/common/silva-138-99-515-806-nb-classifier.qza ) \
|| ( echo "Can't download the SILVA 138 pre-trained classifier." && exit )

# Download the UNITE pre-trained classifier for ITS sequences (99% OTUs from ITS1f/ITS2 region of sequences)
# from https://github.com/colinbrislawn/unite-train/releases (accessed on 19 July 2024)
ls $DIR/ITS/ | ( grep unite_ver9_99_all_29.11.2022-Q2-2023.5.qza > /dev/null && echo "UNITE pre-trained classifier already downloaded" ) \
|| ( echo "Download the UNITE pre-trained classifier" && wget -P $DIR/ITS/ https://github.com/colinbrislawn/unite-train/releases/download/9.0-qiime2-2023.5-demo/unite_ver9_99_all_29.11.2022-Q2-2023.5.qza ) \
|| ( echo "Can't download the UNITE pre-trained classifier." && exit )


# Download SRA files used in this article (16s data)
mkdir $DIR/16s/fastq

PRJNA1210692_16s=(SRR32057061 SRR32057060 SRR32057049 SRR32057030 SRR32057043 \
SRR32057000 SRR32056997 SRR32057026 SRR32057023 SRR32057022 SRR32057059 \
SRR32057058 SRR32057057 SRR32057056 SRR32057055 SRR32057054 SRR32057053 \
SRR32057052 SRR32057051 SRR32057050 SRR32057048 SRR32057047 SRR32057046 \
SRR32057037 SRR32057036 SRR32057035 SRR32057034 SRR32057033 SRR32057032 \
SRR32057031 SRR32057021 SRR32057020 SRR32057019 SRR32057018 SRR32057017 SRR32057016)

all_SRA=(${PRJNA1210692_16s[*]})

for sra in ${all_SRA[*]}
do
ls $DIR/16s/fastq | grep $sra > /dev/null || ( prefetch $sra -p -O $DIR/16s/fastq/ && fasterq-dump $DIR/16s/fastq/$sra/$sra.sra -O $DIR/16s/fastq \
&& gzip $DIR/16s/fastq/$sra"_1.fastq" \
&& gzip $DIR/16s/fastq/$sra"_2.fastq" \
&& rm -R $DIR/16s/fastq/$sra/ )
done

for sra in ${all_SRA[*]}
do
ls $DIR/16s/fastq | grep $sra > /dev/null 
    if [[ $? != 0 ]]
        then
        echo "Can't download $sra. Try run script again." && exit 111
        break
    fi
done

if [[ $? == 111 ]]
    then
    exit
fi

# Download SRA files used in this article (ITS data)

mkdir $DIR/ITS/fastq

PRJNA1210692_ITS=(SRR32057015 SRR32057014 SRR32057045 SRR32057044 SRR32057042 \
SRR32057041 SRR32057040 SRR32057039 SRR32057038 SRR32057005 SRR32057004 \
SRR32057003 SRR32057002 SRR32057001 SRR32056999 SRR32056998 SRR32057013 \
SRR32057012 SRR32057011 SRR32057010 SRR32057009 SRR32057008 SRR32057007 \
SRR32057006 SRR32056996 SRR32056995 SRR32056994 SRR32056993 SRR32056992 \
SRR32056991 SRR32056990 SRR32057029 SRR32057028 SRR32057027 SRR32057025 SRR32057024)

all_SRA=(${PRJNA1210692_ITS[*]})

for sra in ${all_SRA[*]}
do
ls $DIR/ITS/fastq | grep $sra > /dev/null || ( prefetch $sra -p -O $DIR/ITS/fastq/ && fasterq-dump $DIR/ITS/fastq/$sra/$sra.sra -O $DIR/ITS/fastq \
&& gzip $DIR/ITS/fastq/$sra"_1.fastq" \
&& gzip $DIR/ITS/fastq/$sra"_2.fastq" \
&& rm -R $DIR/ITS/fastq/$sra/ )
done

for sra in ${all_SRA[*]}
do
ls $DIR/ITS/fastq | grep $sra > /dev/null 
    if [[ $? != 0 ]]
        then
        echo "Can't download $sra. Try run script again." && exit 111
        break
    fi
done

if [[ $? == 111 ]]
    then
    exit
fi

# Rename SRA for input in qiime2 (16s data)
data_16s_1=($(ls $DIR/16s/fastq | grep _1.fastq.gz))
for oldname in ${data_16s_1[*]}
do
newname=$(echo $oldname | sed "s/_1.fastq.gz/_S1_L001_R1_001.fastq.gz/g")
mv $DIR/16s/fastq/$oldname $DIR/16s/fastq/$newname
done

data_16s_2=($(ls $DIR/16s/fastq | grep _2.fastq.gz))
for oldname in ${data_16s_2[*]}
do
newname=$(echo $oldname | sed "s/_2.fastq.gz/_S1_L001_R2_001.fastq.gz/g")
mv $DIR/16s/fastq/$oldname $DIR/16s/fastq/$newname
done

# Rename SRA for input in qiime2 (ITS data)
data_ITS_1=($(ls $DIR/ITS/fastq | grep _1.fastq.gz))
for oldname in ${data_ITS_1[*]}
do
newname=$(echo $oldname | sed "s/_1.fastq.gz/_S1_L001_R1_001.fastq.gz/g")
mv $DIR/ITS/fastq/$oldname $DIR/ITS/fastq/$newname
done

data_ITS_2=($(ls $DIR/ITS/fastq | grep _2.fastq.gz))
for oldname in ${data_ITS_2[*]}
do
newname=$(echo $oldname | sed "s/_2.fastq.gz/_S1_L001_R2_001.fastq.gz/g")
mv $DIR/ITS/fastq/$oldname $DIR/ITS/fastq/$newname
done