source bash_colors.sh
BASE=$PWD
BLASTP=/home/huerta/blastdb/ncbi-blast-2.3.0+/bin/blastp
TARGET_SPECIES=$@
TARGET_ORTHO="all one2one"
TARGET_TAX_FILTER="NOG auto"
TARGET_GO='non-electronic experimental'
TARGET_METHODS='hmmer diamond'
TARGET_SELF='noself'
CPU=20
REAL=false

function run {
    clr_bold "   $@"
    if [ $REAL == true ]; then
        eval $@ || (clr_red "ERRORS FOUND! Clean up output files" && exit 1)
    fi
    }

for SPECIES in $TARGET_SPECIES; do
    clr_magenta "********** $SPECIES ***********"
    # BLAST ANNOTATION
    PROTEOME_FILE=data/proteomes/$SPECIES.fa
    BLAST_HITS_FILE=blast/$SPECIES/$SPECIES.fa.hits
    BLAST_UNIQUES_FILE=$BLAST_HITS_FILE.unique_ids

    # Run BLASTP
    TASK="BLASTP"
    if [ ! -f  $BLAST_HITS_FILE ]; then
        clr_brown "[RUN] $TASK"
        run "time $BLASTP -db /dag/db/eggnog_blast/eggnog4.proteins.core_periphery.fa -query $PROTEOME_FILE  -evalue 0.001  -num_threads 20 -outfmt 6 -out $BLAST_HITS_FILE -max_target_seqs 1000000"
    else
        clr_green "[DONE] $TASK"
    fi

    # Precompute lists of blastp hits
    TASK="Unique BLAST hits: $BLAST_UNIQUES_FILE"
    if [ ! -f $BLAST_UNIQUES_FILE ]; then
        clr_brown "[RUN] $TASK"
        run "datamash unique 2 < $BLAST_HITS_FILE | tr ',' '\n'  > $BLAST_UNIQUES_FILE"
    else
        clr_green "[DONE] $TASK"
    fi

    for EMAPPER_METHOD in $TARGET_METHODS; do
        for ORTHO_TYPE in $TARGET_ORTHO; do
            for TAX_SCOPE in $TARGET_TAX_FILTER; do
                for GO_METHOD in $TARGET_GO; do
                    for SELF in $TARGET_SELF; do
                        TAG=$SPECIES.$EMAPPER_METHOD.$ORTHO_TYPE.$TAX_SCOPE.$GO_METHOD.$SELF
                        TASK="Blast GO annotation ($GO_METHOD)"
                        if [ ! -f  blast/$SPECIES/$TAG.1e-40.emapper.blast_filtered_annotations ]; then
                            echo blast/$SPECIES/$TAG.1e-40.emapper.blast_filtered_annotations
                            clr_brown "[RUN] $TASK"
                            run "time python blast_annotate.py $BLAST_HITS_FILE $BLAST_UNIQUES_FILE $SPECIES $GO_METHOD $TAG $CPU"
                        else
                            clr_green "[DONE] $TASK"

                        fi
                    done
                done
            done
        done
    done
done
