source common.sh

TARGET_SPECIES=$@
CPU=20
EMAPPER=/home/huerta/eggnog-mapper-latest/emapper.py
BASE=$PWD

TARGET_ORTHO="all one2one"
TARGET_TAX_FILTER="auto NOG"
TARGET_GO='non-electronic experimental'
TARGET_METHODS='diamond hmmer'

REAL=true

for SPECIES in $TARGET_SPECIES; do
    for METHOD in $TARGET_METHODS; do
        PROTEOME_FILE=data/proteomes/$SPECIES.fa
        clr_magenta $SPECIES $METHOD
        BASE_EMAPPER_FILE=$BASE/emapper/$SPECIES/$SPECIES.$METHOD
        if [ "$METHOD" == 'diamond' ]; then
            BASE_EMAPPER_ARGS="-m diamond -d none --cpu $CPU --override"
        elif [ "$METHOD" == 'hmmer' ]; then
            if [ "$SPECIES" == '511145' ]; then
                BASE_EMAPPER_ARGS="-m hmmer -d bact:localhost:51500 --cpu $CPU --override"
            else
                BASE_EMAPPER_ARGS="-m hmmer -d euk:localhost:51400 --cpu $CPU --override"
            fi
        else
            exit 1
        fi

        # GET Seed orthologs excluding self hits
        TASK="emapper seed_orthologs (noself)"
        if [ ! -f $BASE_EMAPPER_FILE.noself.emapper.seed_orthologs ]; then
            clr_brown "[RUN] $TASK"
            run "time $EMAPPER $BASE_EMAPPER_ARGS -i $PROTEOME_FILE --no_annot \
                  -o $BASE_EMAPPER_FILE.noself \
                  --excluded_taxa $SPECIES "
        else
            clr_green "[DONE] $TASK"
        fi

        # GET Seed orthologs
        TASK="emapper seed_orthologs (self)"
        if [ ! -f $BASE_EMAPPER_FILE.self.emapper.seed_orthologs ]; then
            clr_brown "[RUN] $TASK"
            run "time $EMAPPER $BASE_EMAPPER_ARGS -i $PROTEOME_FILE --no_annot \
                     -o $BASE_EMAPPER_FILE.self "
        else
            clr_green "[DONE] $TASK"
        fi

        # Annotate functions based on seed orthologs testing different strategies
        for ORTHO_TYPE in $TARGET_ORTHO; do
            for TAX_SCOPE in $TARGET_TAX_FILTER; do
                for GO_METHOD in $TARGET_GO; do
                    BASE_ANNOTATION_FILE=$BASE_EMAPPER_FILE.$ORTHO_TYPE.$TAX_SCOPE.$GO_METHOD

                    # Excluding orthologs in the self proteome
                    TASK="emapper annotations $BASE_ANNOTATION_FILE (noself)"
                    if [ ! -f $BASE_ANNOTATION_FILE.noself.emapper.annotations ]; then
                        clr_brown "[RUN] $TASK"
                        run "time $EMAPPER $BASE_EMAPPER_ARGS --annotate_hits_table $BASE_EMAPPER_FILE.noself.emapper.seed_orthologs \
                             --tax_scope $TAX_SCOPE \
                             --target_orthologs $ORTHO_TYPE \
                             --go_evidence $GO_METHOD \
                             -o $BASE_ANNOTATION_FILE.noself \
                             --report_orthologs \
                             --excluded_taxa $SPECIES "
                    else
                        clr_green "[DONE] $TASK"
                    fi

                    # Using orthologs from all species (for interpro comparison)
                    TASK="emapper annotations $BASE_ANNOTATION_FILE (self)"
                    if [ ! -f $BASE_ANNOTATION_FILE.self.emapper.annotations ]; then
                        clr_brown "[RUN] $TASK"
                        run "time $EMAPPER $BASE_EMAPPER_ARGS --annotate_hits_table $BASE_EMAPPER_FILE.self.emapper.seed_orthologs \
                             --tax_scope $TAX_SCOPE \
                             --target_orthologs $ORTHO_TYPE \
                             --go_evidence $GO_METHOD \
                             --report_orthologs \
                             -o $BASE_ANNOTATION_FILE.self"
                    else
                        clr_green "[DONE] $TASK"
                    fi

                done
            done

        done
    done
done

