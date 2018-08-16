import "combine_tracks.wdl" as CombineTracks
import "aggregate_combined_tracks.wdl" as Aggregate
workflow IgvConvertFiles {

    Array[File] input_files
    String group_id

    scatter (input_file in input_files) {
        call CombineTracks.IGVConvert as IGVConvertFile {
            input:
                COMMENTCHAR="@",
                INPUT=input_file,
                VALUE=basename(input_file),
                FIELD="SAMPLE",
                OUTPUT = basename(input_file) + ".igv.seg",
                PRE_POST="PRE",
                SEGMENT_MEAN_COL = "MEAN_LOG2_COPY_RATIO"
        }
    }
    call Aggregate.TsvCat as TsvCatFiles {
            input:
                input_files = IGVConvertFile.outFile,
                id = group_id + "_Aggregated"
    }
    output {
            File output_igv_aggregated = TsvCatFiles.aggregated_tsv
    }
}