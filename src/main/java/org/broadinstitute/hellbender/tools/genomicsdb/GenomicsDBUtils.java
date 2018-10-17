package org.broadinstitute.hellbender.tools.genomicsdb;

import com.googlecode.protobuf.format.JsonFormat;
import com.intel.genomicsdb.model.GenomicsDBExportConfiguration;
import com.intel.genomicsdb.model.GenomicsDBVidMapProto;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;

/**
 * Utility class containing various methods for working with GenomicsDB
 * Contains code to modify the GenomicsDB query output format using the Protobuf API
 *
 * References:
 * GenomicsDB Protobuf structs: https://github.com/Intel-HLS/GenomicsDB/blob/master/src/resources/genomicsdb_vid_mapping.proto
 * Protobuf generated Java code guide:
 * https://developers.google.com/protocol-buffers/docs/javatutorial#the-protocol-buffer-api
 * https://developers.google.com/protocol-buffers/docs/reference/java-generated
 */
public class GenomicsDBUtils {

    /**
     *
     * @param reference reference sequence
     * @param workspace path to the GenomicsDB workspace
     * @param callsetJson path to the GenomicsDB callset JSON
     * @param vidmapJson path to the GenomicsDB vidmap JSON
     * @param vcfHeader VCF with header information to use for header lines on export
     * @return a configuration to determine the output format when the GenomicsDB is queried
     */
    public static GenomicsDBExportConfiguration.ExportConfiguration createExportConfiguration(final File reference, final String workspace,
                                                                                               final String callsetJson, final String vidmapJson,
                                                                                               final String vcfHeader) {
        final GenomicsDBExportConfiguration.ExportConfiguration.Builder exportConfigurationBuilder =
                GenomicsDBExportConfiguration.ExportConfiguration.newBuilder()
                        .setWorkspace(workspace)
                        .setReferenceGenome(reference.getAbsolutePath())
                        .setVidMappingFile(vidmapJson)
                        .setCallsetMappingFile(callsetJson)
                        .setVcfHeaderFilename(vcfHeader)
                        .setProduceGTField(false)
                        .setProduceGTWithMinPLValueForSpanningDeletions(false)
                        .setSitesOnlyQuery(false)
                        .setMaxDiploidAltAllelesThatCanBeGenotyped(GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED);
        final Path arrayFolder = Paths.get(workspace, GenomicsDBConstants.DEFAULT_ARRAY_NAME).toAbsolutePath();

        // For the multi-interval support, we create multiple arrays (directories) in a single workspace -
        // one per interval. So, if you wish to import intervals ("chr1", [ 1, 100M ]) and ("chr2", [ 1, 100M ]),
        // you end up with 2 directories named chr1$1$100M and chr2$1$100M. So, the array names depend on the
        // partition bounds.

        // During the read phase, the user only supplies the workspace. The array names are obtained by scanning
        // the entries in the workspace and reading the right arrays. For example, if you wish to read ("chr2",
        // 50, 50M), then only the second array is queried.

        // In the previous version of the tool, the array name was a constant - genomicsdb_array. The new version
        // will be backward compatible with respect to reads. Hence, if a directory named genomicsdb_array is found,
        // the array name is passed to the GenomicsDBFeatureReader otherwise the array names are generated from the
        // directory entries.
        if (Files.exists(arrayFolder)) {
            exportConfigurationBuilder.setArrayName(GenomicsDBConstants.DEFAULT_ARRAY_NAME);
        } else {
            exportConfigurationBuilder.setGenerateArrayNameFromPartitionBounds(true);
        }

        //Parse the vid json and create an in-memory Protobuf structure representing the information in the JSON file
        GenomicsDBVidMapProto.VidMappingPB vidMapPB;
        try {
            vidMapPB = getProtobufVidMappingFromJsonFile(vidmapJson);
        } catch (final IOException e) {
            throw new UserException("Could not open vid json file " + vidmapJson, e);
        }

        //In vidMapPB, fields is a list of GenomicsDBVidMapProto.GenomicsDBFieldInfo objects
        //Each GenomicsDBFieldInfo object contains information about a specific field in the TileDB/GenomicsDB store
        //We iterate over the list and create a field name to list index map
        final HashMap<String, Integer> fieldNameToIndexInVidFieldsList =
                getFieldNameToListIndexInProtobufVidMappingObject(vidMapPB);

        vidMapPB = updateINFOFieldCombineOperation(vidMapPB, fieldNameToIndexInVidFieldsList,
                GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, "element_wise_sum");


        if (vidMapPB != null) {
            //Use rebuilt vidMap in exportConfiguration
            //NOTE: this does NOT update the JSON file, the vidMapPB is a temporary structure that's passed to
            //C++ modules of GenomicsDB for this specific query. Other queries will continue to use the information
            //in the JSON file
            exportConfigurationBuilder.setVidMapping(vidMapPB);
        }

        return exportConfigurationBuilder.build();
    }

    /**
     * Parse the vid json and create an in-memory Protobuf structure representing the
     * information in the JSON file
     *
     * @param vidmapJson vid JSON file
     * @return Protobuf object
     */
    public static GenomicsDBVidMapProto.VidMappingPB getProtobufVidMappingFromJsonFile(final String vidmapJson)
            throws IOException {
        final GenomicsDBVidMapProto.VidMappingPB.Builder vidMapBuilder = GenomicsDBVidMapProto.VidMappingPB.newBuilder();
        try (final Reader reader = Files.newBufferedReader(IOUtils.getPath(vidmapJson))) {
            JsonFormat.merge(reader, vidMapBuilder);
        }
        return vidMapBuilder.build();
    }

    /**
     * In vidMapPB, fields is a list of GenomicsDBVidMapProto.GenomicsDBFieldInfo objects
     * Each GenomicsDBFieldInfo object contains information about a specific field in the TileDB/GenomicsDB store
     * We iterate over the list and create a field name to list index map
     *
     * @param vidMapPB Protobuf vid mapping object
     * @return map from field name to index in vidMapPB.fields list
     */
    public static HashMap<String, Integer> getFieldNameToListIndexInProtobufVidMappingObject(
            final GenomicsDBVidMapProto.VidMappingPB vidMapPB) {
        final HashMap<String, Integer> fieldNameToIndexInVidFieldsList = new HashMap<>();
        for (int fieldIdx = 0; fieldIdx < vidMapPB.getFieldsCount(); ++fieldIdx) {
            fieldNameToIndexInVidFieldsList.put(vidMapPB.getFields(fieldIdx).getName(), fieldIdx);
        }
        return fieldNameToIndexInVidFieldsList;
    }

    /**
     * Update vid Protobuf object with new combine operation for field
     *
     * @param vidMapPB                        input vid object
     * @param fieldNameToIndexInVidFieldsList name to index in list
     * @param fieldName                       INFO field name
     * @param newCombineOperation             combine op ("sum", "median")
     * @return updated vid Protobuf object if field exists, else null
     */
    public static GenomicsDBVidMapProto.VidMappingPB updateINFOFieldCombineOperation(
            final GenomicsDBVidMapProto.VidMappingPB vidMapPB,
            final Map<String, Integer> fieldNameToIndexInVidFieldsList,
            final String fieldName,
            final String newCombineOperation) {
        final int fieldIdx = fieldNameToIndexInVidFieldsList.containsKey(fieldName)
                ? fieldNameToIndexInVidFieldsList.get(fieldName) : -1;
        if (fieldIdx >= 0) {
            //Would need to rebuild vidMapPB - so get top level builder first
            final GenomicsDBVidMapProto.VidMappingPB.Builder updatedVidMapBuilder = vidMapPB.toBuilder();
            //To update the list element corresponding to fieldName, we get the builder for that specific list element
            final GenomicsDBVidMapProto.GenomicsDBFieldInfo.Builder fieldBuilder =
                    updatedVidMapBuilder.getFieldsBuilder(fieldIdx);
            //And update its combine operation
            fieldBuilder.setVCFFieldCombineOperation(newCombineOperation);

            //Rebuild full vidMap
            return updatedVidMapBuilder.build();
        }
        return null;
    }

}
