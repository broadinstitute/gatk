package org.broadinstitute.hellbender.tools.genomicsdb;

import com.googlecode.protobuf.format.JsonFormat;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.genomicsdb.importer.GenomicsDBImporter;
import org.genomicsdb.model.GenomicsDBExportConfiguration;
import org.genomicsdb.model.GenomicsDBVidMapProto;

import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;

/**
 * Utility class containing various methods for working with GenomicsDB
 * Contains code to modify the GenomicsDB import input or query output format using the Protobuf API
 *
 * References:
 * GenomicsDB Protobuf structs: https://github.com/GenomicsDB/GenomicsDB/blob/master/src/resources/genomicsdb_vid_mapping.proto
 * Protobuf generated Java code guide:
 * https://developers.google.com/protocol-buffers/docs/javatutorial#the-protocol-buffer-api
 * https://developers.google.com/protocol-buffers/docs/reference/java-generated
 */
public class GenomicsDBUtils {

    private static final String ELEMENT_WISE_SUM = "element_wise_sum";
    private static final String ELEMENT_WISE_FLOAT_SUM = "element_wise_float_sum";
    private static final String SUM = "sum";
    private static final String HISTOGRAM_SUM = "histogram_sum";
    private static final String STRAND_BIAS_TABLE_COMBINE = "strand_bias_table";
    private static final String GDB_TYPE_FLOAT = "float";
    private static final String GDB_TYPE_INT = "int";

  /**
   * See org.genomicsdb.importer.Constants for hardcoded allele-specific annotation fields to be
   * treated as an array of int/float vectors. All other allele-specific fields that need to be treated
   * as an array of int/float vectors will have to explicitly overridden in this method.
   *
   * @param importer
   */
  public static void updateImportProtobufMapping(GenomicsDBImporter importer) {
      // Example code in case an allele-specific method has to overridden
/*
      GenomicsDBVidMapProto.VidMappingPB vidMapPB = importer.getProtobufVidMapping();
      if (vidMapPB == null) {
          return;
      }

      // In vidMapPB, fields is a list of GenomicsDBVidMapProto.GenomicsDBFieldInfo objects
      // Each GenomicsDBFieldInfo object contains information about a specific field in the
      // GenomicsDB store
      // We iterate over the list and create a field name to list index map
      final HashMap<String, Integer> fieldNameToIndexInVidFieldsList =
              getFieldNameToListIndexInProtobufVidMappingObject(vidMapPB);

      //Fields that need to be treated as array of int/float vectors need to be updated via protobuf
      vidMapPB = updateAlleleSpecificINFOFieldCombineOperation(vidMapPB, fieldNameToIndexInVidFieldsList,
              GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY, ELEMENT_WISE_FLOAT_SUM, true);

      importer.updateProtobufVidMapping(vidMapPB);
*/
  }

    /**
     *
     * @param workspace path to the GenomicsDB workspace
     * @param callsetJson path to the GenomicsDB callset JSON
     * @param vidmapJson path to the GenomicsDB vidmap JSON
     * @param vcfHeader VCF with header information to use for header lines on export
     * @param genomicsDBOptions genotyping parameters to read from a GenomicsDB
     * @return a configuration to determine the output format when the GenomicsDB is queried
     */
    public static GenomicsDBExportConfiguration.ExportConfiguration createExportConfiguration(final String workspace,
                                                                                              final String callsetJson, final String vidmapJson,
                                                                                              final String vcfHeader, final GenomicsDBOptions genomicsDBOptions) {
        final GenomicsDBExportConfiguration.ExportConfiguration.Builder exportConfigurationBuilder =
                GenomicsDBExportConfiguration.ExportConfiguration.newBuilder()
                        .setWorkspace(workspace)
                        .setReferenceGenome(genomicsDBOptions.getReference().toAbsolutePath().toString())
                        .setCallsetMappingFile(callsetJson)
                        .setVcfHeaderFilename(vcfHeader)
                        .setProduceGTField(false)
                        .setProduceGTWithMinPLValueForSpanningDeletions(false)
                        .setSitesOnlyQuery(false)
                        .setMaxDiploidAltAllelesThatCanBeGenotyped(GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED);
        if (genomicsDBOptions != null) {
            exportConfigurationBuilder.setProduceGTField(genomicsDBOptions.doCallGenotypes()).
                    setMaxDiploidAltAllelesThatCanBeGenotyped(genomicsDBOptions.getMaxAlternateAlleles());
        }


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
                GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, ELEMENT_WISE_SUM);

        vidMapPB = updateAlleleSpecificINFOFieldCombineOperation(vidMapPB, fieldNameToIndexInVidFieldsList,
                GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY, ELEMENT_WISE_FLOAT_SUM);

        //Update combine operations for GnarlyGenotyper
        //Note that this MQ format is deprecated, but was used by the prototype version of ReblockGVCF
        vidMapPB = updateINFOFieldCombineOperation(vidMapPB, fieldNameToIndexInVidFieldsList,
                GATKVCFConstants.MAPPING_QUALITY_DEPTH_DEPRECATED, SUM);
        vidMapPB = updateINFOFieldCombineOperation(vidMapPB, fieldNameToIndexInVidFieldsList,
                GATKVCFConstants.RAW_QUAL_APPROX_KEY, SUM);
        vidMapPB = updateINFOFieldCombineOperation(vidMapPB, fieldNameToIndexInVidFieldsList,
                GATKVCFConstants.VARIANT_DEPTH_KEY, SUM);
        vidMapPB = updateINFOFieldCombineOperation(vidMapPB, fieldNameToIndexInVidFieldsList,
                GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY, ELEMENT_WISE_SUM);
        vidMapPB = updateAlleleSpecificINFOFieldCombineOperation(vidMapPB, fieldNameToIndexInVidFieldsList,
                GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY, ELEMENT_WISE_FLOAT_SUM);


        if (vidMapPB != null) {
            //Use rebuilt vidMap in exportConfiguration
            //NOTE: this does NOT update the JSON file, the vidMapPB is a temporary structure that's passed to
            //C++ modules of GenomicsDB for this specific query. Other queries will continue to use the information
            //in the JSON file
            exportConfigurationBuilder.setVidMapping(vidMapPB);
        } else {
            exportConfigurationBuilder.setVidMappingFile(vidmapJson);
        }

        return exportConfigurationBuilder.build();
    }

    public static GenomicsDBExportConfiguration.ExportConfiguration createExportConfiguration(final String workspace,
                                                                                              final String callsetJson, final String vidmapJson,
                                                                                              final String vcfHeader) {
        return createExportConfiguration(workspace, callsetJson, vidmapJson, vcfHeader, null);
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
     * @return updated vid Protobuf object if field exists, else return original vidmap object
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
        return vidMapPB;
    }

    public static GenomicsDBVidMapProto.VidMappingPB updateAlleleSpecificINFOFieldCombineOperation(
            final GenomicsDBVidMapProto.VidMappingPB vidMapPB,
            final Map<String, Integer> fieldNameToIndexInVidFieldsList,
            final String fieldName,
            final String newCombineOperation) {
        return updateAlleleSpecificINFOFieldCombineOperation(vidMapPB, fieldNameToIndexInVidFieldsList, fieldName,
                newCombineOperation, false);
    }

    /**
     * Update vid Protobuf object with a new variable length descriptor, as for allele-specific annotations
     * @param vidMapPB input vid object
     * @param fieldNameToIndexInVidFieldsList name to index in list
     * @param fieldName INFO field name
     * @param newCombineOperation combine op ("histogram_sum", "element_wise_float_sum", "strand_bias_table")
     * @return updated vid Protobuf object if field exists, else null
     */
    public static GenomicsDBVidMapProto.VidMappingPB updateAlleleSpecificINFOFieldCombineOperation(
            final GenomicsDBVidMapProto.VidMappingPB vidMapPB,
            final Map<String, Integer> fieldNameToIndexInVidFieldsList,
            final String fieldName,
            final String newCombineOperation,
            final boolean isGenomicsDBImportOperation)
    {
        int fieldIdx = fieldNameToIndexInVidFieldsList.containsKey(fieldName)
                ? fieldNameToIndexInVidFieldsList.get(fieldName) : -1;
        if(fieldIdx >= 0) {
            //Would need to rebuild vidMapPB - so get top level builder first
            GenomicsDBVidMapProto.VidMappingPB.Builder updatedVidMapBuilder = vidMapPB.toBuilder();
            //To update the list element corresponding to fieldName, we get the builder for that specific list element
            GenomicsDBVidMapProto.GenomicsDBFieldInfo.Builder infoBuilder =
                    updatedVidMapBuilder.getFieldsBuilder(fieldIdx);

            if (isGenomicsDBImportOperation && !infoBuilder.getVCFFieldCombineOperation().equals(newCombineOperation)) {
                GenomicsDBVidMapProto.FieldLengthDescriptorComponentPB.Builder lengthDescriptorComponentBuilder =
                        GenomicsDBVidMapProto.FieldLengthDescriptorComponentPB.newBuilder();
                lengthDescriptorComponentBuilder.setVariableLengthDescriptor("R");

                infoBuilder.clearLength();
                infoBuilder.addLength(lengthDescriptorComponentBuilder.build());
                lengthDescriptorComponentBuilder.setVariableLengthDescriptor("var"); //ignored - can set anything here
                infoBuilder.addLength(lengthDescriptorComponentBuilder.build());

                infoBuilder.clearVcfDelimiter();
                infoBuilder.addVcfDelimiter(AnnotationUtils.ALLELE_SPECIFIC_PRINT_DELIM);
                infoBuilder.addVcfDelimiter(AnnotationUtils.ALLELE_SPECIFIC_REDUCED_DELIM);

                infoBuilder.clearType();
                if (newCombineOperation.equals(HISTOGRAM_SUM)) {
                  infoBuilder.addType(GDB_TYPE_FLOAT);
                  infoBuilder.addType(GDB_TYPE_INT);
                } else {
                  if (newCombineOperation.equals(ELEMENT_WISE_FLOAT_SUM)) {
                    infoBuilder.addType(GDB_TYPE_FLOAT);
                  } else {
                    infoBuilder.addType(GDB_TYPE_INT);
                  }
                }
            }

            if (newCombineOperation.equals(HISTOGRAM_SUM)) {
                infoBuilder.setVCFFieldCombineOperation(HISTOGRAM_SUM);
            } else {
                infoBuilder.setVCFFieldCombineOperation(ELEMENT_WISE_SUM);
            }

            //Rebuild full vidMap
            return updatedVidMapBuilder.build();
        }
        return vidMapPB;
    }

}
