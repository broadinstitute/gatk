package org.broadinstitute.hellbender.tools.variantdb;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StandardAnnotation;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;


@CommandLineProgramProperties(
        summary = "(\"ExtractCohort\") - Filter and extract arrayvariants out of big query.",
        oneLineSummary = "Tool to extract variants out of big query for a subset of samples",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class ArrayExtractCohort extends GATKTool {
    private static final Logger logger = LogManager.getLogger(ExtractCohort.class);
    public static final int DEFAULT_LOCAL_SORT_MAX_RECORDS_IN_RAM = 1000000;
    private VariantContextWriter vcfWriter = null;
    private ArrayExtractCohortEngine engine;

    public enum QueryMode {
        LOCAL_SORT,
        QUERY
    }

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output VCF file to which annotated variants should be written.",
            optional = false
    )
    private String outputVcfPathString = null;

    @Argument(
            fullName = "project-id",
            doc = "ID of the Google Cloud project to use when executing queries",
            optional = false
    )
    private String projectID = null;

    @Argument(
            fullName = "sample-info-table",
            doc = "Fully qualified name of a bigquery table containing a single column `sample` that describes the full list of samples to evoque",
            optional = true
    )
    private String sampleTableName = null;

    @Argument(
            fullName = "probe-info-table",
            doc = "Fully qualified name of a bigquery table containing probe information",
            optional = true
    )
    private String probeTableName = null;

    @Argument(
        fullName = "probe-info-csv",
        doc = "Filepath to CSV export of probe-info table",
        optional = true
)
    private String probeCsvExportFile = null;

    @Argument(
            fullName = "cohort-extract-table",
            doc = "Fully qualified name of the table where the cohort data exists (already subsetted)",
            optional = false
    )
    private String cohortTable = null;

    @Argument(
            fullName = "print-debug-information",
            doc = "If true, print extra debugging output",
            optional = true)
    private boolean printDebugInformation = false;

    @Argument(
            fullName = "local-sort-max-records-in-ram",
            doc = "When doing local sort, store at most this many records in memory at once",
            optional = true
    )
    private int localSortMaxRecordsInRam = DEFAULT_LOCAL_SORT_MAX_RECORDS_IN_RAM;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean useVariantAnnotations() { return true; }

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Arrays.asList(
                StandardAnnotation.class, AS_StandardAnnotation.class
        );
    }

    @Override
    protected void onStartup() {
        super.onStartup();

        //TODO verify what we really need here
        final VariantAnnotatorEngine annotationEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), null, Collections.emptyList(), false, false);

        vcfWriter = createVCFWriter(IOUtils.getPath(outputVcfPathString));

        Map<Integer, String> sampleIdMap = ExtractCohortBQ.getSampleIdMap(sampleTableName, printDebugInformation);

        Collection<String> sampleNames = sampleIdMap.values();
        VCFHeader header = CommonCode.generateRawArrayVcfHeader(new HashSet<>(sampleNames), reference.getSequenceDictionary());

        Map<Long, ProbeInfo> probeIdMap;

        if (probeCsvExportFile == null) {
            probeIdMap = ExtractCohortBQ.getProbeIdMap(probeTableName, printDebugInformation);
        } else {
            probeIdMap = new HashMap<>();
            String line = "";
            try (BufferedReader br = new BufferedReader(new FileReader(probeCsvExportFile))) {
                /// skip the header
                br.readLine();

                while ((line = br.readLine()) != null) {

                    // use comma as separator
                    String[] fields = line.split(",");
                    //ProbeId,Name,GenomeBuild,Chr,Position,Ref,AlleleA,AlleleB,build37Flag
                    //6,ilmnseq_rs9651229_F2BT,37,1,567667,,,,PROBE_SEQUENCE_MISMATCH
                    ProbeInfo p = new ProbeInfo(Long.parseLong(fields[0]),
                                                fields[1], // name
                                                fields[3], // contig
                                                Long.parseLong(fields[4]),   // position
                                                fields[5], // ref
                                                fields[6], // alleleA
                                                fields[7]);// alleleB

                    probeIdMap.put(p.probeId, p);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }



        //ChromosomeEnum.setRefVersion(refVersion);

        engine = new ArrayExtractCohortEngine(
                projectID,
                vcfWriter,
                header,
                annotationEngine,
                reference,
                sampleIdMap,
                probeIdMap,
                cohortTable,
                localSortMaxRecordsInRam,
                false,
                printDebugInformation,
                progressMeter);
        vcfWriter.writeHeader(header);
    }

    @Override
    // maybe think about creating a BigQuery Row walker?
    public void traverse() {
        progressMeter.setRecordsBetweenTimeChecks(100L);
        engine.traverse();
     }

    @Override
    protected void onShutdown() {
        super.onShutdown();

        if ( engine != null ) {
            logger.info(String.format("***Processed %d total sites", engine.getTotalNumberOfSites()));
            logger.info(String.format("***Processed %d total variants", engine.getTotalNumberOfVariants()));
        }

        // Close up our writer if we have to:
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

}
