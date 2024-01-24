package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.FlowBasedProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.groundtruth.SeriesStats;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;
import org.json.JSONObject;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.atomic.AtomicReference;

@CommandLineProgramProperties(
        summary = "SNV Quality Tool (flow space processing)",
        oneLineSummary = "<add description>",
        programGroup = FlowBasedProgramGroup.class
)

@DocumentedFeature
@ExperimentalFeature
public final class ApplySNVQR extends ReadWalker {

    private static final int        VENDOR_QUALITY_CHECK_FLAG = 0x200;

    private static final String     INCLUDE_QC_FAILED_READS_FULL_NAME = "include-qc-failed-reads";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "File to which reads should be written")
    @WorkflowOutput(optionalCompanions={StandardArgumentDefinitions.OUTPUT_INDEX_COMPANION})
    public GATKPath output = null;
    private SAMFileGATKReadWriter outputWriter;


    @Advanced
    @Argument(fullName=INCLUDE_QC_FAILED_READS_FULL_NAME, doc = "include reads with QC failed flag", optional = true)
    public boolean includeQcFailedReads = false;

    @ArgumentCollection
    public FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();

    @ArgumentCollection
    public ApplySNVQRArgumentCollection aqArgs = new ApplySNVQRArgumentCollection();

    @Advanced
    @Hidden
    @Argument(fullName="debug-collect-stats-into", doc = "", optional = true)
    public String debugCollectStatsInto = null;

    // locals
    private FeatureMapper                       mapper;
    private SvnQualReadProcessor                readProcessor;
    private FlowFeatureMapperUtils.Args         utilArgs;
    private JSONObject                          conf;
    private AddFlowSNVQuality                   addFlowSNVQuality = new AddFlowSNVQuality();
    private SeriesStats                         inputQualStats = new SeriesStats();
    private SeriesStats                         outputBQStats = new SeriesStats();
    private SeriesStats                         outputQXStats = new SeriesStats();

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        outputWriter = createSAMWriter(output, true);
        mapper = new SNVMapper(aqArgs);

        // enforce requirement for sorted input if model provided
        if ( aqArgs.model != null ) {
            if (getHeaderForReads().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
                throw new IllegalArgumentException("input file must be coordinated sorted");
            }
        }

        // load configuration
        if ( aqArgs.model != null && aqArgs.conf == null ) {
            throw new CommandLineException.MissingArgument("conf", "Argument 'conf' must be present if model is specified");
        }
        try {
            conf = aqArgs.conf != null ? new JSONObject(new String(Files.readAllBytes(Paths.get(aqArgs.conf)))) : null;
        } catch (IOException e) {
            throw new GATKException("", e);
        }


        // build utils args
        utilArgs = new FlowFeatureMapperUtils.Args();
        utilArgs.header = getHeaderForReads();
        utilArgs.fbArgs = fbargs;
        utilArgs.limitScore = aqArgs.limitScore;
        utilArgs.debugNegatives = aqArgs.debugNegatives;
        utilArgs.keepNegatives = aqArgs.keepNegatives;
        utilArgs.negativeScoreOverride = aqArgs.negativeScoreOverride;

        readProcessor = new SvnQualReadProcessor(aqArgs, conf);
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if ( outputWriter != null ) {
            outputWriter.close();
        }

        try {
            if ( readProcessor != null ) {
                readProcessor.logZeroCounts();;
                if ( aqArgs.statsPathPrefix != null )
                    readProcessor.writeStats(aqArgs.statsPathPrefix);
            }

            if ( debugCollectStatsInto != null )
                printStats(debugCollectStatsInto);
        } catch (IOException e) {
            throw new GATKException("", e);
        }
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        // include dups?
        if ( read.isDuplicate() && !aqArgs.includeDupReads ) {
            return;
        }

        // include supplementary alignments?
        if ( read.isSupplementaryAlignment() && !aqArgs.keepSupplementaryAlignments ) {
            return;
        }

        // include qc-failed reads?
        if ( ((read.getFlags() & VENDOR_QUALITY_CHECK_FLAG) != 0) && !includeQcFailedReads ) {
            return;
        }

        // collect input stats
        if ( debugCollectStatsInto != null )
            collectInputStats(read);

        // add flow SNV
        addFlowSNVQuality.addBaseQuality(read, getHeaderForReads());

        // if more complex, go do it the hard way
        if ( aqArgs.model != null || aqArgs.conf != null || !read.hasAttribute(attrNameForNonCalledBase('A')) ) {

            // find features in read
            final List<MappedFeature> readFeatures = new LinkedList<>();
            final AtomicReference<FlowBasedRead> flowRead = new AtomicReference<>();
            mapper.forEachOnRead(read, referenceContext, fr -> {
                if (logger.isDebugEnabled()) {
                    logger.debug("fr: " + fr);
                }

                // create flow read?
                if (flowRead.get() == null) {
                    final FlowBasedReadUtils.ReadGroupInfo rgInfo = FlowBasedReadUtils.getReadGroupInfo(utilArgs.header, read);
                    flowRead.set(new FlowBasedRead(read, rgInfo.flowOrder, rgInfo.maxClass, utilArgs.fbArgs));
                }
                fr.flowRead = flowRead.get();

                // score the feature
                fr.score = FlowFeatureMapperUtils.scoreFeature(fr, utilArgs);

                // score for validation mode?
                fr.scoreForBase = new double[SequenceUtil.VALID_BASES_UPPER.length];
                for (int baseIndex = 0; baseIndex < fr.scoreForBase.length; baseIndex++) {
                    final byte base = SequenceUtil.VALID_BASES_UPPER[baseIndex];
                    boolean incl = base != fr.readBases[0];
                    if (incl) {
                        fr.scoreForBase[baseIndex] = FlowFeatureMapperUtils.scoreFeature(fr, base, utilArgs);
                    } else {
                        fr.scoreForBase[baseIndex] = Double.NaN;
                    }
                }
                readFeatures.add(fr);
            });

            // if we accumulated features for this read, incorporate them in on the read level
            readProcessor.incorporateReadFeatures(read, readFeatures, getHeaderForReads());
        }

        // collect output stats
        if ( debugCollectStatsInto != null ) {
            collectOutputStats(read);
        }

        // write read to output
        outputWriter.addRead(read);
    }

    private void collectInputStats(GATKRead read) {
        for ( byte q : read.getBaseQualitiesNoCopy() ) {
            inputQualStats.add(q);
        }
    }

    private void collectOutputStats(GATKRead read) {
        if ( read.hasAttribute("BQ") ) {
            for ( byte q : read.getAttributeAsString("BQ").getBytes() ) {
                outputBQStats.add(q - 33);
            }
        }
        final FlowBasedReadUtils.ReadGroupInfo rgInfo = FlowBasedReadUtils.getReadGroupInfo(utilArgs.header, read);
        for ( int i = 0 ; i < 4 ; i++ ) {
            String attrValue = read.getAttributeAsString(attrNameForNonCalledBase(rgInfo.flowOrder.getBytes()[i]));
            for ( byte q : attrValue.getBytes() ) {
                outputQXStats.add(q - 33);
            }
        }
    }

    private void printStats(final String fname) throws IOException {

        inputQualStats.csvWrite(fname + ".inputQual.csv");
        outputBQStats.csvWrite(fname + ".outputBQ.csv");
        outputQXStats.csvWrite(fname + ".outputQX.csv");
    }

    static public String attrNameForNonCalledBase(byte nonCalledBase) {
        return attrNameForNonCalledBase((char)nonCalledBase);
    }

    static public String attrNameForNonCalledBase(char nonCalledBase) {
        return "q" + Character.toLowerCase(nonCalledBase);
    }
}

