/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2014 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.hellbender.tools.walkers.bqsr;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.tribble.Feature;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.HACKRefMetaDataTracker;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilters;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.recalibration.*;
import org.broadinstitute.hellbender.tools.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.baq.BAQ;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.recalibration.EventType;
import org.broadinstitute.hellbender.utils.sam.AlignmentUtils;
import org.broadinstitute.hellbender.utils.sam.ReadUtils;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.function.Predicate;

/**
 * First pass of the base quality score recalibration -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).
 *
 * <p>
 * This walker is designed to work as the first pass in a two-pass processing step. It does a by-locus traversal operating
 * only at sites that are not in dbSNP. We assume that all reference mismatches we see are therefore errors and indicative
 * of poor base quality. This walker generates tables based on various user-specified covariates (such as read group,
 * reported quality score, cycle, and context). Since there is a large amount of data one can then calculate an empirical
 * probability of error given the particular covariates seen at this site, where p(error) = num mismatches / num observations.
 * The output file is a table (of the several covariate values, num observations, num mismatches, empirical quality score).
 * <p>
 * Note: ReadGroupCovariate and QualityScoreCovariate are required covariates and will be added for the user regardless of whether or not they were specified.
 *
 * <p>
 *
 * <h3>Input</h3>
 * <p>
 * The input read data whose base quality scores need to be assessed.
 * <p>
 * A database of known polymorphic sites to skip over.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A GATK Report file with many tables:
 * <ol>
 *     <li>The list of arguments</li>
 *     <li>The quantized qualities table</li>
 *     <li>The recalibration table by read group</li>
 *     <li>The recalibration table by quality score</li>
 *     <li>The recalibration table for all the optional covariates</li>
 * </ol>
 *
 * The GATK Report is intended to be easy to read by humans or computers. Check out the documentation of the GATKReport to learn how to manipulate this table.
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -T BaseRecalibrator \
 *   -I my_reads.bam \
 *   -R resources/Homo_sapiens_assembly18.fasta \
 *   -knownSites bundle/hg18/dbsnp_132.hg18.vcf \
 *   -knownSites another/optional/setOfSitesToMask.vcf \
 *   -o recal_data.table
 * </pre>
 */

@CommandLineProgramProperties(
        usage = "First pass of the Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).",
        usageShort = "Generates recalibration table",
        programGroup = ReadProgramGroup.class
)
public class BaseRecalibrator extends ReadWalker {
    final protected static Logger logger = LogManager.getLogger(BaseRecalibrator.class);

    /**
     * all the command line arguments for BQSR and its covariates
     */
    @ArgumentCollection(doc="all the command line arguments for BQSR and its covariates")
    private final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

    /**
     * When you have nct > 1, BQSR uses nct times more memory to compute its recalibration tables, for efficiency
     * purposes.  If you have many covariates, and therefore are using a lot of memory, you can use this flag
     * to safely access only one table.  There may be some CPU cost, but as long as the table is really big
     * there should be relatively little CPU costs.
     */
    @Argument(fullName = "lowMemoryMode", shortName="lowMemoryMode", doc="Reduce memory usage in multi-threaded code at the expense of threading efficiency", optional = true)
    public boolean lowMemoryMode = false;

    @Argument(fullName = "bqsrBAQGapOpenPenalty", shortName="bqsrBAQGOP", doc="BQSR BAQ gap open penalty (Phred Scaled).  Default value is 40.  30 is perhaps better for whole genome call sets", optional = true)
    public double BAQGOP = BAQ.DEFAULT_GOP;

    /**
     * Enables on-the-fly recalibrate of base qualities, intended primarily for use with BaseRecalibrator and PrintReads
     * (see Best Practices workflow documentation). The covariates tables are produced by the BaseRecalibrator tool.
     * Please be aware that you should only run recalibration with the covariates file created on the same input bam(s).
     */
    @Argument(fullName="BQSR", shortName="BQSR", optional=true, doc="Input covariates table file for on-the-fly base quality score recalibration")
    public File BQSR_RECAL_FILE = null;


    /**
     * This flag tells GATK not to modify quality scores less than this value. Instead they will be written out unmodified in the recalibrated BAM file.
     * In general it's unsafe to change qualities scores below < 6, since base callers use these values to indicate random or bad bases.
     * For example, Illumina writes Q2 bases when the machine has really gone wrong. This would be fine in and of itself,
     * but when you select a subset of these reads based on their ability to align to the reference and their dinucleotide effect,
     * your Q2 bin can be elevated to Q8 or Q10, leading to issues downstream.
     */
    @Argument(fullName = "preserve_qscores_less_than", shortName = "preserveQ", doc = "Don't recalibrate bases with quality scores less than this threshold (with -BQSR)", optional = true)
    public int PRESERVE_QSCORES_LESS_THAN = QualityUtils.MIN_USABLE_Q_SCORE;


    // --------------------------------------------------------------------------------------------------------------
    //
    // quality encoding checking arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * By default the GATK assumes that base quality scores start at Q0 == ASCII 33 according to the SAM specification.
     * However, encoding in some datasets (especially older Illumina ones) starts at Q64. This argument will fix the
     * encodings on the fly (as the data is read in) by subtracting 31 from every quality score. Note that this argument should
     * NEVER be used by default; you should only use it when you have confirmed that the quality scores in your data are
     * not in the correct encoding.
     */
    @Argument(fullName = "fix_misencoded_quality_scores", shortName="fixMisencodedQuals", doc="Fix mis-encoded base quality scores", optional = true)
    public boolean FIX_MISENCODED_QUALS = false;

    /**
     * This flag tells GATK to ignore warnings when encountering base qualities that are too high and that seemingly
     * indicate a problem with the base quality encoding of the BAM file. You should only use this if you really know
     * what you are doing; otherwise you could seriously mess up your data and ruin your analysis.
     */
    @Argument(fullName = "allow_potentially_misencoded_quality_scores", shortName="allowPotentiallyMisencodedQuals", doc="Ignore warnings about base quality score encoding", optional = true)
    public boolean ALLOW_POTENTIALLY_MISENCODED_QUALS = false;

    /**
     * This flag tells GATK to use the original base qualities (that were in the data before BQSR/recalibration) which
     * are stored in the OQ tag, if they are present, rather than use the post-recalibration quality scores. If no OQ
     * tag is present for a read, the standard qual score will be used.
     */
    @Argument(fullName="useOriginalQualities", shortName = "OQ", doc = "Use the base quality scores from the OQ tag", optional = true)
    public Boolean useOriginalBaseQualities = false;

    /**
     * If reads are missing some or all base quality scores, this value will be used for all base quality scores.
     * By default this is set to -1 to disable default base quality assignment.
     */
    //TODO: minValue = 0, maxValue = Byte.MAX_VALUE)
    @Argument(fullName="defaultBaseQualities", shortName = "DBQ", doc = "Assign a default base quality", optional = true)
    public byte defaultBaseQualities = -1;

    /**
     * an object that keeps track of the information necessary for quality score quantization
     */
    private QuantizationInfo quantizationInfo;

    /**
     * list to hold the all the covariate objects that were requested (required + standard + experimental)
     */
    private Covariate[] requestedCovariates;

    private RecalibrationEngine recalibrationEngine;

    private int minimumQToUse;

    private static final String NO_DBSNP_EXCEPTION = "This calculation is critically dependent on being able to skip over known variant sites. Please provide a VCF file containing known sites of genetic variation.";

    private BAQ baq; // BAQ the reads on the fly to generate the alignment uncertainty vector
    private ReferenceDataSource referenceDataSource; // datasource for the reference. We're using a different one from the engine itself to avoid messing with its caches.
    private final static byte NO_BAQ_UNCERTAINTY = (byte)'@';


    private long accumulator;
    /**
     * Parse the -cov arguments and create a list of covariates to be used here
     * Based on the covariates' estimates for initial capacity allocate the data hashmap
     */
    @Override
    public void onTraversalStart() {
        accumulator= 0L;

        baq = new BAQ(BAQGOP); // setup the BAQ object with the provided gap open penalty

        if (RAC.FORCE_PLATFORM != null)
            RAC.DEFAULT_PLATFORM = RAC.FORCE_PLATFORM;

        if (RAC.knownSitesVCF == null && !RAC.RUN_WITHOUT_DBSNP) // Warn the user if no dbSNP file or other variant mask was specified
            throw new UserException.CommandLineException(NO_DBSNP_EXCEPTION);

        RAC.existingRecalibrationReport = BQSR_RECAL_FILE; // if we have a recalibration file, record it so it goes on the report table

        Pair<ArrayList<Covariate>, ArrayList<Covariate>> covariates = RecalUtils.initializeCovariates(RAC); // initialize the required and optional covariates
        ArrayList<Covariate> requiredCovariates = covariates.getLeft();
        ArrayList<Covariate> optionalCovariates = covariates.getRight();

        requestedCovariates = new Covariate[requiredCovariates.size() + optionalCovariates.size()];
        int covariateIndex = 0;
        for (final Covariate covariate : requiredCovariates)
            requestedCovariates[covariateIndex++] = covariate;
        for (final Covariate covariate : optionalCovariates)
            requestedCovariates[covariateIndex++] = covariate;

        logger.info("The covariates being used here: ");
        for (Covariate cov : requestedCovariates) { // list all the covariates being used
            logger.info("\t" + cov.getClass().getSimpleName());
            cov.initialize(RAC); // initialize any covariate member variables using the shared argument collection
        }

        try {
            RAC.RECAL_TABLE = new PrintStream(RAC.RECAL_TABLE_FILE);
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(RAC.RECAL_TABLE_FILE, e);
        }

        initializeRecalibrationEngine();
        minimumQToUse = PRESERVE_QSCORES_LESS_THAN;
        referenceDataSource = REFERENCE_FILE != null ? new ReferenceDataSource(REFERENCE_FILE) : null;
    }

    /**
     * Initialize the recalibration engine
     */
    private void initializeRecalibrationEngine() {
        int numReadGroups = getHeaderForReads().getReadGroups().size();

        recalibrationEngine = new RecalibrationEngine(requestedCovariates, numReadGroups, RAC.RECAL_TABLE_UPDATE_LOG, lowMemoryMode);
    }

    private boolean isLowQualityBase( final SAMRecord read, final int offset ) {
        return read.getBaseQualities()[offset] < minimumQToUse;
    }

    @Override
    public Predicate<SAMRecord> makeReadFilter() {
        return super.makeReadFilter()
                .and(ReadFilters.MAPPING_QUALITY_ZERO.negate())
                .and(ReadFilters.MAPPING_QUALITY_UNAVAIALBLE.negate())
                .and(ReadFilters.UNMAPPED.negate())
                .and(ReadFilters.NOT_PRIMARY_ALIGNMENT.negate())
                .and(ReadFilters.DUPLICATE.negate())
                .and(ReadFilters.FAILS_VENDOR_QUALITY_CHECK.negate());
    }

    private static SAMRecord resetOriginalQualities(SAMRecord read){
        byte[] originalQuals = read.getOriginalBaseQualities();
        if ( originalQuals != null )
            read.setBaseQualities(originalQuals);
        return read;
    }

    private static SAMRecord consolidateCigar(SAMRecord read) {
        // Always consolidate the cigar string into canonical form, collapsing zero-length / repeated cigar elements.
        // Downstream code (like LocusIteratorByState) cannot necessarily handle non-consolidated cigar strings.
        read.setCigar(AlignmentUtils.consolidateCigar(read.getCigar()));
        return read;
    }

    private SAMRecord resetOriginalBaseQualities(SAMRecord read){
        if (! useOriginalBaseQualities) {
            return read;
        }
        byte[] originalQuals = read.getOriginalBaseQualities();
        if ( originalQuals != null )
            read.setBaseQualities(originalQuals);
        return read;
    }

    private SAMRecord setDefaultBaseQualities(SAMRecord read) {
        // if we are using default quals, check if we need them, and add if necessary.
        // 1. we need if reads are lacking or have incomplete quality scores
        // 2. we add if defaultBaseQualities has a positive value
        if (defaultBaseQualities < 0) {
            return read;
        }
        byte reads[] = read.getReadBases();
        byte quals[] = read.getBaseQualities();
        if (quals == null || quals.length < reads.length) {
            byte new_quals[] = new byte[reads.length];
            for (int i = 0; i < reads.length; i++)
                new_quals[i] = defaultBaseQualities;
            read.setBaseQualities(new_quals);
        }
        return read;
    }

    private Function<SAMRecord, SAMRecord> makeReadTransform(){
        Function<SAMRecord, SAMRecord> adaptor = ReadClipper::hardClipAdaptorSequence;

        Function<SAMRecord, SAMRecord> f =
                 adaptor.andThen(ReadClipper::hardClipSoftClippedBases)
                .andThen(BaseRecalibrator::consolidateCigar)
                .andThen(this::setDefaultBaseQualities)
                .andThen(this::resetOriginalBaseQualities);

        return f;
    }
    /**
     * For each read at this locus get the various covariate values and increment that location in the map based on
     * whether or not the base matches the reference at this particular location
     */
    @Override
    public void apply( final SAMRecord originalRead, final ReferenceContext ref) {
        Function<SAMRecord, SAMRecord> transform = makeReadTransform();
        final SAMRecord read = transform.apply(originalRead);

        if( ReadUtils.isEmpty(read) ) { return; } // the whole read was inside the adaptor so skip it

        RecalUtils.parsePlatformForRead(read, RAC);
        if (!RecalUtils.isColorSpaceConsistent(RAC.SOLID_NOCALL_STRATEGY, read)) { // parse the solid color space and check for color no-calls
            return; // skip this read completely
        }

        //temporary HACK for the RMD tracker
        GenomeLocParser glp = new GenomeLocParser(getReferenceDictionary());
        final HACKRefMetaDataTracker metaDataTracker = new HACKRefMetaDataTracker(glp);

        final int[] isSNP = calculateIsSNP(read, ref, originalRead);
        final int[] isInsertion = calculateIsIndel(read, EventType.BASE_INSERTION);
        final int[] isDeletion = calculateIsIndel(read, EventType.BASE_DELETION);
        final int nErrors = nEvents(isSNP, isInsertion, isDeletion);

        // note for efficiency regions we don't compute the BAQ array unless we actually have
        // some error to marginalize over.  For ILMN data ~85% of reads have no error
        final byte[] baqArray = nErrors == 0 ? flatBAQArray(read) : calculateBAQArray(read);

        if( baqArray != null ) { // some reads just can't be BAQ'ed
            final ReadCovariates covariates = RecalUtils.computeCovariates(read, requestedCovariates);
            final boolean[] skip = calculateSkipArray(read, metaDataTracker); // skip known sites of variation as well as low quality and non-regular bases
            final double[] snpErrors = calculateFractionalErrorArray(isSNP, baqArray);
            final double[] insertionErrors = calculateFractionalErrorArray(isInsertion, baqArray);
            final double[] deletionErrors = calculateFractionalErrorArray(isDeletion, baqArray);

            // aggregate all of the info into our info object, and update the data
            final ReadRecalibrationInfo info = new ReadRecalibrationInfo(read, covariates, skip, snpErrors, insertionErrors, deletionErrors);
            recalibrationEngine.updateDataForRead(info);
        }

        accumulator++;
    }

    /**
     * Compute the number of mutational events across all hasEvent vectors
     *
     * Simply the sum of entries in hasEvents
     *
     * @param hasEvents a vector a vectors of 0 (no event) and 1 (has event)
     * @return the total number of events across all hasEvent arrays
     */
    protected static int nEvents(final int[]... hasEvents) {
        int n = 0;
        for ( final int[] hasEvent : hasEvents ) {
            n += MathUtils.sum(hasEvent);
        }
        return n;
    }

    private boolean[] calculateSkipArray( final SAMRecord read, final HACKRefMetaDataTracker metaDataTracker) {
        final byte[] bases = read.getReadBases();
        final boolean[] skip = new boolean[bases.length];
        final boolean[] knownSites = calculateKnownSites(read, metaDataTracker.getValues(RAC.knownSitesVCF, read));
        for( int iii = 0; iii < bases.length; iii++ ) {
            skip[iii] = !BaseUtils.isRegularBase(bases[iii]) || isLowQualityBase(read, iii) || knownSites[iii] || badSolidOffset(read, iii);
        }
        return skip;
    }

    protected boolean badSolidOffset( final SAMRecord read, final int offset ) {
        return ReadUtils.isSOLiDRead(read) && RAC.SOLID_RECAL_MODE != RecalUtils.SOLID_RECAL_MODE.DO_NOTHING && !RecalUtils.isColorSpaceConsistent(read, offset);
    }

    protected boolean[] calculateKnownSites( final SAMRecord read, final List<Feature> features) {
        final int readLength = read.getReadBases().length;
        final boolean[] knownSites = new boolean[readLength];
        Arrays.fill(knownSites, false);
        for( final Feature feat : features ) {
            int featureStartOnRead = ReadUtils.getReadCoordinateForReferenceCoordinate(ReadUtils.getSoftStart(read), read.getCigar(), feat.getStart(), ReadUtils.ClippingTail.LEFT_TAIL, true); // BUGBUG: should I use LEFT_TAIL here?
            if( featureStartOnRead == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
                featureStartOnRead = 0;
            }

            int featureEndOnRead = ReadUtils.getReadCoordinateForReferenceCoordinate(ReadUtils.getSoftStart(read), read.getCigar(), feat.getEnd(), ReadUtils.ClippingTail.LEFT_TAIL, true);
            if( featureEndOnRead == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
                featureEndOnRead = readLength;
            }

            if( featureStartOnRead > readLength ) {
                featureStartOnRead = featureEndOnRead = readLength;
            }

            Arrays.fill(knownSites, Math.max(0, featureStartOnRead), Math.min(readLength, featureEndOnRead + 1), true);
        }
        return knownSites;
    }

    // BUGBUG: can be merged with calculateIsIndel
    protected static int[] calculateIsSNP( final SAMRecord read, final ReferenceContext ref, final SAMRecord originalRead ) {
        final byte[] readBases = read.getReadBases();
        final byte[] refBases = Arrays.copyOfRange(ref.getBases(), read.getAlignmentStart() - originalRead.getAlignmentStart(), ref.getBases().length + read.getAlignmentEnd() - originalRead.getAlignmentEnd());
        final int[] snp = new int[readBases.length];
        int readPos = 0;
        int refPos = 0;
        for ( final CigarElement ce : read.getCigar().getCigarElements() ) {
            final int elementLength = ce.getLength();
            switch (ce.getOperator()) {
                case M:
                case EQ:
                case X:
                    for( int iii = 0; iii < elementLength; iii++ ) {
                        snp[readPos] = ( BaseUtils.basesAreEqual(readBases[readPos], refBases[refPos]) ? 0 : 1 );
                        readPos++;
                        refPos++;
                    }
                    break;
                case D:
                case N:
                    refPos += elementLength;
                    break;
                case I:
                case S: // ReferenceContext doesn't have the soft clipped bases!
                    readPos += elementLength;
                    break;
                case H:
                case P:
                    break;
                default:
                    throw new GATKException("Unsupported cigar operator: " + ce.getOperator());
            }
        }
        return snp;
    }

    protected static int[] calculateIsIndel( final SAMRecord read, final EventType mode ) {
        final int[] indel = new int[read.getReadBases().length];
        int readPos = 0;
        for ( final CigarElement ce : read.getCigar().getCigarElements() ) {
            final int elementLength = ce.getLength();
            switch (ce.getOperator()) {
                case M:
                case EQ:
                case X:
                case S:
                {
                    readPos += elementLength;
                    break;
                }
                case D:
                {
                    final int index = ( read.getReadNegativeStrandFlag() ? readPos : readPos - 1 );
                    updateIndel(indel, index, mode, EventType.BASE_DELETION);
                    break;
                }
                case I:
                {
                    final boolean forwardStrandRead = !read.getReadNegativeStrandFlag();
                    if( forwardStrandRead ) {
                        updateIndel(indel, readPos - 1, mode, EventType.BASE_INSERTION);
                    }
                    readPos += elementLength;
                    if( !forwardStrandRead ) {
                        updateIndel(indel, readPos, mode, EventType.BASE_INSERTION);
                    }
                    break;
                }
                case N:
                case H:
                case P:
                    break;
                default:
                    throw new GATKException("Unsupported cigar operator: " + ce.getOperator());
            }
        }
        return indel;
    }

    private static void updateIndel(final int[] indel, final int index, final EventType mode, final EventType requiredMode) {
        if ( mode == requiredMode && index >= 0 && index < indel.length )
            // protect ourselves from events at the start or end of the read (1D3M or 3M1D)
            indel[index] = 1;
    }

    protected static double[] calculateFractionalErrorArray( final int[] errorArray, final byte[] baqArray ) {
        if(errorArray.length != baqArray.length ) {
            throw new GATKException("Array length mismatch detected. Malformed read?");
        }

        final int BLOCK_START_UNSET = -1;

        final double[] fractionalErrors = new double[baqArray.length];
        Arrays.fill(fractionalErrors, 0.0);
        boolean inBlock = false;
        int blockStartIndex = BLOCK_START_UNSET;
        int iii;
        for( iii = 0; iii < fractionalErrors.length; iii++ ) {
            if( baqArray[iii] == NO_BAQ_UNCERTAINTY ) {
                if( !inBlock ) {
                    fractionalErrors[iii] = (double) errorArray[iii];
                } else {
                    calculateAndStoreErrorsInBlock(iii, blockStartIndex, errorArray, fractionalErrors);
                    inBlock = false; // reset state variables
                    blockStartIndex = BLOCK_START_UNSET; // reset state variables
                }
            } else {
                inBlock = true;
                if( blockStartIndex == BLOCK_START_UNSET ) { blockStartIndex = iii; }
            }
        }
        if( inBlock ) {
            calculateAndStoreErrorsInBlock(iii-1, blockStartIndex, errorArray, fractionalErrors);
        }
        if( fractionalErrors.length != errorArray.length ) {
            throw new GATKException("Output array length mismatch detected. Malformed read?");
        }
        return fractionalErrors;
    }

    private static void calculateAndStoreErrorsInBlock( final int iii,
                                                        final int blockStartIndex,
                                                        final int[] errorArray,
                                                        final double[] fractionalErrors ) {
        int totalErrors = 0;
        for( int jjj = Math.max(0, blockStartIndex - 1); jjj <= iii; jjj++ ) {
            totalErrors += errorArray[jjj];
        }
        for( int jjj = Math.max(0, blockStartIndex - 1); jjj <= iii; jjj++ ) {
            fractionalErrors[jjj] = ((double) totalErrors) / ((double)(iii - Math.max(0, blockStartIndex - 1) + 1));
        }
    }

    /**
     * Create a BAQ style array that indicates no alignment uncertainty
     * @param read the read for which we want a BAQ array
     * @return a BAQ-style non-null byte[] counting NO_BAQ_UNCERTAINTY values
     * // TODO -- could be optimized avoiding this function entirely by using this inline if the calculation code above
     */
    protected  static byte[] flatBAQArray(final SAMRecord read) {
        final byte[] baq = new byte[read.getReadLength()];
        Arrays.fill(baq, NO_BAQ_UNCERTAINTY);
        return baq;
    }

    /**
     * Compute an actual BAQ array for read, based on its quals and the reference sequence
     * @param read the read to BAQ
     * @return a non-null BAQ tag array for read
     */
    private byte[] calculateBAQArray( final SAMRecord read ) {
        baq.baqRead(read,referenceDataSource, BAQ.CalculationMode.RECALCULATE, BAQ.QualityMode.ADD_TAG);
        return BAQ.getBAQTag(read);
    }

    @Override
    public Object onTraversalDone() {
        recalibrationEngine.finalizeData();

        logger.info("Calculating quantized quality scores...");
        quantizeQualityScores();

        logger.info("Writing recalibration report...");
        generateReport();
        logger.info("...done!");

        //logger.info("BaseRecalibrator was able to recalibrate " + result + " reads");
        return accumulator;
    }

    private RecalibrationTables getRecalibrationTable() {
        return recalibrationEngine.getFinalRecalibrationTables();
    }

    /**
     * go through the quality score table and use the # observations and the empirical quality score
     * to build a quality score histogram for quantization. Then use the QuantizeQual algorithm to
     * generate a quantization map (recalibrated_qual -> quantized_qual)
     */
    private void quantizeQualityScores() {
        quantizationInfo = new QuantizationInfo(getRecalibrationTable(), RAC.QUANTIZING_LEVELS);
    }

    private void generateReport() {
        RecalUtils.outputRecalibrationReport(RAC, quantizationInfo, getRecalibrationTable(), requestedCovariates, RAC.SORT_BY_ALL_COLUMNS);
    }
}