package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * TODO
 */
@CommandLineProgramProperties(
        summary = "Apply a score cutoff to filter variants based on a recalibration table",
        oneLineSummary = " Apply a score cutoff to filter variants based on a recalibration table",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public class ApplyVariantAnnotationScores extends MultiVariantWalker {
    public final static double MIN_ACCEPTABLE_LOD_SCORE = -20000.0; // from VariantRecalibratorEngine

    protected static final String LOW_VQSLOD_FILTER_NAME = "LOW_VQSLOD";
    private final double DEFAULT_VQSLOD_CUTOFF = -1.; // all PASS for IsolationForest

    private boolean foundSNPTranches = false;
    private boolean foundINDELTranches = false;

    /////////////////////////////
    // Inputs
    /////////////////////////////
    @Argument(fullName="recal-file", doc="The input recal file used by ApplyVQSR", optional=false)
    private FeatureInput<VariantContext> recal;

    @Argument(fullName="tranches-file", doc="The input tranches file describing where to cut the data", optional=true)
    private GATKPath TRANCHES_FILE;

    /////////////////////////////
    // Outputs
    /////////////////////////////

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output filtered and recalibrated VCF file in which each variant is annotated with its VQSLOD value", optional=false)
    private GATKPath output;

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    @Argument(fullName="truth-sensitivity-filter-level", shortName="ts-filter-level", doc="The truth sensitivity level at which to start filtering", optional=true)
    private Double TS_FILTER_LEVEL = null;

    /**
     *  Filter the input file based on allele-specific recalibration data.  See tool docs for site-level and allele-level filtering details.
     *  Requires a .recal file produced using an allele-specific run of VariantRecalibrator.
     */
    @Argument(fullName="use-allele-specific-annotations", shortName="AS", doc="If specified, the tool will attempt to apply a filter to each allele based on the input tranches and allele-specific .recal file.", optional=true)
    private boolean useASannotations = false;

    @Advanced
    @Argument(fullName="lod-score-cutoff", doc="The VQSLOD score below which to start filtering", optional=true)
    protected Double VQSLOD_CUTOFF = null;

    /**
     * For this to work properly, the --ignore-filter argument should also be applied to the VariantRecalibration command.
     */
    @Argument(fullName="ignore-filter", doc="If specified, the recalibration will be applied to variants marked as filtered by the specified filter name in the input VCF file", optional=true)
    private List<String> IGNORE_INPUT_FILTERS = new ArrayList<>();

    @Argument(fullName="ignore-all-filters", doc="If specified, the variant recalibrator will ignore all input filters. Useful to rerun the VQSR from a filtered output file.", optional=true)
    private boolean IGNORE_ALL_FILTERS = false;

    @Argument(fullName="exclude-filtered", doc="Don't output filtered loci after applying the recalibration", optional=true)
    private boolean EXCLUDE_FILTERED = false;

    @Argument(fullName = "mode", shortName = "mode", doc = "Recalibration mode to employ: 1.) SNP for recalibrating only SNPs (emitting indels untouched in the output VCF); 2.) INDEL for indels; and 3.) BOTH for recalibrating both SNPs and indels simultaneously.", optional=true)
    private VariantTypeMode MODE = VariantTypeMode.SNP;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private VariantContextWriter vcfWriter;
    final private List<TruthSensitivityTranche> tranches = new ArrayList<>();
    final private Set<String> ignoreInputFilterSet = new TreeSet<>();
    final static private String listPrintSeparator = ",";
    final static private String trancheFilterString = "VQSRTranche";
    final static private String arrayParseRegex = "[\\[\\]\\s]";
    final static private String emptyStringValue = "NA";
    final static private String emptyFloatValue = "NaN";


    //---------------------------------------------------------------------------------------------------------------
    //
    // onTraversalStart
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public void onTraversalStart() {
        if( TS_FILTER_LEVEL != null ) {
            try {
                for (final TruthSensitivityTranche t : TruthSensitivityTranche.readTranches(TRANCHES_FILE)) {
                    if (t.getTargetTruthSensitivity() >= TS_FILTER_LEVEL) {
                        tranches.add(t);
                    }
                    logger.info(String.format("Read tranche " + t));
                }
            }
            catch(IOException e ) {
                throw new UserException.CouldNotReadInputFile(TRANCHES_FILE, "Can't read tranches file", e);
            }
            Collections.reverse(tranches); // this algorithm wants the tranches ordered from best (lowest truth sensitivity) to worst (highest truth sensitivity)
        }

        if( IGNORE_INPUT_FILTERS != null ) {
            ignoreInputFilterSet.addAll( IGNORE_INPUT_FILTERS );
        }

        // setup the header fields
        VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();

        final Set<VCFHeaderLine> hInfo = new HashSet<>(inputHeaders);
        hInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.VQS_LOD_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.POSITIVE_LABEL_KEY));
        hInfo.add(GATKVCFHeaderLines.getFilterLine(VCFConstants.PASSES_FILTERS_v4));
        if (useASannotations) {
            hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_FILTER_STATUS_KEY));
            hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_VQS_LOD_KEY));
        }

        checkForPreviousApplyRecalRun(Collections.unmodifiableSet(inputHeaders));

        final TreeSet<String> samples = new TreeSet<>();
        samples.addAll(inputHeader.getGenotypeSamples());

        //generate headers from tranches file
        //TODO: throw away old tranche headers if we're ignoring filters
        //TODO: TS_FILTER_LEVEL and VQSLOD_CUTOFF should use mutex argument declaration
        if( TS_FILTER_LEVEL != null ) {
            // if the user specifies both ts_filter_level and lodCutoff then throw a user error
            if( VQSLOD_CUTOFF != null ) {
                throw new UserException("Arguments --truth-sensitivity-filter-level and --lod-score-cutoff are mutually exclusive. Please only specify one option.");
            }

            if( tranches.size() >= 2 ) {
                for( int iii = 0; iii < tranches.size() - 1; iii++ ) {
                    final TruthSensitivityTranche t = tranches.get(iii);
                    hInfo.add(new VCFFilterHeaderLine(t.getName(), String.format("Truth sensitivity tranche level for " + t.getMode().toString() + " model at VQS Lod: " + t.getMinScore() + " <= x < " + tranches.get(iii+1).getMinScore())));
                }
            }
            if( tranches.size() >= 1 ) {
                hInfo.add(new VCFFilterHeaderLine(tranches.get(0).getName() + "+", String.format("Truth sensitivity tranche level for " + tranches.get(0).getMode().toString() + " model at VQS Lod < " + tranches.get(0).getMinScore())));
            } else {
                throw new UserException("No tranches were found in the file or were above the truth sensitivity filter level " + TS_FILTER_LEVEL);
            }

            logger.info("Keeping all variants in tranche " + tranches.get(tranches.size()-1));
        } else {
            if( VQSLOD_CUTOFF == null ) {
                VQSLOD_CUTOFF = DEFAULT_VQSLOD_CUTOFF;
            }
            hInfo.add(new VCFFilterHeaderLine(LOW_VQSLOD_FILTER_NAME, "VQSLOD < " + VQSLOD_CUTOFF));
            logger.info("Keeping all variants with VQSLOD >= " + VQSLOD_CUTOFF);
        }

        hInfo.addAll(getDefaultToolVCFHeaderLines());
        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter = createVCFWriter(output);
        vcfWriter.writeHeader(vcfHeader);
    }

    private boolean trancheIntervalIsValid(final String sensitivityLimits) {
        final String[] vals = sensitivityLimits.split("to");
        if(vals.length != 2)
            return false;
        try {
            double lowerLimit = Double.parseDouble(vals[0]);
            double upperLimit = Double.parseDouble(vals[1].replace("+",""));    //why does our last tranche end with 100+? Is there anything greater than 100 percent?  Really???
        }
        catch(NumberFormatException e) {
            throw new UserException("Poorly formatted tranche filter name does not contain two sensitivity interval end points.");
        }
        return true;
    }

    /**
     * Check the filter declarations in the input VCF header to see if any ApplyVQSR mode has been run
     * Here we assume that the tranches are named with a specific format: VQSRTranche[SNP|INDEL][lowerLimit]to[upperLimit]
     * @param inputHeaders
     */
    private void checkForPreviousApplyRecalRun(final Set<VCFHeaderLine> inputHeaders) {
        for(final VCFHeaderLine header : inputHeaders) {
            if(header instanceof VCFFilterHeaderLine) {
                final String filterName = ((VCFFilterHeaderLine)header).getID();
                //TODO: clean up these magic numbers
                if(filterName.length() < 12 || !filterName.substring(0, 11).equalsIgnoreCase(trancheFilterString)) {
                    continue;
                }
                if(filterName.charAt(11) == 'S') {
                    //for SNP tranches, get sensitivity limit
                    final String sensitivityLimits = filterName.substring(14);
                    if(trancheIntervalIsValid(sensitivityLimits))
                        foundSNPTranches = true;
                }
                else if(filterName.charAt(11) == 'I') {
                    //for INDEL tranches, get sensitivity limit
                    final String sensitivityLimits = filterName.substring(16);
                    if(trancheIntervalIsValid(sensitivityLimits))
                        foundINDELTranches = true;
                }
            }
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // apply
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext ref, final FeatureContext featureContext) {

        final List<VariantContext> recals =  featureContext.getValues(recal, vc.getStart());
        final boolean evaluateThisVariant = useASannotations || VariantDataManager.checkVariationClass( vc, MODE );

        //vc.isNotFiltered is true for PASS; vc.filtersHaveBeenApplied covers PASS and filters
        final boolean variantIsNotFiltered = IGNORE_ALL_FILTERS || vc.isNotFiltered() ||
                (!ignoreInputFilterSet.isEmpty() && ignoreInputFilterSet.containsAll(vc.getFilters()));

        if( evaluateThisVariant && variantIsNotFiltered) {
            String filterString;
            final VariantContextBuilder builder = new VariantContextBuilder(vc);
            if (!useASannotations) {
                filterString = doSiteSpecificFiltering(vc, recals, builder);
            }
            else {  //allele-specific mode
                filterString = doAlleleSpecificFiltering(vc, recals, builder);
            }

            //for both non-AS and AS modes:
            if( filterString.equals(VCFConstants.PASSES_FILTERS_v4) ) {
                builder.passFilters();
            } else if(filterString.equals(VCFConstants.UNFILTERED)) {
                builder.unfiltered();
            } else {
                builder.filters(filterString);
            }

            final VariantContext outputVC = builder.make();
            if( !EXCLUDE_FILTERED || outputVC.isNotFiltered() ) {
                vcfWriter.add( outputVC );
            }
        } else { // valid VC but not compatible with this mode, so just emit the variant untouched
            vcfWriter.add( vc );
        }
    }

    public double parseFilterLowerLimit(final String trancheFilter) {
        final Pattern pattern = Pattern.compile("VQSRTranche\\S+(\\d+\\.\\d+)to(\\d+\\.\\d+)");
        final Matcher m = pattern.matcher(trancheFilter);
        return m.find() ? Double.parseDouble(m.group(1)) : -1;
    }

    /**
     * Generate the VCF filter string for this record based on the ApplyVQSR modes run so far
     * @param vc the input VariantContext (with at least one ApplyVQSR mode already run)
     * @param bestLod best LOD from the alleles we've seen in this recalibration mode
     * @return the String to use as the VCF filter field
     */
    protected String generateFilterStringFromAlleles(final VariantContext vc, final double bestLod) {
        String filterString = ".";

        final boolean bothModesWereRun = (MODE == VariantTypeMode.SNP && foundINDELTranches) || (MODE == VariantTypeMode.INDEL && foundSNPTranches);
        final boolean onlyOneModeNeeded = !vc.isMixed() && VariantDataManager.checkVariationClass( vc, MODE );

        //if both SNP and INDEL modes have not yet been run (and need to be), leave this variant as unfiltered and add the filters for the alleles in this mode to the INFO field
        if (!bothModesWereRun && !onlyOneModeNeeded) {
            return VCFConstants.UNFILTERED;
        }

        //if both SNP and INDEL modes have been run or the site is not mixed, generate a filter string for this site based on both models
        //pull out the allele filter status from the info field (there may be more than one entry in the list if there were multiple snp/indel alleles assessed in the other mode)
        final String prevFilterStatus = vc.getAttributeAsString(GATKVCFConstants.AS_FILTER_STATUS_KEY, null);

        //if this site hasn't had a filter applied yet
        if (prevFilterStatus != null && !prevFilterStatus.equals(VCFConstants.UNFILTERED)) {
            final String prevAllelesFilterStatusString = vc.getAttributeAsString(GATKVCFConstants.AS_FILTER_STATUS_KEY, null);
            final String[] prevAllelesFilterStatusList = prevAllelesFilterStatusString.split(listPrintSeparator);
            //start with the current best allele filter as the most lenient filter across all modes and all alleles
            String mostLenientFilterName = generateFilterString(bestLod);
            //if the current mode's best allele passes the tranche filter, then let the whole site pass
            if (mostLenientFilterName.equals(VCFConstants.PASSES_FILTERS_v4)) {
                filterString = mostLenientFilterName;
            }
            //if the current mode's best allele does not pass the tranche filter, compare the most lenient filter of this mode with those from the previous mode
            else {
                double mostLenientSensitivityLowerLimit = parseFilterLowerLimit(mostLenientFilterName);
                for (int i = 0; i < prevAllelesFilterStatusList.length; i++) {
                    final String alleleFilterString = prevAllelesFilterStatusList[i].replaceAll(arrayParseRegex, "").trim();
                    //if any allele from the previous mode passed the tranche filter, then let the whole site pass
                    if (alleleFilterString.equals(VCFConstants.PASSES_FILTERS_v4)) { //this allele is PASS
                        mostLenientFilterName = alleleFilterString;
                        break;
                    }
                    //if there's no PASS, then we need to parse the filters to find out how lenient they are
                    else {
                        final double alleleLowerLimit = parseFilterLowerLimit(alleleFilterString);
                        if (alleleLowerLimit == -1)
                            continue;
                        if (alleleLowerLimit < mostLenientSensitivityLowerLimit) {
                            mostLenientSensitivityLowerLimit = alleleLowerLimit;
                            mostLenientFilterName = alleleFilterString;
                        }
                    }

                }
                filterString = mostLenientFilterName;
            }
        }

        //if both modes have been run, but the previous mode didn't apply a filter, use the current mode's best allele VQSLOD filter (shouldn't get run, but just in case)
        else {
            filterString = generateFilterString(bestLod);
        }

        return filterString;
    }

    /**
     * Generate the VCF filter string for this record based on the provided lod score
     * @param lod non-null double
     * @return the String to use as the VCF filter field
     */
    protected String generateFilterString( final double lod ) {
        String filterString = null;
        if( TS_FILTER_LEVEL != null ) {
            for( int i = tranches.size() - 1; i >= 0; i-- ) {
                final TruthSensitivityTranche tranche = tranches.get(i);
                if( lod >= tranche.getMinScore() ) {
                    if( i == tranches.size() - 1 ) {
                        filterString = VCFConstants.PASSES_FILTERS_v4;
                    } else {
                        filterString = tranche.getName();
                    }
                    break;
                }
            }

            if( filterString == null ) {
                filterString = tranches.get(0).getName()+"+";
            }
        } else {
            filterString = ( lod < VQSLOD_CUTOFF ? LOW_VQSLOD_FILTER_NAME : VCFConstants.PASSES_FILTERS_v4 );
        }

        return filterString;
    }

    private VariantContext getMatchingRecalVC(final VariantContext target, final List<VariantContext> recalVCs, final Allele allele) {
        for( final VariantContext recalVC : recalVCs ) {
            if ( target.getEnd() == recalVC.getEnd() ) {
                if (!useASannotations)
                    return recalVC;
                else if (allele.equals(recalVC.getAlternateAllele(0)))
                    return recalVC;
            }
        }
        return null;
    }

    /**
     *
     * @param altIndex current alt allele
     * @param prevLodList lods from previous ApplyVQSR run
     * @param prevASfiltersList AS_filters from previous ApplyVQSR run
     * @param lodString
     * @param AS_filterString
     */
    private void updateAnnotationsWithoutRecalibrating(final int altIndex, final String[] prevLodList, final String[] prevASfiltersList,
                                                       final List<String> lodString, final List<String> AS_filterString) {
        if (foundINDELTranches || foundSNPTranches) {
            if (altIndex < prevLodList.length) {
                lodString.add(prevLodList[altIndex].replaceAll(arrayParseRegex, "").trim());
                AS_filterString.add(prevASfiltersList[altIndex].replaceAll(arrayParseRegex, "").trim());
            }
        } else { //if the other allele type hasn't been processed yet, make sure there are enough entries
            lodString.add(emptyFloatValue);
            AS_filterString.add(emptyStringValue);
        }
    }

    /**
     * Calculate the allele-specific filter status of vc
     * @param vc
     * @param recals
     * @param builder   is modified by adding attributes
     * @return a String with the filter status for this site
     */
    private String doAlleleSpecificFiltering(final VariantContext vc, final List<VariantContext> recals, final VariantContextBuilder builder) {
        double bestLod = MIN_ACCEPTABLE_LOD_SCORE;
        final List<String> lodStrings = new ArrayList<>();
        final List<String> AS_filterStrings = new ArrayList<>();

        String[] prevLodList = null;
        String[] prevASfiltersList = null;

        //get VQSR annotations from previous run of ApplyVQSR, if applicable
        if(foundINDELTranches || foundSNPTranches) {
            final String prevLodString = vc.getAttributeAsString(GATKVCFConstants.AS_VQS_LOD_KEY,"");
            prevLodList = prevLodString.isEmpty()? new String[0] : prevLodString.split(listPrintSeparator);
            final String prevASfilters = vc.getAttributeAsString(GATKVCFConstants.AS_FILTER_STATUS_KEY,"");
            prevASfiltersList = prevASfilters.isEmpty()? new String[0] : prevASfilters.split(listPrintSeparator);
        }

        //for each allele in the current VariantContext
        for (int altIndex = 0; altIndex < vc.getNAlleles()-1; altIndex++) {
            final Allele allele = vc.getAlternateAllele(altIndex);

            //if the current allele is not part of this recalibration mode, add its annotations to the list and go to the next allele
            if (!VariantDataManager.checkVariationClass(vc, allele, MODE)) {
                updateAnnotationsWithoutRecalibrating(altIndex, prevLodList, prevASfiltersList, lodStrings, AS_filterStrings);
                continue;
            }

            //if the current allele does need to have recalibration applied...

            //initialize allele-specific VQSR annotation data with values for spanning deletion
            String alleleLodString = emptyFloatValue;
            String alleleFilterString = emptyStringValue;

            //if it's not a spanning deletion, replace those allele strings with the real values
            if (!GATKVCFConstants.isSpanningDeletion(allele)) {
                VariantContext recalDatum = getMatchingRecalVC(vc, recals, allele);
                if (recalDatum == null) {
                    throw new UserException("Encountered input allele which isn't found in the input recal file. Please make sure VariantRecalibrator and ApplyVQSR were run on the same set of input variants with flag -AS. First seen at: " + vc);
                }

                //compare VQSLODs for all alleles in the current mode for filtering later
                final double lod = recalDatum.getAttributeAsDouble(GATKVCFConstants.VQS_LOD_KEY, MIN_ACCEPTABLE_LOD_SCORE);
                if (lod > bestLod)
                    bestLod = lod;

                alleleLodString = String.format("%.4f", lod);
                alleleFilterString = generateFilterString(lod);

                if(recalDatum != null) {
                    if (recalDatum.hasAttribute(GATKVCFConstants.POSITIVE_LABEL_KEY))
                        builder.attribute(GATKVCFConstants.POSITIVE_LABEL_KEY, true);
                }
            }

            //append per-allele VQSR annotations
            lodStrings.add(alleleLodString);
            AS_filterStrings.add(alleleFilterString);
        }

        // Annotate the new record with its VQSLOD, AS_FilterStatus, and the worst performing annotation
        if(!AS_filterStrings.isEmpty() )
            builder.attribute(GATKVCFConstants.AS_FILTER_STATUS_KEY, AnnotationUtils.encodeStringList(AS_filterStrings));
        if(!lodStrings.isEmpty())
            builder.attribute(GATKVCFConstants.AS_VQS_LOD_KEY, AnnotationUtils.encodeStringList(lodStrings));

        return generateFilterStringFromAlleles(vc, bestLod);
    }

    /**
     * Calculate the filter status for a given VariantContext using the combined data from all alleles at a site
     * @param vc
     * @param recals
     * @param builder   is modified by adding attributes
     * @return a String with the filter status for this site
     */
    private String doSiteSpecificFiltering(final VariantContext vc, final List<VariantContext> recals, final VariantContextBuilder builder) {
        VariantContext recalDatum = getMatchingRecalVC(vc, recals, null);
        if( recalDatum == null ) {
            throw new UserException("Encountered input variant which isn't found in the input recal file. Please make sure VariantRecalibrator and ApplyVQSR were run on the same set of input variants. First seen at: " + vc );
        }

        final String lodString = recalDatum.getAttributeAsString(GATKVCFConstants.VQS_LOD_KEY, null);
        if( lodString == null ) {
            throw new UserException("Encountered a malformed record in the input recal file. There is no lod for the record at: " + vc );
        }
        final double lod;
        try {
            lod = Double.valueOf(lodString);
        } catch (NumberFormatException e) {
            throw new UserException("Encountered a malformed record in the input recal file. The lod is unreadable for the record at: " + vc );
        }

        builder.attribute(GATKVCFConstants.VQS_LOD_KEY, String.format("%.4f", lod));
        if(recalDatum != null) {
            if (recalDatum.hasAttribute(GATKVCFConstants.POSITIVE_LABEL_KEY))
                builder.attribute(GATKVCFConstants.POSITIVE_LABEL_KEY, true);
        }

        return generateFilterString(lod);
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}

