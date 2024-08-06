package org.broadinstitute.hellbender.tools;

import htsjdk.beta.io.bundle.*;
import htsjdk.beta.plugin.registry.HaploidReferenceResolver;
import htsjdk.beta.plugin.variants.VariantsBundle;
import htsjdk.io.HtsPath;
import htsjdk.samtools.util.FileExtensions;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.*;

/**
 * Create a bundle (JSON) file for use with a GATK tool.
 *
 * Since most bundles will contain a primary resource plus at least one secondary resource (typically an index),
 * the tool will attempt to infer standard secondary resources(s) for a given primary resource if no secondary resource
 * is explicitly provided on the command line. Inferred secondary resources are automatically added to the resulting
 * bundle. Secondary resource inference can be suppressed by using the --suppress-resource-resolution argument.
 *
 * Each resource in a bundle must have an associated content type tag. Content types for each resource are either
 * specified on the command line via argument tags, or inferred by the tool. For the primary and secondary resources,
 * when no content type argument tag is provided, the tool will attempt to infer the content type from the file
 * extension. However, the content type for "other" resources (resources that are nether primary nor secondary resources)
 * are NEVER inferred, and must always include a content type argument tag.
 *
 * Bundle output file names must end with the suffix ".json".
 *
 * Common examples:
 *
 * VCF Bundles:
 *
 * 1) Create a resource bundle for a VCF. Let the tool determine the content types, and resolve the secondary resources
 * (which for vcfs is the companion index) automatically by finding a sibling index file. If the sibling file cannot
 * be found, an exception wil lbe thrown:
 *
 *    CreateBundle --primary path/to/my.vcf --output mybundle.json
 *
 * 2) Create a resource bundle for a VCF. Let the tool determine the content types, but suppress resolution of the secondary
 * resources (which for vcfs is the companion index). The resulting bundle will contain only the vcf resource:
 *
 *    CreateBundle --primary path/to/my.vcf --output mybundle.json
 *
 *  3) Create a resource bundle for a VCF. Let the tool determine the content type, but specify the secondary
 *  index resource explicitly (which suppresses secondary resolution). The resulting bundle will contain the vcf
 *  and index resources:
 *
 *     CreateBundle --primary path/to/my.vcf --secondary some/other/path/to/vcd.idx --output mybundle.json
 *
 * Reference bundles: create a bundle using explicitly provided values and content types for the primary and
 * secondary resources:
 *
 *     CreateBundle --primary: path/to/my.fa
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Create a bundle (JSON) file for use with a GATK tool",
        oneLineSummary = "Create a bundle (JSON) file for use with a GATK tool",
        programGroup = OtherProgramGroup.class
)
public class CreateBundle extends CommandLineProgram {
    protected static final Logger logger = LogManager.getLogger(CreateBundle.class);

    @Argument(fullName = StandardArgumentDefinitions.PRIMARY_RESOURCE_LONG_NAME,
            shortName = StandardArgumentDefinitions.PRIMARY_RESOURCE_SHORT_NAME,
            doc="Path to the primary bundle resource (content type will be inferred if no content type tag is specified)")
    GATKPath primaryResource;

    @Argument(fullName = StandardArgumentDefinitions.SECONDARY_RESOURCE_LONG_NAME,
            shortName = StandardArgumentDefinitions.SECONDARY_RESOURCE_SHORT_NAME,
            doc = "Path to a secondary bundle resource for" + StandardArgumentDefinitions.PRIMARY_RESOURCE_LONG_NAME +
                    "(usually an index). The content type will be inferred if no content type tag is specified. If no "+
                    " secondary resource is specified, standard secondary resources for the primary bundle resource " +
                    " will also automatically be inferred unless " + SUPPRESS_SECONDARY_RESOURCE_RESOLUTION_FULL_NAME +
                    " is specified.",
            optional = true)
    GATKPath secondaryResource;

    public static final String OTHER_RESOURCE_FULL_NAME = "other-resource";
    @Argument(fullName = OTHER_RESOURCE_FULL_NAME,
            shortName = OTHER_RESOURCE_FULL_NAME,
            doc = "Path to other bundle resources for " + StandardArgumentDefinitions.PRIMARY_RESOURCE_LONG_NAME
                    + ". The content is not inferred for \"other\" (non-primary and non-secondary) resources, so a" +
                    " content type tag MUST be provided for each \"other\" resource.",
            optional = true)
    List<GATKPath> otherResources;

    public static final String SUPPRESS_SECONDARY_RESOURCE_RESOLUTION_FULL_NAME = "suppress-resource-resolution";
    @Argument(fullName = SUPPRESS_SECONDARY_RESOURCE_RESOLUTION_FULL_NAME,
            doc ="Don't attempt to resolve the primary resource's secondary resources if no secondary resource is explicitly" +
                    " provided (defaults to false) ",
            optional = true)
    boolean suppressResourceResolution = false;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Path the output bundle file (must end with the suffix " + BundleJSON.BUNDLE_EXTENSION + ")")
    GATKPath outputBundlePath;

    private enum BundleType {
        VCF,
        REFERENCE,
        OTHER
    }
    private BundleType outputBundleType;

    @Override
    protected String[] customCommandLineValidation() {
        if (!outputBundlePath.toString().endsWith(BundleJSON.BUNDLE_EXTENSION)){
            return new String[]{ String.format("Output bundle path must end with the suffix %s", BundleJSON.BUNDLE_EXTENSION) };
        }
        outputBundleType = determinePrimaryContentType();
        return super.customCommandLineValidation();
    }

    @Override
    protected Object doWork() {
        try (final BufferedWriter writer = Files.newBufferedWriter(outputBundlePath.toPath(), StandardOpenOption.CREATE)) {
            final Bundle bundle = switch (outputBundleType) {
                case VCF -> createVCFBundle();
                case REFERENCE -> createHaploidReferenceBundle();
                case OTHER -> createOtherBundle();
            };
            writer.write(BundleJSON.toJSON(bundle));
        } catch (final IOException e) {
            throw new RuntimeException(String.format("Failed writing bundle to output %s", outputBundlePath), e);
        }
        return null;
    }

    private BundleType determinePrimaryContentType() {
        BundleType bundleType;

        // determine the type of bundle to create; consult the tag attributes if any, otherwise try to infer from the
        // primary resource's file extension
        final String primaryContentTag = primaryResource.getTag();
        if (primaryContentTag != null && !primaryContentTag.isEmpty()) {
            if (primaryContentTag.equals(BundleResourceType.CT_VARIANT_CONTEXTS)) {
                bundleType = BundleType.VCF;
            } else if (primaryContentTag.equals(BundleResourceType.CT_HAPLOID_REFERENCE)) {
                bundleType = BundleType.REFERENCE;
            } else {
                logger.info(String.format("Primary input content type %s for %s not recognized. A bundle will be created using content types from the provided argument tags.",
                        primaryContentTag,
                        primaryResource));
                bundleType = BundleType.OTHER;
            }
        } else {
            logger.info(String.format("A content type for the primary input was not provided. Attempting to infer the content type from the %s extension.", primaryResource));
            bundleType = inferPrimaryContentType(primaryResource);
        }
        return bundleType;
    }

    private BundleType inferPrimaryContentType(final GATKPath primaryResource) {
        logger.info("Attempting to infer bundle content type from file extension.");
        if (FileExtensions.VCF_LIST.stream().anyMatch(ext -> primaryResource.hasExtension(ext))) {
            return BundleType.VCF;
        } else if (FileExtensions.FASTA.stream().anyMatch(ext -> primaryResource.hasExtension(ext))) {
            return BundleType.REFERENCE;
        } else {
            throw new IllegalArgumentException(String.format("Unable to infer bundle content type from file extension %s. A content type must be provided as part of the argument.", primaryResource));
        }
    }

    private Bundle createVCFBundle() {
        final Collection<BundleResource> bundleResources = new ArrayList<>();

        bundleResources.add(new IOPathResource(primaryResource, BundleResourceType.CT_VARIANT_CONTEXTS));
        if (secondaryResource != null) {
            final String secondaryContentType = secondaryResource.getTag();
            if (secondaryContentType == null) {
                logger.info(String.format("A content type for the secondary resource was not provided. Assuming %s is an index.", secondaryResource));
                bundleResources.add(new IOPathResource(secondaryResource, BundleResourceType.CT_VARIANTS_INDEX));
            } else {
                bundleResources.add(new IOPathResource(secondaryResource, secondaryContentType));
            }
        } else if (!suppressResourceResolution) {
            // secondary resources is null, and resource resolution suppression is off
            final Optional<GATKPath> indexPath = VariantsBundle.resolveIndex(primaryResource, GATKPath::new);
            if (indexPath.isEmpty()) {
                throw new IllegalArgumentException(
                        String.format(
                                "Could not infer an index for %s, you must either specify the index path as a secondary input on the command line or specify the %s argument.",
                                primaryResource.getRawInputString(),
                                SUPPRESS_SECONDARY_RESOURCE_RESOLUTION_FULL_NAME));
            }
            bundleResources.add(new IOPathResource(indexPath.get(), BundleResourceType.CT_VARIANTS_INDEX));
        }
        if (otherResources != null) {
            for (final GATKPath otherInput : otherResources) {
                final String otherContentType = otherInput.getTag();
                if (otherContentType == null) {
                    throw new IllegalArgumentException(
                            String.format(
                                    "A content must be provided for \"other\" input %s.",
                                    otherInput.getRawInputString()));
                } else {
                    bundleResources.add(new IOPathResource(otherInput, otherContentType));
                }
            }
        }
        return new VariantsBundle(bundleResources);
    }

    private Bundle createHaploidReferenceBundle() {
        final Collection<BundleResource> bundleResources = new ArrayList<>();
        if (secondaryResource == null && otherResources == null) {
            // infer dictionary and index
            return HaploidReferenceResolver.referenceBundleFromFastaPath(primaryResource, GATKPath::new);
        }

        bundleResources.add(new IOPathResource(primaryResource, primaryResource.getTag()));
        if (secondaryResource != null) {
            final String secondaryContentType = secondaryResource.getTag();
            if (secondaryContentType == null) {
                throw new IllegalArgumentException(String.format("A content type for the secondary input must be provided."));
            } else {
                bundleResources.add(new IOPathResource(secondaryResource, secondaryContentType));
            }
        }
        if (otherResources != null) {
            for (final GATKPath otherInput : otherResources) {
                final String otherContentType = otherInput.getTag();
                if (otherContentType == null) {
                    throw new IllegalArgumentException(
                            String.format(
                                    "A content type must be provided for \"other\" input %s.",
                                    otherInput.getRawInputString()));
                } else {
                    bundleResources.add(new IOPathResource(otherInput, otherContentType));
                }
            }
        }
        return new Bundle(primaryResource.getTag(), bundleResources);
    }

    private Bundle createOtherBundle() {
        final Collection<BundleResource> bundleResources = new ArrayList<>();
        bundleResources.add(new IOPathResource(primaryResource, primaryResource.getTag()));
        if (secondaryResource != null) {
            final String secondaryContentType = secondaryResource.getTag();
            if (secondaryContentType == null) {
                throw new IllegalArgumentException(String.format("A content type for the secondary input must be provided."));
            } else {
                bundleResources.add(new IOPathResource(secondaryResource, secondaryContentType));
            }
        }
        if (otherResources != null) {
            for (final GATKPath otherInput : otherResources) {
                final String otherContentType = otherInput.getTag();
                if (otherContentType == null) {
                    throw new IllegalArgumentException(
                            String.format(
                                    "A content type must be provided for \"other\" input %s.",
                                    otherInput.getRawInputString()));
                } else {
                    bundleResources.add(new IOPathResource(otherInput, otherContentType));
                }
            }
        }
        return new Bundle(primaryResource.getTag(), bundleResources);
    }
}
