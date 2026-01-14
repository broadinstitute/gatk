package org.broadinstitute.hellbender.tools.reprocessing;

import it.unimi.dsi.util.FrontCodedStringList;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.WorkflowOutput;
import org.broadinstitute.barclay.argparser.WorkflowProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;

@CommandLineProgramProperties(
        summary = "Filters reads by excluding read names from a front-coded binary list (.fcl) or plain text file. " +
                  "Front-coded lists use ~85% less memory than plain text (1.5 GB vs 10.5 GB for 38M names). " +
                  "Reads whose names appear in the exclusion list are discarded; all others are written to output.",
        oneLineSummary = "Filter reads by excluding names from a list",
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
@WorkflowProperties
public class FilterReadsByNameList extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Write filtered output to this BAM/SAM/CRAM file")
    @WorkflowOutput(optionalCompanions = {StandardArgumentDefinitions.OUTPUT_INDEX_COMPANION})
    public GATKPath output;

    @Argument(fullName = "exclude-read-list",
            shortName = "XRL",
            doc = "File containing read names to exclude. Supports .fcl (front-coded binary) or .txt (plain text) formats. " +
                  "Front-coded format dramatically reduces memory usage for large lists (e.g., 1.5 GB vs 10.5 GB for 38M names).")
    public GATKPath excludeReadList;

    private SAMFileGATKReadWriter outputWriter;
    private FrontCodedStringList readNamesToExcludeFCL = null;
    private java.util.Set<String> readNamesToExcludeTXT = null;
    private boolean isFrontCodedList = false;
    
    private long totalReadsProcessed = 0;
    private long readsExcluded = 0;
    private long readsWritten = 0;

    @Override
    public void onTraversalStart() {
        // Determine if input is coordinate-sorted (typical for CRAM)
        final boolean isPreSorted = getHeaderForReads().getSortOrder() != htsjdk.samtools.SAMFileHeader.SortOrder.unsorted;
        outputWriter = createSAMWriter(output, isPreSorted);
        
        // Load exclusion list based on file extension
        final String path = excludeReadList.getRawInputString();
        
        if (path.endsWith(".fcl")) {
            loadFrontCodedList();
            isFrontCodedList = true;
            logger.info("Loaded front-coded read name list from " + path);
            logger.info("List contains " + readNamesToExcludeFCL.size() + " read names");
        } else if (path.endsWith(".txt")) {
            loadTextList();
            isFrontCodedList = false;
            logger.info("Loaded plain text read name list from " + path);
            logger.info("List contains " + readNamesToExcludeTXT.size() + " read names");
        } else {
            throw new UserException("Exclusion list file must have .fcl or .txt extension: " + path);
        }
    }

    private void loadFrontCodedList() {
        try (ObjectInputStream ois = new ObjectInputStream(new FileInputStream(excludeReadList.toPath().toFile()))) {
            readNamesToExcludeFCL = (FrontCodedStringList) ois.readObject();
        } catch (IOException | ClassNotFoundException e) {
            throw new UserException("Failed to load front-coded list from " + excludeReadList.getRawInputString(), e);
        }
    }

    private void loadTextList() {
        readNamesToExcludeTXT = new java.util.HashSet<>();
        try (java.io.BufferedReader reader = new java.io.BufferedReader(
                new java.io.FileReader(excludeReadList.toPath().toFile()))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (!line.isEmpty()) {
                    readNamesToExcludeTXT.add(line);
                }
            }
        } catch (IOException e) {
            throw new UserException("Failed to load text list from " + excludeReadList.getRawInputString(), e);
        }
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        totalReadsProcessed++;
        
        final boolean shouldExclude;
        if (isFrontCodedList) {
            // Binary search: indexOf returns -1 if not found, >= 0 if found
            shouldExclude = readNamesToExcludeFCL.indexOf(read.getName()) >= 0;
        } else {
            shouldExclude = readNamesToExcludeTXT.contains(read.getName());
        }
        
        if (shouldExclude) {
            readsExcluded++;
        } else {
            outputWriter.addRead(read);
            readsWritten++;
        }
        
        // Log progress every 1M reads
        if (totalReadsProcessed % 1_000_000 == 0) {
            logger.info(String.format("Processed %,d reads: excluded %,d, wrote %,d",
                    totalReadsProcessed, readsExcluded, readsWritten));
        }
    }

    @Override
    public void closeTool() {
        if (outputWriter != null) {
            outputWriter.close();
        }
        
        // Final statistics
        logger.info("======================================");
        logger.info("FilterReadsByNameList Complete");
        logger.info(String.format("Total reads processed: %,d", totalReadsProcessed));
        logger.info(String.format("Reads excluded: %,d (%.2f%%)", 
                readsExcluded, 100.0 * readsExcluded / totalReadsProcessed));
        logger.info(String.format("Reads written: %,d (%.2f%%)", 
                readsWritten, 100.0 * readsWritten / totalReadsProcessed));
        logger.info("======================================");
    }
}
