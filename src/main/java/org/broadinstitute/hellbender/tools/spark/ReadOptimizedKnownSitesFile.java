package org.broadinstitute.hellbender.tools.spark;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;

@CommandLineProgramProperties(
        summary = "Converts a VCF of known sites into an optimized format (serialized Java object format)",
        oneLineSummary = "Converts a VCF of known sites into an optimized format (serialized Kryo object format)",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class ReadOptimizedKnownSitesFile extends CommandLineProgram {
    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc="The input kryo file.")
    private String inputFile;

    @Override
    @SuppressWarnings("unchecked")
    protected Object doWork() {
        try {
            IntervalsSkipList<GATKVariant> variants = (IntervalsSkipList<GATKVariant>) deserialize(Files.newInputStream(IOUtils.getPath(inputFile)));
            System.out.println("Variants: " + variants.size());
        } catch (Throwable e) {
            throw new UserException("Problem writing file", e);
        }
        return null;
    }

    private static Object deserialize(InputStream in) throws IOException, ClassNotFoundException {
        Kryo kryo = new Kryo();
        Input input = new Input(in);
        return kryo.readClassAndObject(input);
    }
}
