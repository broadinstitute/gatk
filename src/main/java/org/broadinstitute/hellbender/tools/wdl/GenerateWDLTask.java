package org.broadinstitute.hellbender.tools.wdl;

import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import picard.cmdline.programgroups.OtherProgramGroup;
import scala.Char;

import javax.annotation.Resource;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;


@CommandLineProgramProperties(
        summary="generates WDL task definitions for GATK bundled tools",
        oneLineSummary = "generates WDL task definition for GATK bundled tools",
        programGroup = OtherProgramGroup.class
)

public class GenerateWDLTask extends CommandLineProgram {

    @Argument(shortName="tool", fullName="tool-name")
    public String toolName;

    @Argument(shortName= StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
              fullName=StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              optional = true)
    public File output;

    @Override
    protected Object doWork() {

        final Map<String, Class<?>> toolsByName = Main.findCommandLineProgramClasses(Main.getPackageList(), Collections.emptyList());
        final Class<?> toolClass = toolsByName.get(toolName);
        if (toolClass == null) {
            throw new UserException(Main.getUnknownCommandMessage(new HashSet<>(toolsByName.values()), toolName));
        }
        final Object toolObj;
        try {
            toolObj = toolClass.newInstance();
        } catch (final Exception ex) {
            throw new GATKException("cannot instantiate tool class, probable bug: " + toolClass.getName(), ex);
        }
        final CommandLineArgumentParser parser = parserFor(toolObj);

        final List<NamedArgumentDefinition> namedArgs = parser.getNamedArgumentDefinitions();
        final List<NamedArgumentDefinition> outputs = namedArgs.stream()
                .filter(arg -> {
                    final Field f = arg.getUnderlyingField();
                    final WorkflowResource resource = f.getAnnotation(WorkflowResource.class);
                    return resource != null && resource.output();
                }).collect(Collectors.toList());
        if (outputs.isEmpty()) {
            outputs;
        }
        try (final PrintWriter ow = output == null ? new PrintWriter(System.out) : new PrintWriter(new FileWriter(output))) {
            ow.println("task " + toolName + "{");
            ow.println("   output {");
            if (outputs.isEmpty()) {
                ow.println("     //TODO please add outputs here");
            } else {
                for (final NamedArgumentDefinition arg : outputs) {
                    ow.println("    File" + (arg.isOptional() ? "?" : "") + " " + underscore(arg.getFullName()));
                }
                ow.println("     //TODO add additional outputs here and delete this comment");
            }
            ow.println("   }");
            ow.println("}");
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(output, ex);
        }

    }

    private String underscore(final String name) {
        if (name.isEmpty()) {
            return "_";
        }
        final StringBuilder builder = new StringBuilder(name.length() << 1);
        if (Character.isDigit(name.charAt(0))) {
            builder.append('_');
        }
        for (int i = 0; i < name.length(); i++) {
            final char ch = name.charAt(i);
            if (((ch & 127) == ch) && (ch == '_' || Character.isLetterOrDigit(ch))) {
                builder.append(ch);
            } else {
                builder.append('_');
            }
        }
        return builder.toString();
    }

    private CommandLineArgumentParser parserFor(final Object toolObj) {
        final CommandLineParser parser;
        if (toolObj instanceof CommandLineProgram) {
            parser = ((CommandLineProgram)toolObj).getCommandLineParser();
        } else if (toolObj instanceof picard.cmdline.CommandLineProgram) {
            parser =  ((picard.cmdline.CommandLineProgram)toolObj).getCommandLineParser();
        } else {
            throw new GATKException("unexpected condition, probable bug");
        }

        if (parser instanceof CommandLineArgumentParser) {
            return ((CommandLineArgumentParser)parser);
        } else {
            throw new UserException("Currently we cannot handle Picard tools that use the legacy argument parser");
        }
    }
}
