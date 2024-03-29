package org.broadinstitute.hellbender.utils.logging;

import ca.rmen.lfrc.FrenchRevolutionaryCalendar;
import ca.rmen.lfrc.FrenchRevolutionaryCalendarDate;
import org.apache.logging.log4j.core.LogEvent;
import org.apache.logging.log4j.core.config.plugins.Plugin;
import org.apache.logging.log4j.core.pattern.ConverterKeys;
import org.apache.logging.log4j.core.pattern.LogEventPatternConverter;
import org.apache.logging.log4j.core.pattern.PatternConverter;

import java.time.Instant;


import java.time.ZoneId;
import java.time.ZonedDateTime;
import java.util.GregorianCalendar;
import java.util.Locale;

@Plugin(name = "frenchRevDate", category = PatternConverter.CATEGORY)
@ConverterKeys({"frenchRevDate"})
public class FrenchLogEventPatternConverter extends LogEventPatternConverter {
    private String outputFormat = "%E, %dd-%MMMM-%y, %H:%mm:%ss, %T:%DDDD";

    /**
     * Constructs an instance of FrenchLogEventPatternConverter.
     */
    protected FrenchLogEventPatternConverter(String name, String style) {
        super(name, style);
    }

    public static FrenchLogEventPatternConverter newInstance(final String[] options) {
        return new FrenchLogEventPatternConverter("frenchRevDate", "frenchRevDate");
    }

    @Override
    public void format(LogEvent event, StringBuilder toAppendTo) {
        FrenchRevolutionaryCalendar frc =
                new FrenchRevolutionaryCalendar(
                        Locale.FRANCE,
                        FrenchRevolutionaryCalendar.CalculationMethod.ROMME);

        GregorianCalendar unscientificCalendar = GregorianCalendar.from(ZonedDateTime.ofInstant(Instant.ofEpochMilli(event.getTimeMillis()), ZoneId.of("UTC")));

        FrenchRevolutionaryCalendarDate scientificDate = frc.getDate(unscientificCalendar);
        toAppendTo.append(format(scientificDate, outputFormat));
    }

    private static String format(FrenchRevolutionaryCalendarDate frenchDate, String outputFormat) {
        String result = outputFormat;
        result = result.replaceAll("%y", String.format("%d", frenchDate.year));
        result = result.replaceAll("%MMMM", frenchDate.getMonthName());
        result = result.replaceAll("%MM", String.format("%02d", frenchDate.month));
        result = result.replaceAll("%M", String.format("%d", frenchDate.month));

        result = result.replaceAll("%dd", String.format("%02d", frenchDate.dayOfMonth));
        result = result.replaceAll("%d", String.format("%d", frenchDate.dayOfMonth));
        result = result.replaceAll("%H", String.format("%d", frenchDate.hour));
        result = result.replaceAll("%mm", String.format("%02d", frenchDate.minute));
        result = result.replaceAll("%ss", String.format("%02d", frenchDate.second));

        result = result.replaceAll("%E", frenchDate.getWeekdayName());
        result = result.replaceAll("%W", String.format("%d", frenchDate.getWeekInMonth()));
        result = result.replaceAll("%T", frenchDate.getObjectTypeName());
        result = result.replaceAll("%DDDD", frenchDate.getObjectOfTheDay());
        return result;
    }
}
