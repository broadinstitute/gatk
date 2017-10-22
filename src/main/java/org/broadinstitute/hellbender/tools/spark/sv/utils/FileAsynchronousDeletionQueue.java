package org.broadinstitute.hellbender.tools.spark.sv.utils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.Deque;
import java.util.concurrent.BlockingDeque;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.LinkedBlockingQueue;

/**
 * Created by valentin on 10/18/17.
 */
public class FileAsynchronousDeletionQueue {

    private static final int TRIGGER_SIZE = 100;

    private static final Thread DELETING_THREAD = new Thread() {

        @Override
        public void run() {
            while (!scheduledForDeletion.isEmpty()) {
                final Path next = deathRow.remove();
                try {
                    Files.deleteIfExists(next);
                } catch (final IOException e) {
                    // nothing to do.
                }
            }
        }
    };

    private static BlockingQueue<Path> deathRow = new LinkedBlockingQueue<>();
    private static Deque<Path> scheduledForDeletion = new ArrayDeque<>(100);

    public static void scheduleForDeletionIfExists(final Path path) {
        synchronized (scheduledForDeletion) {
            scheduledForDeletion.add(path);
            if (scheduledForDeletion.size() >= TRIGGER_SIZE) {
                deathRow.addAll(scheduledForDeletion);
                scheduledForDeletion.clear();
                if (!DELETING_THREAD.isAlive()) {
                    DELETING_THREAD.setDaemon(true);
                    DELETING_THREAD.start();
                }
            }
        }
    }
}
