package org.broadinstitute.hellbender.utils;

import java.io.Serializable;
import java.util.function.Consumer;

public interface SerializableConsumer <T> extends Consumer<T>, Serializable {}
