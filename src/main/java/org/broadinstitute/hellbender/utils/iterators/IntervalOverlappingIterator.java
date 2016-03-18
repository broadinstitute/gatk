/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Daniel Gómez-Sánchez
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package org.broadinstitute.hellbender.utils.iterators;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.utils.HasGenomeLocation;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Wraps a {@link htsjdk.samtools.util.Locatable} or {@link org.broadinstitute.hellbender.utils.HasGenomeLocation} with
 * a list of sorted intervals to return only the objects which overlaps with them
 *
 * @author Daniel Gómez-Sánchez (magicDGS)
 */
public class IntervalOverlappingIterator<T extends Locatable & HasGenomeLocation> implements Iterable<T>, Iterator<T> {

	// underlying iterator
	private final Iterator<T> iterator;

	// sorted intervals
	private final Iterator<SimpleInterval> intervals;

	// the current interval to check
	private SimpleInterval currentInterval;

	// the next object to return
	private T next;

	/**
	 * Wraps an iterator to be filter by a sorted list of intervals
	 *
	 * @param iterator		underlying iterator
	 * @param intervals sorted list of intervals to traverse
	 */
	public IntervalOverlappingIterator(Iterator<T> iterator, List<SimpleInterval> intervals) {
		this.iterator = iterator;
		this.intervals = intervals.iterator();
		currentInterval = this.intervals.next();
		advance();
	}

	@Override
	public Iterator<T> iterator() {
		return this;
	}

	@Override
	public boolean hasNext() {
		return next != null;
	}

	@Override
	public T next() {
		if(!hasNext()) {
			throw new NoSuchElementException();
		}
		T toReturn = next;
		advance();
		return toReturn;
	}

	/**
	 * Advance to the next AlignmentContext
	 */
	private void advance() {
		// all sources are finished
		if(!iterator.hasNext() || currentInterval == null) {
			next = null;
		} else {
			next = iterator.next();
			Locatable loc = (next instanceof Locatable) ? next : next.getLocation();
			// if the next AlignmentContext is not in the current interval
			if(!currentInterval.overlaps(loc)) {
				// advance the interval and try with the next one
				currentInterval = (intervals.hasNext()) ? intervals.next() : null;
				advance();
			}
		}
	}
}
