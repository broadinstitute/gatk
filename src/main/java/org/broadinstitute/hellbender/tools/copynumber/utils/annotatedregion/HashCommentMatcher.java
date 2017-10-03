package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import org.apache.commons.collections4.list.UnmodifiableList;
import org.supercsv.comment.CommentMatcher;

import java.util.ArrayList;
import java.util.List;


/**
 * Created by lichtens on 7/14/17.
 */
public class HashCommentMatcher implements CommentMatcher {
    private List<String> comments = new ArrayList<>();

    @Override
    public boolean isComment(String line) {

        final boolean result = line.startsWith("#");
        if (result) {
            comments.add(line);
        }
        return result;
    }

    public HashCommentMatcher() {
    }

    public UnmodifiableList<String> getComments() {
        return new UnmodifiableList<>(comments);
    }
}
