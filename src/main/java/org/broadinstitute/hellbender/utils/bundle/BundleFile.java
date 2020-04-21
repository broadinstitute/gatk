package org.broadinstitute.hellbender.utils.bundle;

class BundleFile {
    private final String path;
    private final String fileType;

    public BundleFile(final String path, final String fileType) {
        this.path = path;
        this.fileType = fileType;
    }

    public String getPath() {
        return path;
    }

    public String getFileType() {
        return fileType;
    }
}
