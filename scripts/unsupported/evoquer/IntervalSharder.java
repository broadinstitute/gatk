public class IntervalSharder {
    public static void main(String[] args) {
        if ( args.length != 3 ) {
            System.out.println("Usage: IntervalSharder contig contigLength shardSize");
            System.exit(1);
        }

        final String contig = args[0];
        final int contigLength = Integer.parseInt(args[1]);
        final int shardSize = Integer.parseInt(args[2]);

        int currentStart = 1 - shardSize;
        int currentEnd = 0;

        while ( currentEnd < contigLength ) {
            currentStart += shardSize;
            currentEnd += shardSize;

            if ( currentEnd > contigLength ) {
                currentEnd = contigLength;
            }

            System.out.printf("%s:%d-%d\n", contig, currentStart, currentEnd);
        } 
    }
}
