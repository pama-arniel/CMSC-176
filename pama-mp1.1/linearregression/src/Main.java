import Jama.Matrix;

public class Main {
    private static final String infile = "HousePricingRelationship.in";
    //private static final String infile = "sample";

    //reads training samples from input file and do their corresponding operations
    public static void main(String[] args) {
        Solution s = new Solution();
        s.load_data(infile);
        //s.print2D(s.getX());
        /*s.print2D( s.to2DArray( s.getX() ) );
        s.print2D( s.generateTheta( s.to2DArray( s.getX() ) ) );*/

        //testing for solveforH
        /*Matrix X = Matrix.constructWithCopy( s.alalTo2DArray( s.getX() ) );
        Matrix theta = Matrix.constructWithCopy( s.generateTheta( s.alalTo2DArray( s.getX() ) ) );
        Matrix H = s.solveforH(X, theta);
        H.print( H.getColumnDimension(), 2);

        Matrix Y = Matrix.constructWithCopy( s.YTo2DArray( s.getY() ) );
        Y.print( Y.getColumnDimension(), 2 );

        //testing for solveforHMinusY
        Matrix hy = s.solveforHMinusY( H, Y );
        hy.print( hy.getColumnDimension(), 2 );

        System.out.println( s.summationOfHMinusYSquared(hy) );*/

        //2nd to the last testing
        /*Matrix X = s.getX();
        Matrix Y = s.getY();
        Matrix theta = s.getTheta();
        System.out.println( X.getColumnDimension() + " == " + theta.getRowDimension() );
        X.print( X.getColumnDimension(), 2);
        theta.print( theta.getColumnDimension(), 2);
        Matrix H = s.solveforH( X, theta );
        H.print( H.getColumnDimension(), 2);
        Y.print( Y.getColumnDimension(), 2);
        Matrix hy = s.solveforHMinusY( H, Y );
        hy.print( hy.getColumnDimension(), 2);
        double sum = s.summationOfHMinusYSquared(hy);
        System.out.println( "Sum: " + sum );*/

        //MAIN TEST
        System.out.println( "COST: " +  s.cost() );
    }
}
