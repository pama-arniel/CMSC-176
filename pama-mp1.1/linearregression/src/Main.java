import Jama.Matrix;

public class Main {
    private static final String infile = "HousePricingRelationship.in";
    //private static final String infile = "sample";

    //reads training samples from input file and do their corresponding operations
    public static void main(String[] args) {
        Solution s = new Solution();
        s.load_data(infile);
        System.out.println( "COST: " +  s.cost() );
    }
}
