import Jama.Matrix;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by Jose Arniel Pama on 8/26/2017.
 * Link to Java Matrices Library, JAMA: http://math.nist.gov/javanumerics/jama/
 */

public class Solution {
    private Matrix X; //the whole input data
    private Matrix Y; //the whole output data
    private Matrix theta; //theta values representing one hypothesis
    private int M; //number of training examples



    public Solution(){
        M = 0;
    }

    public void load_data(String filename){
        BufferedReader br = null;

        try{
            String currline;
            br = new BufferedReader(new FileReader(filename));

            ArrayList<ArrayList<Double>> X = new ArrayList<>();
            ArrayList<Double> Y = new ArrayList<>();

            while( (currline = br.readLine())!= null ){
                ArrayList currentsample = new ArrayList<>(Arrays.asList(currline.split(",")));
                ArrayList<Double> tobeaddedtoX = new ArrayList<>();

                //parse as doubles, since a JAMA Matrix object is represented as double
                double n;
                int len = currentsample.size();
                for(int i = 0; i < (len - 1); i++){
                    n = Double.parseDouble( (String) currentsample.get(i) );
                    tobeaddedtoX.add(n);
                }

                n = Double.parseDouble( (String) currentsample.get(len-1) );
                Y.add(n); //the last item, index len-1, is the Y of every training example

                tobeaddedtoX.add(0, 1.0); //for the X0, which is 1; necessary for matrix operation later
                X.add(tobeaddedtoX);
            }

            M = X.size();
            double[][] x2D = alalTo2DArray(X);
            this.X = Matrix.constructWithCopy( x2D );
            this.theta = Matrix.constructWithCopy( generateTheta( x2D ) );
            this.Y = Matrix.constructWithCopy( YTo2DArray( Y ) );

            System.out.println("===== THE X MATRIX =====");
            this.X.print( this.X.getColumnDimension(), 2);
            System.out.println("===== THE Y MATRIX =====");
            this.Y.print( this.Y.getColumnDimension(), 2);
            System.out.println("===== THE theta MATRIX =====");
            this.theta.print( this.theta.getColumnDimension(), 2);
            System.out.println("NUMBER OF TRAINING EXAMPLES: " + this.M);

        }catch(IOException e){
            e.printStackTrace();
        }finally{
            try{
                if(br != null) br.close();
            }catch(IOException e){
                e.printStackTrace();
            }
        }
    }

    public Matrix getX(){
        return this.X;
    }

    public Matrix getY(){
        return this.Y;
    }

    public Matrix getTheta(){
        return this.theta;
    }

    /////////////// START OF HELPER FUNCTIONS  ///////////////
    public void print2D(ArrayList<ArrayList> alal){
        for(ArrayList al: alal){
            for(Object num: al){
                System.out.print(num.toString()+" ");
            }
            System.out.println();
        }
    }

    public void print2D(double[][] alal){
        for(int i = 0; i < alal.length; i++){
            double[] curr = alal[i];
            for(int j = 0; j < curr.length; j++){
                System.out.print( curr[j] + " ");
            }
            System.out.println();
        }
    }

    //transforms a 2D arraylist to a 2D array
    public double[][] alalTo2DArray(ArrayList<ArrayList<Double>> alal){
        int alallen = alal.size();
        double[][] res = new double[alallen][];
        for(int i = 0; i < alallen; i++){
            ArrayList<Double> row = alal.get(i);

            //perform equivalent 'toArray' operation
            int rowlen = row.size();
            double[] copy = new double[rowlen];

            //manually loop and set individually; simply putting "Integer[] copy = (Integer[]) row.toArray()" does not work
            for(int j = 0; j < rowlen; j++){
                copy[j] = row.get(j);
            }

            res[i] = copy;
        }

        return res;
    }

    //transforms Y to a 2D array
    public double[][] YTo2DArray( ArrayList<Double> Y ){
        int ylen = Y.size();
        double[][] res = new double[ylen][];
        for(int i = 0; i < ylen; i++){
            double[] newrow = new double[1]; //Matrix Y has 1 column always
            newrow[0] = Y.get(i);
            res[i] = newrow;
        }

        return res;
    }

    //this is supposed to generate the theta values randomly, but sets each to 2 at the moment
    public double[][] generateTheta( double[][] X ){
        int xrows = X[0].length;
        double[][] res = new double[xrows][];
        for(int i = 0; i < xrows; i++){
            double[] newrow = new double[1]; //Matrix theta has 1 column always

            //generate random number from 0 to M
            newrow[0] = (int) (0 + (Math.random() * M));
            res[i] = newrow;

            //could be set this way, but would be unnecessary
            /*for(int j = 0; j < newrow.length; j++){
                newrow[j] = 2;
            }*/
        }

        return res;
    }
    /////////////// END OF HELPER FUNCTIONS  ///////////////


    /////////////// START OF COMPUTATION FUNCTIONS  ///////////////
    public Matrix solveforH( Matrix X, Matrix theta ){
        return X.times(theta);
    }

    public Matrix solveforHMinusY( Matrix H, Matrix Y ){
        return H.minus( Y );
    }

    public double summationOfHMinusYSquared( Matrix hminusy ){
        double res = 0;
        double curr;
        int row = hminusy.getRowDimension();
        int col = 1; //it is always 1; could be: "int col = hminusy.getColumnDimension()" but unnecessary
        for(int i = 0; i < row; i++){
            for(int j = 0; j < col; j++){
               curr = hminusy.get(i, j);
               res += curr * curr;
            }
        }
        return res;
    }

    public double cost( Matrix X, Matrix Y, Matrix theta ){
        Matrix H = solveforH(X, theta);
        Matrix hy = solveforHMinusY( H, Y );
        double sumofhysquared = summationOfHMinusYSquared( hy );
        double res = sumofhysquared / (2 * M); //according to the cost function formula: (( h(x) - y )^2 ) / (2m)
        return res;
    }

    public double cost(){
        return cost(this.X, this.Y, this.theta);
    }

    /////////////// END OF COMPUTATION FUNCTIONS  ///////////////

}
