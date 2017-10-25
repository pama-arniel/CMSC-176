import Jama.Matrix;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by Jose Arniel Pama on 10/02/2017.
 * Link to Java Matrices Library, JAMA: http://math.nist.gov/javanumerics/jama/
 */

public class Solution {
    private Matrix X; //the whole input data
    private ArrayList<String> Y; //the whole output data
    private Matrix engineeredX; //the polynomially engineered values of X
    private Matrix multiclassY; //multiclass representation of the response output
    private Matrix theta; //theta values representing one hypothesis
    private int M; //number of training examples; row dimension
    private int N; //number of Xs in one training example; col dimension

    public Solution() {
        M = 0;
    }


    public void load_data(String filename) {
        BufferedReader br = null;

        try {
            String currline;
            br = new BufferedReader(new FileReader(filename));

            ArrayList<ArrayList<Double>> X = new ArrayList<>();
            ArrayList<String> Y = new ArrayList<>();

            currline = br.readLine(); //so that the column labels will be skipped
            while ((currline = br.readLine()) != null) {
                ArrayList currentsample = new ArrayList<>(Arrays.asList(currline.split(",")));
                ArrayList<Double> tobeaddedtoX = new ArrayList<>();

                //parse as doubles, since a JAMA Matrix object is represented as double
                double n;
                int len = currentsample.size();
                for (int i = 0; i < (len - 1); i++) {
                    n = Double.parseDouble((String) currentsample.get(i));
                    tobeaddedtoX.add(n);
                }

                Y.add( (String) currentsample.get(len - 1)); //the last item, index len-1, is the Y of every training example

                tobeaddedtoX.add(0, 1.0); //for the X0, which is 1; necessary for matrix operation later
                X.add(tobeaddedtoX);
            }

            M = X.size();

            int rows = X.size();
            int cols = X.get(0).size();
            //store the contents of double arraylist X to a matrix
            Matrix equivTOX = new Matrix(rows, cols);
            for(int i = 0; i < rows; i++){
                ArrayList<Double> al = X.get(i);
                for(int j = 0; j < cols; j++){
                    equivTOX.set(i, j, al.get(j));
                }
            }

            this.X = equivTOX;
            N = this.X.getColumnDimension();

            this.Y = Y;
            this.theta = new Matrix( N, 1); //row dimension == cols of X matrix

            System.out.println("NUMBER OF TRAINING EXAMPLES: " + this.M);
            System.out.print("\nTHE X MATRIX:");
            this.X.print(N, 2);
            System.out.println("THE OUTPUTS (Y):");

            for(String st: Y){
                System.out.println(st);
            }

            System.out.print("\nINITIAL theta MATRIX: ");
            this.theta.print(this.theta.getColumnDimension(), 2);

            generateMultiClassY();

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (br != null) br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public Matrix scalefeatures(Matrix X){
        int M = X.getRowDimension();
        int N = X.getColumnDimension();
        Matrix resmatrix = new Matrix(M, N);

        for(int col = 0; col < N; col++){

            Matrix Xj = X.getMatrix(0, M - 1, col, col);

            //calculate for the mean
            double currval;
            double total = 0;
            for(int i = 0; i < M; i++){
                currval = Xj.get(i, 0);
                total += currval;
            }

            double currmean = total / M;

            //calculate for the standard deviation
            double xmeansquared = 0; //summation of the (x - mean) squared
            double num;
            for(int i = 0; i < M; i++){
                currval = Xj.get(i, 0);
                num = currval - currmean;
                xmeansquared += (num * num);
            }
            double currstdev = Math.sqrt(xmeansquared / M);

            //scale values according to Z score
            double currzcore;
            for(int i = 0; i < M; i++){
                currval = Xj.get(i, 0);
                currzcore = (currval - currmean) / currstdev;

                if(col==0){ //meaning, it is at the bias column
                    resmatrix.set(i, col, 1.0);
                }else{
                    resmatrix.set(i, col, currzcore);
                }
            }
        }

        System.out.print("\nTHE SCALED X MATRIX:");
        resmatrix.print(N, 2);

        return resmatrix;
    }

    //returns the polynomially engineered values of X
    public Matrix engineerPolynomials(Matrix Xscaled, int degree){
        int someconst = degree + 1; //to be used in the calculations, e.g. calculating the row dimension of the resulting matrix
        int totalcols = Xscaled.getColumnDimension() - 1; //minus 1 for the bias
        int totalrows = (int) Math.pow(someconst, totalcols);
        Matrix exponentmatrix = new Matrix(totalrows, totalcols);

        int exponent = totalcols;
        int interval = 0;

        //generate first the table for exponents to be used in the final engineered polynomials generation
        for(int col = 0; col < totalcols; col++){
            exponent -= 1; //used to know the intervals of each to-be-inputted values which range from 0 to degree
            interval = (int) Math.pow( someconst , exponent ); //refers to how frequent a value (from 0 to degree) alternately appears in a certain column

            int input = 0;
            int num = 0;
            for(int row=0; row < totalrows; row++){
                if(num == interval){
                    if(input < degree){
                        input++;
                    }else{
                        input = 0;
                    }
                    num = 1;
                }else{
                    num++;
                }

                exponentmatrix.set(row, col, input);
            }
        }


        //generate the final polynomially engineered X; use the exponentmatrix for computations
        int scaledXRows = Xscaled.getRowDimension(); //the row dimension of the scaled X
        int expRows = exponentmatrix.getRowDimension(); //the row dimension of the exponents matrix

        Matrix resmatrix = new Matrix(scaledXRows, expRows);
        for(int scaledXrow = 0; scaledXrow < scaledXRows; scaledXrow++){
            for(int expRow = 0; expRow < expRows; expRow++){
                double res = 1.0;
                for(int col = 0; col < totalcols; col++){
                    res *= Math.pow( Xscaled.get(scaledXrow, col+1), exponentmatrix.get(expRow, col));
                }
                resmatrix.set(scaledXrow, expRow, res);
            }
        }

        this.engineeredX = resmatrix;
        System.out.print("\nPolynomially Engineered Values of X: ");
        this.engineeredX.print(this.engineeredX.getColumnDimension(), 2);
        return resmatrix;
    }

    public Matrix getX(){
        return X;
    }

    public void generateMultiClassY(){
        //find first all the possible Y values in Y (e.g. "Iris-setosa", "Iris-versicolor, "Iris-virginica") and store them on a container
        //this should not be hard coded so that this would be usable for other data as well
        ArrayList<String> possibleYVals = new ArrayList<>();
        int ylen = Y.size();
        for(int i = 0; i < ylen; i++){
            String curr = Y.get(i);
            if(!possibleYVals.contains(curr)){
                possibleYVals.add(curr);
            }
        }

        //start generating the multiclass matrix already
        int possibleYValslen = possibleYVals.size(); //number of possible Y values would dictate the number of Y columns
        multiclassY = new Matrix(ylen, possibleYValslen);
        for(int col = 0; col < possibleYValslen; col++){
            String currpossibleYVal = possibleYVals.get(col);

            for(int row = 0; row < ylen; row++){
                String currval = Y.get(row);

                if( currpossibleYVal.equals(currval) ){
                    multiclassY.set(row, col, 1);
                }else{
                    multiclassY.set(row, col, 0);
                }
            }
        }

        System.out.print("\nMultiple Class Y:");
        multiclassY.print(possibleYValslen, 2);
    }

    public void generateTheta(int setval){
        int thetarows = engineeredX.getColumnDimension();
        int thetacols = multiclassY.getColumnDimension();

        this.theta = new Matrix(thetarows, thetacols);

        for(int i=0; i < thetarows; i++){
            for(int j=0; j < thetacols; j++){
                this.theta.set(i, j, setval); //todo try to randomize the values later
            }
        }

        System.out.print("\nMODIFIED theta MATRIX: ");
        this.theta.print(this.theta.getColumnDimension(), 2);
    }

    public void regularizedCost(Matrix engineeredX, Matrix multiclassY, Matrix theta, double lambda){

        int M = engineeredX.getRowDimension(); //M here is for the engineeredX! diff from the orig M above!!
        int thetarows = theta.getRowDimension();
        int thetacols = theta.getColumnDimension();

        ArrayList<Double> constants = new ArrayList<>(); //note: the resulting size of the constants would be equal to the number of theta columns

        //calculate first for the "appended constant"; equal to lambda/2M * summation of the squares of theta (see slides for more details)
        double currtheta = 0.0;
        System.out.println("THE CONSTANT VALUES: ");
        for(int thetacol = 0; thetacol < thetacols; thetacol++){
            double constant = 0.0;
            for(int thetarow = 0; thetarow < thetarows; thetarow++){
                currtheta = theta.get(thetarow,thetacol);
                constant += (currtheta*currtheta);
            }

            constant = constant * (lambda / (2 * M));
            System.out.println("constant: " + constant);
            constants.add(constant);
        }

        System.out.println("Constants array size: " + constants.size() + "\n");

        //calculate now for the regularized cost for each Y column; use the calculated "constant" in the computations
        int ycols = multiclassY.getColumnDimension(); //the number of Y cols in the multiclassY matrix
        int xrows = engineeredX.getRowDimension(); //the number of X rows in the engineeredX matrix
        int xcols = engineeredX.getColumnDimension(); //the number of X cols in the engineeredX matrix

        System.out.println("THE REGULARIZED COSTS: ");
        for(int ycol = 0; ycol < ycols; ycol++){
            double summation = 0.0;

            for(int xrow=0; xrow < xrows; xrow++){
                double yval = multiclassY.get(xrow, ycol);
                double hOfX = calcForH(theta, engineeredX.getMatrix(xrow, xrow, 0, xcols-1));

                double gOfX = calcForG(hOfX);

                double a = yval * Math.log(gOfX);

                double b;
                double b1;

                b1 = Math.log(1 - gOfX);
                if (Double.isInfinite(b1) || Double.isNaN(b1)){
                    b1 = 0.0;
                }

                b = (1 - yval) * b1;
                summation += ( a + b );
            }

            BigDecimal div = new BigDecimal((-summation/M) );
            div = div.setScale(2, BigDecimal.ROUND_HALF_UP);
            double res = div.doubleValue() + constants.get(ycol);
            System.out.print("COST " + (ycol + 1) + " : " + res + "\n");
        }
        System.out.println();
    }

    public double calcForH(Matrix theta, Matrix givenX){
        Matrix resmatrix = givenX.times(theta);
        return resmatrix.get(0,0);
    }

    public double calcForG(double hOfXvalue){
        return (1 / (1 + Math.pow(Math.E, -hOfXvalue)));
    }

    public void regularizedCost(double lambda){
        regularizedCost(engineeredX, multiclassY, theta, lambda);
    }
}