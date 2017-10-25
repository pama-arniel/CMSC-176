import Jama.Matrix;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by Jose Arniel Pama on 9/19/2017.
 * Link to Java Matrices Library, JAMA: http://math.nist.gov/javanumerics/jama/
 */

public class Solution {
    private Matrix X; //the whole input data
    private ArrayList<String> Y; //the whole output data
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

    public void scalefeatures(Matrix X){
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
    }

    public Matrix getX(){
        return X;
    }

}