import java.util.Date;
import mpi.*; // for mpiJava 

public class Heat2D_mpi {
    private static double a = 1.0; // heat speed
    private static double dt = 1.0; // time quantum
    private static double dd = 2.0; // change in system

    public static void main(String[] args) throws MPIException {

        // verify arguments
        if (args.length != 4) {
            System.out.println("usage: " + "java Heat2D size max_time heat_time interval");
            System.exit(-1);
        }

        MPI.Init(args);

        int size = Integer.parseInt(args[0]);
        int max_time = Integer.parseInt(args[1]);
        int heat_time = Integer.parseInt(args[2]);
        int interval = Integer.parseInt(args[3]);
        double r = a * dt / (dd * dd);

        double[] z = null;

        int currentRank = MPI.COMM_WORLD.Rank();
        int mpiSize = MPI.COMM_WORLD.Size();
        int xSize = size / mpiSize;        
        
        // Initializing printing array for rank 0
        double[] printingArray = null; 
        if (currentRank == 0) {
            printingArray = new double[size * size];
        }

        // later used for boundary exchange
        double[] leftBoundary = null; 
        double[] rightBoundary = null; 

        // boundary container init
        if (mpiSize > 1) {
            if (currentRank == 0) {
                rightBoundary = new double[size]; 
            } else if (currentRank == mpiSize - 1) {
                leftBoundary = new double[size]; 
            } else {
                leftBoundary = new double[size]; 
                rightBoundary = new double[size]; 
            }
        }

        // 1. divide and allocate 2D spaces in 1D for each processes
        if (currentRank != mpiSize - 1) { // if not the last rank
            z = new double[2 * xSize * size];
        } else { // last Rank (this is for when using 3 processes
            // (100 - ((3 - 1) * 33)) = 100 - 66; 34 * 100 * 2
            // (100 - ((4 - 1) * 25)) = 100 - 75; 25 * 100 * 2
            xSize = (size - ((mpiSize - 1) * xSize));
            z = new double[2 * xSize * size];
        }

        // initialize double[] with 0s
        initArray(z);

        // start a timer
        Date startTime = null;
        if (currentRank == 0) {
            startTime = new Date();
        }

        // simulate heat diffusion
        for (int t = 0; t < max_time; t++) {

            int p = t % 2;

            // 2. padding (top, down, left, right)
            //      2.1 padding left and right 
            // if mpiSize == 1 left and right should be applied on both sides
            // other wise, first and last method require paddings of left and right 
            if (currentRank == 0) { 
                leftPadding(z, p, xSize, size);          
            } 
            
            if (currentRank == mpiSize - 1) {
                rightPadding(z, p, xSize, size);            
            } 

            //      2.2 padding top and down 
            //          all should do the top down padding 
            topDownPadding(z, p, xSize, size);

            //      2.3 heating the middle 
            if (t < heat_time) {
                heating(z, p, xSize, size, currentRank);
            }
            
            // 3. exchange boundary data
            // even rank send & odd rank recv first 
            if (mpiSize > 1) { // 2 <= mpiSize 
                boundaryComm(z, rightBoundary, leftBoundary, currentRank, mpiSize, p, xSize, size); 
            }



            // 4. print out the result
            // System.out.println("interval = " + interval); 
            if (interval != 0 && (t % interval == 0 || t == max_time - 1)) { 
                if (currentRank == 0) { // master process

                    // put master calc value in printing array
                    for (int x = 0; x < xSize; x++) {
                        for (int y = 0; y < size; y++) {
                            printingArray[size * x + y] = z[p * xSize * size + x * size + y];
                        }
                    }

                    // receive calculated values from slaves
                    for (int i = 1; i < mpiSize; i++) {
                        int offset = i * xSize * size;
                        if (i != mpiSize - 1) {     // not from the last process
                            MPI.COMM_WORLD.Recv(printingArray, offset, xSize * size, MPI.DOUBLE, i, 0);
                        } else {                    // from the last process
                            MPI.COMM_WORLD.Recv(printingArray, offset, (size - ((mpiSize - 1) * xSize)) * size, MPI.DOUBLE,
                                    i, 0);
                        }
                    }

                    printStatus(printingArray, size, t);
                } else { // other processes
                     // send array
                    MPI.COMM_WORLD.Send(z, p * xSize * size, z.length / 2, MPI.DOUBLE, 0, 0);
                }
            } 
            // 5. Forward Euler method computation 
            heatComputation(z, rightBoundary, leftBoundary, p, xSize, size, r);

        }
        MPI.Finalize();
        
        // finishe the timer
        
        if (currentRank == 0) {
            Date endTime = new Date();
            System.out.println("Elapsed time = " +  (endTime.getTime() - startTime.getTime()));
        }
    }

    public static void heatComputation(double[] z, double[] rightBoundary, double[] leftBoundary, int p, int xSize, int size, double r) {
        int p2 = (p + 1) % 2;
        int previousPlane = p * xSize * size;
        int currentPlane = p2 * xSize * size; 

        double current;
        double east;
        double west;
        double north;
        double south;

        // Compute with left exchanged boundary 
        if (leftBoundary != null) {
            for (int y = 1; y < size - 1; y++) {
                // x = 0
                current = z[previousPlane + y]; 
                east = z[previousPlane + 1 * size + y];
                west = leftBoundary[y];
                north = z[previousPlane + y - 1];
                south = z[previousPlane + y + 1];
                z[currentPlane + y] = current + r * (east - 2 * current + west) + r * (south - 2 * current + north);
            }
        }

        // Compute with right exchanged boundary 
        if (rightBoundary != null) {
            for (int y = 1; y < size - 1; y++) {
                // x = 0
                current = z[previousPlane + (xSize - 1) * size + y]; 
                east = rightBoundary[y];
                west =  z[previousPlane + (xSize - 2) * size + y]; 
                north = z[previousPlane + (xSize - 1) * size + y - 1];
                south = z[previousPlane + (xSize - 1) * size + y + 1];
                z[currentPlane + (xSize - 1) * size + y] = current + r * (east - 2 * current + west) + r * (south - 2 * current + north);
            }
        }

        // heat computation in the middle 
        for (int x = 1; x < xSize - 1; x++) {
            for (int y = 1; y < size - 1; y++) {
                current = z[previousPlane + x * size + y];
                east = z[previousPlane + (x + 1) * size + y];
                west = z[previousPlane + (x - 1) * size + y];
                north = z[previousPlane + x * size + y - 1];
                south = z[previousPlane + x * size + y + 1];
                z[currentPlane + x * size + y] = current + r * (east - 2 * current + west) + r * (south - 2 * current + north);
            }
        }

    }

    public static void boundaryComm(double[] z, double[] rightBoundary, double[] leftBoundary, 
        int currentRank, int mpiSize, int p, int xSize, int size) throws MPIException { 
        
        int leftOffset = p * xSize * size; 
        int rightOffset = p * xSize * size + size * (xSize - 1); 

        if (currentRank == 0) { // left most node

            // send boundary 
            MPI.COMM_WORLD.Send(z, rightOffset, size, MPI.DOUBLE, currentRank + 1, 0); 
            // System.out.println(currentRank + "sent (even)"); 

            // recv boundary 
            MPI.COMM_WORLD.Recv(rightBoundary, 0, size, MPI.DOUBLE, currentRank + 1, 0); 
            // System.out.println(currentRank + "recved (even)"); 

        } else if (currentRank == mpiSize - 1) { // right most node 

            // send boundary 
            MPI.COMM_WORLD.Send(z, leftOffset, size, MPI.DOUBLE, currentRank - 1, 0); 
            // System.out.println(currentRank + "sent (even)");

            // recv boundary 
            MPI.COMM_WORLD.Recv(leftBoundary, 0, size, MPI.DOUBLE, currentRank - 1, 0);
            // System.out.println(currentRank + "recved (Even)");

        } else {
            if (currentRank % 2 == 0) { // even node send first 

                // send boundary 
                MPI.COMM_WORLD.Send(z, leftOffset, size, MPI.DOUBLE, currentRank - 1, 0); 
                // System.out.println(currentRank + "sent left (even) ");

                MPI.COMM_WORLD.Send(z, rightOffset, size, MPI.DOUBLE, currentRank + 1, 0); 
                // System.out.println(currentRank + "sent right (even) ");

                // recv boundary 
                MPI.COMM_WORLD.Recv(leftBoundary, 0, size, MPI.DOUBLE, currentRank - 1, 0);
                // System.out.println(currentRank + "recved left (even)");
                MPI.COMM_WORLD.Recv(rightBoundary, 0, size, MPI.DOUBLE, currentRank + 1, 0);
                // System.out.println(currentRank + "recved right (even)");

            } else { // odd node recv first 

                // recv boundary 
                MPI.COMM_WORLD.Recv(leftBoundary, 0, size, MPI.DOUBLE, currentRank - 1, 0);
                // System.out.println(currentRank + "recved left (odd)");
                MPI.COMM_WORLD.Recv(rightBoundary, 0, size, MPI.DOUBLE, currentRank + 1, 0);
                // System.out.println(currentRank + "recved right (odd)");

                // send boundary 
                MPI.COMM_WORLD.Send(z, leftOffset, size, MPI.DOUBLE, currentRank - 1, 0); 
                // System.out.println(currentRank + "sent left (odd)");
                MPI.COMM_WORLD.Send(z, rightOffset, size, MPI.DOUBLE, currentRank + 1, 0); 
                // System.out.println(currentRank + "sent right (odd)");

            }

        }
    
    }

    public static void heating(double[] array, int p, int xSize, int size, int currentRank) {
        int heatingStartIndex = size / 3;   // inclusive
        int heatingEndIndex = size / 3 * 2; // exclusive

        for (int i = 0; i < xSize; i++) {
            // x index of current stripe index in Whole array  
            int xInTotal = currentRank * xSize + i; 
            if (xInTotal >= heatingStartIndex && xInTotal < heatingEndIndex) {
                // z[p][x][y] = 19.0 
                array[p * xSize * size + i * size] = 19.0; 
            }

            // else ignore 
        }
    }

    public static void leftPadding(double[] array, int p, int xSize, int size) {
       
        for (int y = 0; y < size; y++) {
            array[p * xSize * size + y] = array[p * xSize * size + 1 * size + y]; 
        }        
    }

    public static void rightPadding(double[] array, int p, int xSize, int size) {
        int x = xSize - 1; 
        for (int y = 0; y < size; y++) {
            array[p * xSize * size + x * size + y] = array[p * xSize * size + (x - 1) * size + y]; 
        }
    }

    public static void topDownPadding(double[] array, int p, int xSize, int size) {
        for (int x = 0; x < xSize; x++) {
            array[p * xSize * size + x * size] = array[p * xSize * size + x * size + 1];
            array[p * xSize * size + x * size + (size - 1)] = array[p * xSize * size + x * size + (size - 2)];
        }
    }

    public static void initArray(double[] array) {
        for (int i = 0; i < array.length; i++) {
            array[i] = 0.0;
        }
    }

    public static void printStatus(double[] printingArray, int size, int t) {
        System.out.println("time = " + t);
        for (int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {
                System.out.print((int) (Math.floor(printingArray[x * size + y] / 2)) + " ");
            }
            System.out.println();
        }
        System.out.println();
    }
}
