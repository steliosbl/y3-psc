#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

// You can compile this file with
// g++ -O3 assignment-code.cpp -o assignment-code
// or with the Makefile  and run it with
// ./assignment-code

// Results will be added to the paraview-output directory. In it you will find
// a result.pvd file that you can open with ParaView. To see the points you will
// need to look a the properties of result.pvd and select the representation
// "Point Gaussian". Pressing play will play your time steps.

class NBodySimulation {

  double t;
  double tFinal;
  double tPlot;
  double tPlotDelta;

  int NumberOfBodies;

  /**
   * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
   * each pointer represents one molecule/particle/body.
   */
  double** x;

  /**
   * Equivalent to x storing the velocities.
   */
  double** v;

  /**
   * One mass entry per molecule/particle.
   */
  double*  mass;

  /**
   * Global time step size used.
   */
  double timeStepSize;

  /**
   * Maximum velocity of all particles.
   */
  double maxV;

  /**
   * Minimum distance between two elements.
   */
  double minDx;

  /**
   * Stream for video output file.
   */
  std::ofstream videoFile;

  /**
   * Output counters.
   */
  int snapshotCounter;
  int timeStepCounter;


public:
  /**
   * Constructor.
   */
  NBodySimulation () :
    t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
    x(nullptr), v(nullptr), mass(nullptr),
    timeStepSize(0), maxV(0), minDx(0), videoFile(nullptr),
    snapshotCounter(0), timeStepCounter(0) {};

  /**
   * Destructor.
   */
  ~NBodySimulation () {
    if (x != nullptr) {
      for (int i=0; i<NumberOfBodies; i++)
        delete [] x[i];
      delete [] x;
    }
    if (v != nullptr) {
      for (int i=0; i<NumberOfBodies; i++)
        delete [] v[i];
      delete [] v;
    }
    if (mass != nullptr) {
      delete [] mass;
    }
  }

  /**
   * Set up scenario from the command line.
   *
   * If you need additional helper data structures, you can initialise them
   * here. Alternatively, you can introduce a totally new function to initialise
   * additional data fields and call this new function from main after setUp().
   * Either way is fine.
   *
   * This operation's semantics is not to be changed in the assignment.
   */
  void setUp(int argc, char** argv) {
    NumberOfBodies = (argc-4) / 7;

    x    = new double*[NumberOfBodies];
    v    = new double*[NumberOfBodies];
    mass = new double [NumberOfBodies];

    int readArgument = 1;

    tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
    tFinal       = std::stof(argv[readArgument]); readArgument++;
    timeStepSize = std::stof(argv[readArgument]); readArgument++;

    for (int i=0; i<NumberOfBodies; i++) {
      x[i] = new double[3];
      v[i] = new double[3];

      x[i][0] = std::stof(argv[readArgument]); readArgument++;
      x[i][1] = std::stof(argv[readArgument]); readArgument++;
      x[i][2] = std::stof(argv[readArgument]); readArgument++;

      v[i][0] = std::stof(argv[readArgument]); readArgument++;
      v[i][1] = std::stof(argv[readArgument]); readArgument++;
      v[i][2] = std::stof(argv[readArgument]); readArgument++;

      mass[i] = std::stof(argv[readArgument]); readArgument++;

      if (mass[i]<=0.0 ) {
        std::cerr << "invalid mass for body " << i << std::endl;
        exit(-2);
      }
    }

    std::cout << "created setup with " << NumberOfBodies << " bodies"
              << std::endl;

    if (tPlotDelta<=0.0) {
      std::cout << "plotting switched off" << std::endl;
      tPlot = tFinal + 1.0;
    }
    else {
      std::cout << "plot initial setup plus every " << tPlotDelta
                << " time units" << std::endl;
      tPlot = 0.0;
    }
  }

  /**
   * This is where the force is calculated
   * Currently: gravity, ie (x_i-x_j) * m_i * m_j/r^3
   **/
  double force_calculation(int i, int j, int direction){
    // Euclidean distance
    const double distance = sqrt(
                                   (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
                                   (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
                                   (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
                                   );
    const double distance3 = distance * distance * distance;
    minDx = std::min( minDx,distance );

    return (x[i][direction]-x[j][direction]) * mass[i]*mass[j] / distance3;
  }

  /**
   * This is where the timestepping scheme and force updates are implemented
   */
  void updateBody() {

    timeStepCounter++;
    maxV   = 0.0;
    minDx  = std::numeric_limits<double>::max();

    // force0 = force along x direction
    // force1 = force along y direction
    // force2 = force along z direction
    double* force0 = new double[NumberOfBodies];
    double* force1 = new double[NumberOfBodies];
    double* force2 = new double[NumberOfBodies];

    force0[0] = 0.0;
    force1[0] = 0.0;
    force2[0] = 0.0;

    for (int i=1; i<NumberOfBodies; i++) {
      // x,y,z forces acting on particle 0.
      force0[0] += force_calculation(i,0,0);
      force1[0] += force_calculation(i,0,1);
      force2[0] += force_calculation(i,0,2);
    }

    x[0][0] = x[0][0] + timeStepSize * v[0][0];
    x[0][1] = x[0][1] + timeStepSize * v[0][1];
    x[0][2] = x[0][2] + timeStepSize * v[0][2];

    v[0][0] = v[0][0] + timeStepSize * force0[0] / mass[0];
    v[0][1] = v[0][1] + timeStepSize * force1[0] / mass[0];
    v[0][2] = v[0][2] + timeStepSize * force2[0] / mass[0];

    maxV = std::sqrt( v[0][0]*v[0][0] + v[0][1]*v[0][1] + v[0][2]*v[0][2] );

    t += timeStepSize;

    delete[] force0;
    delete[] force1;
    delete[] force2;
  }

  /**
   * Check if simulation has been completed.
   */
  bool hasReachedEnd() {
    return t > tFinal;
  }

  /**
   * This operation is not to be changed in the assignment.
   */
  void openParaviewVideoFile() {
    videoFile.open("paraview-output/result.pvd");
    videoFile << "<?xml version=\"1.0\"?>" << std::endl
              << "<VTKFile type=\"Collection\""
                 " version=\"0.1\""
                 " byte_order=\"LittleEndian\""
                 " compressor=\"vtkZLibDataCompressor\">" << std::endl
              << "<Collection>";
  }

  /**
   * This operation is not to be changed in the assignment.
   */
  void closeParaviewVideoFile() {
    videoFile << "</Collection>"
              << "</VTKFile>" << std::endl;
    videoFile.close();
  }

  /**
   * The file format is documented at
   * http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
   *
   * This operation is not to be changed in the assignment.
   */
  void printParaviewSnapshot() {
    static int counter = -1;
    counter++;
    std::stringstream filename, filename_nofolder;
    filename << "paraview-output/result-" << counter <<  ".vtp";
    filename_nofolder << "result-" << counter <<  ".vtp";
    std::ofstream out( filename.str().c_str() );
    out << "<VTKFile type=\"PolyData\" >" << std::endl
        << "<PolyData>" << std::endl
        << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
        << "  <Points>" << std::endl
        << "   <DataArray type=\"Float64\""
                  " NumberOfComponents=\"3\""
                  " format=\"ascii\">";

    for (int i=0; i<NumberOfBodies; i++) {
      out << x[i][0]
          << " "
          << x[i][1]
          << " "
          << x[i][2]
          << " ";
    }

    out << "   </DataArray>" << std::endl
        << "  </Points>" << std::endl
        << " </Piece>" << std::endl
        << "</PolyData>" << std::endl
        << "</VTKFile>"  << std::endl;

    out.close();

    videoFile << "<DataSet timestep=\"" << counter
              << "\" group=\"\" part=\"0\" file=\"" << filename_nofolder.str()
              << "\"/>" << std::endl;
  }

  /**
   * This operations is not to be changed in the assignment.
   */
  void printSnapshotSummary() {
    std::cout << "plot next snapshot"
              << ",\t time step=" << timeStepCounter
              << ",\t t="         << t
              << ",\t dt="        << timeStepSize
              << ",\t v_max="     << maxV
              << ",\t dx_min="    << minDx
              << std::endl;
  }

  /**
   * This operations is not to be changed in the assignment.
   */
  void takeSnapshot() {
    if (t >= tPlot) {
      printParaviewSnapshot();
      printSnapshotSummary();
      tPlot += tPlotDelta;
    }
  }

  /**
   * This operations is not to be changed in the assignment.
   */
  void printSummary() {
    std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
    std::cout << "Position of first remaining object: "
              << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;
  }
};


/**
 * Main routine.
 *
 * No major changes in assignment. You can add a few initialisation
 * or stuff if you feel the need to do so. But keep in mind that you
 * may not alter what the program plots to the terminal.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " << std::string(argv[0])
              << " plot-time final-time dt objects" << std::endl
              << " Details:" << std::endl
              << " ----------------------------------" << std::endl
              << "  plot-time:        interval after how many time units to plot."
                 " Use 0 to switch off plotting" << std::endl
              << "  final-time:      simulated time (greater 0)" << std::endl
              << "  dt:              time step size (greater 0)" << std::endl
              << "  objects:         any number of bodies, specified by position, velocity, mass" << std::endl
              << std::endl
              << "Examples of arguments:" << std::endl
              << "+ One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "    0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0" << std::endl
              << "+ One body spiralling around the other" << std::endl
              << "    0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0" << std::endl
              << "+ Three-body setup from first lecture" << std::endl
              << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0" << std::endl
              << "+ Five-body setup" << std::endl
              << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0" << std::endl
              << std::endl;

    return -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries"
                 " (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments"
                 " (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  // Code that initialises and runs the simulation.
  NBodySimulation nbs;
  nbs.setUp(argc,argv);
  nbs.openParaviewVideoFile();
  nbs.takeSnapshot();

  while (!nbs.hasReachedEnd()) {
    nbs.updateBody();
    nbs.takeSnapshot();
  }

  nbs.printSummary();
  nbs.closeParaviewVideoFile();

  return 0;
}
