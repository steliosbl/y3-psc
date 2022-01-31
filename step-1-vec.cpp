#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <cstring>

// You can compile this file with
// g++ -O3 assignment-code.cpp -o assignment-code
// or with the Makefile  and run it with
// ./assignment-code

// Results will be added to the paraview-output directory. In it you will find
// a result.pvd file that you can open with ParaView. To see the points you will
// need to look a the properties of result.pvd and select the representation
// "Point Gaussian". Pressing play will play your time steps.

struct vec3
{
  double x, y, z;
  vec3()
  {
  }

  vec3(double x, double y, double z)
  {
    this->x = x;
    this->y = y;
    this->z = z;
  }

  inline vec3 &operator+=(const vec3 &b)
  {
    this->x += b.x;
    this->y += b.y;
    this->z += b.z;
    return *this;
  }

  inline vec3 &operator-=(const vec3 &b)
  {
    this->x -= b.x;
    this->y -= b.y;
    this->z -= b.z;
    return *this;
  }

  inline vec3 &operator*=(const vec3 &b)
  {
    this->x *= b.x;
    this->y *= b.y;
    this->z *= b.z;
    return *this;
  }

  inline vec3 &operator/=(const vec3 &b)
  {
    this->x /= b.x;
    this->y /= b.y;
    this->z /= b.z;
    return *this;
  }

  inline vec3 &operator*=(const double &b)
  {
    this->x *= b;
    this->y *= b;
    this->z *= b;
    return *this;
  }

  inline vec3 &operator/=(const double &b)
  {
    this->x /= b;
    this->y /= b;
    this->z /= b;
    return *this;
  }
};

inline vec3 operator+(vec3 a, const vec3 &b)
{
  a += b;
  return a;
}

inline vec3 operator-(vec3 a, const vec3 &b)
{
  a -= b;
  return a;
}

inline vec3 operator*(vec3 a, const vec3 &b)
{
  a *= b;
  return a;
}

inline vec3 operator/(vec3 a, const vec3 &b)
{
  a /= b;
  return a;
}

inline vec3 operator*(vec3 a, const double &b)
{
  a *= b;
  return a;
}

inline vec3 operator/(vec3 a, const double &b)
{
  a /= b;
  return a;
}

inline double magnitude(const vec3 &a)
{
  return sqrt(
      a.x * a.x +
      a.y * a.y +
      a.z * a.z);
}

class NBodySimulation
{

  double t;
  double tFinal;
  double tPlot;
  double tPlotDelta;

  int NumberOfBodies;   // Total number of particles (*at a given time*)
  double C;             // Collision detection constant
  vec3 *x;              // Position components
  vec3 *v;              // Velocity components
  vec3 *f;              // Force components
  int *collisions;      // Collisions to track
  double *mass;         // Masses of particles
  double timeStepSize;  // Global time step size
  double timeStepSizeInitial;  // Global time step size
  double maxV;          // Maximum velocity of all particles
  double minDx;         // Minimum distance between two elements
  double min_mass;

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
  NBodySimulation() : t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
                      x(nullptr), v(nullptr), f(nullptr), mass(nullptr), collisions(nullptr),
                      timeStepSize(0), timeStepSizeInitial(0), maxV(0), minDx(0), videoFile(),
                      snapshotCounter(0), timeStepCounter(0){};

  /**
   * Destructor.
   */
  ~NBodySimulation()
  {
    if (x != nullptr)
    {
      delete[] x;
    }
    if (v != nullptr)
    {
      delete[] v;
    }
    if (f != nullptr)
    {
      delete[] f;
    }
    if (mass != nullptr)
    {
      delete[] mass;
    }
    if (collisions != nullptr)
    {
      delete[] collisions;
    }
  }

  inline void zero_forces()
  {
    for (int i = 0; i < NumberOfBodies; i++)
    {
      f[i] = {0.0, 0.0, 0.0};
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
  void setUp(int argc, char **argv)
  {
    NumberOfBodies = (argc - 4) / 7;
    C = 0.01 / NumberOfBodies;

    x = new vec3[NumberOfBodies];
    v = new vec3[NumberOfBodies];
    f = new vec3[NumberOfBodies];
    mass = new double[NumberOfBodies];
    collisions = new int[NumberOfBodies];
    zero_forces();
    min_mass = std::numeric_limits<double>::max();

    int readArgument = 1;

    tPlotDelta = std::stof(argv[readArgument]);
    readArgument++;
    tFinal = std::stof(argv[readArgument]);
    readArgument++;
    timeStepSizeInitial = timeStepSize = std::stof(argv[readArgument]);
    readArgument++;

    for (int i = 0; i < NumberOfBodies; i++)
    {
      x[i] = vec3(std::stof(argv[readArgument]), std::stof(argv[readArgument + 1]), std::stof(argv[readArgument + 2]));
      readArgument += 3;

      v[i] = vec3(std::stof(argv[readArgument]), std::stof(argv[readArgument + 1]), std::stof(argv[readArgument + 2]));
      readArgument += 3;

      mass[i] = std::stof(argv[readArgument]);
      min_mass = std::min(min_mass, mass[i]);
      readArgument++;

      if (mass[i] <= 0.0)
      {
        std::cerr << "invalid mass for body " << i << std::endl;
        exit(-2);
      }
    }

    std::cout << "created setup with " << NumberOfBodies << " bodies"
              << std::endl;

    if (tPlotDelta <= 0.0)
    {
      std::cout << "plotting switched off" << std::endl;
      tPlot = tFinal + 1.0;
    }
    else
    {
      std::cout << "plot initial setup plus every " << tPlotDelta
                << " time units" << std::endl;
      tPlot = 0.0;
    }
  }

  void join_particles(int i, int n_collisions)
  {
    vec3 x_new = x[i] * mass[i]; // Accumulate the numerator of the mass-weighed mean of position
    vec3 v_new = v[i] * mass[i]; // Same for velocity
    double mass_new = mass[i];   // Accumulate total mass

    // Iterate the bodies to join with this one in reverse order
    // We know that they will be in increasing order of index
    for (int c = 0; c < n_collisions; c++)
    {
      int j = collisions[n_collisions - c - 1];

      // Accumulate numerator
      x_new += x[j] * mass[j];
      v_new += v[j] * mass[j];

      // Accumulate combined mass (the denominator of the mass-weighted mean)
      mass_new += mass[j];

      // Replace body j with the last element in global arrays
      // That way we effectively shrink their size by 1
      x[j] = x[NumberOfBodies - c - 1];
      v[j] = v[NumberOfBodies - c - 1];
      mass[j] = mass[NumberOfBodies - c - 1];
    }

    // Now divide by the denominator and assign
    x[i] = x_new / mass_new;
    v[i] = v_new / mass_new;
    mass[i] = mass_new;
  }

  void collision_detection()
  {
    // Iterate upper triangle of all particles
    for (int i = 0; i < NumberOfBodies; i++)
    {
      // We will detect all collisions with this particle
      int n_collisions = 0; // And count them here
      for (int j = i + 1; j < NumberOfBodies; j++)
      {
        if (magnitude(x[j] - x[i]) <= C * (mass[i] + mass[j]))
        {
          collisions[n_collisions] = j;
          n_collisions += 1;
        }
      }

      // Run a second loop to join all particles with this one
      join_particles(i, n_collisions);
      NumberOfBodies -= n_collisions;
    }
  }

  inline void force_update()
  {
    // Iterate upper triangle of all particles
    for (int i = 0; i < NumberOfBodies; i++)
    {
      for (int j = i + 1; j < NumberOfBodies; j++)
      {
        // Compute distance and track max
        double distance = magnitude(x[j] - x[i]);
        minDx = std::min(minDx, distance);

        // Vector of new force. Normally we would divide by m to get acceleration
        // Instead we skip the mass component of the force and multiply by the other one
        // That way we avoid multiplying and then dividing in the next step
        vec3 temp = (x[j] - x[i]) / (distance * distance * distance);
        f[i] += temp * mass[j];
        f[j] -= temp * mass[i];
      }
    }
  }

  /**
   * This is where the timestepping scheme and force updates are implemented
   */
  void updateBody()
  {
    timeStepCounter++;
    maxV = 0.0;
    minDx = std::numeric_limits<double>::max();

    for (int i = 0; i < NumberOfBodies; i++)
    {
      // Step 1
      // Compute half the next Euler time-step for velocity
      v[i] += f[i] * (timeStepSize / 2);

      // Step 2
      // Update positions
      x[i] += v[i] * timeStepSize;
    }

    // Zero out forces
    zero_forces();

    // Step 3
    // Detect collisions and fuse bodies
    collision_detection();

    // Step 4
    // Calculate new forces from the new positions
    force_update();

    // Step 5
    // Update the velocities by full time step
    for (int i = 0; i < NumberOfBodies; i++)
    {
      v[i] += f[i] * (timeStepSize / 2);
      maxV = std::max(magnitude(v[i]), maxV);
    }

    t += timeStepSize;
  }

  /**
   * Check if simulation has been completed.
   */
  bool hasReachedEnd()
  {
    return t > tFinal;
  }

  /**
   * This operation is not to be changed in the assignment.
   */
  void openParaviewVideoFile()
  {
    videoFile.open("paraview-output/result.pvd");
    videoFile << "<?xml version=\"1.0\"?>" << std::endl
              << "<VTKFile type=\"Collection\""
                 " version=\"0.1\""
                 " byte_order=\"LittleEndian\""
                 " compressor=\"vtkZLibDataCompressor\">"
              << std::endl
              << "<Collection>";
  }

  /**
   * This operation is not to be changed in the assignment.
   */
  void closeParaviewVideoFile()
  {
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
  void printParaviewSnapshot()
  {
    static int counter = -1;
    counter++;
    std::stringstream filename, filename_nofolder;
    filename << "paraview-output/result-" << counter << ".vtp";
    filename_nofolder << "result-" << counter << ".vtp";
    std::ofstream out(filename.str().c_str());
    out << "<VTKFile type=\"PolyData\" >" << std::endl
        << "<PolyData>" << std::endl
        << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
        << "  <Points>" << std::endl
        << "   <DataArray type=\"Float64\""
           " NumberOfComponents=\"3\""
           " format=\"ascii\">";

    for (int i = 0; i < NumberOfBodies; i++)
    {
      out << x[i].x
          << " "
          << x[i].y
          << " "
          << x[i].z
          << " ";
    }

    out << "   </DataArray>" << std::endl
        << "  </Points>" << std::endl
        << " </Piece>" << std::endl
        << "</PolyData>" << std::endl
        << "</VTKFile>" << std::endl;

    out.close();

    videoFile << "<DataSet timestep=\"" << counter
              << "\" group=\"\" part=\"0\" file=\"" << filename_nofolder.str()
              << "\"/>" << std::endl;
  }

  /**
   * This operations is not to be changed in the assignment.
   */
  void printSnapshotSummary()
  {
    std::cout << "plot next snapshot"
              << ",\t time step=" << timeStepCounter
              << ",\t t=" << t
              << ",\t dt=" << timeStepSize
              << ",\t v_max=" << maxV
              << ",\t dx_min=" << minDx
              << std::endl;
  }

  /**
   * This operations is not to be changed in the assignment.
   */
  void takeSnapshot()
  {
    if (t >= tPlot)
    {
      printParaviewSnapshot();
      printSnapshotSummary();
      tPlot += tPlotDelta;
    }
  }

  /**
   * This operations is not to be changed in the assignment.
   */
  void printSummary()
  {
    std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
    std::cout << "Position of first remaining object: "
              << x[0].x << ", " << x[0].y << ", " << x[0].z << std::endl;
  }
};

/**
 * Main routine.
 *
 * No major changes in assignment. You can add a few initialisation
 * or stuff if you feel the need to do so. But keep in mind that you
 * may not alter what the program plots to the terminal.
 */
int main(int argc, char **argv)
{
  if (argc == 1)
  {
    std::cerr << "usage: " << std::string(argv[0])
              << " plot-time final-time dt objects" << std::endl
              << " Details:" << std::endl
              << " ----------------------------------" << std::endl
              << "  plot-time:        interval after how many time units to plot."
                 " Use 0 to switch off plotting"
              << std::endl
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
  else if ((argc - 4) % 7 != 0)
  {
    std::cerr << "error in arguments: each planet is given by seven entries"
                 " (position, velocity, mass)"
              << std::endl;
    std::cerr << "got " << argc << " arguments"
                                   " (three of them are reserved)"
              << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  // Code that initialises and runs the simulation.
  NBodySimulation nbs;
  nbs.setUp(argc, argv);
  nbs.openParaviewVideoFile();
  nbs.takeSnapshot();

  while (!nbs.hasReachedEnd())
  {
    nbs.updateBody();
    nbs.takeSnapshot();
  }

  nbs.printSummary();
  nbs.closeParaviewVideoFile();

  return 0;
}
