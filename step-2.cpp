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

// Constants
#define RCUT 2.5
#define DOMAIN 10.0
#define ALPHA 1.0

class NBodySimulation
{

  double t;
  double tFinal;
  double tPlot;
  double tPlotDelta;

  int NumberOfBodies; // Total number of particles (*at a given time*)
  int NC = DOMAIN / (RCUT / ALPHA);
  int N_CELLS = NC * NC * NC;
  double CELL_SIDE_LEN = RCUT / ALPHA;
  double RCUT_2 = RCUT * RCUT;

  // Position components
  double *xx;
  double *xy;
  double *xz;

  // Velocity components
  double *vx;
  double *vy;
  double *vz;

  // Force components
  double *fx;
  double *fy;
  double *fz;

  // Linked-cell algorithm cells & linked-list
  int *cells;
  int *next_particle;

  int *neighbour_indices;
  int n_neighbours;

  double *velocities;
  double *mass;            // Masses of particles
  double timeStepSize;     // Global time step size
  double timeStepSizeHalf; // Global time step size halved
  double maxV;             // Maximum velocity of all particles
  double minDx;            // Minimum distance between two elements

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
  NBodySimulation() : t(0), tFinal(0), tPlot(0), tPlotDelta(0),
                      xx(nullptr), xy(nullptr), xz(nullptr),
                      vx(nullptr), vy(nullptr), vz(nullptr),
                      fx(nullptr), fy(nullptr), fz(nullptr),
                      mass(nullptr), neighbour_indices(nullptr),
                      velocities(nullptr), timeStepSize(0),
                      maxV(0), minDx(0), videoFile(), NumberOfBodies(0),
                      snapshotCounter(0), timeStepCounter(0){};

  /**
   * Destructor.
   */
  ~NBodySimulation()
  {
    if (xx != nullptr)
      delete[] xx;
    if (xy != nullptr)
      delete[] xy;
    if (xz != nullptr)
      delete[] xz;
    if (vx != nullptr)
      delete[] vx;
    if (vy != nullptr)
      delete[] vy;
    if (vz != nullptr)
      delete[] vz;
    if (fx != nullptr)
      delete[] fx;
    if (fy != nullptr)
      delete[] fy;
    if (fz != nullptr)
      delete[] fz;
    if (velocities != nullptr)
      delete[] velocities;
    if (mass != nullptr)
      delete[] mass;
    if (neighbour_indices != nullptr)
      delete[] neighbour_indices;
  }

  inline double distance_squared(int i, int j)
  {
    return (xx[i] - xx[j]) * (xx[i] - xx[j]) +
           (xy[i] - xy[j]) * (xy[i] - xy[j]) +
           (xz[i] - xz[j]) * (xz[i] - xz[j]);
  }

  inline double magnitude_squared(double a, double b, double c)
  {
    return a * a + b * b + c * c;
  }

  inline void zero_forces()
  {
    std::fill(fx, fx + NumberOfBodies, 0);
    std::fill(fy, fy + NumberOfBodies, 0);
    std::fill(fz, fz + NumberOfBodies, 0);
  }

  void init_neighbour_indices()
  {
    // int n_max = (2*ALPHA+1)*(2*ALPHA+1)*(2*ALPHA+1);
    int n_max = (ALPHA + 1) * (ALPHA + 1) * (ALPHA + 1);
    int *results = new int[n_max];
    int n_results = 0;

    for (int x = 1; x <= ALPHA + 1; x++)
    {
      for (int y = 1; y <= ALPHA + 1; y++)
      {
        for (int z = 1; z <= ALPHA + 1; z++)
        {
          int distance = (x - 1) * (x - 1) + (y - 1) * (y - 1) + (z - 1) * (z - 1);
          if (distance < ALPHA * ALPHA)
          {
            results[n_results] = scalar_cell_index(x, y, z);
            n_results++;
          }
        }
      }
    }

    for (int i = 1; i <= ALPHA; i++)
    {
      results[n_results] = scalar_cell_index(0, 0, i);
      results[n_results + 1] = scalar_cell_index(0, i, 0);
      results[n_results + 2] = scalar_cell_index(i, 0, 0);
      n_results += 3;
    }

    results[n_results] = 0;
    n_results += 1;

    std::sort(results, results + n_results);

    neighbour_indices = new int[n_results];
    n_neighbours = n_results;
    std::copy(results, results + n_results, neighbour_indices);

    delete[] results;
  }

  inline int scalar_cell_index(int x, int y, int z)
  {
    return x + NC * y + NC * NC * z;
  }

  inline int scalar_cell_index(int i)
  {
    double x = (xx[i] + DOMAIN / 2) / CELL_SIDE_LEN;
    double y = (xy[i] + DOMAIN / 2) / CELL_SIDE_LEN;
    double z = (xz[i] + DOMAIN / 2) / CELL_SIDE_LEN;
    return scalar_cell_index(x, y, z);
  }

  void sort_particles()
  {
    // Iterate all cells
    for (int c = 0; c < N_CELLS; c++)
    {

      // Take the head of the linked list of that cell
      int *prev = &cells[c];
      int i = cells[c];

      // As long as we have particles in this cell
      while (i != -1)
      {

        // Get the new cell index for the particle
        int p_c = scalar_cell_index(i);
        // If it is different from the current cell
        if (p_c != c && p_c >= 0 && p_c < N_CELLS)
        {

          // Make the next of the particle be the new head of the cell
          *prev = next_particle[i];

          // Make the next of the particle be the head of the new cell
          next_particle[i] = cells[p_c];

          // Make the head of the new cell be the particle
          cells[p_c] = i;
        }
        else
        {
          prev = &next_particle[i];
        }

        // Move to the next particle in the linked list

        i = *prev;
      }
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

    // Initialise data structures
    cells = new int[N_CELLS];
    next_particle = new int[NumberOfBodies];
    xx = new double[NumberOfBodies];
    xy = new double[NumberOfBodies];
    xz = new double[NumberOfBodies];
    vx = new double[NumberOfBodies];
    vy = new double[NumberOfBodies];
    vz = new double[NumberOfBodies];
    fx = new double[NumberOfBodies];
    fy = new double[NumberOfBodies];
    fz = new double[NumberOfBodies];
    mass = new double[NumberOfBodies];
    velocities = new double[NumberOfBodies];

    std::fill(cells, cells + N_CELLS, -1);

    // Zero out the forces
    zero_forces();

    // Get the array of relative (scalar) indices for the
    // upper-triangle neighbours of a cell
    init_neighbour_indices();

    int readArgument = 1;

    tPlotDelta = std::stof(argv[readArgument]);
    readArgument++;
    tFinal = std::stof(argv[readArgument]);
    readArgument++;
    timeStepSize = std::stof(argv[readArgument]);
    timeStepSizeHalf = timeStepSize / 2;
    readArgument++;

    for (int i = 0; i < NumberOfBodies; i++)
    {
      xx[i] = std::stof(argv[readArgument]);
      xy[i] = std::stof(argv[readArgument + 1]);
      xz[i] = std::stof(argv[readArgument + 2]);
      readArgument += 3;

      vx[i] = std::stof(argv[readArgument]);
      vy[i] = std::stof(argv[readArgument + 1]);
      vz[i] = std::stof(argv[readArgument + 2]);
      readArgument += 3;

      next_particle[i] = i + 1;

      mass[i] = std::stof(argv[readArgument]);
      readArgument++;

      if (mass[i] <= 0.0)
      {
        std::cerr << "invalid mass for body " << i << std::endl;
        exit(-2);
      }
    }

    // So far we have assigned each particle's next to be the following particle we read
    // For the last one, make it an null pointer
    next_particle[NumberOfBodies - 1] = -1;

    // Make the head of the 0th cell be the 0th particle which then links to the rest
    cells[0] = 0;

    // Sort the particles into actual cells
    sort_particles();

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

  inline void force_update_cell_pair(int cell_i, int cell_j)
  {
    double min_dx = minDx;
    double f_new_x, f_new_y, f_new_z;
    double f_x, f_y, f_z;
    double dx, dy, dz;

    // Get the head of the cell_i's linked list
    int i = cells[cell_i];

    // While there are particles to go through
    while (i != -1)
    {
      int j = (cell_i == cell_j) ? next_particle[i] : cells[cell_j];
      f_new_x = fx[i];
      f_new_y = fy[i];
      f_new_z = fz[i];

      // While there are particles to go through
      while (j != -1)
      {
        // If the two particles are distinct
        if (j != i)
        {
          // Compute distance and track max
          dx = xx[i] - xx[j];
          dy = xy[i] - xy[j];
          dz = xz[i] - xz[j];

          double d2 = magnitude_squared(dx, dy, dz);
          min_dx = std::min(d2, min_dx);

          // Check whether the particles are in range
          if (d2 <= RCUT_2)
          {
            // d^-4 = 1/(d^2*d^2)
            double d4 = 1 / (d2 * d2);

            // d^-8 = d^-4 * d^-4
            double d8 = d4 * d4;

            // LJ = d^-14 - d^-8 = d^-8 (d^-8*d^2 - 1)
            double lj = d8 * (d8 * d2 - 1);

            f_x = dx * lj;
            f_y = dy * lj;
            f_z = dz * lj;

            f_new_x += f_x;
            f_new_y += f_y;
            f_new_z += f_z;

            fx[j] -= f_x;
            fy[j] -= f_y;
            fz[j] -= f_z;
          }
        }
        j = next_particle[j];
      }

      // Assign final force (actually acceleration) value
      fx[i] = f_new_x;
      fy[i] = f_new_y;
      fz[i] = f_new_z;
      i = next_particle[i];
    }

    minDx = min_dx < minDx ? min_dx : minDx;
  }

  inline void force_update_all()
  {
    // Iterate all cells from bottom left corner to top right
    for (int cell_i = 0; cell_i < N_CELLS; cell_i++)
    {

      // Iterate all the upper-triangle of
      // neighbours (within range) of the current cell
      for (int n = 0; n < n_neighbours; n++)
      {

        // The scalar index of cell-j is the index of cell-a
        // offset by (each) neighbour index
        int cell_j = neighbour_indices[n] + cell_i;

        // Check we have not gone over
        if (cell_j < N_CELLS)
        {
          force_update_cell_pair(cell_i, cell_j);
        }
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
      vx[i] += fx[i] * timeStepSizeHalf;
      vy[i] += fy[i] * timeStepSizeHalf;
      vz[i] += fz[i] * timeStepSizeHalf;

      // Step 2
      // Update positions
      xx[i] += vx[i] * timeStepSize;
      xy[i] += vy[i] * timeStepSize;
      xz[i] += vz[i] * timeStepSize;
    }

    // Step 2.5
    // Update cells
    sort_particles();

    // Step 3
    // Zero out old forces
    zero_forces();

    // Step 4
    // Calculate new forces from the new positions
    force_update_all();

    // Step 5
    // Update the velocities by full time step
    for (int i = 0; i < NumberOfBodies; i++)
    {
      vx[i] += (fx[i] / mass[i]) * timeStepSizeHalf;
      vy[i] += (fy[i] / mass[i]) * timeStepSizeHalf;
      vz[i] += (fz[i] / mass[i]) * timeStepSizeHalf;
      velocities[i] = magnitude_squared(vx[i], vy[i], vz[i]);
    }

    t += timeStepSize;
    if (t >= tPlot)
    {
      maxV = std::sqrt(*std::max_element(velocities, velocities + NumberOfBodies));
      minDx = std::sqrt(minDx);
    }
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
      out << xx[i]
          << " "
          << xy[i]
          << " "
          << xz[i]
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
              << xx[0] << ", " << xy[0] << ", " << xz[0] << std::endl;
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
