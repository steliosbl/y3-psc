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

#define vec3_iterate(idx, start, stop)             \
  for (idx.x = start.x; idx.x < stop.x; idx.x++)   \
    for (idx.y = start.y; idx.y < stop.y; idx.y++) \
      for (idx.z = start.z; idx.z < stop.z; idx.z++)

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

  vec3(double x)
  {
    this->x = x;
    this->y = x;
    this->z = x;
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

inline double magnitude_squared(const vec3 &a)
{
  return a.x * a.x +
         a.y * a.y +
         a.z * a.z;
}

struct Particle
{
  struct vec3 x;
  struct vec3 v;
  struct vec3 f;
  struct Particle *next;

  Particle() {}

  Particle(vec3 x, vec3 v, vec3 f, Particle *next)
  {
    this->x = x;
    this->v = v;
    this->f = f;
    this->next = next;
  }
};

inline void insert_particle(Particle **root, Particle *i)
{
  i->next = *root;
  *root = i;
}

inline void remove_particle(Particle **i)
{
  *i = (*i)->next;
}

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

  Particle *particles;
  Particle **cells;

  vec3 *f; // Force components

  int *neighbour_indices;
  int n_neighbours;

  double *mass;               // Masses of particles
  double timeStepSize;        // Global time step size
  double timeStepSizeInitial; // Global time step size
  double maxV;                // Maximum velocity of all particles
  double minDx;               // Minimum distance between two elements
  vec3 momentum;

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
                      cells(nullptr), particles(nullptr), f(nullptr), mass(nullptr), neighbour_indices(nullptr),
                      timeStepSize(0), timeStepSizeInitial(0), maxV(0), minDx(0), videoFile(),
                      snapshotCounter(0), timeStepCounter(0){};

  /**
   * Destructor.
   */
  ~NBodySimulation()
  {
    if (cells != nullptr)
    {
      delete[] cells;
    }
    if (particles != nullptr)
    {
      delete[] particles;
    }
    if (f != nullptr)
    {
      delete[] f;
    }
    if (mass != nullptr)
    {
      delete[] mass;
    }
    if (neighbour_indices != nullptr)
    {
      delete[] neighbour_indices;
    }
  }

  void init_neighbour_indices()
  {
    // int n_max = (2*ALPHA+1)*(2*ALPHA+1)*(2*ALPHA+1);
    int n_max = (ALPHA + 1) * (ALPHA + 1) * (ALPHA + 1);
    int *results = new int[n_max];
    int n_results = 0;
    vec3 corner;

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

  inline void zero_forces()
  {
    for (int i = 0; i < NumberOfBodies; i++)
    {
      f[i] = {0.0, 0.0, 0.0};
    }
  }

  inline int scalar_cell_index(int x, int y, int z)
  {
    return x + NC * y + NC * NC * z;
  }

  inline int scalar_cell_index(vec3 i)
  {
    return scalar_cell_index(i.x, i.y, i.z);
  }

  inline int scalar_cell_index(Particle *i)
  {
    return scalar_cell_index((i->x + DOMAIN / 2) / CELL_SIDE_LEN);
  }

  void sort_particles()
  {
    // Iterate all cells
    for (int c = 0; c < N_CELLS; c++)
    {

      // Take the head of the linked list of that cell
      Particle **prev = &cells[c];
      Particle *i = *prev;

      // As long as we have particles in this cell
      while (i != nullptr)
      {

        // Get the new cell index for the particle
        int p_c = scalar_cell_index(i);
        // If it is different from the current cell
        if (p_c != c && p_c >= 0 && p_c < N_CELLS)
        {

          // Make the next of the particle be the new head of the cell
          *prev = i->next;

          // Make the next of the particle be the head of the new cell
          i->next = cells[p_c];

          // Make the head of the new cell be the particle
          cells[p_c] = i;
        }
        else
        {
          prev = &i->next;
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
    cells = new Particle *[N_CELLS];
    particles = new Particle[NumberOfBodies];
    f = new vec3[NumberOfBodies];
    mass = new double[NumberOfBodies];

    // Init cells as empty pointers
    for (int i = 0; i < N_CELLS; i++)
    {
      cells[i] = nullptr;
    }

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
    timeStepSizeInitial = timeStepSize = std::stof(argv[readArgument]);
    readArgument++;

    for (int i = 0; i < NumberOfBodies; i++)
    {
      particles[i] = Particle(
          vec3(
              std::stof(argv[readArgument]),
              std::stof(argv[readArgument + 1]),
              std::stof(argv[readArgument + 2])),
          vec3(
              std::stof(argv[readArgument + 3]),
              std::stof(argv[readArgument + 4]),
              std::stof(argv[readArgument + 5])),
          vec3(0.0), &particles[i + 1]);
      readArgument += 6;

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
    particles[NumberOfBodies - 1].next = nullptr;

    // Make the head of the 0th cell be the 0th particle which then links to the rest
    cells[0] = &particles[0];

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
    // Get the head of the cell_i's linked list
    Particle *i = cells[cell_i];

    // While there are particles to go through
    while (i != nullptr)
    {
      Particle *j = (cell_i == cell_j) ? i->next : cells[cell_j];

      // While there are particles to go through
      while (j != nullptr)
      {
        // If the two particles are distinct
        if (j != i)
        {
          // Get the distance between the particles
          vec3 distance_vec = i->x - j->x;
          double d2 = magnitude_squared(distance_vec);
          minDx = std::min(d2, minDx);

          // Check whether the particles are in range
          if (d2 <= RCUT_2)
          {
            // d^6 = d^2^3
            double d6 = 1 / (d2 * d2 * d2);
            // LJ = (d^12-d^6)/d = d^6 * d *(d^6 -1)
            double lj = d6 * d6 - d6;

            vec3 force = lj * distance_vec;

            i->f += force;
            j->f -= force;
          }
        }
        j = j->next;
      }
      i = i->next;
    }
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
    momentum = vec3(0.0);

    for (int i = 0; i < NumberOfBodies; i++)
    {
      Particle *p = &particles[i];
      // Step 1
      // Compute half the next Euler time-step for velocity
      p->v += p->f * (timeStepSize / 2);
      p->x += p->v * timeStepSize;
    }

    sort_particles();

    // Step 3
    // Zero out old forces
    zero_forces();

    // Step 4
    // Calculate new forces from the new positions
    force_update_all();

    for (int i = 0; i < NumberOfBodies; i++)
    {
      Particle *p = &particles[i];
      p->f *= mass[i];
      // Step 1
      // Compute half the next Euler time-step for velocity
      p->v += p->f * (timeStepSize / 2);
      maxV = std::max(magnitude_squared(p->v), maxV);
      momentum += p->v * mass[i];
    }

    maxV = std::sqrt(maxV);
    minDx = std::sqrt(minDx);

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
      out << particles[i].x.x
          << " "
          << particles[i].x.y
          << " "
          << particles[i].x.z
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
              << ",\t momentum=" << magnitude(momentum)
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
              << particles[0].x.x << ", " << particles[0].x.y << ", " << particles[0].x.z << std::endl;
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
