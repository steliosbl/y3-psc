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

class NBodySimulation
{

    double t;
    double tFinal;
    double tPlot;
    double tPlotDelta;

    int NumberOfBodies; // Total number of particles (*at a given time*)
    double C;           // Collision detection constant (squared)

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

    double *velocities;
    double *distances;
    int *collisions;         // Collisions to track
    double *mass;            // Masses of particles
    double timeStepSize;     // Global time step size
    double timeStepSizeHalf; // Global time step size halved
    double maxV;             // Maximum velocity of all particles
    double minDx;            // Minimum distance between two elements
    double max_mass;

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
                        NumberOfBodies(0), max_mass(1.0),
                        xx(nullptr), xy(nullptr), xz(nullptr),
                        vx(nullptr), vy(nullptr), vz(nullptr),
                        fx(nullptr), fy(nullptr), fz(nullptr),
                        velocities(nullptr), distances(nullptr),
                        mass(nullptr), collisions(nullptr),
                        timeStepSize(0), timeStepSizeHalf(0),
                        maxV(0), minDx(0), videoFile(),
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
        if (distances != nullptr)
            delete[] distances;
        if (mass != nullptr)
            delete[] mass;
        if (collisions != nullptr)
            delete[] collisions;
    }

    inline void zero_forces()
    {
        std::fill(fx, fx + NumberOfBodies, 0);
        std::fill(fy, fy + NumberOfBodies, 0);
        std::fill(fz, fz + NumberOfBodies, 0);
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
        C = 1.0 / (NumberOfBodies * 100);
        max_mass = 0.0;

        xx = new double[NumberOfBodies];
        xy = new double[NumberOfBodies];
        xz = new double[NumberOfBodies];
        vx = new double[NumberOfBodies];
        vy = new double[NumberOfBodies];
        vz = new double[NumberOfBodies];
        fx = new double[NumberOfBodies];
        fy = new double[NumberOfBodies];
        fz = new double[NumberOfBodies];
        velocities = new double[NumberOfBodies];
        distances = new double[NumberOfBodies];
        mass = new double[NumberOfBodies];
        collisions = new int[NumberOfBodies];
        zero_forces();

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

            mass[i] = std::stof(argv[readArgument]);
            max_mass = std::max(max_mass, mass[i]);
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

    inline void join_particles(int i, int n_collisions)
    {
        // Accumulate total mass
        double mass_new = mass[i];

        // Accumulate the numerator of the mass-weighed mean of position
        double x_new_x = xx[i] * mass_new;
        double x_new_y = xy[i] * mass_new;
        double x_new_z = xz[i] * mass_new;

        // Same for velocity
        double v_new_x = vx[i] * mass_new;
        double v_new_y = vy[i] * mass_new;
        double v_new_z = vz[i] * mass_new;

// Iterate the bodies to join with this one in reverse order
// We know that they will be in increasing order of index
#pragma omp simd reduction(+ \
                           : x_new_x, x_new_y, x_new_z, v_new_x, v_new_y, v_new_z, mass_new)
        for (int c = 0; c < n_collisions; c++)
        {
            int j = collisions[n_collisions - c - 1];

            // Accumulate numerator
            x_new_x += xx[j] * mass[j];
            x_new_y += xy[j] * mass[j];
            x_new_z += xz[j] * mass[j];
            v_new_x += vx[j] * mass[j];
            v_new_y += vy[j] * mass[j];
            v_new_z += vz[j] * mass[j];

            // Accumulate combined mass (the denominator of the mass-weighted mean)
            mass_new += mass[j];

            // Replace body j with the last element in global arrays
            // That way we effectively shrink their size by 1
            int old_idx = NumberOfBodies - c - 1;
            if (old_idx != i)
            {
                xx[j] = xx[old_idx];
                xy[j] = xy[old_idx];
                xz[j] = xz[old_idx];
                vx[j] = vx[old_idx];
                vy[j] = vy[old_idx];
                vz[j] = vz[old_idx];
                fx[j] = fx[old_idx];
                fy[j] = fy[old_idx];
                fz[j] = fz[old_idx];
                mass[j] = mass[old_idx];
            }
        }

        // Now divide by the denominator and assign
        mass[i] = mass_new;
        max_mass = std::max(max_mass, mass_new);

        mass_new = 1.0 / mass_new;
        xx[i] = x_new_x * mass_new;
        xy[i] = x_new_y * mass_new;
        xz[i] = x_new_z * mass_new;
        vx[i] = v_new_x * mass_new;
        vy[i] = v_new_y * mass_new;
        vz[i] = v_new_z * mass_new;
    }

    inline void collision_detection(int i)
    {
        // We will detect all collisions with this particle
        int n_collisions = 0; // And count them here
        for (int j = i + 1; j < NumberOfBodies; j++)
        {
            if (distance_squared(i, j) <= C * (mass[i] + mass[j]))
            {
                collisions[n_collisions] = j;
                n_collisions += 1;
            }
        }

        // Run a second loop to join all particles with this one
        join_particles(i, n_collisions);
        NumberOfBodies -= n_collisions;

        // Re-calc the forces for the resulting particle
        fx[i] = 0.0;
        fy[i] = 0.0;
        fz[i] = 0.0;
        // force_update_single(i);
    }

    inline void force_update()
    {
        double f_new_x, f_new_y, f_new_z;
        double f_x, f_y, f_z;
        double dx, dy, dz;

        for (int i = 0; i < NumberOfBodies; i++)
        {
            f_new_x = fx[i];
            f_new_y = fy[i];
            f_new_z = fz[i];
#pragma omp simd reduction(+ \
                           : f_new_x, f_new_y, f_new_z)
            for (int j = i + 1; j < NumberOfBodies; j++)
            {
                // Compute distance and track max
                dx = xx[j] - xx[i];
                dy = xy[j] - xy[i];
                dz = xz[j] - xz[i];

                double distance2 = magnitude_squared(dx, dy, dz);
                double distance = std::sqrt(distance2);
                double denom = 1 / (distance2 * distance);
                distances[j] = distance;

                // Compute new acceleration.
                // Normally we would divide by m to get acceleration
                // Instead we skip the mass component of the force and multiply by the other one
                // That way we avoid multiplying and then dividing in the next step
                f_x = dx * denom;
                f_y = dy * denom;
                f_z = dz * denom;

                f_new_x += f_x * mass[j];
                f_new_y += f_y * mass[j];
                f_new_z += f_z * mass[j];

                fx[j] -= f_x * mass[i];
                fy[j] -= f_y * mass[i];
                fz[j] -= f_z * mass[i];
            }

            // Assign final force (actually acceleration) value
            fx[i] = f_new_x;
            fy[i] = f_new_y;
            fz[i] = f_new_z;

            // Find min Dx
            double min_dx = std::numeric_limits<double>::max();
#pragma omp simd reduction(min \
                           : min_dx)
            for (int j = i + 1; j < NumberOfBodies; j++)
            {
                min_dx = min_dx < distances[j] ? min_dx : distances[j];
            }

            minDx = minDx < min_dx ? minDx : min_dx;
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

#pragma omp simd
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

        // Step 3
        // Zero out old forces
        zero_forces();

        // Step 4
        // Calculate new forces from the new positions
        // Iterate upper triangle of all particles
        force_update();

// Step 5
// Update the velocities by full time step
#pragma omp simd
        for (int i = 0; i < NumberOfBodies; i++)
        {
            vx[i] += fx[i] * timeStepSizeHalf;
            vy[i] += fy[i] * timeStepSizeHalf;
            vz[i] += fz[i] * timeStepSizeHalf;
            velocities[i] = magnitude_squared(vx[i], vy[i], vz[i]);
        }

        // Step 6
        // Collisions
        // Largest possible collision radius is 2C * the maximum mass of any current particle
        // So if the smallest distance between any two particles is leq this, we should check
        if (minDx <= C * max_mass * 2)
        {
            // Iterate upper triangle of all particles
            for (int i = 0; i < NumberOfBodies; i++)
            {
                collision_detection(i);
            }
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
        maxV = std::sqrt(*std::max_element(velocities, velocities + NumberOfBodies));
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
