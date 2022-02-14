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
#include <benchmark/benchmark.h>
#include <x86intrin.h>

struct vec3
{
  // double x, y, z, n;
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

inline double magnitude_squared(const vec3 &a)
{
  return a.x * a.x +
         a.y * a.y +
         a.z * a.z;
}

static void Magnitude_ThreeList(benchmark::State &state)
{
  int NumberOfBodies = 1337 * 2;
  double *xx = new double[NumberOfBodies];
  double *xy = new double[NumberOfBodies];
  double *xz = new double[NumberOfBodies];

  int i = 34;
  int j = 364;
  for (auto _ : state)
  {
    double result = (xx[i] - xx[j]) * (xx[i] - xx[j]) +
                    (xy[i] - xy[j]) * (xy[i] - xy[j]) +
                    (xz[i] - xz[j]) * (xz[i] - xz[j]);
    benchmark::DoNotOptimize(result);
  }
}

static void Magnitude_TripleList_NoLoop(benchmark::State &state)
{
  double *x = new double[6];
  int i = 0;
  int j = 1;
  for (auto _ : state)
  {
    double result = (x[3 * i] - x[3 * j]) * (x[3 * i] - x[3 * j]) +
                    (x[3 * i + 1] - x[3 * j + 1]) * (x[3 * i + 1] - x[3 * j + 1]) +
                    (x[3 * i + 2] - x[3 * j + 2]) * (x[3 * i + 2] - x[3 * j + 2]);
    benchmark::DoNotOptimize(result);
  }
}

static void Magnitude_TripleList_Loop(benchmark::State &state)
{
  double *x = new double[6];
  int i = 0;
  int j = 1;
  for (auto _ : state)
  {
    double result;
#pragma omp simd reduction(+ \
                           : result)
    for (int d = 0; d < 3; d++)
    {
      result += (x[3 * i + d] - x[3 * j + d]) * (x[3 * i + d] - x[3 * j + d]);
    }
    benchmark::DoNotOptimize(result);
  }
}

static void Magnitude_Vec3(benchmark::State &state)
{
  vec3 *x = new vec3[2];
  int i = 0;
  int j = 1;
  for (auto _ : state)
  {
    double result = magnitude_squared(x[i] - x[j]);
    benchmark::DoNotOptimize(result);
  }
}

static void Magnitude_NestedList_NoLoop(benchmark::State &state)
{
  double **x = new double *[2];
  x[0] = new double[3];
  x[1] = new double[3];

  int i = 0;
  int j = 1;

  for (auto _ : state)
  {
    double result = (x[i][0] - x[j][0]) * (x[i][0] - x[j][0]) +
                    (x[i][1] - x[j][1]) * (x[i][1] - x[j][1]) +
                    (x[i][2] - x[j][2]) * (x[i][2] - x[j][2]);
    benchmark::DoNotOptimize(result);
  }
}

static void Magnitude_NestedList_Loop(benchmark::State &state)
{
  double **x = new double *[2];
  x[0] = new double[3];
  x[1] = new double[3];

  int i = 0;
  int j = 1;

  for (auto _ : state)
  {
    double result = 0;

#pragma omp simd reduction(+ \
                           : result)
    for (int d = 0; d < 3; d++)
    {
      result += (x[i][d] - x[j][d]) * (x[i][d] - x[j][d]);
    }

    benchmark::DoNotOptimize(result);
  }
}

// BENCHMARK(Magnitude_TripleList_Loop);
// BENCHMARK(Magnitude_TripleList_NoLoop);
// BENCHMARK(Magnitude_Vec3);
// BENCHMARK(Magnitude_NestedList_Loop);
// BENCHMARK(Magnitude_NestedList_NoLoop);
// BENCHMARK(Magnitude_ThreeList);

inline double triplelist_magnitude(double *x, int i3, int j3)
{
  return (x[i3] - x[j3]) * (x[i3] - x[j3]) + (x[i3 + 1] - x[j3 + 1]) * (x[i3 + 1] - x[j3 + 1]) + (x[i3 + 2] - x[j3 + 2]) * (x[i3 + 2] - x[j3 + 2]);
}

static void CollisionDetection_TripleList(benchmark::State &state)
{
  int NumberOfBodies = 1337 * 2;
  double *x = new double[NumberOfBodies * 3];
  double *v = new double[NumberOfBodies * 3];
  double *mass = new double[NumberOfBodies];
  double C = 1;

  for (auto _ : state)
  {

    for (int i = 0; i < NumberOfBodies; i++)
    {
#pragma openmp simd
      for (int j = i + 1; j < NumberOfBodies; j++)
      {
        if (triplelist_magnitude(x, i * 3, j * 3) <= (mass[i] * mass[j]))
        {
          v[i] = v[j];
        }
        else
        {
          v[j] = v[i];
        }
      }

      benchmark::DoNotOptimize(v);
    }
  }
}

static void CollisionDetection_TripleList_Aux(benchmark::State &state)
{
  int NumberOfBodies = 1337 * 2;
  double *x = new double[NumberOfBodies * 3];
  double *v = new double[NumberOfBodies * 3];
  double *mass = new double[NumberOfBodies];
  double C = 1;

  for (auto _ : state)
  {
    bool *collisions = new bool[NumberOfBodies];

    for (int i = 0; i < NumberOfBodies; i++)
    {

#pragma openmp simd
      for (int j = i + 1; j < NumberOfBodies; j++)
      {
        int c = triplelist_magnitude(x, i * 3, j * 3) <= (mass[i] * mass[j]);
        v[i] = c * v[j] + (1 - c) * v[i];
      }

      benchmark::DoNotOptimize(v);
    }
  }
}

static void CollisionDetection_Vec3(benchmark::State &state)
{
  int NumberOfBodies = 1337 * 2;
  vec3 *x = new vec3[NumberOfBodies];
  double *v = new double[NumberOfBodies];
  double *mass = new double[NumberOfBodies];

  for (auto _ : state)
  {

    for (int i = 0; i < NumberOfBodies; i++)
    {
      for (int j = i + 1; j < NumberOfBodies; j++)
      {
        if (magnitude_squared(x[j] - x[i]) <= (mass[i] + mass[j]))
        {
          v[i] = v[j];
        }
        else
        {
          v[j] = v[i];
        }
      }
    }
  }
  benchmark::DoNotOptimize(v);
}

// BENCHMARK(CollisionDetection_TripleList);
// BENCHMARK(CollisionDetection_TripleList_Aux);
// BENCHMARK(CollisionDetection_Vec3);

inline double magnitude_squared(double x, double y, double z)
{
  return x * x + y * y + z * z;
}

static void Forces_TripleList_NoMP(benchmark::State &state)
{
  int NumberOfBodies = 1337 * 2;
  double *x = new double[NumberOfBodies * 3];
  double *f = new double[NumberOfBodies * 3];
  double *mass = new double[NumberOfBodies * 3];
  int i3 = 60;

  for (auto _ : state)
  {
    for (int i = 0; i < NumberOfBodies; i++)
    {
      double f_new_x = f[i3];
      double f_new_y = f[i3 + 1];
      double f_new_z = f[i3 + 2];
      double f_x, f_y, f_z;
      for (int j3 = i3 + 3; j3 < NumberOfBodies * 3; j3 += 3)
      {
        // Compute distance components
        double dx = x[j3] - x[i3];
        double dy = x[j3 + 1] - x[i3 + 1];
        double dz = x[j3 + 2] - x[i3 + 2];

        // Compute distance and track max
        double distance2 = magnitude_squared(dx, dy, dz);
        double distance3 = 1.0 / (distance2 * std::sqrt(distance2));

        f_x = dx * distance3;
        f_y = dy * distance3;
        f_z = dz * distance3;

        f_new_x += f_x * mass[j3];
        f_new_y += f_y * mass[j3];
        f_new_z += f_z * mass[j3];

        f[j3] -= f_x * mass[i3];
        f[j3 + 1] -= f_y * mass[i3];
        f[j3 + 2] -= f_z * mass[i3];
      }

      f[i3] = f_new_x;
      f[i3 + 1] = f_new_y;
      f[i3 + 2] = f_new_z;
    }
  }
  benchmark::DoNotOptimize(f);
}

static void Forces_TripleList_Reduction(benchmark::State &state)
{
  int NumberOfBodies = 1337 * 2;
  double *x = new double[NumberOfBodies * 3];
  double *f = new double[NumberOfBodies * 3];
  double *mass = new double[NumberOfBodies * 3];
  int i3 = 60;

  for (auto _ : state)
  {
    for (int i = 0; i < NumberOfBodies; i++)
    {
      double f_new_x = f[i3];
      double f_new_y = f[i3 + 1];
      double f_new_z = f[i3 + 2];
      double f_x, f_y, f_z;

#pragma omp simd reduction(+ \
                           : f_new_x, f_new_y, f_new_z)
      for (int j3 = i3 + 3; j3 < NumberOfBodies * 3; j3 += 3)
      {
        // Compute distance components
        double dx = x[j3] - x[i3];
        double dy = x[j3 + 1] - x[i3 + 1];
        double dz = x[j3 + 2] - x[i3 + 2];

        // Compute distance and track max
        double distance2 = magnitude_squared(dx, dy, dz);
        double distance3 = 1.0 / (distance2 * std::sqrt(distance2));

        f_x = dx * distance3;
        f_y = dy * distance3;
        f_z = dz * distance3;

        f_new_x += f_x * mass[j3];
        f_new_y += f_y * mass[j3];
        f_new_z += f_z * mass[j3];

        f[j3] -= f_x * mass[i3];
        f[j3 + 1] -= f_y * mass[i3];
        f[j3 + 2] -= f_z * mass[i3];
      }

      f[i3] = f_new_x;
      f[i3 + 1] = f_new_y;
      f[i3 + 2] = f_new_z;
    }
  }
  benchmark::DoNotOptimize(f);
}

static void Forces_Vec3(benchmark::State &state)
{
  int NumberOfBodies = 1337 * 2;
  vec3 *x = new vec3[NumberOfBodies];
  vec3 *f = new vec3[NumberOfBodies];
  double *mass = new double[NumberOfBodies];
  int i = 60;

  for (auto _ : state)
  {
    for (int i = 0; i < NumberOfBodies; i++)
    {

      vec3 f_new = f[i];

      for (int j = i + 1; j < NumberOfBodies; j++)
      {
        vec3 distance_vec = x[j] - x[i];
        double distance2 = magnitude_squared(distance_vec);
        double distance3 = 1.0 / (distance2 * std::sqrt(distance2));

        vec3 temp = (distance_vec)*distance3;
        f_new += temp * mass[j];
        f[j] -= temp * mass[i];
      }

      f[i] = f_new;
    }
  }
  benchmark::DoNotOptimize(f);
}

static void Forces_ThreeLists_NoMP(benchmark::State &state)
{
  int NumberOfBodies = 1337 * 2;
  double *xx = new double[NumberOfBodies];
  double *xy = new double[NumberOfBodies];
  double *xz = new double[NumberOfBodies];
  double *fx = new double[NumberOfBodies];
  double *fy = new double[NumberOfBodies];
  double *fz = new double[NumberOfBodies];
  double *mass = new double[NumberOfBodies];

  for (auto _ : state)
  {
    for (int i = 0; i < NumberOfBodies; i++)
    {
      double f_new_x = fx[i];
      double f_new_y = fy[i];
      double f_new_z = fz[i];
      double f_x, f_y, f_z;

      for (int j = i + 1; j < NumberOfBodies; j++)
      {
        double dx = xx[j] - xx[i];
        double dy = xy[j] - xy[i];
        double dz = xz[j] - xz[i];

        // Compute distance and track max
        double distance2 = magnitude_squared(dx, dy, dz);
        double distance3 = 1.0 / (distance2 * std::sqrt(distance2));

        f_x = dx * distance3;
        f_y = dy * distance3;
        f_z = dz * distance3;

        f_new_x += f_x * mass[j];
        f_new_y += f_y * mass[j];
        f_new_z += f_z * mass[j];

        fx[j] -= f_x * mass[i];
        fy[j] -= f_y * mass[i];
        fz[j] -= f_z * mass[i];
      }

      fx[i] = f_new_x;
      fy[i] = f_new_y;
      fz[i] = f_new_z;
    }
  }
  benchmark::DoNotOptimize(fx);
  benchmark::DoNotOptimize(fy);
  benchmark::DoNotOptimize(fz);
}

static void Forces_ThreeLists_Reduction(benchmark::State &state)
{
  int NumberOfBodies = 1337 * 2;
  double *xx = new double[NumberOfBodies];
  double *xy = new double[NumberOfBodies];
  double *xz = new double[NumberOfBodies];
  double *fx = new double[NumberOfBodies];
  double *fy = new double[NumberOfBodies];
  double *fz = new double[NumberOfBodies];
  double *mass = new double[NumberOfBodies];

  for (auto _ : state)
  {
    for (int i = 0; i < NumberOfBodies; i++)
    {
      double f_new_x = fx[i];
      double f_new_y = fy[i];
      double f_new_z = fz[i];
      double f_x, f_y, f_z;
#pragma omp simd reduction(+ \
                           : f_new_x, f_new_y, f_new_z)
      for (int j = i + 1; j < NumberOfBodies; j++)
      {
        double dx = xx[j] - xx[i];
        double dy = xy[j] - xy[i];
        double dz = xz[j] - xz[i];

        // Compute distance and track max
        double distance2 = magnitude_squared(dx, dy, dz);
        double distance3 = 1 / (distance2 * std::sqrt(distance2));

        f_x = dx * distance3;
        f_y = dy * distance3;
        f_z = dz * distance3;

        f_new_x += f_x * mass[j];
        f_new_y += f_y * mass[j];
        f_new_z += f_z * mass[j];

        fx[j] -= f_x * mass[i];
        fy[j] -= f_y * mass[i];
        fz[j] -= f_z * mass[i];
      }

      fx[i] = f_new_x;
      fy[i] = f_new_y;
      fz[i] = f_new_z;
    }
  }
  benchmark::DoNotOptimize(fx);
  benchmark::DoNotOptimize(fy);
  benchmark::DoNotOptimize(fz);
}

static void Forces_ThreeLogicalLists_NoMP(benchmark::State &state)
{
  int NumberOfBodies = 1337 * 2;
  double *x = new double[NumberOfBodies * 3];
  double *f = new double[NumberOfBodies * 3];
  double *mass = new double[NumberOfBodies];
  int i = 60;

  for (auto _ : state)
  {
    for (int i = 0; i < NumberOfBodies; i++)
    {
      double f_new_x = f[i];
      double f_new_y = f[i + NumberOfBodies];
      double f_new_z = f[i + 2 * NumberOfBodies];
      double f_x, f_y, f_z;

      for (int j = i + 1; j < NumberOfBodies; j++)
      {
        double dx = x[j] - x[i];
        double dy = x[j + NumberOfBodies] - x[i + NumberOfBodies];
        double dz = x[j + 2 * NumberOfBodies] - x[i + 2 * NumberOfBodies];

        // Compute distance and track max
        double distance2 = magnitude_squared(dx, dy, dz);
        double distance3 = 1.0 / (distance2 * std::sqrt(distance2));

        f_x = dx * distance3;
        f_y = dy * distance3;
        f_z = dz * distance3;

        f_new_x += f_x * mass[j];
        f_new_y += f_y * mass[j];
        f_new_z += f_z * mass[j];

        f[j] -= f_x * mass[i];
        f[j + NumberOfBodies] -= f_y * mass[i];
        f[j + 2 * NumberOfBodies] -= f_z * mass[i];
      }

      f[i] = f_new_x;
      f[i + NumberOfBodies] = f_new_y;
      f[i + 2 * NumberOfBodies] = f_new_z;
    }
  }
  benchmark::DoNotOptimize(f);
}

static void Forces_ThreeLogicalLists_Reduction(benchmark::State &state)
{
  int NumberOfBodies = 1337 * 2;
  double *x = new double[NumberOfBodies * 3];
  double *f = new double[NumberOfBodies * 3];
  double *mass = new double[NumberOfBodies];
  int i = 60;

  for (auto _ : state)
  {
    for (int i = 0; i < NumberOfBodies; i++)
    {
      double f_new_x = f[i];
      double f_new_y = f[i + NumberOfBodies];
      double f_new_z = f[i + 2 * NumberOfBodies];
      double f_x, f_y, f_z;

#pragma omp simd reduction(+ \
                           : f_new_x, f_new_y, f_new_z)
      for (int j = i + 1; j < NumberOfBodies; j++)
      {
        double dx = x[j] - x[i];
        double dy = x[j + NumberOfBodies] - x[i + NumberOfBodies];
        double dz = x[j + 2 * NumberOfBodies] - x[i + 2 * NumberOfBodies];

        // Compute distance and track max
        double distance2 = magnitude_squared(dx, dy, dz);
        double distance3 = 1.0 / (distance2 * std::sqrt(distance2));

        f_x = dx * distance3;
        f_y = dy * distance3;
        f_z = dz * distance3;

        f_new_x += f_x * mass[j];
        f_new_y += f_y * mass[j];
        f_new_z += f_z * mass[j];

        f[j] -= f_x * mass[i];
        f[j + NumberOfBodies] -= f_y * mass[i];
        f[j + 2 * NumberOfBodies] -= f_z * mass[i];
      }

      f[i] = f_new_x;
      f[i + NumberOfBodies] = f_new_y;
      f[i + 2 * NumberOfBodies] = f_new_z;
    }
  }
  benchmark::DoNotOptimize(f);
}

// BENCHMARK(Forces_TripleList_NoMP);
// BENCHMARK(Forces_TripleList_Reduction);
// BENCHMARK(Forces_ThreeLists_NoMP);
// BENCHMARK(Forces_ThreeLists_Reduction);
// BENCHMARK(Forces_ThreeLogicalLists_NoMP);
// BENCHMARK(Forces_ThreeLogicalLists_Reduction);
// BENCHMARK(Forces_Vec3);

static void TimeStep_ThreeLists_OneLoop(benchmark::State &state)
{
  int NumberOfBodies = 1337 * 2;
  double *xx = new double[NumberOfBodies];
  double *xy = new double[NumberOfBodies];
  double *xz = new double[NumberOfBodies];
  double *fx = new double[NumberOfBodies];
  double *fy = new double[NumberOfBodies];
  double *fz = new double[NumberOfBodies];
  double *vx = new double[NumberOfBodies];
  double *vy = new double[NumberOfBodies];
  double *vz = new double[NumberOfBodies];
  double *mass = new double[NumberOfBodies];

  double timeStep = 0.534;
  for (auto _ : state)
  {
#pragma omp simd
    for (int i = 0; i < NumberOfBodies; i++)
    {
      vx[i] += fx[i] * timeStep;
      vy[i] += fy[i] * timeStep;
      vz[i] += fz[i] * timeStep;

      xx[i] += vx[i] * timeStep;
      xy[i] += vy[i] * timeStep;
      xz[i] += vz[i] * timeStep;
    }
  }

  benchmark::DoNotOptimize(vx);
  benchmark::DoNotOptimize(vy);
  benchmark::DoNotOptimize(vz);
  benchmark::DoNotOptimize(xx);
  benchmark::DoNotOptimize(xy);
  benchmark::DoNotOptimize(xz);
}

static void TimeStep_ThreeLists_TwoLoop(benchmark::State &state)
{
  int NumberOfBodies = 1337 * 2;
  double *xx = new double[NumberOfBodies];
  double *xy = new double[NumberOfBodies];
  double *xz = new double[NumberOfBodies];
  double *fx = new double[NumberOfBodies];
  double *fy = new double[NumberOfBodies];
  double *fz = new double[NumberOfBodies];
  double *vx = new double[NumberOfBodies];
  double *vy = new double[NumberOfBodies];
  double *vz = new double[NumberOfBodies];
  double *mass = new double[NumberOfBodies];

  double timeStep = 0.534;
  for (auto _ : state)
  {
#pragma omp simd
    for (int i = 0; i < NumberOfBodies; i++)
    {
      vx[i] += fx[i] * timeStep;
      vy[i] += fy[i] * timeStep;
      vz[i] += fz[i] * timeStep;
    }
#pragma omp simd
    for (int i = 0; i < NumberOfBodies; i++)
    {
      xx[i] += vx[i] * timeStep;
      xy[i] += vy[i] * timeStep;
      xz[i] += vz[i] * timeStep;
    }
  }

  benchmark::DoNotOptimize(vx);
  benchmark::DoNotOptimize(vy);
  benchmark::DoNotOptimize(vz);
  benchmark::DoNotOptimize(xx);
  benchmark::DoNotOptimize(xy);
  benchmark::DoNotOptimize(xz);
}

static void Max_Builtin(benchmark::State &state)
{
  double *x = new double[9999];
  double m;
  for (auto _ : state)
  {
    m = *std::max_element(x, x + 9999);
    benchmark::DoNotOptimize(m);
  }
}

static void Max_MP(benchmark::State &state)
{
  double *x = new double[9999];
  for (auto _ : state)
  {
    double maxval;
#pragma omp simd reduction(max \
                           : maxval)
    for (int i = 0; i < 9999; i++)
    {
      maxval = maxval > x[i] ? maxval : x[i];
    }
    benchmark::DoNotOptimize(maxval);
  }
}

static void Fill_Builtin(benchmark::State &state)
{
  double *x = new double[9999];
  double *y = new double[9999];
  double *z = new double[9999];
  for (auto _ : state)
  {
    std::fill_n(x, 9999, 0.0);
    std::fill_n(y, 9999, 0.0);
    std::fill_n(z, 9999, 0.0);
    benchmark::DoNotOptimize(x);
    benchmark::DoNotOptimize(y);
    benchmark::DoNotOptimize(z);
  }
}

static void Fill_MP(benchmark::State &state)
{
  double *x = new double[9999];
  double *y = new double[9999];
  double *z = new double[9999];
  for (auto _ : state)
  {
#pragma omp simd
    for (int i = 0; i < 9999; i++)
    {
      x[i] = 0.0;
      y[i] = 0.0;
      z[i] = 0.0;
    }
    benchmark::DoNotOptimize(x);
    benchmark::DoNotOptimize(y);
    benchmark::DoNotOptimize(z);
  }
}

BENCHMARK(Max_Builtin);
BENCHMARK(Max_MP);
BENCHMARK(Fill_Builtin);
BENCHMARK(Fill_MP);

// BENCHMARK(TimeStep_ThreeLists_OneLoop);
// BENCHMARK(TimeStep_ThreeLists_TwoLoop);

BENCHMARK_MAIN();