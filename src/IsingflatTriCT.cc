#include "CtIsing.h"
#include <getopt.h>   // for getopt_long, optarg
#include <cstdlib>    // for atoi, atol, std::stod
#include <iostream>   // for std::cout, std::endl (if printing)
#include <numeric>
#include "statistics.h"
#include <chrono>
#include <fstream>
using namespace std;



//To Do: Checking that the statistics is respected when cutting ladder
//       Add wrappers to prevent contamination from nodes in different lists. 
int main(int argc, char *argv[])
{
    //default parameters
    unsigned int seed = 1234u;
    int Nx = 6;
    int Ny = 6; 
    double k1 = 1.0; 
    double k2 = 1.0; 
    double k3 = 1.0; 
    double invT = 10.0; 
    double t_coup = 1;
    int    n_iter= 1000; 
    int n_therm = 1000;
    int n_skip = 10;
    std::string data_dir = "/Users/jinyunlin/QFE-Research/QFEgeometry/data/IsingCT";



    const struct option long_options[] = {
    {"seed", required_argument, 0, 'S'},
    {"n_iter", required_argument, 0, 'h'},
    {"n_therm", required_argument, 0, 't'},
    {"n_skip", required_argument, 0, 's'},
    { "Nx", required_argument, 0, 'x'},
    { "Ny", required_argument, 0, 'y'},
    { "k1", required_argument, 0, 'a' },
    { "k2", required_argument, 0, 'b' },
    { "k3", required_argument, 0, 'c' },
    {"t_coupling", required_argument, 0, 'f'},
    {"invT", required_argument, 0, 'T'},
    {"data_dir", required_argument, 0, 'd'},
    {0, 0, 0, 0}};
    const char* short_options = "S:h:t:s:x:y:a:b:c:f:T:d:";
 while (true) {

    int o = 0;
    int c = getopt_long(argc, argv, short_options, long_options, &o);
    if (c == -1) break;

    switch (c) {
      case 'S': seed = atol(optarg); break;
      case 'h': n_iter = atoi(optarg); break;
      case 't': n_therm = atoi(optarg); break;
      case 's': n_skip = atoi(optarg); break;
      case 'x': Nx = atoi(optarg); break;
      case 'y': Ny = atoi(optarg); break;
      case 'a': k1 = std::stod(optarg); break;
      case 'b': k2 = std::stod(optarg); break;
      case 'c': k3 = std::stod(optarg); break;
      case 'f': t_coup = std::stod(optarg); break;
      case 'T': invT = std::stod(optarg); break;
      case 'd': data_dir = optarg; break;
      default: break;
    }
  }
std::cout << "Running simulation with parameters:\n";
std::cout << "  Lattice size: " << Nx << " x " << Ny << "\n";
std::cout << "  invT: " << invT << ", t_coup: " << t_coup << "\n";
std::cout << "  k1: " << k1 << ", k2: " << k2 << ", k3: " << k3 << "\n";
std::cout << "  Seed: " << seed << ", Iterations: " << n_iter 
          << ", Therm: " << n_therm << ", Skip: " << n_skip << "\n";

    //Compute Critical Coupling 
    // double l_ratio = k2/k1; 
    // double K1 = 0.5 * asinh(l_ratio);
    // double K2 = 0.5 * asinh(1.0 / l_ratio);

    QfeMeasReal mag;     // magnetization
    QfeMeasReal mag_2;   // magnetization^2
    QfeMeasReal mag_4;   // magnetization^4
    QfeLattice lattice;
    QfeMeasReal clustersize;
    lattice.SeedRng(seed);
    lattice.InitTriangle(Nx, Ny, k1, k2, k3);
    ContinuousTimeLattice field(&lattice, invT, t_coup);
    CtIsing Ising(&field); 
    Ising.HotStart(); 

    int n_threads = omp_get_max_threads();
    std::vector<QfeRng> thread_rngs(n_threads);

    for (int i = 0; i < n_threads; ++i) {
        thread_rngs[i] = QfeRng(seed+ i * 982451653);  // Prime-offset to reduce correlation
    }
    std::vector<double> res(3, 0.0);
    for (int step = 0; step < n_iter + n_therm; ++step) 
    {
        //printf("iteration: %d\n", step);
        Ising.CutLadder(thread_rngs);
        res = Ising.BuildRung(thread_rngs);
        int cluster_num = Ising.CountClusters();
        double av_cluster = Ising.n_sites*invT/(cluster_num);
        Ising.UpdateSpins();
        if (step >= n_therm && (step - n_therm) % n_skip == 0) {
        //printf("entering loop----\n");
        printf("iteration: %d\n", step);
        double m = Ising.computeMeanSpin();
        double m_sq = m * m;
        printf("m: %.12f m^2: %.12f clustersize: %.12f cluster number: %d p_rate: %.12f p_av: %.12f overlap: %.12f\n", m, m_sq, av_cluster, cluster_num, res[0], res[1], res[2]);
        mag.Measure(fabs(m)); 
        mag_2.Measure(m_sq); 
        mag_4.Measure(m_sq*m_sq); 
        clustersize.Measure(av_cluster);
        }
    }

    double m_mean  = mag.Mean();
    double m_err   = mag.Error();
    double m2_mean = mag_2.Mean();
    double m2_err  = mag_2.Error();
    double m4_mean = mag_4.Mean();
    double m4_err  = mag_4.Error();
    double cluster_mean  = clustersize.Mean();
    double cluster_err  = clustersize.Error();
    double U4 = 1-m4_mean/(3*m2_mean*m2_mean);

  printf("spin: %.12e %.12e %.4f %.4f\n", \
      m_mean, m_err, \
      mag.AutocorrFront(), mag.AutocorrBack());
  printf("m^2: %.12e %.12e %.4f %.4f\n", \
      m2_mean, m2_err, \
      mag_2.AutocorrFront(), mag_2.AutocorrBack());
  printf("m^4: %.12e %.12e %.4f %.4f\n", \
      m4_mean, m4_err, \
      mag_4.AutocorrFront(), mag_4.AutocorrBack());
    printf("cluster: %.12e %.12e %.4f %.4f\n", \
      cluster_mean, cluster_err, \
      clustersize.AutocorrFront(), clustersize.AutocorrBack());
  printf("binder cumulant: %.12f\n", U4);   

    //Ising.ExportWorldlinestoCSV("/Users/jinyunlin/QFE-Research/QFEgeometry/data/worldline/spin_segments_output.csv");

    //Ising.CheckRungClusterConsistency();
    // // field.TimeLadders[0].insertAtEnd(8.0);
    // for (int i=0; i<Ising.n_sites; i++)
    // {
    //     std::cout<< i <<std::endl;
    //     field.TimeLadders[i].display();
    // }
       char path[200];
        sprintf(path, "%s/IsingCT.dat", data_dir.c_str());
        FILE* file = fopen(path, "a");
        assert(file!= nullptr);

    fprintf(file, "%d %d %.12f %.12f", Nx, Ny, invT, t_coup);  // [0]
    fprintf(file, " %.12e %.12e %.12e", k1, k2, k3);
    fprintf(file, " %.12e %.12e", m_mean, m_err);  // [11,12]
    fprintf(file, " %.12e %.12e", m2_mean, m2_err); 
    fprintf(file, " %.12e %.12e", m4_mean, m4_err);  
    fprintf(file, " %.12e %.12e", cluster_mean, cluster_err);
    fprintf(file, " %.12e \n", U4);  
    fclose(file);
    return 0; 
}
