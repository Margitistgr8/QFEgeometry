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
    double k1 = 1.0; 
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
    { "k1", required_argument, 0, 'a' },
    {"t_coupling", required_argument, 0, 'b'},
    {"invT", required_argument, 0, 'T'},
    {"data_dir", required_argument, 0, 'd'},
    {0, 0, 0, 0}};
    const char* short_options = "S:h:t:s:x:a:b:T:d:";
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
      case 'a': k1=  std::stod(optarg); break;
      case 'b': t_coup = std::stod(optarg); break;
      case 'T': invT = std::stod(optarg); break;
      case 'd': data_dir = optarg; break;
      default: break;
    }
  }
std::cout << "Running simulation with parameters:\n";
std::cout << "  Lattice size: " << Nx << "\n";
std::cout << "  invT: " << invT << ", t_coup: " << t_coup << "\n";
std::cout << "  k1: " << k1 << "\n";
std::cout << "  Seed: " << seed << ", Iterations: " << n_iter 
          << ", Therm: " << n_therm << ", Skip: " << n_skip << "\n";

    QfeMeasReal mag;     // magnetization
    QfeMeasReal mag_2;   // magnetization^2
    QfeMeasReal mag_4;   // magnetization^4
    QfeLattice lattice;
    QfeMeasReal overlap;
    QfeMeasReal clustersize;
    QfeMeasReal accept_p; 
    std::vector<double> all_overlaps;

    lattice.SeedRng(seed);
    lattice.InitChain(Nx, k1);


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
        //std::vector<double> trial_overlaps = Ising.CollectOverlapHistogramVector();
        //all_overlaps.insert(all_overlaps.end(), trial_overlaps.begin(), trial_overlaps.end());  
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
        accept_p.Measure(res[1]);
        overlap.Measure(res[2]);
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
    printf("acceptance probability: %.12e %.12e %.4f %.4f\n", \
      accept_p.Mean(), accept_p.Error(), \
      accept_p.AutocorrFront(), accept_p.AutocorrBack());
    printf("overlap: %.12e %.12e %.4f %.4f\n", \
      overlap.Mean(), overlap.Error(), \
      overlap.AutocorrFront(), overlap.AutocorrBack());
  printf("binder cumulant: %.12f\n", U4);   


    char worldline[200];
    sprintf(worldline, "%s/1Dspin_segments_output_%d_%.12f_%.12f_%.12f.csv", data_dir.c_str(), Nx, invT, t_coup, k1);

    Ising.ExportWorldlinestoCSV(worldline);

#if 0
    char hist[200];
    sprintf(hist, "%s/1Dspin_segments_histogram_%d_%.12f_%.12f_%.12f.csv", data_dir.c_str(), Nx, invT, t_coup, k1);
    FILE* histfile = fopen(hist, "a");
    assert(histfile!= nullptr);
    for (double val : all_overlaps)
    {fprintf(histfile, "%.12f\n", val);}
#endif
    char path[200];
    sprintf(path, "%s/Ising1DCT.dat", data_dir.c_str());
    FILE* file = fopen(path, "a");
    assert(file!= nullptr);

    fprintf(file, "%d %.12f %.12f %.12e", Nx, invT, t_coup, k1);
    fprintf(file, " %.12e %.12e", m_mean, m_err);  
    fprintf(file, " %.12e %.12e", m2_mean, m2_err); 
    fprintf(file, " %.12e %.12e", m4_mean, m4_err);  
    fprintf(file, " %.12e %.12e", cluster_mean, cluster_err);
    fprintf(file, " %.12e %.12e", accept_p.Mean(), accept_p.Error());
    fprintf(file, " %.12e %.12e", overlap.Mean(), overlap.Error());
    fprintf(file, " %.12e \n", U4);  
    fclose(file);
    return 0; 
}
