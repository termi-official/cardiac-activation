// Author: Francesco Regazzoni - MOX, Politecnico di Milano
// Email:  francesco.regazzoni@polimi.it
// Date:   2020

#include <chrono>
#include <fstream>
#include <iomanip>

#include "sarcomere.hpp"

std::map<std::string, std::vector<double>>
sarcomere::solve(const std::function<double(const double &time)> &Ca,
                 const std::function<double(const double &time)> &SL,
                 const double &Tmax, const double &time_step) {
  double calcium;
  double sarcomere_length;
  double dSL_dt;

  std::map<std::string, std::vector<double>> results;

  double dt_write = 0.01; // [s]
  double last_time_write = 0.0; // [s]

  double time = 0.0; // [s]
  std::vector<double> state(initial_state);
  std::vector<double> prevstate(initial_state);
  std::vector<double> state_at_t(initial_state);

  std::cout << model_name << " model. Computing... " << std::flush;
  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();
  // Time loop
  while (time <= Tmax + 1e-10) {
    // Evaluate inputs
    calcium = Ca(time);
    sarcomere_length = SL(time);
    dSL_dt = (SL(time) - SL(time - time_step)) / time_step;

    // Update state
    for(int i = 0; i < 20; i++) {
      prevstate[i] = state[i];
    }
    solve_time_step(state, calcium, sarcomere_length, dSL_dt, time_step);

    // Store the results
    if (time + time_step >= last_time_write) {
      // Update next write time
      while (time + time_step >= last_time_write) {
        double tprev = time; 
        double tnext = time + time_step;
        double tnow  = last_time_write;
        // Interpolate state linearly
        for(int i = 0; i < 20; i++) {
          state_at_t[i] = prevstate[i] + (tnow-tprev)*(state[i]-prevstate[i])/(tnext-tprev);
        }
        calcium = Ca(tnow);
        sarcomere_length = SL(tnow);
        dSL_dt = (SL(tnow) - SL(tnow - time_step)) / time_step;
        results["time"].push_back(tnow);
        results["Ca"].push_back(calcium);
        results["SL"].push_back(sarcomere_length);
        results["dSL_dt"].push_back(dSL_dt);
        results["Ta"].push_back(get_active_tension(state_at_t, sarcomere_length));
        results["As"].push_back(get_active_stiffness(state_at_t, sarcomere_length));
        for(int i=0;i<20;i++)
          results["S" + std::__cxx11::to_string(i)].push_back(state_at_t[i]);

        last_time_write += dt_write;
      }
  }

    // Increase time
    time += time_step;
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  double duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end - begin)
          .count() *
      1e-6;
  std::cout << "done. Time elapsed: " << duration << " s" << std::endl;

  return results;
}

void sarcomere::write_csv(std::map<std::string, std::vector<double>> results,
                          std::string file_name) {
  std::ofstream csvfile;
  csvfile.open(file_name);
  csvfile << "t,Ca,SL,dSL_dt,Ta,As";
  for(int i=0;i<20;i++)
          csvfile << ",S"+ std::__cxx11::to_string(i);
  csvfile << std::endl;
  csvfile << std::scientific << std::setprecision(8);
  for (unsigned int i = 0; i < results["time"].size(); ++i) {
    csvfile << results["time"][i] << "," << results["Ca"][i] << ","
            << results["SL"][i] << "," << results["dSL_dt"][i] << ","
            << results["Ta"][i] << "," << results["As"][i];
    for(int j=0;j<20;j++)
          csvfile << "," << results["S" + std::__cxx11::to_string(j)][i];
    csvfile << std::endl;
  }

  csvfile.close();
}
